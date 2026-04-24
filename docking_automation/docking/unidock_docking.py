from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path
from typing import List, Optional, Tuple

from docking_automation.converters.molecule_converter import MoleculeConverter
from docking_automation.docking.screening_tool import ScreeningTool
from docking_automation.docking.docking_parameters import DockingParameters, UniDockParameters
from docking_automation.docking.docking_result import DockingResult
from docking_automation.docking.preprocessed_compound_set import PreprocessedCompoundSet
from docking_automation.docking.preprocessed_protein import PreprocessedProtein
from docking_automation.molecule.compound_set import CompoundSet
from docking_automation.molecule.protein import Protein


class UniDockDocking(ScreeningTool):
    """Uni-Dock v1.1.0 をサブプロセス呼び出しで利用する DockingTool 実装。

    1 run = 1 protein × バッチ ligand（--ligand_index で列挙）。
    AutoDockVina と同一の前処理（obabel + Meeko）を使用する。
    """

    UNIDOCK_BINARY: str = "unidock"

    def __init__(self, binary_path: Optional[str] = None) -> None:
        self.binary = binary_path or self.UNIDOCK_BINARY
        self.converter = MoleculeConverter()

    def _preprocess_protein(self, protein: Protein) -> PreprocessedProtein:
        temp_dir = Path(tempfile.mkdtemp())
        pdbqt_path = temp_dir / f"{protein.id}.pdbqt"
        self.converter.protein_to_pdbqt(protein, pdbqt_path)
        return PreprocessedProtein(file_path=pdbqt_path)

    def _preprocess_compound_set(self, compound_set: CompoundSet) -> PreprocessedCompoundSet:
        temp_dir = Path(tempfile.mkdtemp())
        pdbqt_paths = self.converter.compound_to_pdbqt(compound_set, temp_dir)
        compound_hash_cache = getattr(compound_set, "_CompoundSet__compound_hash_cache").copy()
        return PreprocessedCompoundSet(file_paths=pdbqt_paths, compound_hash_cache=compound_hash_cache)

    def dock(self, parameters: DockingParameters) -> List[DockingResult]:
        common = parameters.common
        specific = parameters.specific

        if not isinstance(specific, UniDockParameters):
            raise ValueError("specific は UniDockParameters のインスタンスである必要があります。")

        protein = common.protein
        compound_set = common.compound_set
        grid_box = common.grid_box

        all_indices = list(range(len(compound_set.file_paths)))
        results, failed_indices = self._run_batch(
            protein, compound_set, all_indices, grid_box, specific, specific.search_mode
        )

        if specific.rescue_mode and failed_indices:
            rescue_results, _ = self._run_batch(
                protein, compound_set, failed_indices, grid_box, specific, specific.rescue_search_mode
            )
            results.extend(rescue_results)

        return results

    def _run_batch(
        self,
        protein: PreprocessedProtein,
        compound_set: PreprocessedCompoundSet,
        ligand_indices: List[int],
        grid_box,
        params: UniDockParameters,
        search_mode: str,
    ) -> Tuple[List[DockingResult], List[int]]:
        ligand_paths = [compound_set.file_paths[i] for i in ligand_indices]

        if protein.file_path is None:
            raise ValueError("PreprocessedProtein.file_path が未設定のため UniDock を実行できません")

        with tempfile.TemporaryDirectory() as tmpdir_str:
            tmpdir = Path(tmpdir_str)

            ligand_index = tmpdir / "ligands.txt"
            ligand_index.write_text("\n".join(str(f) for f in ligand_paths))

            out_dir = tmpdir / "output"
            out_dir.mkdir()

            cmd = self._build_cli_command(
                protein.file_path, ligand_index, grid_box, out_dir, params, search_mode=search_mode
            )
            proc_result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)

            return self._parse_output(protein, compound_set, ligand_indices, out_dir, proc_result, params)

    def _build_cli_command(
        self,
        receptor: Path,
        ligand_index: Path,
        grid_box,
        out_dir: Path,
        params: UniDockParameters,
        search_mode: Optional[str] = None,
    ) -> List[str]:
        center = grid_box.center
        size = grid_box.size
        _search_mode = search_mode if search_mode is not None else params.search_mode
        return [
            self.binary,
            "--receptor", str(receptor),
            "--ligand_index", str(ligand_index),
            "--center_x", str(center[0]),
            "--center_y", str(center[1]),
            "--center_z", str(center[2]),
            "--size_x", str(size[0]),
            "--size_y", str(size[1]),
            "--size_z", str(size[2]),
            "--scoring", params.scoring,
            "--search_mode", _search_mode,
            "--num_modes", str(params.num_modes),
            "--seed", str(params.seed),
            "--dir", str(out_dir),
            "--verbosity", str(params.verbosity),
        ]

    def _parse_output(
        self,
        protein: PreprocessedProtein,
        compound_set: PreprocessedCompoundSet,
        ligand_indices: List[int],
        out_dir: Path,
        proc_result: subprocess.CompletedProcess,
        params: UniDockParameters,
    ) -> Tuple[List[DockingResult], List[int]]:
        results = []
        failed_indices = []
        for orig_idx in ligand_indices:
            ligand_path = compound_set.file_paths[orig_idx]
            stem = ligand_path.stem
            out_pdbqt = out_dir / f"{stem}_out.pdbqt"

            if not out_pdbqt.exists():
                failed_indices.append(orig_idx)
                continue

            score = self._extract_score(out_pdbqt)
            if score is not None:
                if score >= params.score_threshold_max or score <= params.score_threshold_min:
                    score = None
            if score is None:
                failed_indices.append(orig_idx)
                continue

            sdf_path = out_pdbqt.with_suffix(".sdf")
            try:
                self.converter.pdbqt_to_sdf(out_pdbqt, sdf_path)
                result_path = sdf_path
            except Exception:
                result_path = out_pdbqt

            compound_hash = compound_set.get_compound_hash(orig_idx)
            compound_set_id = stem.rsplit("_", 1)[0] if "_" in stem else stem

            assert protein.file_path is not None  # 上位メソッドでガード済み
            results.append(DockingResult(
                result_path=result_path,
                protein_id=protein.file_path.stem,
                compound_set_id=compound_set_id,
                compound_index=orig_idx,
                docking_score=score,
                protein_content_hash=protein.content_hash,
                compound_content_hash=compound_hash,
                metadata={"tool": "Uni-Dock", "source": "unidock"},
            ))

        return results, failed_indices

    def _extract_score(self, pdbqt_path: Path) -> Optional[float]:
        for line in pdbqt_path.read_text().splitlines():
            if "VINA RESULT" in line:
                parts = line.split()
                for i, p in enumerate(parts):
                    if p.rstrip(":") == "RESULT" and i + 1 < len(parts):
                        try:
                            return float(parts[i + 1])
                        except ValueError:
                            pass
        return None

    # --- ScreeningTool (旧 CacheableReceptorDocking) 準拠 ---
    # Uni-Dock v1 CLI は --write_maps / --maps で Vina 互換の map cache を
    # 生成・読込できる。v2 のような根本的な高速化ではないが、受容体毎の
    # map 計算 (~数秒) を N×M で繰り返さずに済む。
    # v2 への移行が推奨だが、Phase 3 production path の互換性保持のため対応。

    def prepare_receptor_cache(
        self,
        protein: Protein,
        grid_box,
        out_cache: Path,
        force: bool = False,
    ) -> Path:
        """Uni-Dock v1 CLI `--write_maps` で map 群を生成してキャッシュする。

        Args:
            protein: 受容体。PDBQT 化されていない場合は内部で変換。
            grid_box: ドッキングボックス。
            out_cache: 出力パス prefix。`{out_cache}.C.map` 等が生成される。
            force: True なら既存 .map 群を無視して再計算。

        Returns:
            out_cache (= .map 群の prefix)。dock_with_cache にそのまま渡せる。

        Note:
            v1 の map 計算は内部で GPU を使うため、ホストに CUDA デバイスが必須。
        """
        out_cache = Path(out_cache)
        out_cache.parent.mkdir(parents=True, exist_ok=True)
        existing = list(out_cache.parent.glob(f'{out_cache.name}.*.map'))
        if existing and not force:
            return out_cache

        # PDBQT 化 (未処理ならここで)
        preprocessed = self._preprocess_protein(protein)
        if preprocessed.file_path is None:
            raise ValueError('PreprocessedProtein.file_path が未設定のため map を生成できません')

        center = grid_box.center
        size = grid_box.size

        # Uni-Dock v1 は --write_maps に ligand と組み合わせて起動する必要あり
        # (ligand 無しだと docking 工程がないためエラー)。ダミー ligand で呼ぶ。
        with tempfile.TemporaryDirectory() as tmp_str:
            tmp = Path(tmp_str)
            dummy_ligand = self._write_dummy_pdbqt_ligand(tmp / 'dummy.pdbqt')
            ligand_index = tmp / 'ligands.txt'
            ligand_index.write_text(str(dummy_ligand))
            cmd = [
                self.binary,
                '--receptor', str(preprocessed.file_path),
                '--ligand_index', str(ligand_index),
                '--center_x', str(center[0]),
                '--center_y', str(center[1]),
                '--center_z', str(center[2]),
                '--size_x', str(size[0]),
                '--size_y', str(size[1]),
                '--size_z', str(size[2]),
                '--write_maps', str(out_cache),
                '--dir', str(tmp / 'out'),
                '--num_modes', '1',
                '--exhaustiveness', '8',
                '--scoring', 'vina',
            ]
            subprocess.run(cmd, capture_output=True, text=True, timeout=600)

        return out_cache

    def dock_with_cache(
        self,
        cache: Path,
        ligand_paths: List[Path],
        grid_box,
        protein_content_hash: str,
        compound_content_hashes: Optional[List[str]] = None,
        exhaustiveness: int = 8,
        num_modes: int = 1,
        scoring: str = 'vina',
    ) -> List[DockingResult]:
        """キャッシュ済み map を `--maps` で読み込み、ligand_paths を一括 docking。

        Args:
            cache: `prepare_receptor_cache` が返した prefix。
            ligand_paths: PDBQT 形式 ligand リスト。
            grid_box: grid (center と size は map と一致必須)。
            protein_content_hash: DockingResult.protein_content_hash に格納。
            compound_content_hashes: None の場合はファイル stem。
            exhaustiveness, num_modes, scoring: Uni-Dock CLI オプション。

        Returns:
            各 ligand の best pose DockingResult リスト。
        """
        center = grid_box.center
        size = grid_box.size
        results: List[DockingResult] = []
        with tempfile.TemporaryDirectory() as tmp_str:
            tmp = Path(tmp_str)
            out_dir = tmp / 'out'
            out_dir.mkdir()
            ligand_index = tmp / 'ligands.txt'
            ligand_index.write_text('\n'.join(str(p) for p in ligand_paths))

            cmd = [
                self.binary,
                '--maps', str(cache),
                '--ligand_index', str(ligand_index),
                '--center_x', str(center[0]),
                '--center_y', str(center[1]),
                '--center_z', str(center[2]),
                '--size_x', str(size[0]),
                '--size_y', str(size[1]),
                '--size_z', str(size[2]),
                '--dir', str(out_dir),
                '--num_modes', str(num_modes),
                '--exhaustiveness', str(exhaustiveness),
                '--scoring', scoring,
            ]
            subprocess.run(cmd, capture_output=True, text=True, timeout=600)

            for idx, lig_path in enumerate(ligand_paths):
                out_pdbqt = out_dir / f'{lig_path.stem}_out.pdbqt'
                if not out_pdbqt.exists():
                    continue
                score = self._extract_score(out_pdbqt)
                if score is None:
                    continue
                ch = (
                    compound_content_hashes[idx]
                    if compound_content_hashes is not None
                    else lig_path.stem
                )
                results.append(DockingResult(
                    result_path=out_pdbqt,
                    protein_id='',
                    compound_set_id=lig_path.parent.name,
                    compound_index=idx,
                    docking_score=score,
                    protein_content_hash=protein_content_hash,
                    compound_content_hash=ch,
                    metadata={'tool': 'Uni-Dock', 'source': 'unidock_cached_maps'},
                ))
        return results

    @staticmethod
    def _write_dummy_pdbqt_ligand(path: Path) -> Path:
        """cache prep 起動用の極小 PDBQT ligand (atom 1 個)。

        Uni-Dock v1 は --write_maps でも ligand 引数が必須のため、軽量なダミーを
        渡す。map 生成自体は ligand に依存しないため任意の PDBQT で OK。
        """
        path.write_text(
            'REMARK  Dummy ligand for map generation\n'
            'ROOT\n'
            'ATOM      1  C   UNL     1       0.000   0.000   0.000  1.00  0.00     0.000 C\n'
            'ENDROOT\n'
            'TORSDOF 0\n'
        )
        return path
