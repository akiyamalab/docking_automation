"""Uni-Dock 2 (v0.6.x) ラッパー: receptor JSON キャッシュ前提の高速 docking。

Uni-Dock 2 は docking 1 回あたり `analyze_receptor_topology` で ~5 分の前処理オーバーヘッド
が発生するが、`engine_checkpoint=true` で保存した `ud2_engine_inputs.json` を
`UnidockProtocolRunner(receptor_file_name=...json, ...)` に渡すと 196× (約 1.6 s) に短縮される。

本クラスはこの事実を前提とし、以下 2 ステージに分離:
  1. `prepare_receptor_cache(protein, grid_box, out_json)` — 前処理を 1 回実行 → JSON 保存
  2. `dock_with_cache(cache_json, compound_set, grid_box, ...)` — JSON からキャッシュ済み docking

通常運用ではステージ 1 を事前に全受容体に対し `OMP_NUM_THREADS=1` で並列実行し
(`scripts/prepare_unidock2_caches.py` 参照)、ステージ 2 をスクリーニング時に使用する。

Requirements:
- conda env: `unidock2` (pip では配布されていない)
- `from unidock_processing.unidocktools.unidock_protocol_runner import UnidockProtocolRunner`
"""
from __future__ import annotations

import json
import os
import shutil
import tempfile
from pathlib import Path
from typing import List, Optional, Tuple

from docking_automation.docking.docking_result import DockingResult
from docking_automation.docking.grid_box import GridBox
from docking_automation.molecule.protein import Protein


class UniDock2Docking:
    """Uni-Dock 2 の 2 ステージ (cache prep + cached docking) インターフェース。

    Uni-Dock 2 の特徴:
    - docking 1 回の 99.8% は `analyze_receptor_topology` が占める (~5 分/receptor)
    - このステップの結果は receptor + grid_box に依存せず、JSON 化可能
    - キャッシュ後の kernel 実行は 10 ligand で 1.6 秒程度

    Note:
        本クラスは既存の `DockingToolABC` は継承しない。Uni-Dock 2 の
        設計は v1 と根本的に異なり (SDF 入力 / DMS 受容体 / 階層 config YAML)、
        共通 ABC に押し込めると各メソッドの意味が歪むため、独立クラスとする。
    """

    def __init__(self, cache_dir: Optional[Path] = None) -> None:
        """
        Args:
            cache_dir: receptor JSON cache の格納ディレクトリ。
                未指定時は各メソッド呼び出し時に指定必須。
        """
        self.cache_dir = Path(cache_dir) if cache_dir else None

    def screen_against_repo(
        self,
        cache_json: Path,
        ligand_paths: List[Path],
        grid_box: GridBox,
        protein_content_hash: str,
        compound_content_hashes: List[str],
        compound_set_id: str,
        repo,
        timeout_sec: float = 600.0,
        max_retries: int = 2,
    ) -> List[DockingResult]:
        """content_hash ベースで HDF5 repo から再利用しつつ、未計算ペアだけ cached docking。

        既存の HDF5DockingResultRepository と組み合わせて Phase 4 規模のスクリーニングで
        冪等性 (途中中断からの再開) を担保する。

        Args:
            cache_json: 受容体 JSON cache (prepare_receptor_cache の戻り値)。
            ligand_paths: SDF パスのリスト (3D 化済)。
            grid_box: GridBox。
            protein_content_hash: Protein.content_hash。
            compound_content_hashes: ligand_paths と同順序の content_hash リスト。
            compound_set_id: CompoundSet 識別子。
            repo: HDF5DockingResultRepository。
            timeout_sec/max_retries: dock_with_cache_robust に渡す。

        Returns:
            再利用 + 新規計算 を合わせた全 DockingResult。

        Side effect:
            repo.save() を新規計算結果に対して呼ぶ。既に保存済みのペアは再保存しない。
        """
        n = len(ligand_paths)
        assert len(compound_content_hashes) == n, 'ligand paths と hashes は同じ長さ必要'

        reused_results: List[DockingResult] = []
        new_indices: List[int] = []
        for i, ch in enumerate(compound_content_hashes):
            if repo._exists(protein_content_hash, ch):
                try:
                    existing = repo.load_by_hashes(protein_content_hash, ch)
                    if existing is not None:
                        reused_results.append(existing)
                        continue
                except Exception:
                    pass
            new_indices.append(i)

        if not new_indices:
            return reused_results

        new_ligands = [ligand_paths[i] for i in new_indices]
        new_hashes = [compound_content_hashes[i] for i in new_indices]
        new_results = self.dock_with_cache_robust(
            cache_json=cache_json,
            ligand_sdf_list=new_ligands,
            grid_box=grid_box,
            protein_content_hash=protein_content_hash,
            compound_content_hashes=new_hashes,
            timeout_sec=timeout_sec,
            max_retries=max_retries,
        )

        # compound_set_id, compound_index を再配置して save
        # dock_with_cache_robust は compound_index に ligand_sdf_list の絶対 index を入れるので、
        # ここで元の CompoundSet の index に戻す。
        for r in new_results:
            r.compound_set_id = compound_set_id
            lig_local_idx = r.compound_index
            if 0 <= lig_local_idx < len(new_indices):
                r.compound_index = new_indices[lig_local_idx]
            repo.save(r)

        return reused_results + new_results

    def cache_path_for(self, protein: Protein) -> Path:
        """content_hash ベースの cache file path を返す。

        Phase 4 のような大規模運用では protein_content_hash をキーにすれば
        同一内容の receptor が複数 UniProt ID で現れても 1 回の prep で済む。
        """
        if self.cache_dir is None:
            raise ValueError("cache_dir が未設定です。__init__ で指定してください。")
        return self.cache_dir / f"{protein.content_hash}.json"

    def prepare_receptor_cache(
        self,
        protein: Protein,
        grid_box: GridBox,
        out_json: Path,
        force: bool = False,
    ) -> Path:
        """受容体を 1 回だけ Uni-Dock 2 で前処理し、結果を JSON キャッシュに保存する。

        ~5 分の重い処理なので、複数 receptor の場合は
        `scripts/prepare_unidock2_caches.py` を使って並列実行すべし。

        Args:
            protein: 前処理する受容体 (PDB ベース)。
            grid_box: ドッキング box。cache 生成のトリガーとして 10 ligand dock を 1 回実行するため必要。
            out_json: 出力 JSON パス。既存かつ force=False なら何もしない。
            force: True なら既存キャッシュを上書き。

        Returns:
            保存された JSON cache のパス (= out_json)。

        Side effect:
            `OMP_NUM_THREADS` 等を 1 に設定する (複数 proc 並列化のため)。
        """
        out_json = Path(out_json)
        if out_json.exists() and not force:
            return out_json

        self._enforce_omp_single_thread()

        # Uni-Dock 2 の docking 呼び出しは必ず ligand が必要なため、最小 1 ligand を用意。
        # dummy SDF (ligand 内容は cache には含まれない ← receptor 部分のみ抽出するため OK)
        from unidock_processing.unidocktools.unidock_protocol_runner import (
            UnidockProtocolRunner,
        )

        with tempfile.TemporaryDirectory(prefix='ud2_prep_') as tmp_str:
            tmp_dir = Path(tmp_str)
            dummy_sdf = self._write_dummy_sdf(tmp_dir / 'dummy.sdf')

            workdir = tmp_dir / 'wd'
            workdir.mkdir()
            runner = UnidockProtocolRunner(
                receptor_file_name=str(self._receptor_path_or_dms(protein, tmp_dir)),
                ligand_sdf_file_name_list=[str(dummy_sdf)],
                target_center=tuple(grid_box.center),
                working_dir_name=str(workdir),
                docking_pose_sdf_file_name=str(tmp_dir / 'dummy_pose.sdf'),
                engine_checkpoint=True,
            )
            runner.run_unidock_protocol()

            # ud2_engine_inputs.json を探す → receptor キーだけ取り出す
            ckpt = next(workdir.rglob('ud2_engine_inputs.json'), None)
            if ckpt is None:
                raise RuntimeError(
                    f'ud2_engine_inputs.json not generated under {workdir}. '
                    'engine_checkpoint pipeline failed.'
                )
            with ckpt.open() as f:
                data = json.load(f)
            if 'receptor' not in data:
                raise RuntimeError("cache JSON is missing 'receptor' key.")

            out_json.parent.mkdir(parents=True, exist_ok=True)
            with out_json.open('w') as f:
                json.dump({'receptor': data['receptor']}, f)

        return out_json

    def dock_with_cache_robust(
        self,
        cache_json: Path,
        ligand_sdf_list: List[Path],
        grid_box: GridBox,
        protein_content_hash: str,
        compound_content_hashes: Optional[List[str]] = None,
        working_dir: Optional[Path] = None,
        docking_pose_sdf: Optional[Path] = None,
        timeout_sec: float = 600.0,
        max_retries: int = 2,
    ) -> List[DockingResult]:
        """subprocess + timeout + retry ラッパー。内部デッドロック対策。

        `UnidockProtocolRunner` は内部で pathos 経由の multiprocessing プールを
        使うが、N=16 等の高並列時に稀に futex デッドロックが観察された。
        本メソッドは dock_with_cache を独立プロセスで実行し、timeout 超過時に
        process group ごと kill して retry する。

        Args:
            timeout_sec: 1 回あたりの上限時間。超過すると kill + retry。
            max_retries: タイムアウト時の再試行回数 (0 なら 1 回だけ実行)。
            他の引数は dock_with_cache と同じ。

        Raises:
            TimeoutError: retry 含め全試行でタイムアウト。
        """
        import json as _json
        import shlex
        import signal
        import subprocess
        import sys

        if working_dir is None:
            working_dir_ctx = tempfile.TemporaryDirectory(prefix='ud2_robust_')
            working_dir = Path(working_dir_ctx.name)
        else:
            working_dir_ctx = None
            working_dir = Path(working_dir)
            working_dir.mkdir(parents=True, exist_ok=True)

        if docking_pose_sdf is None:
            docking_pose_sdf = working_dir / 'pose.sdf'

        try:
            args_file = working_dir / '_worker_args.json'
            args_file.write_text(_json.dumps({
                'cache_json': str(cache_json),
                'ligand_sdf_list': [str(p) for p in ligand_sdf_list],
                'grid_center': list(grid_box.center),
                'grid_size': list(grid_box.size),
                'working_dir': str(working_dir),
                'docking_pose_sdf': str(docking_pose_sdf),
            }))

            # `-m` だと package __init__ が openbabel を早期 load するため、環境によっては
            # libstdc++ と msys 拡張の ABI が衝突する。worker はスタンドアロン設計なので
            # 直接ファイルパスで呼び出す。
            worker = Path(__file__).parent / '_unidock2_worker.py'
            cmd = [sys.executable, str(worker), str(args_file)]

            last_exc: Optional[BaseException] = None
            for attempt in range(max_retries + 1):
                try:
                    subprocess.run(
                        cmd,
                        timeout=timeout_sec,
                        start_new_session=True,  # os.killpg() で確実に全子孫を終了させるため
                        check=True,
                        capture_output=True,
                    )
                    break
                except subprocess.TimeoutExpired as e:
                    last_exc = e
                    self._kill_session_from(e.cmd)
                    continue
                except subprocess.CalledProcessError as e:
                    last_exc = e
                    # 異常終了は内部エラー。retry 対象に含める。
                    continue
            else:
                raise TimeoutError(
                    f'dock_with_cache_robust: all {max_retries + 1} attempts failed '
                    f'(last: {type(last_exc).__name__})'
                ) from last_exc

            return self._parse_pose_sdf(
                Path(docking_pose_sdf),
                ligand_sdf_list,
                protein_content_hash,
                compound_content_hashes,
            )
        finally:
            if working_dir_ctx is not None:
                working_dir_ctx.cleanup()

    @staticmethod
    def _kill_session_from(cmd) -> None:
        """subprocess.TimeoutExpired 時に session 単位で強制終了。

        start_new_session=True で起動した child は独立 session。
        孫プロセス (pathos プール等) も同じ session に属するので、
        session leader に SIGKILL → setsid 配下全部が終了する。
        """
        import os
        import signal
        import subprocess as sp
        # TimeoutExpired は Popen を自動で kill するが、孫プロセスは別 PID group にある
        # かもしれないので pgrep で残存を拾って SIGKILL する。
        try:
            pids = sp.check_output(['pgrep', '-f', '_unidock2_worker']).decode().split()
            for pid in pids:
                try:
                    os.kill(int(pid), signal.SIGKILL)
                except (ProcessLookupError, ValueError):
                    pass
        except sp.CalledProcessError:
            pass

    def dock_with_cache(
        self,
        cache_json: Path,
        ligand_sdf_list: List[Path],
        grid_box: GridBox,
        protein_content_hash: str,
        compound_content_hashes: Optional[List[str]] = None,
        working_dir: Optional[Path] = None,
        docking_pose_sdf: Optional[Path] = None,
    ) -> List[DockingResult]:
        """キャッシュ済み receptor JSON を用いて複数 ligand を一括 docking。

        Args:
            cache_json: `prepare_receptor_cache` で生成した JSON。
            ligand_sdf_list: 既に 3D 化済み SDF のパスリスト。
            grid_box: docking box (center と size)。
            protein_content_hash: 結果の `DockingResult.protein_content_hash` に使う。
            compound_content_hashes: None の場合は SDF stem を使用。
            working_dir: Uni-Dock 2 が中間生成物を置くディレクトリ。None → tempdir (削除)。
            docking_pose_sdf: 統合ポーズ SDF 出力先。None → tempdir 内の `pose.sdf`。

        Returns:
            各 ligand の best pose 1 つだけを含む DockingResult リスト
            (Uni-Dock 2 は 1 ligand あたり num_pose 個返すが、ここでは best 1 件のみ採用)。
        """
        self._enforce_omp_single_thread()

        from unidock_processing.unidocktools.unidock_protocol_runner import (
            UnidockProtocolRunner,
        )

        tmp_ctx = None
        if working_dir is None:
            tmp_ctx = tempfile.TemporaryDirectory(prefix='ud2_cached_')
            working_dir = Path(tmp_ctx.name)
        else:
            working_dir = Path(working_dir)
            working_dir.mkdir(parents=True, exist_ok=True)

        if docking_pose_sdf is None:
            docking_pose_sdf = working_dir / 'pose.sdf'

        try:
            runner = UnidockProtocolRunner(
                receptor_file_name=str(cache_json),  # .json 拡張子で analyze_receptor_topology を bypass
                ligand_sdf_file_name_list=[str(p) for p in ligand_sdf_list],
                target_center=tuple(grid_box.center),
                working_dir_name=str(working_dir),
                docking_pose_sdf_file_name=str(docking_pose_sdf),
            )
            runner.run_unidock_protocol()

            return self._parse_pose_sdf(
                docking_pose_sdf,
                ligand_sdf_list,
                protein_content_hash,
                compound_content_hashes,
            )
        finally:
            if tmp_ctx is not None:
                tmp_ctx.cleanup()

    # ----- internals -----
    @staticmethod
    def _enforce_omp_single_thread() -> None:
        """複数 proc 並列化時の OpenMP オーバーサブスクリプションを防ぐ。

        OMP=24 (default) だと 2 proc 並列時に 48 threads が 24 cores を奪い合い、
        計測上は非対称な tail latency が発生する。OMP=1 に固定すると並列効率 ~99%。
        """
        for k in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS'):
            os.environ.setdefault(k, '1')

    @staticmethod
    def _receptor_path_or_dms(protein: Protein, tmp_dir: Path) -> Path:
        """DMS があればそれを使い、なければ PDB をそのまま返す (Uni-Dock 2 は両方対応)。"""
        # Uni-Dock 2 は PDB 入力時に内部で DMS 化するので、PDB のままでも動作する。
        # 既に DMS 化済みなら短縮化のため優先。
        dms_candidate = Path(str(protein.path).rsplit('.', 1)[0] + '.dms')
        if dms_candidate.exists():
            return dms_candidate
        return Path(protein.path)

    @staticmethod
    def _write_dummy_sdf(path: Path) -> Path:
        """ベンゼン 1 分子の SDF を書き出す。cache prep トリガー用。"""
        from rdkit import Chem
        from rdkit.Chem import AllChem

        m = Chem.MolFromSmiles('c1ccccc1')
        m = Chem.AddHs(m)
        AllChem.EmbedMolecule(m, randomSeed=0)
        AllChem.MMFFOptimizeMolecule(m)
        m.SetProp('_Name', 'benzene_dummy')
        w = Chem.SDWriter(str(path))
        w.write(m)
        w.close()
        return path

    @staticmethod
    def _parse_pose_sdf(
        pose_sdf: Path,
        ligand_sdf_list: List[Path],
        protein_content_hash: str,
        compound_content_hashes: Optional[List[str]],
    ) -> List[DockingResult]:
        """Uni-Dock 2 の統合 pose SDF を ligand 毎の DockingResult に分解する。

        Uni-Dock 2 の pose SDF は 1 ligand × N pose のフラット列。
        `ud2_molecule_name` (例: `MOL_0_unidock2_pose_0`) で ligand を識別する。
        同一 ligand の最良スコア (pose_0) のみ採用。
        """
        from rdkit import Chem

        if not pose_sdf.exists():
            return []

        sup = Chem.SDMolSupplier(str(pose_sdf), removeHs=False)
        best_by_mol: dict = {}
        for mol in sup:
            if mol is None:
                continue
            name = mol.GetProp('ud2_molecule_name') if mol.HasProp('ud2_molecule_name') else ''
            if not name:
                continue
            # name 例: "MOL_3_unidock2_pose_0"
            # pose_0 が best (ranked)、他は破棄
            parts = name.rsplit('_', 1)
            if len(parts) != 2 or not parts[1].isdigit():
                continue
            mol_key, pose_idx = parts[0], int(parts[1])
            if pose_idx != 0:
                continue
            score = float(mol.GetProp('vina_binding_free_energy')) if mol.HasProp('vina_binding_free_energy') else None
            # ligand index を MOL_<n>_unidock2_pose から抽出
            # "MOL_3_unidock2" → 3
            try:
                lig_idx = int(mol_key.split('_')[1])
            except (ValueError, IndexError):
                continue
            best_by_mol[lig_idx] = (score, mol)

        results: List[DockingResult] = []
        for lig_idx in sorted(best_by_mol.keys()):
            if lig_idx >= len(ligand_sdf_list):
                continue
            score, mol = best_by_mol[lig_idx]
            lig_path = ligand_sdf_list[lig_idx]
            compound_hash = (
                compound_content_hashes[lig_idx] if compound_content_hashes is not None
                else lig_path.stem
            )

            # 単一ポーズ SDF に書き出し DockingResult.result_path にする
            tmp_sdf = Path(tempfile.mkstemp(suffix='.sdf')[1])
            w = Chem.SDWriter(str(tmp_sdf))
            w.write(mol)
            w.close()

            results.append(
                DockingResult(
                    result_path=tmp_sdf,
                    protein_id='',
                    compound_set_id=lig_path.parent.name,
                    compound_index=lig_idx,
                    docking_score=score,
                    protein_content_hash=protein_content_hash,
                    compound_content_hash=compound_hash,
                    compoundset_content_hash=compound_hash,
                    metadata={'tool': 'Uni-Dock 2', 'source': 'unidock2_cached'},
                )
            )

        return results
