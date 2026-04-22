from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path
from typing import List, Optional

from docking_automation.converters.molecule_converter import MoleculeConverter
from docking_automation.docking.docking import DockingToolABC
from docking_automation.docking.docking_parameters import DockingParameters, UniDockParameters
from docking_automation.docking.docking_result import DockingResult
from docking_automation.docking.preprocessed_compound_set import PreprocessedCompoundSet
from docking_automation.docking.preprocessed_protein import PreprocessedProtein
from docking_automation.molecule.compound_set import CompoundSet
from docking_automation.molecule.protein import Protein


class UniDockDocking(DockingToolABC):
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

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            ligand_index = tmpdir / "ligands.txt"
            ligand_index.write_text("\n".join(str(f) for f in compound_set.file_paths))

            out_dir = tmpdir / "output"
            out_dir.mkdir()

            cmd = self._build_cli_command(
                protein.file_path, ligand_index, grid_box, out_dir, specific
            )
            proc_result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)

            return self._parse_output(protein, compound_set, out_dir, proc_result, specific)

    def _build_cli_command(
        self,
        receptor: Path,
        ligand_index: Path,
        grid_box,
        out_dir: Path,
        params: UniDockParameters,
    ) -> List[str]:
        center = grid_box.center
        size = grid_box.size
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
            "--search_mode", params.search_mode,
            "--num_modes", str(params.num_modes),
            "--seed", str(params.seed),
            "--dir", str(out_dir),
            "--verbosity", str(params.verbosity),
        ]

    def _parse_output(
        self,
        protein: PreprocessedProtein,
        compound_set: PreprocessedCompoundSet,
        out_dir: Path,
        proc_result: subprocess.CompletedProcess,
        params: UniDockParameters,
    ) -> List[DockingResult]:
        results = []
        for idx, ligand_path in enumerate(compound_set.file_paths):
            stem = ligand_path.stem
            out_pdbqt = out_dir / f"{stem}_out.pdbqt"

            if not out_pdbqt.exists():
                continue

            score = self._extract_score(out_pdbqt)
            if score is not None:
                if score >= params.score_threshold_max or score <= params.score_threshold_min:
                    score = None
            if score is None:
                continue

            sdf_path = out_pdbqt.with_suffix(".sdf")
            try:
                self.converter.pdbqt_to_sdf(out_pdbqt, sdf_path)
                result_path = sdf_path
            except Exception:
                result_path = out_pdbqt

            compound_hash = compound_set.get_compound_hash(idx)
            compound_set_id = stem.rsplit("_", 1)[0] if "_" in stem else stem

            results.append(DockingResult(
                result_path=result_path,
                protein_id=protein.file_path.stem,
                compound_set_id=compound_set_id,
                compound_index=idx,
                docking_score=score,
                protein_content_hash=protein.content_hash,
                compound_content_hash=compound_hash,
                metadata={"tool": "Uni-Dock", "source": "unidock"},
            ))

        return results

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
