"""化合物前処理パイプライン（決定論的）

仕様書: context/docking_automation_compound_pipeline.md

pipeline:
    SMILES → Dimorphite-DL 2.0.2 (protonation, pH 7.4)
           → RDKit TautomerEnumerator.Canonicalize()
           → RDKit ETKDGv3 (seed=42)
           → MMFF94 最適化 (maxIters=2000)
           → Chem.AddHs(addCoords=True)
           → SDF 書き出し
"""
from __future__ import annotations

import hashlib
import random
from io import StringIO
from pathlib import Path
from typing import Optional

import numpy as np

random.seed(42)
np.random.seed(42)


def preprocess_smiles(
    smiles: str,
    seed: int = 42,
    ph_min: float = 6.4,
    ph_max: float = 8.4,
    max_variants: int = 1,
) -> Optional[str]:
    """SMILES → 前処理済み SDF ブロック（文字列）を返す。

    失敗時は None を返す（例外を投げない）。

    Args:
        smiles: 入力 SMILES 文字列
        seed: ETKDGv3 の randomSeed（デフォルト 42。殿ご裁可必須値）
        ph_min / ph_max: Dimorphite-DL の pH 範囲
        max_variants: Dimorphite-DL の max_variants（1 = 主要マイクロ種のみ、必須）

    Returns:
        SDF ブロック文字列（\\n 終端）、失敗時は None
    """
    try:
        from dimorphite_dl import protonate_smiles
        from rdkit import Chem
        from rdkit.Chem import AllChem, SDWriter
        from rdkit.Chem.MolStandardize import rdMolStandardize

        # Step 1: Dimorphite-DL protonation
        protonated = protonate_smiles(
            smiles_input=smiles,
            ph_min=ph_min,
            ph_max=ph_max,
            precision=1.0,
            max_variants=max_variants,
            label_identifiers=False,
            label_states=False,
            validate_output=True,
        )
        if not protonated:
            return None
        protonated_smiles = protonated[0]

        # Step 2: Parse and tautomer canonicalization
        mol = Chem.MolFromSmiles(protonated_smiles)
        if mol is None:
            return None
        te = rdMolStandardize.TautomerEnumerator()
        mol = te.Canonicalize(mol)
        if mol is None:
            return None

        # Step 3: Add hydrogens before 3D generation
        mol = Chem.AddHs(mol)

        # Step 4: ETKDGv3 3D conformer generation
        params = AllChem.ETKDGv3()
        params.randomSeed = seed
        result = AllChem.EmbedMolecule(mol, params)
        if result == -1:
            return None

        # Step 5: MMFF94 optimization
        AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)

        # Step 6: Add explicit H coordinates (ensure coords=True)
        mol = Chem.AddHs(mol, addCoords=True)

        # Step 7: SDF 書き出し
        buf = StringIO()
        writer = SDWriter(buf)
        writer.write(mol)
        writer.close()
        return buf.getvalue()

    except Exception:
        return None


def preprocess_and_hash(smiles: str, seed: int = 42) -> str:
    """preprocess_smiles で SDF を生成し SHA-256 hash を返す。

    テスト用ユーティリティ関数。
    """
    sdf = preprocess_smiles(smiles, seed=seed)
    if sdf is None:
        raise ValueError(f"前処理失敗: {smiles}")
    return hashlib.sha256(sdf.encode()).hexdigest()
