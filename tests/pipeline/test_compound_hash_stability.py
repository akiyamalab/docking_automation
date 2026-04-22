import hashlib
import json
from pathlib import Path

from docking_automation.compound_pipeline.preprocess import preprocess_smiles, preprocess_and_hash

ASPIRIN = "CC(=O)OC1=CC=CC=C1C(=O)O"
ETHANOL = "CCO"
IBUPROFEN = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"


class TestPreprocessDeterminism:
    def test_same_smiles_same_sdf(self):
        """同一 SMILES → 同一 SDF ブロック（セッション内決定論）"""
        sdf1 = preprocess_smiles(ASPIRIN)
        sdf2 = preprocess_smiles(ASPIRIN)
        assert sdf1 == sdf2, "Same-input determinism broken"

    def test_same_smiles_same_hash(self):
        """同一 SMILES → 同一 hash（セッション内）"""
        h1 = preprocess_and_hash(ETHANOL)
        h2 = preprocess_and_hash(ETHANOL)
        assert h1 == h2, "Hash determinism broken"

    def test_different_smiles_different_hash(self):
        """異なる SMILES → 異なる hash"""
        h1 = preprocess_and_hash(ASPIRIN)
        h2 = preprocess_and_hash(IBUPROFEN)
        assert h1 != h2

    def test_seed_stability(self):
        """seed=42 で繰り返し呼んでも同じ結果"""
        h1 = preprocess_and_hash(ASPIRIN, seed=42)
        h2 = preprocess_and_hash(ASPIRIN, seed=42)
        assert h1 == h2


class TestGoldenHash:
    """ゴールデンハッシュテスト: 初回実行で golden_hashes.json を生成し、
    以降はその値と比較する。"""

    GOLDEN_PATH = "tests/pipeline/golden_hashes.json"
    SMILES_LIST = [ASPIRIN, ETHANOL, IBUPROFEN]

    def test_golden_hash_stability(self):
        """golden_hashes.json の値と一致すること。
        ファイルが無ければ現在の結果で初期化し PASS（初回のみ）。
        """
        golden_path = Path(self.GOLDEN_PATH)
        hashes = {s: preprocess_and_hash(s) for s in self.SMILES_LIST}
        if not golden_path.exists():
            golden_path.parent.mkdir(parents=True, exist_ok=True)
            with open(golden_path, "w") as f:
                json.dump(hashes, f, indent=2)
            return  # 初回はスキップ（golden を生成）
        with open(golden_path) as f:
            golden = json.load(f)
        for smiles, expected in golden.items():
            actual = preprocess_and_hash(smiles)
            assert actual == expected, f"Hash drift: {smiles}"

    def test_preprocess_returns_valid_sdf(self):
        """前処理結果が有効な SDF フォーマットであること"""
        sdf = preprocess_smiles(ASPIRIN)
        assert sdf is not None
        assert "$$$$" in sdf  # SDF 末尾区切り
