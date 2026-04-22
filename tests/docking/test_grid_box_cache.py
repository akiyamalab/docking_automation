import pytest

from docking_automation.docking.grid_box import GridBox
from docking_automation.docking.grid_box_cache import GridBoxCache, GridBoxCacheEntry
from docking_automation.molecule.protein import Protein


PDB_CONTENT = "ATOM      1  N   ALA A   1      10.000  10.000  10.000  1.00  0.00           N"


@pytest.fixture
def make_protein(tmp_path):
    def _make(name: str, content: str = PDB_CONTENT) -> Protein:
        p = tmp_path / f"{name}.pdb"
        p.write_text(content)
        return Protein(path=p, id=name)
    return _make


@pytest.fixture
def sample_grid_box():
    return GridBox(center=(40.0, 27.0, 41.0), size=(22.0, 22.0, 22.0))


@pytest.fixture
def mock_protein_set(make_protein):
    class MockProteinSet:
        def __init__(self, proteins):
            self._proteins = proteins
        def __iter__(self):
            return iter(self._proteins)

    p1 = make_protein("prot_A")
    p2 = make_protein("prot_B")
    return MockProteinSet([p1, p2])


class TestGridBoxCacheBasicOps:
    def test_put_and_get_hit(self, tmp_path, make_protein, sample_grid_box):
        cache = GridBoxCache(path=tmp_path / "cache.json")
        protein = make_protein("protein1")
        cache.put(protein, sample_grid_box)
        result = cache.get(protein)
        assert result is not None
        assert list(result.center) == pytest.approx(list(sample_grid_box.center))
        assert list(result.size) == pytest.approx(list(sample_grid_box.size))

    def test_get_miss_returns_none(self, tmp_path, make_protein):
        cache = GridBoxCache(path=tmp_path / "cache.json")
        protein = make_protein("unknown")
        assert cache.get(protein) is None

    def test_has_with_hash_check(self, tmp_path, make_protein, sample_grid_box):
        cache = GridBoxCache(path=tmp_path / "cache.json")
        protein = make_protein("protein1")
        cache.put(protein, sample_grid_box)
        assert cache.has(protein) is True

        # 同じIDで中身が違うファイルを作成 → hash不一致でFalse
        other = make_protein("protein1", content=PDB_CONTENT + "\nEND")
        assert cache.has(other) is False

    def test_invalidate(self, tmp_path, make_protein, sample_grid_box):
        cache = GridBoxCache(path=tmp_path / "cache.json")
        protein = make_protein("protein1")
        cache.put(protein, sample_grid_box)
        assert cache.has(protein) is True
        cache.invalidate("protein1")
        assert cache.has(protein) is False


class TestGridBoxCachePersistence:
    def test_save_and_reload(self, tmp_path, make_protein, sample_grid_box):
        cache_path = tmp_path / "cache.json"
        cache = GridBoxCache(path=cache_path)
        protein = make_protein("protein1")
        cache.put(protein, sample_grid_box)
        cache.save()

        reloaded = GridBoxCache.from_file(cache_path)
        result = reloaded.get(protein)
        assert result is not None
        assert list(result.center) == pytest.approx(list(sample_grid_box.center))
        assert list(result.size) == pytest.approx(list(sample_grid_box.size))

    def test_from_file_empty_when_not_exists(self, tmp_path):
        cache = GridBoxCache.from_file(tmp_path / "nonexistent.json")
        assert len(cache) == 0

    def test_atomic_save(self, tmp_path, make_protein, sample_grid_box):
        cache_path = tmp_path / "cache.json"
        cache = GridBoxCache(path=cache_path)
        protein = make_protein("protein1")
        cache.put(protein, sample_grid_box)
        cache.save(atomic=True)

        # tmpファイルが残っていないことを確認
        tmp_files = list(tmp_path.glob("*.tmp"))
        assert len(tmp_files) == 0
        assert cache_path.exists()


class TestGridBoxCacheMissingIds:
    def test_missing_ids(self, tmp_path, make_protein, sample_grid_box, mock_protein_set):
        cache = GridBoxCache(path=tmp_path / "cache.json")
        proteins = list(mock_protein_set)

        # prot_A のみキャッシュ済み
        cache.put(proteins[0], sample_grid_box)

        missing = cache.missing_ids(mock_protein_set)
        assert "prot_B" in missing
        assert "prot_A" not in missing

    def test_missing_ids_hash_mismatch(self, tmp_path, make_protein, sample_grid_box):
        cache = GridBoxCache(path=tmp_path / "cache.json")
        original = make_protein("prot_X")
        cache.put(original, sample_grid_box)

        # 同一IDで内容が違うProteinを使ったProteinSet
        modified = make_protein("prot_X", content=PDB_CONTENT + "\nEND")

        class OneItemSet:
            def __iter__(self):
                return iter([modified])

        missing = cache.missing_ids(OneItemSet())
        assert "prot_X" in missing


class TestGridBoxCacheLen:
    def test_len(self, tmp_path, make_protein, sample_grid_box):
        cache = GridBoxCache(path=tmp_path / "cache.json")
        assert len(cache) == 0
        protein = make_protein("protein1")
        cache.put(protein, sample_grid_box)
        assert len(cache) == 1

    def test_contains(self, tmp_path, make_protein, sample_grid_box):
        cache = GridBoxCache(path=tmp_path / "cache.json")
        protein = make_protein("protein1")
        assert "protein1" not in cache
        cache.put(protein, sample_grid_box)
        assert "protein1" in cache


class TestGridBoxCacheEntries:
    def test_entries_iterator(self, tmp_path, make_protein, sample_grid_box):
        cache = GridBoxCache(path=tmp_path / "cache.json")
        p1 = make_protein("prot_A")
        p2 = make_protein("prot_B")
        cache.put(p1, sample_grid_box)
        cache.put(p2, sample_grid_box)

        ids = {e.protein_id for e in cache.entries()}
        assert ids == {"prot_A", "prot_B"}


class TestGridBoxCacheMissingPolicy:
    def test_get_with_policy_skip_when_missing(self, tmp_path):
        cache = GridBoxCache(path=tmp_path / "cache.json")
        result = cache.get_with_policy("unknown_protein", missing_policy="skip")
        assert result is None

    def test_get_with_policy_error_when_missing(self, tmp_path):
        cache = GridBoxCache(path=tmp_path / "cache.json")
        with pytest.raises(KeyError, match="unknown_protein"):
            cache.get_with_policy("unknown_protein", missing_policy="error")

    def test_get_with_policy_fallback_centroid(self, tmp_path):
        cache = GridBoxCache(path=tmp_path / "cache.json")
        result = cache.get_with_policy(
            "unknown_protein",
            missing_policy="fallback_centroid",
            fallback_center=[0, 0, 0],
            fallback_size=[20, 20, 20],
        )
        assert result is not None
        assert list(result.center) == pytest.approx([0, 0, 0])
        assert list(result.size) == pytest.approx([20, 20, 20])

    def test_get_with_policy_returns_cached_regardless_policy(
        self, tmp_path, make_protein, sample_grid_box
    ):
        cache = GridBoxCache(path=tmp_path / "cache.json")
        protein = make_protein("prot_A")
        cache.put(protein, sample_grid_box)

        for policy in ("skip", "error", "fallback_centroid"):
            result = cache.get_with_policy("prot_A", missing_policy=policy)
            assert result is not None
            assert list(result.center) == pytest.approx(list(sample_grid_box.center))

    def test_get_with_policy_fallback_missing_args_raises(self, tmp_path):
        cache = GridBoxCache(path=tmp_path / "cache.json")
        with pytest.raises(ValueError):
            cache.get_with_policy(
                "unknown_protein",
                missing_policy="fallback_centroid",
                fallback_center=None,
                fallback_size=None,
            )
