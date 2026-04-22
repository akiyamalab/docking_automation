from pathlib import Path

import pytest

from docking_automation.molecule.protein import Protein
from docking_automation.molecule.protein_set import ProteinSet

PDB_ALA_N = "ATOM      1  N   ALA A   1      10.000  10.000  10.000  1.00  0.00           N"
PDB_ALA_CA = (
    "ATOM      1  N   ALA A   1      10.000  10.000  10.000  1.00  0.00           N\n"
    "ATOM      2  CA  ALA A   1      11.000  10.000  10.000  1.00  0.00           C"
)
PDB_GLY_N = "ATOM      1  N   GLY A   1      20.000  20.000  20.000  1.00  0.00           N"


def _make_pdb(tmp_path: Path, name: str, content: str = PDB_ALA_N) -> Path:
    p = tmp_path / name
    p.write_text(content)
    return p


def _make_protein(tmp_path: Path, name: str, content: str = PDB_ALA_N) -> Protein:
    path = _make_pdb(tmp_path, name, content)
    return Protein(path)


@pytest.fixture
def protein_set(tmp_path):
    p1 = _make_protein(tmp_path, "prot_a.pdb", PDB_ALA_N)
    p2 = _make_protein(tmp_path, "prot_b.pdb", PDB_GLY_N)
    p3 = _make_protein(tmp_path, "prot_c.pdb", PDB_ALA_CA)
    return ProteinSet([p1, p2, p3])


class TestProteinSetCreation:
    def test_create_from_proteins(self, tmp_path):
        p1 = _make_protein(tmp_path, "a.pdb")
        p2 = _make_protein(tmp_path, "b.pdb")
        ps = ProteinSet([p1, p2])
        assert len(ps) == 2
        assert list(ps.protein_ids) == ["a", "b"]
        proteins = list(ps)
        assert proteins[0].id == "a"
        assert proteins[1].id == "b"

    def test_duplicate_protein_id_raises(self, tmp_path):
        p1 = _make_protein(tmp_path, "dup.pdb")
        p2 = Protein(p1.path, id=p1.id)
        with pytest.raises(ValueError, match="重複"):
            ProteinSet([p1, p2])

    def test_auto_generate_id_when_omitted(self, tmp_path):
        p = _make_protein(tmp_path, "x.pdb")
        ps = ProteinSet([p])
        assert ps.id is not None
        assert len(ps.id) > 0

    def test_explicit_id_preserved(self, tmp_path):
        p = _make_protein(tmp_path, "x.pdb")
        ps = ProteinSet([p], id="my_set")
        assert ps.id == "my_set"


class TestProteinSetCollection:
    def test_getitem(self, protein_set):
        p = protein_set["prot_a"]
        assert p.id == "prot_a"

    def test_getitem_missing_raises(self, protein_set):
        with pytest.raises(KeyError):
            _ = protein_set["nonexistent"]

    def test_contains(self, protein_set):
        assert "prot_a" in protein_set
        assert "nonexistent" not in protein_set

    def test_subset(self, protein_set):
        sub = protein_set.subset(["prot_a", "prot_c"])
        assert len(sub) == 2
        assert "prot_b" not in sub
        assert "prot_a" in sub
        assert "prot_c" in sub

    def test_subset_missing_raises(self, protein_set):
        with pytest.raises(KeyError):
            protein_set.subset(["prot_a", "nonexistent"])

    def test_filter(self, protein_set):
        filtered = protein_set.filter(lambda p: p.id != "prot_b")
        assert len(filtered) == 2
        assert "prot_b" not in filtered


class TestProteinSetHashes:
    def test_content_hashes_returns_all(self, protein_set):
        hashes = protein_set.content_hashes()
        assert set(hashes.keys()) == set(protein_set.protein_ids)

    def test_content_hashes_cached(self, protein_set):
        h1 = protein_set.content_hashes()
        h2 = protein_set.content_hashes()
        assert h1 == h2

    def test_same_content_same_hash(self, tmp_path):
        p1 = _make_protein(tmp_path, "x1.pdb", PDB_ALA_N)
        p2 = _make_protein(tmp_path, "x2.pdb", PDB_ALA_N)
        ps = ProteinSet([p1, p2])
        hashes = ps.content_hashes()
        assert hashes["x1"] == hashes["x2"]


class TestProteinSetFactory:
    def test_from_directory(self, tmp_path):
        _make_pdb(tmp_path, "p1.pdb")
        _make_pdb(tmp_path, "p2.pdb", PDB_GLY_N)
        _make_pdb(tmp_path, "ignore.txt")

        ps = ProteinSet.from_directory(tmp_path)
        assert len(ps) == 2
        assert "p1" in ps
        assert "p2" in ps

    def test_from_directory_recursive(self, tmp_path):
        sub = tmp_path / "sub"
        sub.mkdir()
        _make_pdb(tmp_path, "top.pdb")
        _make_pdb(sub, "nested.pdb")

        ps = ProteinSet.from_directory(tmp_path, recursive=True)
        assert len(ps) == 2

    def test_from_directory_custom_id(self, tmp_path):
        _make_pdb(tmp_path, "a.pdb")
        ps = ProteinSet.from_directory(tmp_path, id="custom")
        assert ps.id == "custom"

    def test_from_afdb_mouse_filters_large(self, tmp_path):
        pdb_dir = tmp_path / "pdb"
        pdb_dir.mkdir()

        small_content = PDB_ALA_CA  # 1 CA atom
        large_content = "\n".join(
            f"ATOM  {i+1:4d}  CA  ALA A{i+1:4d}      0.000   0.000   0.000  1.00  0.00           C"
            for i in range(5)
        )

        (pdb_dir / "small.pdb").write_text(small_content)
        (pdb_dir / "large.pdb").write_text(large_content)

        ps = ProteinSet.from_afdb_mouse(root=pdb_dir, max_residues=3, id="test_afdb")
        assert len(ps) == 1
        assert "small" in ps
        assert "large" not in ps

    def test_from_afdb_mouse_limit(self, tmp_path):
        pdb_dir = tmp_path / "pdb"
        pdb_dir.mkdir()
        for i in range(5):
            (pdb_dir / f"p{i}.pdb").write_text(PDB_ALA_N)

        ps = ProteinSet.from_afdb_mouse(root=pdb_dir, limit=3, max_residues=None)
        assert len(ps) == 3

    def test_from_afdb_mouse_empty_dir(self, tmp_path):
        pdb_dir = tmp_path / "pdb"
        pdb_dir.mkdir()
        ps = ProteinSet.from_afdb_mouse(root=pdb_dir, id="empty")
        assert len(ps) == 0
