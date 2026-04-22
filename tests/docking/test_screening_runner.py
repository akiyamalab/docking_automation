from __future__ import annotations

import hashlib
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from docking_automation.docking.screening_runner import ScreeningResult, ScreeningRunner


# ── helpers ──────────────────────────────────────────────────────────────────

def _make_protein(protein_id: str, content_hash: str) -> MagicMock:
    p = MagicMock()
    p.id = protein_id
    p.content_hash = content_hash
    return p


def _make_protein_set(n: int) -> MagicMock:
    proteins = [_make_protein(f"protein_{i}", f"hash_p{i}") for i in range(n)]
    ps = MagicMock()
    ps.__iter__ = MagicMock(side_effect=lambda: iter(proteins))
    ps.__getitem__ = MagicMock(side_effect=lambda pid: next(p for p in proteins if p.id == pid))
    ps.content_hashes = MagicMock(
        return_value={p.id: p.content_hash for p in proteins}
    )
    return ps


def _make_compound_set(n: int) -> MagicMock:
    cs = MagicMock()
    cs.get_compound_count = MagicMock(return_value=n)
    cs.get_compound_hash = MagicMock(
        side_effect=lambda i: hashlib.sha256(f"compound_{i}".encode()).hexdigest()
    )
    return cs


def _make_grid_box_cache(proteins_with_box: list[str] | None = None) -> MagicMock:
    """proteins_with_box=None の場合は全タンパク質に GridBox あり。"""
    cache = MagicMock()
    if proteins_with_box is None:
        cache.get = MagicMock(return_value=MagicMock())
    else:
        def _get(protein):
            return MagicMock() if protein.id in proteins_with_box else None
        cache.get = MagicMock(side_effect=_get)
    return cache


def _make_repo(existing_pairs: set[tuple[str, str]] | None = None) -> MagicMock:
    repo = MagicMock()
    if existing_pairs is None:
        repo._exists = MagicMock(return_value=False)
    else:
        repo._exists = MagicMock(
            side_effect=lambda p_hash, c_hash: (p_hash, c_hash) in existing_pairs
        )
    return repo


def _make_runner(
    n_proteins: int = 2,
    n_compounds: int = 3,
    tmp_path: Path | None = None,
    grid_box_cache=None,
    grid_box_missing_policy: str = "skip",
) -> ScreeningRunner:
    hdf5_path = (tmp_path or Path("/tmp")) / "test.h5"
    log_path = (tmp_path or Path("/tmp")) / "log.jsonl"
    return ScreeningRunner(
        protein_set=_make_protein_set(n_proteins),
        compound_set=_make_compound_set(n_compounds),
        grid_box_cache=grid_box_cache or _make_grid_box_cache(),
        hdf5_path=hdf5_path,
        log_path=log_path,
        grid_box_missing_policy=grid_box_missing_policy,
    )


# ── tests ─────────────────────────────────────────────────────────────────────

def test_screening_result_is_frozen():
    result = ScreeningResult(
        total_pairs=12,
        new_pairs=10,
        reused_pairs=2,
        failed_pairs=0,
        elapsed_sec=1.23,
        hdf5_path=Path("/tmp/test.h5"),
        log_path=Path("/tmp/log.jsonl"),
    )
    with pytest.raises(Exception):
        result.total_pairs = 999  # type: ignore[misc]


def test_enumerate_pairs_count(tmp_path):
    runner = _make_runner(n_proteins=3, n_compounds=4, tmp_path=tmp_path)
    pairs = list(runner._enumerate_pairs())
    assert len(pairs) == 12
    protein_ids = {pid for pid, _ in pairs}
    assert len(protein_ids) == 3
    compound_indices = {ci for _, ci in pairs}
    assert compound_indices == {0, 1, 2, 3}


def test_filter_unprocessed_skips_existing(tmp_path):
    n_proteins = 2
    n_compounds = 3
    runner = _make_runner(n_proteins=n_proteins, n_compounds=n_compounds, tmp_path=tmp_path)

    protein_hashes = runner.protein_set.content_hashes()
    c_hash_0 = runner.compound_set.get_compound_hash(0)
    existing = {(protein_hashes["protein_0"], c_hash_0)}

    repo = _make_repo(existing_pairs=existing)
    unprocessed = runner._filter_unprocessed(repo)

    total = n_proteins * n_compounds
    assert len(unprocessed) == total - 1
    for protein_id, compound_index in unprocessed:
        p_hash = protein_hashes[protein_id]
        c_hash = runner.compound_set.get_compound_hash(compound_index)
        assert (p_hash, c_hash) not in existing


def test_run_resume_true_skips_done(tmp_path):
    n_proteins = 2
    n_compounds = 3
    runner = _make_runner(n_proteins=n_proteins, n_compounds=n_compounds, tmp_path=tmp_path)

    protein_hashes = runner.protein_set.content_hashes()
    all_existing = {
        (protein_hashes[f"protein_{i}"], runner.compound_set.get_compound_hash(j))
        for i in range(n_proteins)
        for j in range(n_compounds)
    }
    mock_repo = _make_repo(existing_pairs=all_existing)

    result = runner.run(resume=True, _repo=mock_repo)

    assert result.new_pairs == 0
    assert result.reused_pairs == n_proteins * n_compounds
    assert result.total_pairs == n_proteins * n_compounds
    assert result.failed_pairs == 0


def test_grid_box_missing_skip(tmp_path):
    n_proteins = 1
    n_compounds = 2
    cache = _make_grid_box_cache(proteins_with_box=[])  # GridBox なし
    runner = _make_runner(
        n_proteins=n_proteins,
        n_compounds=n_compounds,
        tmp_path=tmp_path,
        grid_box_cache=cache,
        grid_box_missing_policy="skip",
    )

    mock_repo = _make_repo(existing_pairs=set())

    result = runner.run(resume=True, _repo=mock_repo)

    assert result.new_pairs == 0
    assert result.failed_pairs == 0
    assert result.total_pairs == n_proteins * n_compounds

    log_lines = runner.log_path.read_text().strip().splitlines()
    assert len(log_lines) == n_proteins * n_compounds
    for line in log_lines:
        import json
        entry = json.loads(line)
        assert entry["error"] == "grid_box_missing"
