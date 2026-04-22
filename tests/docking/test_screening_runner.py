from __future__ import annotations

import gzip
import hashlib
import json
import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from docking_automation.docking.screening_runner import ScreeningResult, ScreeningRunner


# ── helpers ──────────────────────────────────────────────────────────────────

def _make_protein(protein_id: str, content_hash: str, tmp_path: Path | None = None) -> MagicMock:
    p = MagicMock()
    p.id = protein_id
    p.content_hash = content_hash
    p.path = (tmp_path or Path("/tmp")) / f"{protein_id}.pdb"
    return p


def _make_protein_set(n: int, tmp_path: Path | None = None) -> MagicMock:
    proteins = [_make_protein(f"protein_{i}", f"hash_p{i}", tmp_path) for i in range(n)]
    ps = MagicMock()
    ps.__iter__ = MagicMock(side_effect=lambda: iter(proteins))
    ps.__getitem__ = MagicMock(side_effect=lambda pid: next(p for p in proteins if p.id == pid))
    ps.content_hashes = MagicMock(
        return_value={p.id: p.content_hash for p in proteins}
    )
    return ps


def _make_compound_set(n: int, tmp_path: Path | None = None) -> MagicMock:
    cs = MagicMock()
    cs.get_compound_count = MagicMock(return_value=n)
    cs.get_compound_hash = MagicMock(
        side_effect=lambda i: hashlib.sha256(f"compound_{i}".encode()).hexdigest()
    )
    cs.path = (tmp_path or Path("/tmp")) / "compounds.sdf"
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


def _make_grid_box(center=(0.0, 0.0, 0.0), size=(20.0, 20.0, 20.0)) -> MagicMock:
    gb = MagicMock()
    gb.center = center
    gb.size = size
    return gb


def _make_runner(
    n_proteins: int = 2,
    n_compounds: int = 3,
    tmp_path: Path | None = None,
    grid_box_cache=None,
    grid_box_missing_policy: str = "skip",
    dask_n_workers: int = 4,
    _dock_fn=None,
    _cluster_kwargs: dict | None = None,
) -> ScreeningRunner:
    hdf5_path = (tmp_path or Path("/tmp")) / "test.h5"
    log_path = (tmp_path or Path("/tmp")) / "log.jsonl"
    return ScreeningRunner(
        protein_set=_make_protein_set(n_proteins, tmp_path),
        compound_set=_make_compound_set(n_compounds, tmp_path),
        grid_box_cache=grid_box_cache or _make_grid_box_cache(),
        hdf5_path=hdf5_path,
        log_path=log_path,
        grid_box_missing_policy=grid_box_missing_policy,
        dask_n_workers=dask_n_workers,
        _dock_fn=_dock_fn,
        _cluster_kwargs=_cluster_kwargs,
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


# ── Dask統合テスト用フェイクdock関数 (モジュールレベル必須: Daskがcloudpickleでシリアライズするため) ──

def _fake_dock_one_protein(
    protein_path,
    protein_id,
    protein_content_hash,
    compound_sdf_path,
    compound_indices,
    compound_hashes,
    grid_center,
    grid_size,
    exhaustiveness=1,
    top_n_poses=1,
):
    """テスト用: ファイルアクセスなしに即座にフェイク結果を返す。"""
    import gzip

    return [
        {
            "protein_id": protein_id,
            "compound_index": idx,
            "protein_content_hash": protein_content_hash,
            "compound_content_hash": compound_hashes.get(idx, f"chash_{idx}"),
            "score": round(-7.0 - idx * 0.1, 3),
            "pose_blob": gzip.compress(b"fake_sdf"),
            "elapsed_sec": 0.001,
            "error": None,
        }
        for idx in compound_indices
    ]


@pytest.mark.slow
def test_run_with_dask_small(tmp_path):
    """2タンパク質×2化合物のDask実行が完了する統合テスト。"""
    n_proteins = 2
    n_compounds = 2

    grid_box = _make_grid_box()
    cache = MagicMock()
    cache.get = MagicMock(return_value=grid_box)

    runner = _make_runner(
        n_proteins=n_proteins,
        n_compounds=n_compounds,
        tmp_path=tmp_path,
        grid_box_cache=cache,
        dask_n_workers=2,
        _dock_fn=_fake_dock_one_protein,
        _cluster_kwargs={"processes": False},
    )
    mock_repo = _make_repo(existing_pairs=set())

    result = runner.run(resume=True, _repo=mock_repo)

    total = n_proteins * n_compounds
    assert result.total_pairs == total
    assert result.new_pairs == total
    assert result.reused_pairs == 0
    assert result.failed_pairs == 0
    assert result.elapsed_sec > 0


@pytest.mark.slow
def test_jsonl_output_format(tmp_path):
    """JSONL ファイルが正しいフォーマットで出力される。"""
    n_proteins = 2
    n_compounds = 2

    grid_box = _make_grid_box()
    cache = MagicMock()
    cache.get = MagicMock(return_value=grid_box)

    runner = _make_runner(
        n_proteins=n_proteins,
        n_compounds=n_compounds,
        tmp_path=tmp_path,
        grid_box_cache=cache,
        dask_n_workers=2,
        _dock_fn=_fake_dock_one_protein,
        _cluster_kwargs={"processes": False},
    )
    mock_repo = _make_repo(existing_pairs=set())

    runner.run(resume=True, _repo=mock_repo)

    log_text = runner.log_path.read_text().strip()
    assert log_text, "JSONL ファイルが空"

    lines = log_text.splitlines()
    assert len(lines) == n_proteins * n_compounds, (
        f"期待行数: {n_proteins * n_compounds}, 実際: {len(lines)}"
    )

    for line in lines:
        entry = json.loads(line)
        assert "protein_id" in entry
        assert "compound_index" in entry
        assert "status" in entry
        assert "score" in entry
        assert "elapsed_sec" in entry
        assert entry["status"] in ("new", "reused", "failed")
        assert isinstance(entry["compound_index"], int)
        assert isinstance(entry["elapsed_sec"], (int, float))
