from __future__ import annotations

import gzip
import hashlib
import json
import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from docking_automation.docking.autodockvina_docking import AutoDockVina
from docking_automation.docking.screening_runner import (
    ScreeningResult,
    ScreeningRunner,
    _make_docking_tool,
)
from docking_automation.docking.unidock_docking import UniDockDocking


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
    backend="vina",
    search_mode="balance",
    rescue_mode=False,
    rescue_search_mode="detail",
    **kwargs,
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


# ── backend 切替テスト ───────────────────────────────────────────────────────

def test_backend_default_is_vina(tmp_path):
    """ScreeningRunner() のデフォルト backend が "vina" であること。"""
    runner = _make_runner(tmp_path=tmp_path)
    assert runner.backend == "vina"


def test_backend_unidock_accepted(tmp_path):
    """ScreeningRunner(backend="unidock") が例外なく生成できること。"""
    runner = _make_runner(tmp_path=tmp_path)
    runner2 = ScreeningRunner(
        protein_set=runner.protein_set,
        compound_set=runner.compound_set,
        grid_box_cache=runner.grid_box_cache,
        hdf5_path=tmp_path / "test2.h5",
        log_path=tmp_path / "log2.jsonl",
        backend="unidock",
    )
    assert runner2.backend == "unidock"


def test_make_docking_tool_vina():
    """_make_docking_tool("vina") が AutoDockVina インスタンスを返すこと。"""
    tool = _make_docking_tool("vina")
    assert isinstance(tool, AutoDockVina)


def test_make_docking_tool_unidock():
    """_make_docking_tool("unidock") が UniDockDocking インスタンスを返すこと。"""
    tool = _make_docking_tool("unidock")
    assert isinstance(tool, UniDockDocking)


def test_make_docking_tool_unknown():
    """_make_docking_tool("unknown_backend") が ValueError を raise すること。"""
    with pytest.raises(ValueError, match="Unknown backend"):
        _make_docking_tool("unknown_backend")


def test_backend_propagated_to_dock_one_protein(tmp_path):
    """run() 実行時に client.submit が backend="unidock" 引数付きで呼ばれること。"""
    grid_box = _make_grid_box()
    cache = MagicMock()
    cache.get = MagicMock(return_value=grid_box)

    protein_set = _make_protein_set(1, tmp_path)
    compound_set = _make_compound_set(1, tmp_path)

    runner = ScreeningRunner(
        protein_set=protein_set,
        compound_set=compound_set,
        grid_box_cache=cache,
        hdf5_path=tmp_path / "test.h5",
        log_path=tmp_path / "log.jsonl",
        backend="unidock",
        _dock_fn=_fake_dock_one_protein,
    )

    mock_future = MagicMock()
    mock_client = MagicMock()
    mock_client.submit.return_value = mock_future

    runner._submit_to_dask(mock_client, {"protein_0": [0]})

    assert mock_client.submit.called, "client.submit が呼ばれていない"
    call_kwargs = mock_client.submit.call_args[1]
    assert call_kwargs.get("backend") == "unidock", (
        f"backend kwarg が 'unidock' でない: {call_kwargs}"
    )


# ── extra_padding テスト ────────────────────────────────────────────────────────

def test_extra_padding_increases_box_size(tmp_path):
    """extra_padding=5.0 でgrid_sizeが各次元 +10.0 (= 5.0*2) 拡大されること。"""
    base_size = (20.0, 22.0, 24.0)
    grid_box = _make_grid_box(center=(0.0, 0.0, 0.0), size=base_size)
    cache = MagicMock()
    cache.get = MagicMock(return_value=grid_box)

    protein_set = _make_protein_set(1, tmp_path)
    compound_set = _make_compound_set(1, tmp_path)

    runner = ScreeningRunner(
        protein_set=protein_set,
        compound_set=compound_set,
        grid_box_cache=cache,
        hdf5_path=tmp_path / "test.h5",
        log_path=tmp_path / "log.jsonl",
        extra_padding=5.0,
        _dock_fn=_fake_dock_one_protein,
    )

    mock_future = MagicMock()
    mock_client = MagicMock()
    mock_client.submit.return_value = mock_future

    runner._submit_to_dask(mock_client, {"protein_0": [0]})

    assert mock_client.submit.called
    call_args = mock_client.submit.call_args[0]
    # positional args: dock_fn, protein_path, protein_id, p_hash, sdf_path,
    #                  compound_indices, compound_hashes, grid_center, grid_size, ...
    passed_size = call_args[8]
    expected = [s + 10.0 for s in base_size]
    assert passed_size == pytest.approx(expected), (
        f"期待サイズ: {expected}, 実際: {passed_size}"
    )




def test_extra_padding_zero_unchanged(tmp_path):
    """extra_padding=0.0 でgrid_sizeが変化しないこと（後方互換）。"""
    base_size = (20.0, 22.0, 24.0)
    grid_box = _make_grid_box(center=(0.0, 0.0, 0.0), size=base_size)
    cache = MagicMock()
    cache.get = MagicMock(return_value=grid_box)

    protein_set = _make_protein_set(1, tmp_path)
    compound_set = _make_compound_set(1, tmp_path)

    runner = ScreeningRunner(
        protein_set=protein_set,
        compound_set=compound_set,
        grid_box_cache=cache,
        hdf5_path=tmp_path / "test.h5",
        log_path=tmp_path / "log.jsonl",
        extra_padding=0.0,
        _dock_fn=_fake_dock_one_protein,
    )

    mock_future = MagicMock()
    mock_client = MagicMock()
    mock_client.submit.return_value = mock_future

    runner._submit_to_dask(mock_client, {"protein_0": [0]})

    assert mock_client.submit.called
    call_args = mock_client.submit.call_args[0]
    # positional args: dock_fn, protein_path, protein_id, p_hash, sdf_path,
    #                  compound_indices, compound_hashes, grid_center, grid_size, ...
    passed_size = call_args[8]
    expected = list(base_size)
    assert passed_size == pytest.approx(expected), (
        f"期待サイズ: {expected}, 実際: {passed_size}"
    )


# --- ScreeningTool 統合テスト (Step 1 刷新) ---

def test_make_docking_tool_unidock2():
    """unidock2 backend で UniDock2Docking が返ること。"""
    from docking_automation.docking.unidock2_docking import UniDock2Docking
    from docking_automation.docking.screening_runner import _make_docking_tool
    tool = _make_docking_tool('unidock2')
    assert isinstance(tool, UniDock2Docking)


def test_make_docking_tool_returns_screening_tool():
    """全 backend が ScreeningTool 派生を返すこと。"""
    from docking_automation.docking.screening_tool import ScreeningTool
    from docking_automation.docking.screening_runner import _make_docking_tool
    for backend in ('vina', 'unidock', 'unidock2'):
        assert isinstance(_make_docking_tool(backend), ScreeningTool), (
            f'{backend} backend should return ScreeningTool instance'
        )


# --- cache_dir 統合テスト ---

def test_screening_runner_init_cache_params(tmp_path):
    """cache_dir / force_cache / cache_n_workers が正しく格納されること。"""
    ps = MagicMock()
    cs = MagicMock()
    gbc = MagicMock()
    cache_dir = tmp_path / "cache"
    runner = ScreeningRunner(
        protein_set=ps,
        compound_set=cs,
        grid_box_cache=gbc,
        hdf5_path=tmp_path / "out.h5",
        cache_dir=cache_dir,
        force_cache=True,
        cache_n_workers=12,
    )
    assert runner.cache_dir == cache_dir
    assert runner.force_cache is True
    assert runner.cache_n_workers == 12


def test_screening_runner_init_cache_defaults(tmp_path):
    """cache_dir がデフォルトで None であること。"""
    runner = ScreeningRunner(
        protein_set=MagicMock(),
        compound_set=MagicMock(),
        grid_box_cache=MagicMock(),
        hdf5_path=tmp_path / "out.h5",
    )
    assert runner.cache_dir is None
    assert runner.force_cache is False
    assert runner.cache_n_workers is None


def test_prepare_receptor_caches_skips_non_unidock2(tmp_path):
    """backend=vina のとき prepare_receptor_caches は空 dict を返す。"""
    runner = ScreeningRunner(
        protein_set=MagicMock(),
        compound_set=MagicMock(),
        grid_box_cache=MagicMock(),
        hdf5_path=tmp_path / "out.h5",
        backend="vina",
        cache_dir=tmp_path / "cache",
    )
    result = runner.prepare_receptor_caches({"prot_A": [0, 1]})
    assert result == {}


def test_prepare_receptor_caches_skips_when_no_cache_dir(tmp_path):
    """cache_dir=None のとき空 dict を返す。"""
    runner = ScreeningRunner(
        protein_set=MagicMock(),
        compound_set=MagicMock(),
        grid_box_cache=MagicMock(),
        hdf5_path=tmp_path / "out.h5",
        backend="unidock2",
        cache_dir=None,
    )
    result = runner.prepare_receptor_caches({"prot_A": [0, 1]})
    assert result == {}


def test_prepare_receptor_caches_skips_existing(tmp_path):
    """既存 cache がある場合はスキップし、Dask クラスタを作らない。"""
    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()

    protein = _make_protein("pA", "hash_aaa", tmp_path)
    # protein.path ファイルも存在させる
    protein.path.write_text("ATOM dummy")

    # 既存 cache ファイルを作成
    (cache_dir / "hash_aaa.json").write_text("{}")

    ps = MagicMock()
    ps.__getitem__ = MagicMock(return_value=protein)
    ps.content_hashes.return_value = {"pA": "hash_aaa"}

    gbc = MagicMock()
    gbc.get.return_value = MagicMock(center=(0, 0, 0), size=(10, 10, 10))

    runner = ScreeningRunner(
        protein_set=ps,
        compound_set=MagicMock(),
        grid_box_cache=gbc,
        hdf5_path=tmp_path / "out.h5",
        backend="unidock2",
        cache_dir=cache_dir,
        force_cache=False,
    )
    # force=False なので全部スキップ → Dask cluster を作らない → 空 dict
    result = runner.prepare_receptor_caches({"pA": [0, 1]})
    assert result == {}


@patch("docking_automation.docking.screening_runner._make_docking_tool")
def test_dock_one_protein_via_tool_uses_cache_dir(mock_make_tool, tmp_path):
    """cache_dir 指定時、既存 cache があれば prepare_receptor_cache を呼ばない。"""
    from docking_automation.docking.screening_runner import dock_one_protein_via_tool

    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()

    # dummy cache を作成
    chash = "dummyhash123"
    (cache_dir / f"{chash}.json").write_text('{"receptor": []}')

    # dummy compound SDF
    sdf = tmp_path / "compounds.sdf"
    sdf.write_text(
        "benzene\n\n\n  1  0  0  0  0  0  0  0  0  0  1 V2000\n"
        "    0.0    0.0    0.0 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "M  END\n$$$$\n"
    )

    # dummy protein PDB
    pdb = tmp_path / "pA.pdb"
    pdb.write_text("ATOM      1  N   ALA A   1       0.0   0.0   0.0  1.00  0.00\nEND\n")

    mock_tool = MagicMock()
    mock_tool.dock_with_cache.return_value = []
    mock_make_tool.return_value = mock_tool

    dock_one_protein_via_tool(
        protein_path=str(pdb),
        protein_id="pA",
        protein_content_hash=chash,
        compound_sdf_path=str(sdf),
        compound_indices=[0],
        compound_hashes={0: "chash_0"},
        grid_center=[0.0, 0.0, 0.0],
        grid_size=[10.0, 10.0, 10.0],
        backend="unidock2",
        cache_dir=str(cache_dir),
    )

    # prepare_receptor_cache は呼ばれない (cache が既にある)
    mock_tool.prepare_receptor_cache.assert_not_called()
    # dock_with_cache は cache path で呼ばれる
    mock_tool.dock_with_cache.assert_called_once()
    call_kwargs = mock_tool.dock_with_cache.call_args
    assert str(call_kwargs.kwargs.get("cache", call_kwargs[1].get("cache", ""))) == str(
        cache_dir / f"{chash}.json"
    ) or str(call_kwargs[0][0] if call_kwargs[0] else "") == str(
        cache_dir / f"{chash}.json"
    )


@patch("docking_automation.docking.screening_runner._make_docking_tool")
def test_dock_one_protein_via_tool_fallback_no_cache(mock_make_tool, tmp_path):
    """cache_dir なしの場合、prepare_receptor_cache が呼ばれること。"""
    from docking_automation.docking.screening_runner import dock_one_protein_via_tool

    sdf = tmp_path / "compounds.sdf"
    sdf.write_text(
        "benzene\n\n\n  1  0  0  0  0  0  0  0  0  0  1 V2000\n"
        "    0.0    0.0    0.0 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "M  END\n$$$$\n"
    )
    pdb = tmp_path / "pA.pdb"
    pdb.write_text("ATOM      1  N   ALA A   1       0.0   0.0   0.0  1.00  0.00\nEND\n")

    mock_tool = MagicMock()
    mock_tool.dock_with_cache.return_value = []
    mock_make_tool.return_value = mock_tool

    dock_one_protein_via_tool(
        protein_path=str(pdb),
        protein_id="pA",
        protein_content_hash="somehash",
        compound_sdf_path=str(sdf),
        compound_indices=[0],
        compound_hashes={0: "chash_0"},
        grid_center=[0.0, 0.0, 0.0],
        grid_size=[10.0, 10.0, 10.0],
        backend="unidock2",
    )

    # cache_dir なしなので prepare_receptor_cache が呼ばれる
    mock_tool.prepare_receptor_cache.assert_called_once()


def test_submit_to_dask_passes_cache_dir(tmp_path):
    """_submit_to_dask が cache_dir を kwargs に含めること。"""
    protein = _make_protein("pA", "hash_aaa", tmp_path)
    ps = MagicMock()
    ps.__getitem__ = MagicMock(return_value=protein)
    ps.content_hashes.return_value = {"pA": "hash_aaa"}

    cs = MagicMock()
    cs.path = tmp_path / "compounds.sdf"
    cs.get_compound_hash.return_value = "chash_0"

    gbc = MagicMock()
    gbc.get.return_value = MagicMock(center=(0, 0, 0), size=(10, 10, 10))

    cache_dir = tmp_path / "cache"
    runner = ScreeningRunner(
        protein_set=ps,
        compound_set=cs,
        grid_box_cache=gbc,
        hdf5_path=tmp_path / "out.h5",
        cache_dir=cache_dir,
        _dock_fn=_fake_dock_one_protein,
    )

    mock_client = MagicMock()
    runner._submit_to_dask(mock_client, {"pA": [0]})

    # submit の kwargs に cache_dir が含まれること
    call_kwargs = mock_client.submit.call_args.kwargs
    assert call_kwargs.get("cache_dir") == str(cache_dir)
