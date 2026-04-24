"""Uni-Dock 2 wrapper (UniDock2Docking) のテスト。

v1 の UniDockDocking と違い、Uni-Dock 2 は SDF ligand + DMS/JSON receptor で
動作する新エンジン。本テストは以下の 2 層に分ける:

- unit: GPU + unidock2 conda env 不要 (pose パーサ / omp 設定 / cache path)
- integration: GPU + unidock2 env 必須 (pytest.mark.skipif で動的 skip)

CI (Dockerfile が unidock2 env を含む) で GPU 付きホストなら全テスト実行。
GPU 無し CI (現状の GitHub Actions 標準 runner) では integration は skip。
"""
from __future__ import annotations

import json
import os
import shutil
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from docking_automation.docking.grid_box import GridBox
from docking_automation.docking.unidock2_docking import UniDock2Docking

UNIDOCK2_ENV = Path('/opt/miniforge/envs/unidock2')
SEED_LIGAND = Path('/tmp/ud2_test/ligands_sdf/lig_000.sdf')
SEED_CACHE = Path('/tmp/ud2_test/rec0_cache.json')


def _has_gpu() -> bool:
    try:
        import subprocess
        r = subprocess.run(['nvidia-smi'], capture_output=True, timeout=5)
        return r.returncode == 0
    except Exception:
        return False


def _has_unidock_processing() -> bool:
    """現在の interpreter で unidock_processing が import できるか。

    base env で pytest を走らせている場合は False になり、integration tests が
    skip される (base env には UD2 が無い)。unidock2 env で pytest を走らせる
    場合は True になり integration tests 実行。
    """
    try:
        import importlib.util
        return importlib.util.find_spec('unidock_processing') is not None
    except Exception:
        return False


# --- unit tests (GPU 不要) ---

def test_cache_path_for_uses_content_hash(tmp_path):
    """cache_path は protein.content_hash ベース (UniProt ID が重複する
    同内容受容体を 1 cache に寄せるため)。"""
    tool = UniDock2Docking(cache_dir=tmp_path)
    mock = MagicMock()
    mock.content_hash = 'deadbeef1234'
    assert tool.cache_path_for(mock) == tmp_path / 'deadbeef1234.json'


def test_cache_path_requires_cache_dir():
    tool = UniDock2Docking()
    mock = MagicMock()
    mock.content_hash = 'abc'
    with pytest.raises(ValueError):
        tool.cache_path_for(mock)


def test_enforce_omp_single_thread_sets_env(monkeypatch):
    """OMP=1 固定は v2 並列化時の非対称 tail latency 回避に必須。"""
    monkeypatch.delenv('OMP_NUM_THREADS', raising=False)
    monkeypatch.delenv('MKL_NUM_THREADS', raising=False)
    monkeypatch.delenv('OPENBLAS_NUM_THREADS', raising=False)
    UniDock2Docking._enforce_omp_single_thread()
    assert os.environ['OMP_NUM_THREADS'] == '1'
    assert os.environ['MKL_NUM_THREADS'] == '1'
    assert os.environ['OPENBLAS_NUM_THREADS'] == '1'


def test_enforce_omp_does_not_override_user_setting(monkeypatch):
    """既に設定済みの OMP 値は尊重 (ユーザが意図して 2 等にした場合)。"""
    monkeypatch.setenv('OMP_NUM_THREADS', '2')
    UniDock2Docking._enforce_omp_single_thread()
    assert os.environ['OMP_NUM_THREADS'] == '2'  # 上書きしない


def test_receptor_path_uses_dms_if_exists(tmp_path):
    """DMS が同階層に存在すれば優先 (analyze_receptor_topology を 2.2s に
    短縮できる) されるが、無ければ PDB フォールバック。"""
    pdb = tmp_path / 'rec.pdb'
    pdb.write_text('dummy')
    # DMS 無し → PDB
    out = UniDock2Docking._receptor_path_or_dms(MagicMock(path=pdb), tmp_path)
    assert out == pdb

    dms = tmp_path / 'rec.dms'
    dms.write_text('dummy')
    # DMS 有り → DMS
    out = UniDock2Docking._receptor_path_or_dms(MagicMock(path=pdb), tmp_path)
    assert out == dms


def test_parse_pose_sdf_empty_returns_empty(tmp_path):
    """pose SDF が無い/空のケース (docking 失敗時) は空リスト。"""
    empty = tmp_path / 'nonexistent.sdf'
    results = UniDock2Docking._parse_pose_sdf(empty, [], 'hash', None)
    assert results == []


def test_parse_pose_sdf_takes_only_best_pose(tmp_path):
    """各 ligand について pose_0 (best) のみ採用、他の pose は破棄。

    Uni-Dock 2 の pose SDF は 1 ligand × N pose のフラット列で、
    `ud2_molecule_name` 属性 (例: `MOL_3_unidock2_pose_0`) で判別する。
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem

    # 2 ligand × 2 pose = 4 mol の pose SDF を合成
    sdf_path = tmp_path / 'pose.sdf'
    w = Chem.SDWriter(str(sdf_path))
    for lig_idx in range(2):
        for pose_idx in range(2):
            m = Chem.MolFromSmiles('c1ccccc1')
            m = Chem.AddHs(m)
            AllChem.EmbedMolecule(m, randomSeed=42)
            m.SetProp('ud2_molecule_name', f'MOL_{lig_idx}_unidock2_pose_{pose_idx}')
            # pose_0 は良いスコア、pose_1 は悪いスコア
            m.SetProp('vina_binding_free_energy', str(-8.0 - lig_idx + pose_idx * 2))
            w.write(m)
    w.close()

    lig_paths = [Path(f'/dummy/lig_{i}.sdf') for i in range(2)]
    results = UniDock2Docking._parse_pose_sdf(sdf_path, lig_paths, 'phash', None)
    assert len(results) == 2
    # lig_0 の best pose は -8.0、lig_1 は -9.0
    by_idx = {r.compound_index: r for r in results}
    assert by_idx[0].docking_score == pytest.approx(-8.0)
    assert by_idx[1].docking_score == pytest.approx(-9.0)
    # protein_content_hash が伝搬している
    assert all(r.protein_content_hash == 'phash' for r in results)


# --- integration tests (GPU + unidock2 env 必須) ---

_integration = pytest.mark.skipif(
    not _has_gpu()
    or not UNIDOCK2_ENV.exists()
    or not SEED_CACHE.exists()
    or not _has_unidock_processing(),
    reason='GPU or unidock2 env or seed cache or unidock_processing missing',
)


@_integration
def test_dock_with_cache_returns_valid_scores(tmp_path):
    """キャッシュ済み receptor JSON + 1 ligand で実 docking (GPU)。

    v1 で破綻した lig_019 でも v2 では正常値が返る (2026-04-24 検証)。
    """
    tool = UniDock2Docking()
    lig019 = Path('/tmp/ud2_test/ligands_sdf/lig_019.sdf')
    if not lig019.exists():
        pytest.skip('seed ligand lig_019.sdf not present')

    gb = GridBox(center=[-2.0023, -1.3722, -1.7703], size=[30.0, 30.0, 30.0])
    results = tool.dock_with_cache(
        cache_json=SEED_CACHE,
        ligand_sdf_list=[lig019],
        grid_box=gb,
        protein_content_hash='integration_rec0',
        working_dir=tmp_path,
    )
    assert len(results) == 1
    # v2 は v1 の FLT_MAX 問題を回避。正常レンジ [-15, +3] 内であること
    assert -15.0 < results[0].docking_score < 3.0


@_integration
def test_dock_with_cache_robust_timeout_propagates(tmp_path):
    """timeout が短すぎる場合は TimeoutError が上がることを確認。

    v2 docking は cached でも最低 1 秒程度かかる。
    `timeout_sec=0.05` で確実にタイムアウトさせ、retry 枯渇で例外。
    """
    tool = UniDock2Docking()
    gb = GridBox(center=[-2.0023, -1.3722, -1.7703], size=[30.0, 30.0, 30.0])
    with pytest.raises(TimeoutError):
        tool.dock_with_cache_robust(
            cache_json=SEED_CACHE,
            ligand_sdf_list=[SEED_LIGAND],
            grid_box=gb,
            protein_content_hash='integration_rec0',
            working_dir=tmp_path,
            timeout_sec=0.05,
            max_retries=0,
        )
