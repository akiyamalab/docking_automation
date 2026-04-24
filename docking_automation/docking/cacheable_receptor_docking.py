"""受容体レベル前処理のキャッシュ API の型定義。

AutoDock Vina / Uni-Dock v1 / Uni-Dock v2 はいずれも「受容体 + grid ボックスごとに
重い前処理があり、以降の多数 ligand 処理では再利用可能」という構造を持つ。
本プロトコルはその共通インターフェースを宣言する (必須継承ではない型ヒント用 Protocol)。

具体的な重い処理:
- AutoDock Vina (CPU): `compute_vina_maps` で atom-type 別 grid potential map を計算
  (receptor サイズに依存、典型的に 5〜30 秒)。`.C.map`, `.N.map` 等の複数ファイルで保存。
- Uni-Dock 2 (GPU): `analyze_receptor_topology` (msys.LoadDMS + residue mol リスト構築)
  で ~5 分。`ud2_engine_inputs.json` の `receptor` キーに保存。
- Uni-Dock 1 (GPU): CLI の `--write_cache` / `--read_cache` でバイナリ cache 保存。

Phase 4 スケールの N×M スクリーニングでは N 受容体のキャッシュを 1 回だけ作成、
M ligand を `dock_with_cache` で高速処理する。24 core ノード上で OMP=1 並列化
すれば N 件の prep が ~N/24 × 受容体あたり処理時間で完了する。
"""
from __future__ import annotations

from pathlib import Path
from typing import List, Optional, Protocol, runtime_checkable

from docking_automation.docking.docking_result import DockingResult
from docking_automation.docking.grid_box import GridBox
from docking_automation.molecule.protein import Protein


@runtime_checkable
class CacheableReceptorDocking(Protocol):
    """受容体+grid の前処理結果をファイルキャッシュし、多数 ligand で再利用するドッキングツール。

    本 Protocol は `isinstance()` で runtime check 可能 (`@runtime_checkable`)。
    """

    def prepare_receptor_cache(
        self,
        protein: Protein,
        grid_box: GridBox,
        out_cache: Path,
        force: bool = False,
    ) -> Path:
        """受容体と grid box 情報からキャッシュを生成してファイルに保存。

        Args:
            protein: 受容体 (Protein ドメインオブジェクト)。
            grid_box: ドッキングボックス (center + size)。
            out_cache: キャッシュ出力先パス。ツールにより単一ファイルか、拡張子付き
                ファイル群の prefix かは異なる (Vina は prefix, v2 は単一 JSON)。
            force: 既存キャッシュがあっても上書きするか (既定 False)。

        Returns:
            キャッシュファイルの代表パス (複数ファイルの場合は prefix または代表)。
        """
        ...

    def dock_with_cache(
        self,
        cache: Path,
        ligand_paths: List[Path],
        grid_box: GridBox,
        protein_content_hash: str,
        compound_content_hashes: Optional[List[str]] = None,
    ) -> List[DockingResult]:
        """キャッシュを利用して多数 ligand をドッキング。

        ligand_paths の形式はツール依存:
        - AutoDock Vina: PDBQT
        - Uni-Dock 2: SDF (3D 化済み)
        - Uni-Dock 1: PDBQT

        Args:
            cache: `prepare_receptor_cache` が返したパス。
            ligand_paths: ligand ファイル群。形式はツール依存。
            grid_box: ドッキングボックス (cache と一致する必要あり)。
            protein_content_hash: 結果 `DockingResult.protein_content_hash` に格納。
            compound_content_hashes: None の場合はファイル stem を使用。

        Returns:
            各 ligand の best pose の DockingResult (1 ligand 1 件)。
        """
        ...
