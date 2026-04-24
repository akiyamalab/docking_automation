"""ドッキング結果スコアの解析・可視化 (matplotlib)。

入力ソース 2 種類に対応:
  1. HDF5 リポジトリ (output/phase4_hdf5/*.h5) — Phase 4 本番形式
  2. Uni-Dock 生出力ディレクトリ (--dir で書かれた *_out.pdbqt 群) — ベンチ等暫定データ

使用例:
    python3 examples/analyze_scores.py --hdf5-dir output/phase4_hdf5
    python3 examples/analyze_scores.py --pdbqt-root /tmp/mps_bench/grid
"""
from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path
from typing import Iterable, List, Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

VINA_RESULT_RE = re.compile(r"REMARK VINA RESULT:\s+(-?\d+\.\d+)")


def parse_pdbqt_score(pdbqt: Path) -> Optional[float]:
    try:
        with pdbqt.open() as f:
            for line in f:
                m = VINA_RESULT_RE.match(line)
                if m:
                    return float(m.group(1))
                if not line.startswith("REMARK") and not line.startswith("MODEL"):
                    return None
    except OSError:
        return None
    return None


def load_from_pdbqt_tree(root: Path) -> pd.DataFrame:
    rows: List[dict] = []
    for pdbqt in root.rglob("*_out.pdbqt"):
        score = parse_pdbqt_score(pdbqt)
        if score is None:
            continue
        lig = pdbqt.stem.replace("_out", "")
        rel = pdbqt.relative_to(root)
        group = rel.parts[0] if len(rel.parts) > 1 else ""
        run = rel.parts[1] if len(rel.parts) > 2 else ""
        rows.append(
            {
                "group": group,
                "run": run,
                "compound_id": lig,
                "score": score,
            }
        )
    return pd.DataFrame(rows)


def load_from_hdf5(hdf5_dir: Path) -> pd.DataFrame:
    import h5py

    rows: List[dict] = []
    for h5 in sorted(hdf5_dir.glob("*.h5")):
        with h5py.File(h5, "r") as f:
            if "results" not in f:
                continue
            for ph in f["results"]:
                for ch in f[f"results/{ph}"]:
                    g = f[f"results/{ph}/{ch}"]
                    rows.append(
                        {
                            "shard": h5.name,
                            "protein_id": g.attrs.get("protein_id", ""),
                            "protein_hash": ph,
                            "compound_hash": ch,
                            "compound_index": int(g.attrs.get("compound_index", 0)),
                            "score": float(g["docking_score"][()]),
                        }
                    )
    return pd.DataFrame(rows)


def plot_score_histogram(df: pd.DataFrame, out_png: Path) -> None:
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(df["score"], bins=60, color="steelblue", edgecolor="white")
    ax.set_xlabel("Docking score [kcal/mol]")
    ax.set_ylabel("# pairs")
    ax.set_title(f"Score distribution (N={len(df)})")
    ax.axvline(df["score"].median(), color="red", linestyle="--",
               label=f"median={df['score'].median():.2f}")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_png, dpi=120)
    plt.close(fig)


def plot_score_heatmap(df: pd.DataFrame, out_png: Path,
                       index_col: str, columns_col: str,
                       max_rows: int = 40, max_cols: int = 40) -> None:
    pivot = df.pivot_table(index=index_col, columns=columns_col,
                           values="score", aggfunc="min")
    if pivot.shape[0] > max_rows:
        pivot = pivot.loc[pivot.min(axis=1).nsmallest(max_rows).index]
    if pivot.shape[1] > max_cols:
        pivot = pivot.loc[:, pivot.min(axis=0).nsmallest(max_cols).index]

    fig, ax = plt.subplots(figsize=(min(12, 0.3 * pivot.shape[1] + 3),
                                    min(10, 0.3 * pivot.shape[0] + 2)))
    im = ax.imshow(pivot.values, aspect="auto", cmap="viridis_r")
    ax.set_xticks(range(pivot.shape[1]))
    ax.set_xticklabels(pivot.columns, rotation=90, fontsize=7)
    ax.set_yticks(range(pivot.shape[0]))
    ax.set_yticklabels(pivot.index, fontsize=7)
    ax.set_xlabel(columns_col)
    ax.set_ylabel(index_col)
    ax.set_title(f"Docking score heatmap (top {pivot.shape[0]}×{pivot.shape[1]})")
    fig.colorbar(im, ax=ax, label="score [kcal/mol]")
    fig.tight_layout()
    fig.savefig(out_png, dpi=120)
    plt.close(fig)


def plot_top_receptors(df: pd.DataFrame, out_png: Path,
                       group_col: str, top_n: int = 30) -> None:
    top1 = (df.groupby(group_col)["score"].min()
            .sort_values().head(top_n))
    fig, ax = plt.subplots(figsize=(6, 0.25 * len(top1) + 1))
    ax.barh(range(len(top1)), top1.values, color="steelblue")
    ax.set_yticks(range(len(top1)))
    ax.set_yticklabels(top1.index, fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel("best score [kcal/mol]")
    ax.set_title(f"Top {top_n} by {group_col}")
    fig.tight_layout()
    fig.savefig(out_png, dpi=120)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--hdf5-dir", type=Path, help="HDF5 シャードディレクトリ")
    src.add_argument("--pdbqt-root", type=Path, help="Uni-Dock 生出力ルート")
    ap.add_argument("--out-dir", type=Path,
                    default=Path(__file__).parent / "output" / "analysis")
    ap.add_argument("--min-score", type=float, default=-20.0,
                    help="下限 (計算破綻値 除外用、既定 -20 kcal/mol)")
    ap.add_argument("--max-score", type=float, default=5.0,
                    help="上限 (計算破綻値 除外用、既定 5 kcal/mol)")
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    if args.hdf5_dir:
        df = load_from_hdf5(args.hdf5_dir)
        group_col = "protein_id"
        columns_col = "compound_index"
    else:
        df = load_from_pdbqt_tree(args.pdbqt_root)
        group_col = "group"
        columns_col = "compound_id"

    if df.empty:
        print("[ERROR] スコアが 1 件も読み込めませんでした。")
        return

    n_raw = len(df)
    df = df[(df["score"] >= args.min_score) & (df["score"] <= args.max_score)]
    n_filt = len(df)
    if n_filt < n_raw:
        print(f"計算破綻スコア除外: {n_raw} -> {n_filt} pairs "
              f"(range [{args.min_score}, {args.max_score}])")

    csv_path = args.out_dir / "scores.csv"
    df.to_csv(csv_path, index=False)

    print(f"総ペア数: {len(df)}")
    print(f"スコア統計:\n{df['score'].describe()}")
    print(f"CSV: {csv_path}")

    plot_score_histogram(df, args.out_dir / "score_histogram.png")
    plot_score_heatmap(df, args.out_dir / "score_heatmap.png",
                       index_col=group_col, columns_col=columns_col)
    plot_top_receptors(df, args.out_dir / "top_receptors.png",
                       group_col=group_col)
    print(f"図版: {args.out_dir}/*.png")


if __name__ == "__main__":
    main()
