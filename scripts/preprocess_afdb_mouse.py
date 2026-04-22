#!/usr/bin/env python3
"""AFDB マウス (UP000000589) v6 bulk DL + 前処理スクリプト

Usage:
    python scripts/preprocess_afdb_mouse.py [--resume] [--skip-download] [--limit N]
"""
import argparse
import gzip
import random
import shutil
import subprocess
import sys
import tarfile
from pathlib import Path

RAW_URL = "https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000000589_10090_MOUSE_v6.tar"
RAW_DIR = Path("data/afdb/raw")
EXTRACT_DIR = Path("data/afdb/extracted")
PDB_DIR = Path("data/afdb/pdb")
TAR_PATH = RAW_DIR / "UP000000589_10090_MOUSE_v6.tar"
SKIPPED_LARGE_PATH = Path("data/afdb/skipped_large.txt")
MAX_RESIDUES_DEFAULT = 2000
MIN_DISK_GB = 15


def preflight(required_gb: float = MIN_DISK_GB) -> None:
    """ディスク容量確認。不足時は RuntimeError を raise。"""
    stat = shutil.disk_usage(".")
    available_gb = stat.free / (1024 ** 3)
    print(f"[preflight] 利用可能ディスク容量: {available_gb:.1f} GB (必要: {required_gb} GB)")
    if available_gb < required_gb:
        raise RuntimeError(
            f"ディスク容量不足: {available_gb:.1f} GB < {required_gb} GB"
        )


def download(resume: bool = True) -> None:
    """wget --continue で tar を DL。サイズ検証あり。"""
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    cmd = ["wget", "--tries=3", "--timeout=60", "-O", str(TAR_PATH), RAW_URL]
    if resume:
        cmd.insert(1, "--continue")
    print(f"[download] 実行コマンド: {' '.join(cmd)}")
    result = subprocess.run(cmd)
    if result.returncode != 0:
        raise RuntimeError(f"wget 失敗 (returncode={result.returncode})")
    size_gb = TAR_PATH.stat().st_size / (1024 ** 3)
    print(f"[download] 完了。ファイルサイズ: {size_gb:.2f} GB")


def validate_tar() -> int:
    """tar -tf でエントリ数カウント。破損時は RuntimeError。"""
    print(f"[validate_tar] tar 整合性チェック中: {TAR_PATH}")
    result = subprocess.run(
        ["tar", "-tf", str(TAR_PATH)],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        raise RuntimeError(f"tar 検証失敗: {result.stderr}")
    entries = [line for line in result.stdout.splitlines() if line.strip()]
    count = len(entries)
    print(f"[validate_tar] エントリ数: {count} (期待値: ~64,800)")
    return count


def extract(skip_existing: bool = True) -> None:
    """tar 展開。pdb.gz のみ extracted/ へ。Python tarfile使用（NFS互換）。"""
    EXTRACT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"[extract] pdb.gz を {EXTRACT_DIR} へ展開中...")

    with tarfile.open(TAR_PATH) as tf:
        members = [m for m in tf.getmembers() if m.name.endswith(".pdb.gz")]
        print(f"[extract] pdb.gz エントリ数: {len(members)}")

        extracted_count = 0
        skipped_count = 0
        for member in members:
            out_path = EXTRACT_DIR / Path(member.name).name
            if skip_existing and out_path.exists():
                skipped_count += 1
                continue
            f = tf.extractfile(member)
            if f is None:
                continue
            out_path.write_bytes(f.read())
            extracted_count += 1
            if extracted_count % 1000 == 0:
                print(f"[extract] 進捗: {extracted_count}/{len(members)} 件完了")

    extracted = list(EXTRACT_DIR.rglob("*.pdb.gz"))
    print(f"[extract] 展開完了。ファイル数: {len(extracted)} (skipped={skipped_count})")


def decompress_pdb(limit: int | None = None) -> None:
    """extracted/*.pdb.gz を pdb/*.pdb に展開。limit で件数制限可。"""
    PDB_DIR.mkdir(parents=True, exist_ok=True)
    gz_files = sorted(EXTRACT_DIR.rglob("*.pdb.gz"))
    if limit is not None:
        gz_files = gz_files[:limit]
    print(f"[decompress_pdb] {len(gz_files)} 件を展開中...")

    for gz_path in gz_files:
        stem = gz_path.name[:-3]  # .pdb.gz -> .pdb
        out_path = PDB_DIR / stem
        if out_path.exists():
            continue
        with gzip.open(gz_path, "rb") as f_in, open(out_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

    count = len(list(PDB_DIR.glob("*.pdb")))
    print(f"[decompress_pdb] 完了。PDB ファイル数: {count}")


def validate(sample_n: int = 10) -> dict:
    """PDBファイル数・ランダムサンプルのパース確認。結果を dict で返す。"""
    pdb_files = list(PDB_DIR.glob("*.pdb"))
    total = len(pdb_files)
    print(f"[validate] PDB ファイル総数: {total}")

    sample = random.sample(pdb_files, min(sample_n, total)) if pdb_files else []
    parse_ok = 0
    parse_fail = []

    for path in sample:
        try:
            with open(path) as f:
                lines = f.readlines()
            atom_lines = [l for l in lines if l.startswith("ATOM") or l.startswith("HETATM")]
            if atom_lines:
                parse_ok += 1
            else:
                parse_fail.append(path.name)
        except Exception as e:
            parse_fail.append(f"{path.name}:{e}")

    result = {
        "pdb_count": total,
        "sample_n": len(sample),
        "parse_ok": parse_ok,
        "parse_fail": parse_fail,
    }
    print(f"[validate] サンプルパース結果: {result}")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description="AFDB マウス v6 bulk DL + 前処理")
    parser.add_argument("--resume", action="store_true", help="wget --continue で再開")
    parser.add_argument("--skip-download", action="store_true", help="DL をスキップ")
    parser.add_argument("--limit", type=int, default=None, help="展開する PDB 数の上限（テスト用）")
    args = parser.parse_args()

    preflight()
    if not args.skip_download:
        download(resume=args.resume)
    validate_tar()
    extract()
    decompress_pdb(limit=args.limit)
    result = validate()
    print(f"[DONE] {result}")


if __name__ == "__main__":
    main()
