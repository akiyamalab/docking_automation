#!/usr/bin/env python
"""AFDB PDB files を UniProt ID 先頭2文字でシャーディングする."""
import argparse
import shutil
from pathlib import Path


def shard_pdb_files(src_dir: Path, dst_dir: Path, dry_run: bool = False) -> dict:
    """src_dir 内の AF-*.pdb* を UniProt ID 先頭2文字のバケットに振り分ける.

    src_dir == dst_dir の場合はその場でサブディレクトリへ移動（in-place sharding）。
    Returns: {"moved": int, "buckets": set}
    """
    moved = 0
    buckets: set = set()

    pdb_files = sorted(src_dir.glob("AF-*.pdb*"))
    if not pdb_files:
        print(f"[shard] No AF-*.pdb* files found in {src_dir}")
        return {"moved": 0, "buckets": set()}

    for pdb_file in pdb_files:
        parts = pdb_file.name.split("-")
        if len(parts) < 2:
            print(f"[shard] SKIP (unexpected name): {pdb_file.name}")
            continue
        bucket = parts[1][:2]
        dest = dst_dir / bucket / pdb_file.name
        if dest == pdb_file:
            continue
        buckets.add(bucket)
        if dry_run:
            print(f"[dry-run] {pdb_file} -> {dest}")
        else:
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(pdb_file), str(dest))
        moved += 1

    if not dry_run:
        print(f"Sharding complete: {dst_dir} ({moved} files, {len(buckets)} buckets)")
    else:
        print(f"[dry-run] Would move {moved} files into {len(buckets)} buckets under {dst_dir}")
    return {"moved": moved, "buckets": buckets}


def main() -> None:
    parser = argparse.ArgumentParser(description="Shard AFDB PDB files by UniProt ID prefix")
    parser.add_argument("src_dir", type=Path, help="Source directory containing AF-*.pdb* files")
    parser.add_argument(
        "dst_dir",
        type=Path,
        nargs="?",
        help="Destination root (default: same as src_dir for in-place sharding)",
    )
    parser.add_argument("--dry-run", action="store_true", help="Print moves without executing")
    args = parser.parse_args()

    src = args.src_dir.resolve()
    dst = args.dst_dir.resolve() if args.dst_dir else src

    result = shard_pdb_files(src, dst, dry_run=args.dry_run)
    print(f"Result: {result['moved']} files moved into {len(result['buckets'])} buckets")


if __name__ == "__main__":
    main()
