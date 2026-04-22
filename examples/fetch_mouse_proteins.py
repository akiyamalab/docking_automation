"""
Fetch 10 mouse (Mus musculus) protein structures from AlphaFold DB.
- 3 size buckets: small(100-150), medium(200-300), large(350-500)
- Excluded: AF-A0A0A2JW93, AF-A0A0H2UNG0 (rule violation: >500 residues)
"""

import json
import time
from pathlib import Path

import requests

WORK_DIR = Path(__file__).parent
OUTPUT_DIR = WORK_DIR / "input" / "afdb_mouse"
UNIPROT_API = "https://rest.uniprot.org/uniprotkb/search"
AFDB_API = "https://alphafold.ebi.ac.uk/api/prediction/{}"
SLEEP_SEC = 0.5

EXCLUDED_IDS = {"AF-A0A0A2JW93", "AF-A0A0H2UNG0"}

BUCKETS = [
    {"name": "small",  "query": "organism_id:10090 AND reviewed:true AND length:[100 TO 150]", "target": 4},
    {"name": "medium", "query": "organism_id:10090 AND reviewed:true AND length:[200 TO 300]", "target": 3},
    {"name": "large",  "query": "organism_id:10090 AND reviewed:true AND length:[350 TO 500]", "target": 3},
]


def fetch_uniprot_ids_for_bucket(query: str, n: int) -> list[str]:
    """Fetch reviewed mouse UniProt IDs for a given query."""
    params = {
        "query": query,
        "format": "list",
        "size": n * 5,
        "sort": "length asc",
    }
    resp = requests.get(UNIPROT_API, params=params, timeout=30)
    resp.raise_for_status()
    ids = [line.strip() for line in resp.text.splitlines() if line.strip()]
    return ids


def fetch_pdb_url(uniprot_id: str) -> str | None:
    """Return the PDB download URL from AFDB, or None if not available."""
    url = AFDB_API.format(uniprot_id)
    try:
        resp = requests.get(url, timeout=30)
        if resp.status_code == 404:
            return None
        resp.raise_for_status()
        data = resp.json()
        if isinstance(data, list) and data:
            return data[0].get("pdbUrl")
        return None
    except requests.exceptions.RequestException as e:
        print(f"[AFDB] {uniprot_id}: request error — {e}")
        return None


def download_pdb(pdb_url: str, dest: Path) -> bool:
    """Download PDB file to dest. Returns True on success."""
    try:
        resp = requests.get(pdb_url, timeout=60)
        resp.raise_for_status()
        dest.write_bytes(resp.content)
        return True
    except requests.exceptions.RequestException as e:
        print(f"[DL] {pdb_url}: {e}")
        return False


def fetch_bucket(bucket: dict, already_collected: set[str]) -> list[str]:
    """Fetch proteins for one bucket. Returns list of downloaded PDB paths."""
    name = bucket["name"]
    target = bucket["target"]
    print(f"\n[Bucket:{name}] target={target}")

    uniprot_ids = fetch_uniprot_ids_for_bucket(bucket["query"], target)
    print(f"[UniProt:{name}] {len(uniprot_ids)} IDs retrieved")

    paths: list[str] = []
    for uid in uniprot_ids:
        if len(paths) >= target:
            break

        # 除外IDチェック（UniProt IDベース）
        af_id = f"AF-{uid}"
        if af_id in EXCLUDED_IDS:
            print(f"[skip] {uid}: excluded")
            continue
        if uid in already_collected:
            print(f"[skip] {uid}: already collected")
            continue

        time.sleep(SLEEP_SEC)
        pdb_url = fetch_pdb_url(uid)
        if pdb_url is None:
            print(f"[skip] {uid}: no AFDB entry")
            continue

        filename = pdb_url.split("/")[-1]
        # ファイル名ベースでも除外チェック
        stem = filename.split("-F1-")[0] if "-F1-" in filename else filename
        if stem in EXCLUDED_IDS:
            print(f"[skip] {uid}: excluded (by filename)")
            continue

        dest = OUTPUT_DIR / filename
        time.sleep(SLEEP_SEC)
        if dest.exists():
            print(f"[exist] {uid} → {filename} (skip download)")
            paths.append(str(dest))
            already_collected.add(uid)
        elif download_pdb(pdb_url, dest):
            paths.append(str(dest))
            already_collected.add(uid)
            print(f"[ok]   {uid} → {filename}")
        else:
            print(f"[fail] {uid}: download failed")

    print(f"[Bucket:{name}] collected {len(paths)}/{target}")
    return paths


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    already_collected: set[str] = set()
    all_paths: list[str] = []

    for bucket in BUCKETS:
        paths = fetch_bucket(bucket, already_collected)
        all_paths.extend(paths)

    list_path = OUTPUT_DIR / "protein_list.json"
    list_path.write_text(json.dumps(all_paths, indent=2, ensure_ascii=False))

    print(f"\n=== Result ===")
    print(f"Total: {len(all_paths)} proteins")
    for p in all_paths:
        print(f"  {Path(p).name}")
    print(f"protein_list.json → {list_path}")

    if len(all_paths) < 10:
        print(f"[WARN] Expected 10 proteins, got {len(all_paths)}")


if __name__ == "__main__":
    main()
