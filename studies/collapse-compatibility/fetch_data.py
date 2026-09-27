"""Download immutable source objects into a cache outside the repository."""
import argparse
import json
from pathlib import Path
import urllib.request
from io_data import HERE,digest,verify_sources


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--data-dir",type=Path,required=True)
    p.add_argument("--verify-only",action="store_true")
    args = p.parse_args()
    directory = args.data_dir.resolve()
    if directory.is_relative_to(HERE.parents[1]): raise ValueError("Use external cache")
    if not args.verify_only:
        directory.mkdir(parents=True,exist_ok=True)
        for item in json.loads((HERE/"manifest.json").read_text())["files"]:
            target = directory/item["name"]
            if target.exists(): continue
            with urllib.request.urlopen(item["url"],timeout=60) as response:
                raw = response.read(item["bytes"]+1024 if item.get("normalization") else item["bytes"]+1)
            if item.get("normalization") == "horizons_generation_timestamp":
                from motion_sources import normalize_horizons
                raw = normalize_horizons(raw)
            if len(raw)!=item["bytes"] or digest(raw)!=item["sha256"]: raise ValueError("Download integrity failure")
            target.write_bytes(raw)
    verify_sources(directory)
    print("All external source objects match pinned hashes")


if __name__=="__main__": main()
