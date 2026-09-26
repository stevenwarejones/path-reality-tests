#!/usr/bin/env python3
"""Verify bundled or user-downloaded Dryad files against the pinned manifest.

Verification is the main job. Automatic download is a best-effort fallback
that Dryad may block with 401/403; do not bypass its access controls.
Never accepts an HTML error page as a workbook; never changes the manifest.
No credentials, access-control workarounds, or automatic ontology conclusions.
"""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
import urllib.error
import urllib.request

ROOT = Path(__file__).resolve().parent


def verify(path, entry):
    data = path.read_bytes()
    if len(data) != entry["size"] or hashlib.sha256(data).hexdigest() != entry["sha256"]:
        raise ValueError(f"Size or SHA-256 mismatch: {path.name}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify-only", action="store_true")
    parser.add_argument("--directory", type=Path, default=ROOT / "raw" / "dataset")
    args = parser.parse_args()
    args.directory.mkdir(parents=True, exist_ok=True)
    records = []
    blocked = False
    entries = json.loads((ROOT / "manifest.json").read_text())["files"]
    expected = {e["name"] for e in entries}
    unexpected = sorted(p.name for p in args.directory.iterdir() if p.name not in expected)
    if unexpected:
        print("Unexpected dataset entries:", ", ".join(unexpected))
        return 2
    for entry in entries:
        dest = args.directory / entry["name"]
        record = {"file": entry["name"], "url": entry["url"]}
        try:
            if not dest.exists():
                if args.verify_only or blocked:
                    raise FileNotFoundError("Missing original; download disabled or access blocked")
                with urllib.request.urlopen(entry["url"], timeout=30) as response:
                    data = response.read()
                if len(data) != entry["size"] or hashlib.sha256(data).hexdigest() != entry["sha256"]:
                    raise ValueError("Downloaded bytes do not match pinned manifest")
                temp = dest.with_suffix(dest.suffix + ".part")
                temp.write_bytes(data)
                temp.replace(dest)
            verify(dest, entry)
            record["status"] = "verified"
        except (OSError, ValueError) as error:
            record.update(status="unavailable", error=str(error))
            if isinstance(error, ValueError):
                blocked = True  # never fetch further files after an integrity mismatch
            if isinstance(error, urllib.error.HTTPError) and error.code in (401, 403, 405, 429):
                blocked = True  # stop repeated requests after a control/rate-limit response
        records.append(record)
        print(entry["name"], record["status"])
    out = ROOT / "results" / "access-log.json"
    out.parent.mkdir(exist_ok=True)
    out.write_text(json.dumps({"utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
                              "files": records}, indent=2) + "\n")
    return 0 if all(r["status"] == "verified" for r in records) else 2


if __name__ == "__main__":
    raise SystemExit(main())
