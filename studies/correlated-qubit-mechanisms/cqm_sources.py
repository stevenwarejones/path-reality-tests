#!/usr/bin/env python3
"""Acquire pinned external data; never execute the authors' notebooks/pickles."""
import argparse
import hashlib
import json
from pathlib import Path
import struct
import urllib.request
import zipfile
import zlib

from cqm_matched import TABLE_URL, TABLE_SHA256, extract_table

HERE = Path(__file__).resolve().parent


def digest(data):
    return hashlib.sha256(data).hexdigest()


def get(url, start=None, end=None):
    headers = {} if start is None else {"Range": f"bytes={start}-{end}"}
    with urllib.request.urlopen(urllib.request.Request(url, headers=headers), timeout=120) as r:
        if start is not None:
            if r.status != 206 or not r.headers.get("Content-Range", "").startswith(f"bytes {start}-{end}/"):
                raise ValueError("server did not honor the exact byte range")
        return r.read()


def decode_member(block, member):
    fields = struct.unpack_from("<4s5H3L2H", block)
    if fields[0] != b"PK\x03\x04" or fields[2] & 1:
        raise ValueError("invalid/encrypted ZIP local header")
    name_end = 30 + fields[-2]
    if block[30:name_end].decode() != member["name"]:
        raise ValueError("wrong ZIP member")
    pos = name_end + fields[-1]
    compressed = block[pos:pos + member["compressed_size"]]
    if fields[3] == 8:
        data = zlib.decompress(compressed, -15)
    elif fields[3] == 0:
        data = compressed
    else:
        raise ValueError("unsupported ZIP compression")
    if len(data) != member["size"] or digest(data) != member["sha256"]:
        raise ValueError("member integrity failure")
    if not fields[2] & 8 and zlib.crc32(data) != fields[6]:
        raise ValueError("ZIP CRC mismatch")
    return data


def acquire(cache, verify_only=False, full_mechanical=False, mechanical_prediction=False, mechanical_interventions=False, mechanical_transfer=False):
    if mechanical_transfer:
        mechanical_prediction = mechanical_interventions = True
    cache = cache.resolve()
    if cache.is_relative_to(HERE.parents[1]):
        raise ValueError("source cache must be outside the git repository")
    cache.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((HERE / "manifest.json").read_text())
    table_path = cache / "gamma-table2.html"
    if not table_path.exists() and not verify_only:
        table_path.write_bytes(extract_table(get(TABLE_URL).decode()))
    if digest(table_path.read_bytes()) != TABLE_SHA256:
        raise ValueError("paper Table 2 integrity failure")
    print("Verified gamma paper Table 2")
    for source in manifest["sources"]:
        if source["id"] != "mechanical":
            path = cache / source["cache_name"]
            if not path.exists() and not verify_only:
                path.write_bytes(get(source["url"]))
            data = path.read_bytes()
            if len(data) != source["bytes"] or digest(data) != source["sha256"]:
                raise ValueError(f"archive integrity failure: {path}")
            if hashlib.md5(data).hexdigest() != source["md5"]:
                raise ValueError("repository MD5 mismatch")
            with zipfile.ZipFile(path) as z:
                if {i.filename for i in z.infolist() if not i.is_dir()} != {m["name"] for m in source["members"]}:
                    raise ValueError("archive schema changed")
                for m in source["members"]:
                    if digest(z.read(m["name"])) != m["sha256"]:
                        raise ValueError("member hash mismatch")
        else:
            full = cache / source["cache_name"]
            if full_mechanical:
                if not full.exists() and not verify_only:
                    # Stream the optional 3.1 GB archive instead of buffering it.
                    with urllib.request.urlopen(source["url"], timeout=120) as r, full.open("wb") as f:
                        while chunk := r.read(8 * 1024 * 1024):
                            f.write(chunk)
                md5 = hashlib.md5()
                with full.open("rb") as f:
                    while chunk := f.read(8 * 1024 * 1024):
                        md5.update(chunk)
                if full.stat().st_size != source["bytes"] or md5.hexdigest() != source["md5"]:
                    raise ValueError("full mechanical archive integrity failure")
            tail = cache / "mechanical_tail.bin"
            cd = source["central_directory"]
            if not tail.exists() and not verify_only:
                tail.write_bytes(get(source["url"], source["bytes"] - cd["tail_bytes"], source["bytes"] - 1))
            if digest(tail.read_bytes()) != cd["tail_sha256"]:
                raise ValueError("central-directory tail hash mismatch")
            selected = source["members"] + (source.get("prediction_members", []) if mechanical_prediction else [])
            if mechanical_interventions:
                selected += source.get("intervention_members", [])
            if mechanical_transfer:
                selected += source.get("transfer_members", [])
            for m in selected:
                path = cache / "mechanical" / m["name"]
                if not path.exists() and not verify_only:
                    path.parent.mkdir(parents=True, exist_ok=True)
                    if full_mechanical:
                        with zipfile.ZipFile(full) as z:
                            path.write_bytes(z.read(m["name"]))
                    else:
                        start = m["offset"]
                        end = min(source["bytes"] - 1, start + 30 + len(m["name"].encode()) + m["compressed_size"] + 1024)
                        path.write_bytes(decode_member(get(source["url"], start, end), m))
                if len(path.read_bytes()) != m["size"] or digest(path.read_bytes()) != m["sha256"]:
                    raise ValueError(f"member integrity failure: {path}")
        total_members = len(selected) if source["id"] == "mechanical" else len(source["members"])
        print(f"Verified {source['id']}: {total_members} members")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--verify-only", action="store_true")
    parser.add_argument("--full-mechanical", action="store_true")
    parser.add_argument("--mechanical-prediction", action="store_true",
                        help="also acquire run-level dwell/phase records and intervention clocks (~68 MB compressed)")
    parser.add_argument("--mechanical-interventions", action="store_true",
                        help="also acquire calibrated 1 ms control records (~33 MB compressed)")
    parser.add_argument("--mechanical-transfer", action="store_true",
                        help="include matched PT forcing records and the prediction/control sources")
    args = parser.parse_args()
    acquire(args.cache, args.verify_only, args.full_mechanical, args.mechanical_prediction, args.mechanical_interventions, args.mechanical_transfer)
