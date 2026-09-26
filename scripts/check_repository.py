#!/usr/bin/env python3
"""Check the tracked data allowlist and repository-local Markdown links."""
import hashlib
import json
from pathlib import Path
import re
import subprocess

ROOT = Path(__file__).resolve().parents[1]
STUDY = "studies/wen-2026-propagator"


def main():
    paths = subprocess.check_output(["git", "ls-files", "-z"], cwd=ROOT).decode().split("\0")
    tracked = {p for p in paths if p}
    entries = json.loads((ROOT / STUDY / "manifest.json").read_text())["files"]
    expected = {f"{STUDY}/raw/dataset/{e['name']}" for e in entries}
    originals = {p for p in tracked if "/raw/" in p and not p.endswith("/raw/README.md")}
    if originals != expected:
        raise SystemExit(f"Tracked originals differ from manifest: {sorted(originals ^ expected)}")
    problems = []
    if any(p.startswith("studies/nist-bell-causal-audit/fixtures/") for p in tracked):
        problems.append("NIST source extracts require explicit authorization; use external pinned inputs")
    for path in sorted(tracked):
        if path in expected:
            # The author's README is preserved byte-for-byte, not rewritten for local links.
            if (ROOT / path).is_symlink():
                problems.append(f"Original must be a regular file: {path}")
            continue
        if path.lower().endswith((".pdf", ".zip", ".xlsx")) or "/results/workbooks/" in path:
            problems.append(f"Unapproved source/extracted file: {path}")
        if path.endswith(".md"):
            text = (ROOT / path).read_text()
            for link in re.findall(r"\]\(([^)]+)\)", text):
                if "://" in link or link.startswith(("#", "mailto:")):
                    continue
                target = ((ROOT / path).parent / link.split("#")[0]).resolve()
                if not target.is_relative_to(ROOT) or not target.exists():
                    problems.append(f"Broken/nonlocal Markdown link: {path}: {link}")
    if problems:
        raise SystemExit("\n".join(problems))
    print("Tracked data allowlist and local documentation links passed")


if __name__ == "__main__":
    main()
