"""Pinned external sources; fail closed on mismatches."""
import argparse
import hashlib
import json
from pathlib import Path
import urllib.request
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def manifest():
    return json.loads((HERE/'manifest.json').read_text())


def verify(path, entry):
    data = Path(path).read_bytes()
    if len(data) != entry['bytes'] or hashlib.sha256(data).hexdigest() != entry['sha256']:
        raise ValueError(f'Source integrity failure: {Path(path).name}')
    return Path(path)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cache', type=Path, required=True)
    p.add_argument('--verify-only', action='store_true')
    p.add_argument('--portfolio', action='store_true', help='Also fetch discovery portfolio sources')
    a = p.parse_args()
    cache = a.cache.resolve()
    if cache.is_relative_to(ROOT):
        p.error('Source cache must be outside repository')
    cache.mkdir(parents=True, exist_ok=True)
    entries = manifest()['files']
    if a.portfolio:
        entries += json.loads((HERE/'portfolio-manifest.json').read_text())['files']
    for e in entries:
        path = cache/e.get('cache_name', e['name'])
        path.parent.mkdir(parents=True, exist_ok=True)
        if not path.exists() and not a.verify_only:
            temp = path.with_suffix('.download')
            try:
                with urllib.request.urlopen(e['url'], timeout=90) as response:
                    temp.write_bytes(response.read())
                verify(temp, e)
                temp.replace(path)
            finally:
                temp.unlink(missing_ok=True)
        verify(path, e)
        print('Verified', path.name)
