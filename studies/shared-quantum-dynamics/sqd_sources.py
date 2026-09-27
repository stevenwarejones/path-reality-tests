"""Strict parser of deposited gate strings and count denominators; no source execution."""
import hashlib
import json
from pathlib import Path
import re
import urllib.request
import numpy as np

HERE = Path(__file__).resolve().parent
GATES = {'Gi': 0, 'Gx': 1, 'Gy': 2}


def digest(b):
    return hashlib.sha256(b).hexdigest()


def parse_word(text):
    """Only the documented literal/parenthesized-power grammar, never eval."""
    if text == '{}':
        return ()
    pos = 0
    out = []
    token = re.compile(r'G[ixy]|\((?:G[ixy])+\)\^[1-9][0-9]*')
    while pos < len(text):
        m = token.match(text, pos)
        if not m:
            raise ValueError('invalid circuit syntax')
        term = m.group()
        if term.startswith('('):
            body, exponent = term[1:].split(')^')
            n = int(exponent)
            seq = tuple(GATES[g] for g in re.findall(r'G[ixy]', body))
            if n * len(seq) > 20000:
                raise ValueError('oversized sequence')
            out.extend(seq * n)
        else:
            out.append(GATES[term])
        pos = m.end()
    if not out or len(out) > 20000:
        raise ValueError('invalid sequence length')
    return tuple(out)


def read_counts(path, family):
    lines = Path(path).read_text().splitlines()
    if lines[0] != '## Columns = plus count, count total':
        raise ValueError('unrecognized count convention')
    rows = []
    for i, line in enumerate(lines[1:]):
        fields = line.split()
        if len(fields) != 5 or fields[3:] != ['0', '0']:
            raise ValueError('unsupported trailing columns')
        word = parse_word(fields[0])
        k, n = map(float, fields[1:3])
        if not (np.isfinite(k) and np.isfinite(n) and k.is_integer() and n.is_integer() and 0 <= k <= n and n > 0):
            raise ValueError('invalid counts/denominator')
        rows.append({'family': family, 'row': i, 'word': word, 'expression': fields[0],
                     'word_sha256': digest(bytes(word)), 'k': int(k), 'n': int(n), 'length': len(word)})
    return rows


def external_cache(path):
    p = Path(path).resolve()
    repo = HERE.parents[1]
    if p == repo or repo in p.parents:
        raise ValueError('source cache must be outside repository')
    return p


def sources(cache, download=False):
    cache = external_cache(cache)
    entries = json.loads((HERE / 'manifest.json').read_text())['files']
    rows = []
    for e in entries:
        p = cache / e['cache_name']
        if download and not p.exists():
            cache.mkdir(parents=True, exist_ok=True)
            b = urllib.request.urlopen(e['url'], timeout=90).read()
            if digest(b) != e['sha256']:
                raise ValueError('download checksum mismatch')
            p.write_bytes(b)
        if digest(p.read_bytes()) != e['sha256']:
            raise ValueError('source checksum mismatch: ' + p.name)
        if e['role'] in ('GST', 'RB'):
            rows.extend(read_counts(p, e['role']))
    return rows


def split(rows):
    protocol = json.loads((HERE / 'protocol.json').read_text())
    return np.array([r['length'] <= protocol['training'][r['family'] + '_max_primitive_length'] for r in rows])


if __name__ == '__main__':
    import argparse
    a = argparse.ArgumentParser()
    a.add_argument('--cache', required=True)
    a.add_argument('--download', action='store_true')
    args = a.parse_args()
    rows = sources(args.cache, args.download)
    for fam in ('GST', 'RB'):
        rr = [r for r in rows if r['family'] == fam]
        print(fam, len(rr), 'rows;', sum(r['n'] for r in rr), 'shots')
