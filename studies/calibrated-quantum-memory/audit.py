#!/usr/bin/env python3
"""Hash-checked NMN inventory and conditional past-marginal diagnostic.

No process reconstruction, parameter optimization, or memory claim is made.
The regrouping is an inferred outcome-dependent label permutation; see README.
"""
import argparse
from collections import defaultdict
import hashlib
import itertools
import json
import math
from pathlib import Path
import urllib.request

HERE = Path(__file__).resolve().parent
STATES = ('xp', 'xm', 'yp', 'ym', 'zp', 'zm')
AXES = ('x', 'y', 'z')
OUTCOMES = ('00', '01', '10', '11')
LABELS = {','.join(t) for t in itertools.product(STATES, AXES, STATES, AXES)}


def opposite(state):
    if state not in STATES:
        raise ValueError('invalid state')
    return state[0] + ('m' if state[1] == 'p' else 'p')


def validate_run(rows):
    if set(rows) != LABELS:
        raise ValueError('expected the full 324-circuit label grid')
    for counts in rows.values():
        if set(counts) != set(OUTCOMES):
            raise ValueError('missing or unexpected outcome')
        if any(type(n) is not int or n < 0 for n in counts.values()):
            raise ValueError('counts must be nonnegative integers')
        if sum(counts.values()) == 0:
            raise ValueError('empty row')


def regroup(rows):
    """Involution: b=1 cells move between opposite preparation labels."""
    validate_run(rows)
    result = {}
    for label, counts in rows.items():
        a, b, p, z = label.split(',')
        partner = rows[','.join((a, b, opposite(p), z))]
        result[label] = {o: (counts if o[0] == '0' else partner)[o]
                         for o in OUTCOMES}
    return result


def past_bound(rows, alpha, comparisons):
    """Uniform bound on distance to any shared past marginal.

    Each per-row Bernoulli interval has error <= alpha/comparisons. Independence
    between rows is unnecessary; independence within rows is an assumption.
    For all contexts in a group, |p_i-q|<=epsilon requires
    epsilon >= (max lower_i-min upper_i)/2.
    """
    if not 0 < alpha < 1 or comparisons < len(rows):
        raise ValueError('invalid simultaneous budget')
    groups = defaultdict(list)
    for label, counts in sorted(rows.items()):
        n = sum(counts.values())
        if n <= 0:
            raise ValueError('empty physical circuit')
        phat = (counts['00'] + counts['01']) / n
        rad = math.sqrt(math.log(2 * comparisons / alpha) / (2 * n))
        lo, hi = max(0., phat-rad), min(1., phat+rad)
        groups[tuple(label.split(',')[:2])].append((lo, hi, phat, label, n))
    records = []
    for group, values in sorted(groups.items()):
        upper_row = max(values, key=lambda r: r[0])
        lower_row = min(values, key=lambda r: r[1])
        records.append({
            'past_labels': list(group),
            'epsilon_lower': max(0., (upper_row[0]-lower_row[1])/2),
            'empirical_range': max(v[2] for v in values)-min(v[2] for v in values),
            'high_row': upper_row[3], 'low_row': lower_row[3],
            'high_frequency': upper_row[2], 'low_frequency': lower_row[2],
            'high_interval_lower': upper_row[0], 'low_interval_upper': lower_row[1],
        })
    # Deterministic tie handling also makes the chosen strongest contrast reproducible.
    return sorted(records, key=lambda r: (r['epsilon_lower'], r['empirical_range']),
                  reverse=True)


def verify_file(path, entry):
    data = path.read_bytes()
    if len(data) != entry['bytes'] or hashlib.sha256(data).hexdigest() != entry['sha256']:
        raise ValueError(f'source integrity mismatch: {path.name}')
    return data


def analyze(source_dir):
    manifest = json.loads((HERE/'sources.json').read_text())
    loaded = []
    for entry in manifest['nmn_files']:
        data = verify_file(source_dir/entry['name'], entry)
        if entry['name'].endswith('.json'):
            loaded.append((entry['name'], json.loads(data)))
    total_rows = sum(len(rows) for _, runs in loaded for rows in runs.values())
    output = {'kind': 'measured_conditional_acquisition_diagnostic',
              'alpha': .05, 'simultaneous_rows': total_rows,
              'mapping_status': 'inferred_not_verified_against_original_job_records',
              'quantum_memory_certified': False, 'runs': []}
    for filename, runs in loaded:
        for run, rows in sorted(runs.items()):
            physical = regroup(rows)
            if regroup(physical) != rows:
                raise AssertionError('regrouping must be invertible')
            totals = [sum(c.values()) for c in rows.values()]
            physical_totals = sorted({sum(c.values()) for c in physical.values()})
            expected = 9216 if filename == 'NMN_lab_rslts.json' else 8000
            if physical_totals != [expected]:
                raise ValueError('regrouped denominators differ from pinned audit')
            if sum(totals) != sum(sum(c.values()) for c in physical.values()):
                raise AssertionError('counts were lost or duplicated')
            bounds = past_bound(physical, .05, total_rows)
            output['runs'].append({
                'source': filename, 'delay_key': run, 'circuits': len(rows),
                'archived_row_total_min': min(totals),
                'archived_row_total_max': max(totals),
                'total_counts': sum(totals), 'regrouped_shots': expected,
                'simultaneous_bernoulli_radius': math.sqrt(
                    math.log(2*total_rows/.05)/(2*expected)),
                'strongest_past_marginal_contrast': bounds[0],
                'positive_groups': sum(b['epsilon_lower'] > 0 for b in bounds),
            })
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-dir', required=True, type=Path)
    parser.add_argument('--download', action='store_true')
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    if args.source_dir.resolve().is_relative_to(HERE.parents[1]):
        raise ValueError('keep source downloads outside the repository')
    if args.download:
        args.source_dir.mkdir(parents=True, exist_ok=True)
        manifest = json.loads((HERE/'sources.json').read_text())
        for entry in manifest['nmn_files']:
            path = args.source_dir/entry['name']
            if not path.exists():
                with urllib.request.urlopen(entry['url'], timeout=60) as response:
                    data = response.read()
                if hashlib.sha256(data).hexdigest() != entry['sha256']:
                    raise ValueError('download hash mismatch')
                path.write_bytes(data)
            verify_file(path, entry)
    result = json.dumps(analyze(args.source_dir), indent=2, sort_keys=True)+'\n'
    target = HERE/'results/nmn-audit.json'
    if args.check:
        if target.read_text() != result:
            raise SystemExit('audit snapshot differs')
        print('Pinned NMN counts, inverse regrouping and conditional bounds reproduce')
    else:
        target.parent.mkdir(exist_ok=True)
        target.write_text(result)
        print(target)


if __name__ == '__main__':
    main()
