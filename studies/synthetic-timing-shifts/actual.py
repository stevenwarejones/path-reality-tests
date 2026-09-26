#!/usr/bin/env python3
"""Conditional, memory-robust bounds for fixed actual-setting record contrasts."""
import argparse
import json
from math import expm1, isfinite, log
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
LAMBDAS = (.005, .01, .02, .04, .08, .16, .32, .64, 1., 2.)
LABELS = 2*4*4
STARTS = 2
ALPHA = .01
EPSILONS = (0., .0001, .001, .01, .05, .1)
LEAKAGES = (0., 1e-6, 1e-5, 5e-5)
FEATURES = ('any', 'early', 'core', 'late')


def count_interval(k, u, n):
    if any(type(v) is not int for v in (k, u, n)) or not 0 <= k <= u <= n or n == 0:
        raise ValueError('invalid binary count completion')
    budget = log(2*len(LAMBDAS)*LABELS*STARTS/ALPHA)
    lower = max([0.] + [(v*k-budget)/expm1(v) for v in LAMBDAS])
    upper = min([float(n)] + [(-v*u-budget)/expm1(-v) for v in LAMBDAS])
    return lower, upper


def effect_interval(arms, n, epsilon):
    if not isfinite(epsilon) or not 0 <= epsilon < .25 or len(arms) != 2:
        raise ValueError('invalid assignment sensitivity')
    intervals = [count_interval(int(k), int(u), n) for k, u in arms]
    bounds = [(max(0., l/(n*(.25+epsilon))), min(1., u/(n*(.25-epsilon))))
              for l, u in intervals]
    if any(l > u for l, u in bounds):
        return None  # incompatible confidence event/model; never call an empty set a detection
    return [max(-1., bounds[1][0]-bounds[0][1]), min(1., bounds[1][1]-bounds[0][0])]


def tv_interval(contrasts, leakage=0.):
    """Any-event absolute contrast equals the omitted no-event contrast."""
    if not isfinite(leakage) or not 0 <= leakage <= 1 or len(contrasts) != 4:
        raise ValueError('invalid TV inputs')
    if any(c is None or c[0] > c[1] for c in contrasts):
        return None
    nearest = [max(0., lo, -hi) for lo, hi in contrasts]
    farthest = [max(abs(lo), abs(hi)) for lo, hi in contrasts]
    lo = max(max(nearest), .5*sum(nearest))
    hi = min(1., .5*sum(farthest))
    if lo > hi+1e-14:
        return None
    lo = min(lo, hi)  # permit floating-point roundoff at a singleton interval
    return [max(0., lo-leakage), min(1., hi+leakage)]


def analyze(directory=HERE):
    results = []
    sources = {}
    for side in ('bob', 'alice'):
        source = json.loads((directory/f'{side}-actual-counts.json').read_text())
        sources[side] = {k: source[k] for k in ('rows', 'uncertain_rows', 'hdf_sha256', 'archive_sha256',
            'overlap_verified', 'selected_events_including_uncertain', 'multiclick_rows')}
        known = np.array(source['known_event_counts'], dtype=np.int64)
        unknown = np.array(source['compatible_uncertainty'], dtype=np.int64)
        exposures = np.array(source['trusted_exposures'], dtype=np.int64)
        assert known.shape == (2, 2, 2, 4) and unknown.shape == (2, 2, 2, 2)
        assert np.array_equal(known[..., 0], known[..., 1:].sum(axis=-1))
        assert np.all(known[..., 0] <= exposures)
        assert int(exposures.sum())+source['uncertain_rows'] == source['rows']
        for period, halves in [('full', [0, 1]), ('first_half', [0]), ('second_half', [1])]:
            n = sum(source['half_rows'][h] for h in halves)
            k, u = known[halves].sum(axis=0), unknown[halves].sum(axis=0)
            for local in (0, 1):
                for envelope_index, envelope in enumerate(('unrestricted', 'event_supported')):
                    for epsilon in EPSILONS:
                        contrasts = [effect_interval([(k[x, local, f], k[x, local, f]+u[x, local, envelope_index])
                                      for x in (0, 1)], n, epsilon) for f in range(4)]
                        tv = tv_interval(contrasts)
                        results.append(dict(receiver=side, period=period, rows=n, local_setting=local,
                            envelope=envelope, assumed_assignment_epsilon=epsilon,
                            known_event_counts=k[:, local].tolist(),
                            compatible_uncertain_rows=u[:, local, envelope_index].tolist(),
                            feature_contrasts=dict(zip(FEATURES, contrasts)), tv_interval=tv,
                            incompatible_model_or_confidence_event=tv is None,
                            nuisance_tv_intervals=[dict(assumed_sum_arm_tv=B, interval=tv_interval(contrasts, B))
                                                   for B in LEAKAGES]))
    return dict(schema_version=1, analysis_kind='retrospective_actual_remote_setting_conditional_bounds',
        alpha=ALPHA, simultaneous_labels=LABELS, deterministic_starts=STARTS,
        log_threshold=log(2*len(LAMBDAS)*LABELS*STARTS/ALPHA), sources=sources,
        calibration_certified=False, spacelike_causal_interpretation_certified=False,
        independent_run_replication=False, comparisons=results)


def interval_text(interval):
    return 'incompatible model/event' if interval is None else f'[{interval[0]*1e6:.1f}, {interval[1]*1e6:.1f}]'


def report(result):
    comparisons = result['comparisons']
    lines = ['# Actual-setting timing feasibility', '',
        'Retrospective exploration of the pinned run4 overlap. These are 99% simultaneous',
        '**assumption-conditional** confidence bounds, allowing arbitrary temporal memory.',
        'Assignment, record completeness, pairing and physical nuisance calibration remain unverified.',
        'No result below proves or disproves a law of physics.', '',
        '## Reconciliation and population', '',
        '| Receiver | Archived rows | Uncertain rows retained | Selected detector events | Multiclick rows |',
        '|---|---:|---:|---:|---:|']
    for side, source in result['sources'].items():
        lines.append(f'| {side} | {source["rows"]:,} | {source["uncertain_rows"]:,} | {source["selected_events_including_uncertain"]:,} | {source["multiclick_rows"]:,} |')
    lines += ['', 'Every local raw setting, legacy click word and timestamp excursion reconciles exactly.',
        'No-event rows remain in the denominator. The fixed outcome is no event / early / core / late',
        'for the first detector record in baseline pulse bits 4, 5, 6, without a phase-radius cut.', '',
        '## Full-record bounds under ideal assignment (epsilon = 0)', '',
        'TV is total variation between the two **time-averaged four-category record laws**.',
        'Units are parts per million of all archived trials, not of detected photons.',
        'These do not bound fine-scale timing shifts, instantaneous effects or effects that cancel over time.', '',
        '| Receiver | Local setting | Unrestricted TV interval (ppm) | Event-supported TV interval (ppm) |',
        '|---|---:|---:|---:|']
    for side in ('bob', 'alice'):
        for y in (0, 1):
            selected = [r for r in comparisons if r['receiver'] == side and r['local_setting'] == y
                        and r['period'] == 'full' and r['assumed_assignment_epsilon'] == 0]
            lines.append(f'| {side} | {y} | {interval_text(selected[0]["tv_interval"])} | {interval_text(selected[1]["tv_interval"])} |')
    lines += ['', 'The event-supported envelope assumes complete receiver detector records; the unrestricted',
        'envelope allows an event in every compatible uncertain row. Neither covers missing physical trials.', '',
        '## Fixed chronological cross-check', '',
        '| Receiver | Local setting | Period | Event-supported TV interval, epsilon=0 (ppm) |',
        '|---|---:|---|---:|']
    for r in comparisons:
        if r['period'] != 'full' and r['envelope'] == 'event_supported' and r['assumed_assignment_epsilon'] == 0:
            lines.append(f'| {r["receiver"]} | {r["local_setting"]} | {r["period"]} | {interval_text(r["tv_interval"])} |')
    lines += ['', 'Alice and Bob share this run and apparatus. This is a second-receiver cross-check,',
        '**not an independent-run replication**. The halves also share apparatus and are not independent replications.', '',
        '## Assignment sensitivity', '',
        '| Assumed per-history epsilon | Largest full-record event-supported TV upper bound (ppm) |',
        '|---:|---:|']
    for epsilon in EPSILONS:
        selected = [r for r in comparisons if r['period'] == 'full' and r['envelope'] == 'event_supported'
                    and r['assumed_assignment_epsilon'] == epsilon]
        upper = max(r['tv_interval'][1] for r in selected if r['tv_interval'] is not None)
        lines.append(f'| {epsilon:g} | {upper*1e6:.1f} |')
    ideal = [r for r in comparisons if r['assumed_assignment_epsilon'] == 0]
    incompatible = sum(r['tv_interval'] is None for r in ideal)
    nonzero = sum(r['tv_interval'] is not None and r['tv_interval'][0] > 0 for r in ideal)
    lines += ['', f'At epsilon=0: {nonzero} of {len(ideal)} fixed receiver/setting/period/envelope TV intervals',
        f'have a positive lower bound; {incompatible} have incompatible confidence/model constraints.',
        'These are overlapping comparisons with joint error control, not independent tests.',
        'Observed setting frequencies cannot certify the required per-history bound.', '',
        'For a sum-of-arm nuisance TV budget B, subtract B from the lower TV bound (floor zero)',
        'and add B to the upper bound (cap one). The JSON includes B=0, 1, 10 and 50 ppm.',
        'Those are sensitivity assumptions, not measured apparatus tolerances.', '',
        '## Feasibility decision and next gate', '',
        '- **Proceed with conditional archive methods:** both receivers reconstruct and the same fixed',
        '  analysis runs with no-clicks, temporal-memory bounds and missing-record sensitivity.',
        '- **Do not advance to a physics claim yet:** no certified per-history assignment bound,',
        '  causal nuisance budget or independent-run replication has been supplied.',
        '- Next obtain calibration evidence and an independently selected run, freeze this protocol,',
        '  then rerun it without tuning gates or categories. Decide the scientifically useful effect',
        '  threshold before examining that run; compare it with the calibrated upper bounds.', '',
        'All signed feature intervals and sensitivity comparisons are in `actual-results.json`.',
        'See `actual-protocol.md` for the derivation, estimand and limitations.', '']
    return '\n'.join(lines)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--directory', type=Path, default=HERE)
    p.add_argument('--check', action='store_true')
    a = p.parse_args()
    result = analyze(a.directory)
    for name, text in [('actual-results.json', json.dumps(result, sort_keys=True, separators=(',', ':'))+'\n'),
                       ('actual-report.md', report(result))]:
        path = a.directory/name
        if a.check:
            if path.read_text() != text:
                raise SystemExit(f'{name} differs')
        else:
            path.write_text(text)
    print('actual-setting results '+('verified' if a.check else 'written'))


if __name__ == '__main__':
    main()
