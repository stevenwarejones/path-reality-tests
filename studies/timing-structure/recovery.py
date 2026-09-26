#!/usr/bin/env python3
"""Frozen synthetic recovery and false-interpretation controls, not physical power."""
import argparse
import json
from math import log
from pathlib import Path

import numpy as np
from scipy.stats import beta
from common import HERE
from inference import contrast, difference, tv, cancellation, future_result

SEED = 20260928
REPETITIONS = 400


def estimate(k, n):
    return dict(count=k, repetitions=n, rate=k/n, interval95=[
        0. if k == 0 else float(beta.ppf(.025, k, n-k+1)),
        1. if k == n else float(beta.ppf(.975, k+1, n-k))])


def category_trial(rng, mode, amplitude, n=107_109_468, state0=None):
    # Fixed balanced state populations. Each row receives independent fair X,Y.
    # Sparse rate is a declared toy value, not a fit to the archive or analog model.
    state0 = n//2 if state0 is None else state0
    counts = np.zeros((2, 2, 2, 4), dtype=np.int64)
    for g in (0, 1):
        probabilities = []
        for x in (0, 1):
            for y in (0, 1):
                p = .0003*(.5 if y == 0 else 1.5)
                if mode == 'drift_null':
                    p *= (.6 if g == 0 else 1.4)
                orientation = (-1 if g == 0 else state0/(n-state0)) if mode == 'reversing' else 1
                shift = amplitude*(2*x-1)*orientation
                probabilities.extend([.25*(1-p), .25*p*(.25-shift), .25*p*.5, .25*p*(.25+shift)])
        raw = rng.multinomial(state0 if g == 0 else n-state0, probabilities).reshape(2, 2, 4)
        counts[g, ..., 0] = raw[..., 1:].sum(axis=-1)
        counts[g, ..., 1:] = raw[..., 1:]
    pooled_found = interaction_found = hidden_found = False
    for y in (0, 1):
        records = []
        for raw in (counts.sum(axis=0), counts[0], counts[1]):
            records.append([contrast([(raw[x, y, f], raw[x, y, f]) for x in (0, 1)], n, 0)
                            for f in range(4)])
        pooled_found |= any(c[0] > 0 or c[1] < 0 for c in records[0])
        interaction_found |= any((c := difference(records[2][f], records[1][f]))[0] > 0 or c[1] < 0 for f in range(4))
        hidden_found |= cancellation(*[tv(v) for v in records])['hidden_tv'][0] > 0
    return [bool(pooled_found), bool(interaction_found), bool(hidden_found)]


def future_trial(rng, mode, n=200000):
    if mode == 'persistent_delayed':
        switches = rng.random(n+2) < .1
        x = (np.cumsum(switches) % 2).astype('i1')
        # Stationary initial label; no deterministic warm-start bias.
        x ^= rng.integers(0, 2)
    else:
        x = rng.integers(0, 2, n+2)
    y = rng.integers(0, 2, n)
    # Ordinary outcome memory: delayed dependence on the preceding remote bit.
    cause = x[:n] if mode != 'one_row_misalignment' else x[1:n+1]
    probability = .01*(1+.8*(2*cause-1))
    probability *= np.where(np.arange(n) < n//2, .6, 1.4)  # nonstationary outcomes
    detection = rng.random(n) < probability
    future = x[2:n+2] if mode != 'one_row_misalignment' else x[1:n+1]
    mask = detection & (y == 0)
    k0, k1 = int((mask & (future == 0)).sum()), int((mask & (future == 1)).sum())
    return [future_result(k0, k1, 0, eta)['reject_freshness'] for eta in (0., .4)]


def run(repetitions=REPETITIONS, seed=SEED):
    if repetitions < 1:
        raise ValueError('repetitions must be positive')
    rng = np.random.default_rng(seed)
    category = []
    for mode, amplitude in [('null', 0.), ('drift_null', 0.), ('fixed', .2), ('reversing', .1), ('reversing', .2), ('reversing', .25)]:
        hits = np.zeros(3, dtype=int)
        for _ in range(repetitions):
            hits += category_trial(rng, mode, amplitude)
        category.append(dict(mode=mode, amplitude=amplitude, outcomes={name: estimate(int(k), repetitions)
            for name, k in zip(('pooled_detected', 'interaction_detected', 'hidden_tv_resolved'), hits)}))
    occupancy = []
    # Descriptive marginal occupancy only: no remote-label association is fitted.
    for side in ('alice', 'bob'):
        source = json.loads((HERE/f'{side}-counts.json').read_text())
        for condition, g0, g1 in [('clock', 1, 2), ('recovery', 3, 4)]:
            m0, m1 = [sum(source['group_membership'][g]) for g in (g0, g1)]
            n = 107109468
            state0 = round(n*m0/(m0+m1))
            hits = np.zeros(3, dtype=int)
            for _ in range(repetitions):
                hits += category_trial(rng, 'reversing', .25, n, state0)
            occupancy.append(dict(receiver=side, condition=condition, synthetic_state0_rows=state0,
                rare_state_amplitude=.25, common_state_amplitude=.25*state0/(n-state0),
                outcomes={name: estimate(int(k), repetitions) for name, k in
                          zip(('pooled_detected', 'interaction_detected', 'hidden_tv_resolved'), hits)}))
    temporal = []
    for mode in ('fresh_with_outcome_memory', 'persistent_delayed', 'one_row_misalignment'):
        hits = np.zeros(2, dtype=int)
        for _ in range(repetitions):
            hits += future_trial(rng, mode)
        temporal.append(dict(mode=mode, outcomes={str(eta): estimate(int(k), repetitions)
                        for eta, k in zip((0., .4), hits)}))
    return dict(seed=seed, repetitions=repetitions, category_trials=107109468, future_trials=200000,
                calibration='theoretical bounds; simulations check implementation and power in declared toy models',
                category=category, occupancy_matched=occupancy, temporal=temporal)


def report(r):
    lines = ['# Frozen recovery experiments', '',
        'These are independent synthetic trials, not new physical experiments and not',
        'estimates of power against a calibrated picosecond shift. Rates and 95% exact',
        'binomial intervals are in recovery-results.json. No seeds or amplitudes were',
        'changed after seeing the results. Statistical validity comes from the proofs.', '',
        f'Seed {r["seed"]}; {r["repetitions"]} repetitions per case.', '',
        '## Sparse categorical responses', '',
        '107,109,468 trials; two balanced fixed local states; fair independent settings.',
        'Local-setting event probabilities .00015/.00045; conditional early/core/late',
        'probabilities .25/.5/.25. Amplitude moves mass between early and late in',
        'opposite directions for the two remote settings. Reversing also flips by state.',
        'The same conservative inference constants as the actual study are used.', '',
        '| Model | Amplitude | Pooled detected | State interaction detected | Hidden TV resolved |',
        '|---|---:|---:|---:|---:|']
    for v in r['category']:
        a = v['outcomes']
        lines.append(f'| {v["mode"]} | {v["amplitude"]} | {a["pooled_detected"]["count"]} | {a["interaction_detected"]["count"]} | {a["hidden_tv_resolved"]["count"]} |')
    lines += ['', '## Occupancy-matched cancellation stress test', '',
        'Additional toy designs use each observed marginal condition occupancy, with',
        'the same declared sparse event law. The rare-state amplitude is .25 and the',
        'common-state amplitude is scaled by N_rare/N_common, forcing exact cancellation',
        'in expectation. This panel uses descriptive occupancies, not fitted remote effects.',
        'It is not a calibrated physical power curve. Strong rare-state effects can be',
        'hard to detect when their weighted contribution is small.', '',
        '| Receiver | Condition | Rare-state rows | Interaction detected | Hidden TV resolved |',
        '|---|---|---:|---:|---:|']
    for v in r['occupancy_matched']:
        lines.append(f'| {v["receiver"]} | {v["condition"]} | {v["synthetic_state0_rows"]:,} | {v["outcomes"]["interaction_detected"]["count"]} | {v["outcomes"]["hidden_tv_resolved"]["count"]} |')
    lines += ['', '## Temporal controls', '',
        '200,000 trials, a drifting .01 baseline event rate, and ordinary outcome',
        'dependence on an earlier setting. One fixed future-lag/feature/local-setting',
        'test is evaluated using the full 64-test family threshold.', '',
        '| Model | Freshness rejection at eta=0 | Rejection allowing eta=.4 |',
        '|---|---:|---:|']
    for v in r['temporal']:
        lines.append(f'| {v["mode"]} | {v["outcomes"]["0.0"]["count"]} | {v["outcomes"]["0.4"]["count"]} |')
    lines += ['', 'Persistent settings have stay probability .9, so eta=.4 is necessary.',
        'Their future association is produced by ordinary memory, with no influence',
        'from the future. The misalignment case deliberately pairs the current causal',
        'setting as if it were the next row: it is an indexing/order failure, not',
        'retrocausality. The large eta=.4 row is a sensitivity illustration, not an',
        'assertion that it certifies the misaligned process. Exact quantization and',
        'coarse-category invisibility witnesses are verified separately in results.json.', '']
    return '\n'.join(lines)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--check', action='store_true')
    a = p.parse_args()
    r = run()
    for name, content in [('recovery-results.json', json.dumps(r, indent=2, sort_keys=True)+'\n'), ('recovery-report.md', report(r))]:
        if a.check:
            if (HERE/name).read_text() != content:
                raise SystemExit(name+' differs')
        else:
            (HERE/name).write_text(content)
    print('timing-structure recovery '+('verified' if a.check else 'written'))


if __name__ == '__main__':
    main()
