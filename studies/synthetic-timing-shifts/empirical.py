#!/usr/bin/env python3
"""Artificial assignments and shifts on a fixed empirical LOCAL background.

No actual remote settings or outcomes are read. Repeated random assignments
measure conditional injection recovery, not experimental significance.
"""
import argparse
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location('timing_synthetic_core', HERE / 'design.py')
core = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = core
SPEC.loader.exec_module(core)


def load_background(path):
    with np.load(path, allow_pickle=False) as z:
        rows, first = np.unique(z['event_row'], return_index=True)
        times = z['event_phase'][first]
        state = z['event_state'][first].astype(np.int64)
        n = int(z['source_rows'])
        exposures = z['exposures'].copy()
        provenance = json.loads(str(z['provenance']))
        if len(z['event_row']) != provenance['selected_events']:
            raise ValueError('cache event-count/provenance mismatch')
    if (not len(rows) or np.any(rows < 0) or np.any(rows >= n)
            or np.any(~np.isfinite(times)) or np.any(~np.isin(state, [0, 1]))
            or exposures.shape != (32, 2) or np.any(exposures < 0)):
        raise ValueError('invalid background arrays')
    block = rows*32//n
    fold = (block >= 8).astype(np.int64)
    regime = block % 2
    strata = 2*regime + state
    population = np.zeros((2, 4), dtype=np.int64)
    for b in range(32):
        for y in (0, 1):
            population[int(b >= 8), 2*(b % 2)+y] += exposures[b, y]
    event_counts = np.bincount(4*fold+strata, minlength=8).reshape(2, 4)
    if np.any(event_counts > population) or population.sum() != provenance['retained_trials']:
        raise ValueError('cache population conservation failure')
    return dict(times=times, fold=fold, regime=regime, strata=strata,
                population=population, event_counts=event_counts, provenance=provenance)


def assign(background, rng):
    # One artificial bit per event-bearing TRIAL, never one per photon.
    x = rng.integers(0, 2, len(background['times']), dtype=np.int64)
    group = 4*background['fold'] + background['strata']
    n1_events = np.bincount(group, weights=x, minlength=8).reshape(2, 4).astype(np.int64)
    n1_empty = rng.binomial(background['population']-background['event_counts'], .5)
    n1 = n1_events+n1_empty
    exposures = np.stack((background['population']-n1, n1), axis=2)
    return x, exposures


def histogram(background, x, shift, mode, latent=None):
    direction = np.ones(len(x), dtype=np.int64)
    if mode in ('reversing', 'digital_reversing', 'redigitized_reversing', 'centered_latent_reversing'):
        direction = 1-2*background['regime']
    elif mode == 'transfer_failure':
        direction = np.where(background['fold'] == 0, 1, 1-2*background['regime'])
    elif mode not in ('fixed', 'digital_fixed', 'redigitized_fixed', 'local_only'):
        raise ValueError('unknown injection mode')
    # The copied local timestamps all receive the same displacement per trial.
    # The first-event feature therefore shifts identically; no photon is added,
    # discarded, reranked or moved to another baseline pulse/row by this model.
    displacement = shift*direction*(2*x-1)/2
    if mode == 'local_only':
        displacement = shift*(1-2*background['regime'])
    elif mode.startswith('digital_'):
        if shift/2 != int(shift/2):
            raise ValueError('symmetric digital arm shifts must each be whole bins')
        # Symmetric integer offsets to copied detector time tags, preserving
        # the same reversal/cancellation structure as the continuous injection.
    elif mode.startswith('redigitized_'):
        if (latent is None or len(latent) != len(x) or np.any(~np.isfinite(latent))
                or np.any(latent < -.5) or np.any(latent >= .5)):
            raise ValueError('redigitization requires latent offsets in [-.5, .5]')
        # Assumed latent position within the originally recorded integer tag.
        # Rounding after a physical shift produces integer increments only.
        displacement = np.floor(latent + displacement + .5)
    elif mode == 'centered_latent_reversing':
        # Adversarial latent model: all original true times at tag-bin centers.
        # Sub-half-bin arm displacements then change no recorded time tag.
        displacement = np.floor(displacement + .5)
    times = background['times'] + displacement
    bins = np.searchsorted(core.EDGES[1:-1], times, side='right')
    index = (((background['fold']*4+background['strata'])*2+x)*18+bins)
    return np.bincount(index, minlength=2*4*2*18).reshape(2, 4, 2, 18)


def learn(counts, exposures):
    differences = (counts[:, 1]/np.maximum(1, exposures[:, 1, None]) -
                   counts[:, 0]/np.maximum(1, exposures[:, 0, None]))
    differences[np.any(exposures == 0, axis=1)] = 0
    return {'click_count': np.ones((4, 18), dtype=np.int64),
            'pooled_timing': np.tile(core.TIMING_SIGN, (4, 1)),
            'trained_timing': (np.sign(differences @ core.TIMING_SIGN)[:, None]
                               * core.TIMING_SIGN).astype(np.int64),
            'trained_histogram': np.sign(differences).astype(np.int64)}


def test_scores(counts, features):
    pvalues = {}
    for name, weights in features.items():
        positive = int(counts[:, 1][weights > 0].sum()+counts[:, 0][weights < 0].sum())
        negative = int(counts[:, 0][weights > 0].sum()+counts[:, 1][weights < 0].sum())
        pvalues[name] = core.exact_score_pvalue(positive, negative)
    return pvalues


def run(background, repeats=400, seed=20260927):
    cases = [('null', 0., 'fixed'), ('local_clock_control', 1., 'local_only')]
    cases += [(f'{mode}_{s:g}', s, mode) for mode in ('fixed', 'reversing')
              for s in (.125, .25, .5, 1.)]
    cases += [('regime_transfer_failure', 1., 'transfer_failure')]
    cases += [(f'{mode}_{s:g}', s, mode) for mode in ('digital_fixed', 'digital_reversing')
              for s in (2., 4.)]
    cases += [(f'{mode}_{s:g}', s, mode) for mode in ('redigitized_fixed', 'redigitized_reversing')
              for s in (.125, .25, .5, 1.)]
    cases += [('centered_latent_reversing_0.25', .25, 'centered_latent_reversing')]
    if repeats <= 0:
        raise ValueError('positive repeat count required')
    totals = {name: dict.fromkeys(core.TESTS, 0) for name, _, _ in cases}
    family = dict.fromkeys(totals, 0)
    # Same assignments across injections for paired comparisons. Different
    # repetitions are independent artificial relabelings of one fixed background.
    rng = np.random.Generator(np.random.PCG64(seed))
    for _ in range(repeats):
        x, exposures = assign(background, rng)
        latent = rng.uniform(-.5, .5, len(x))
        for name, shift, mode in cases:
            counts = histogram(background, x, shift, mode, latent)
            features = learn(counts[0], exposures[0])
            pvalues = test_scores(counts[1], features)
            decisions = {test: p <= .01/len(core.TESTS) for test, p in pvalues.items()}
            for test, decision in decisions.items():
                totals[name][test] += int(decision)
            family[name] += int(any(decisions.values()))
    rows = []
    for name, shift, mode in cases:
        rows.append(dict(name=name, shift_bins=shift, mode=mode, repeats=repeats,
            rejections=totals[name], family_rejections=family[name],
            intervals_95={test: core.binomial_interval(k, repeats) for test, k in totals[name].items()},
            family_interval_95=core.binomial_interval(family[name], repeats)))
    return dict(schema_version=1, data_kind='empirical_background_artificial_assignments',
                seed=seed, source=background['provenance'], family_alpha=.01,
                actual_remote_labels_used=False, experimental_significance_enabled=False,
                causal_claim_enabled=False, spacelike_claim_enabled=False, rows=rows)


def report(result):
    s = result['source']
    lines = ['# Empirical-background injection recovery', '',
             '**This is an injection study on one fixed local record, not a test of',
             'the actual remote settings and not experimental evidence of signaling.**', '',
             'See [protocol and limitations](empirical-protocol.md). Fresh independent',
             'artificial setting bits are assigned to trials. All repetitions reuse the',
             'same measured background; their intervals concern conditional recovery.', '',
             f"Receiver: {s['side']}. Verified overlap: {s['overlap_verified']:,} rows.",
             f"Retained: {s['retained_trials']:,}; exclusions before artificial assignment: "
             f"{sum(s['exclusions'].values()):,}.",
             f"Selected baseline-pulse events: {s['selected_events']:,} in {s['event_trials']:,} trials; "
             f"{s['multiclick_trials']:,} multiclick trials.",
             f"No-event trials retained: {s['no_click_trials']:,}. No phase-radius cut.", '',
             'All selected events remain in the external cache. The tested feature is',
             'the first selected event per trial in original record order. All timing',
             'bins and overflow tails survive injection. Other events in that trial',
             'do not become independent observations.', '',
             '## Actual local timing background', '',
             'Residual phases are in 78.125 ps tag-bin units relative to the configured',
             'local peak. Quantiles are descriptive; they do not establish clock accuracy.', '',
             '| Segment | Local setting | Trials | Event trials | 1% | 25% | Median | 75% | 99% |',
             '|---|---:|---:|---:|---:|---:|---:|---:|---:|']
    for row in s['background_by_segment']:
        q = ' | '.join(f'{v:.3f}' for v in row['first_event_phase_quantiles_bins'])
        lines.append(f"| {row['segment']} | {row['receiver_setting']} | {row['trials']:,} | "
                     f"{row['event_trials']:,} | {q} |")
    lines += ['', '## Conditional detection frequency', '',
              f"{result['rows'][0]['repeats']} artificial relabelings per injection; identical labels are reused",
              'across injection sizes within each repetition. Training uses the first',
              'quarter; all four tests use the remaining three quarters. Each rejects',
              'at 0.01/4. Trained features use local-setting × block-parity strata.', '',
              '| Injection | Arm shift (ps) | Count | Pooled timing | Trained timing | Trained histogram | Any test (95% MC interval) |',
              '|---|---:|---:|---:|---:|---:|---|']
    for row in result['rows']:
        rates = [100*row['rejections'][test]/row['repeats'] for test in core.TESTS]
        lo, hi = [100*p for p in row['family_interval_95']]
        shift = '0 (local-only shift)' if row['mode'] == 'local_only' else f"{row['shift_bins']*78.125:.2f}"
        lines.append(f"| {row['name']} | {shift} | "+' | '.join(f'{p:.2f}%' for p in rates)+
                     f" | {100*row['family_rejections']/row['repeats']:.2f}% ({lo:.2f}–{hi:.2f}%) |")
    lines += ['', 'The local-clock control adds a regime-dependent 78.125 ps offset to',
              'both artificial assignment arms equally. It is a null control.',
              'The transfer-failure injection changes the relation between regime and',
              'shift direction after training. Unlike the ideal Gaussian simulation,',
              'real regime weights need not balance, so pooled cancellation need not',
              'be exact and this case need not be completely invisible.', '',
              'The fixed/reversing rows use symmetric floating-point perturbations of',
              'decoded phases. Digital rows instead add symmetric whole-tag offsets,',
              'which are exact integer offsets to copied detector time tags. Redigitized',
              'rows assume uniform latent positions inside each original tag bin and',
              'round after the shift; this is an explicit, uncalibrated rounding model.',
              'The centered-latent counterexample places true arrivals at tag-bin centers;',
              'its small analog shift changes no stored tag. Sub-bin power therefore',
              'depends on an unobserved latent assumption. None establishes',
              'instrumental picosecond resolution or arbitrary-memory real-label',
              'significance. Null rates validate artificial randomization',
              'on this fixed background only. Full per-test Monte Carlo intervals and',
              'source hashes are in [empirical-results.json](empirical-results.json).', '']
    return '\n'.join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--background', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, default=HERE)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = run(load_background(args.background))
    outputs = {'empirical-results.json': json.dumps(result, indent=2, sort_keys=True)+'\n',
               'empirical-report.md': report(result)}
    if args.check:
        for name, content in outputs.items():
            if (args.output_dir/name).read_text() != content:
                raise SystemExit(f'Stale empirical-background artifact: {name}')
        print('Reproduced all empirical-background injections and report')
    else:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        for name, content in outputs.items():
            (args.output_dir/name).write_text(content)
        print('Wrote empirical-background injection results')


if __name__ == '__main__':
    main()
