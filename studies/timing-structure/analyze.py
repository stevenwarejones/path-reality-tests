#!/usr/bin/env python3
"""Reproduce every declared conditional bound, lag diagnostic and exact witness."""
import argparse
from fractions import Fraction
import json
from pathlib import Path

import numpy as np
from common import HERE, GROUPS, LAGS, FEATURES, ENVELOPES, artifact_json
from inference import (EPSILONS, ETAS, contrast, difference, tv, cancellation,
                       descriptive_score, future_result, COUNT_ALPHA, FUTURE_ALPHA,
                       COUNT_LABELS, FUTURE_TESTS, minimum_imbalance)
from witnesses import all_witnesses, quantization


def analyze(directory=HERE):
    strata, hidden, lag_scores, future, sources, anchors = [], [], [], [], {}, []
    for side in ('alice', 'bob'):
        source = json.loads((directory/f'{side}-counts.json').read_text())
        sources[side] = dict(provenance=source['provenance'], interior_rows=source['interior_rows'],
                            half_rows=source['half_rows'], group_membership=source['group_membership'],
                            group_unknown=source['group_unknown'])
        gc = np.array(source['groups']['counts'], dtype=np.int64)
        gu = np.array(source['groups']['unknown'], dtype=np.int64)
        lc = np.array(source['lags']['counts'], dtype=np.int64)
        lu = np.array(source['lags']['unknown'], dtype=np.int64)
        fu = np.array(source['lags']['future_unknown'], dtype=np.int64)
        if gc.shape != (5, 2, 2, 2, 4) or lc.shape != (9, 2, 2, 2, 4):
            raise ValueError('aggregate shape mismatch')
        if not np.array_equal(gc[..., 0], gc[..., 1:].sum(axis=-1)):
            raise ValueError('category count mismatch')
        for period, halves in [('full', [0, 1]), ('first_half', [0]), ('second_half', [1])]:
            n = sum(source['half_rows'][h] for h in halves)
            k, u = gc[:, halves].sum(axis=1), gu[:, halves].sum(axis=1)
            for y in (0, 1):
                for ei, envelope in enumerate(ENVELOPES):
                    for epsilon in EPSILONS:
                        records = []
                        for g, group in enumerate(GROUPS):
                            intervals = [contrast([(k[g, x, y, f], k[g, x, y, f]+u[g, x, y, ei])
                                                    for x in (0, 1)], n, epsilon) for f in range(4)]
                            record = dict(receiver=side, period=period, local_setting=y, group=group,
                                envelope=envelope, epsilon=epsilon, rows=n,
                                feature_contrasts=dict(zip(FEATURES, intervals)), weighted_tv=tv(intervals))
                            strata.append(record)
                            records.append(record)
                        for condition, a, b in [('clock', 1, 2), ('recovery', 3, 4)]:
                            hidden.append(dict(receiver=side, period=period, local_setting=y,
                                condition=condition, envelope=envelope, epsilon=epsilon,
                                interactions={f: difference(records[b]['feature_contrasts'][f],
                                                            records[a]['feature_contrasts'][f]) for f in FEATURES},
                                **cancellation(records[0]['weighted_tv'], records[a]['weighted_tv'], records[b]['weighted_tv'])))
            k, u = lc[:, halves].sum(axis=1), lu[:, halves].sum(axis=1)
            for j, lag in enumerate(LAGS):
                for y in (0, 1):
                    for ei, envelope in enumerate(ENVELOPES):
                        for f, feature in enumerate(FEATURES):
                            lag_scores.append(dict(receiver=side, period=period, lag=lag, local_setting=y,
                                feature=feature, envelope=envelope, kind='descriptive_not_confidence',
                                **descriptive_score(int(k[j, 0, y, f]), int(k[j, 1, y, f]),
                                                    int(u[j, 0, y, ei]), int(u[j, 1, y, ei]), n)))
        k, unknown = lc.sum(axis=1), fu.sum(axis=1)
        for j, lag in enumerate(LAGS):
            if lag <= 0:
                continue
            for y in (0, 1):
                for f, feature in enumerate(FEATURES):
                    for ei, envelope in enumerate(ENVELOPES):
                        for eta in ETAS:
                            counts = [int(k[j, x, y, f]) for x in (0, 1)]
                            uncertain = int(unknown[j, y, ei])
                            future.append(dict(receiver=side, lag=lag, local_setting=y, feature=feature,
                                envelope=envelope, known_event_counts=counts, compatible_unknown_rows=uncertain,
                                minimum_rejecting_known_imbalance=minimum_imbalance(sum(counts), uncertain, eta),
                                **future_result(*counts, uncertain, eta)))
        pooled_counts = gc[0].sum(axis=(0, 1))  # local setting, feature
        exposure = np.array(source['groups']['exposures'])[0].sum(axis=(0, 1))
        for y in (0, 1):
            p = Fraction(int(pooled_counts[y, 0]), int(exposure[y]))
            anchors.append(dict(receiver=side, local_setting=y, trusted_event_frequency=str(p),
                interpretation='hypothetical latent-law witness anchored to a descriptive frequency, not an inferred effect',
                witness=quantization(p)))
    return dict(schema_version=1, count_alpha=COUNT_ALPHA, future_alpha=FUTURE_ALPHA,
        count_labels=COUNT_LABELS, future_tests=FUTURE_TESTS, sources=sources,
        calibration_certified=False, independent_run_replication=False,
        strata=strata, cancellation=hidden, lag_scores=lag_scores, future_tests_results=future,
        witnesses=all_witnesses(), empirical_witness_anchors=anchors)


def interval(v):
    return 'incompatible' if v is None else f'[{v[0]*1e6:.1f}, {v[1]*1e6:.1f}]'


def report(r):
    lines = ['# Timing structure and what the archive can hide', '',
        'One existing run; both receivers. Retrospective, assumption-conditional analysis.',
        'The declared family separates current-setting effects, lag diagnostics and exact',
        'information-loss witnesses. No physics violation or RNG certification follows automatically.', '',
        '## Source reconciliation', '',
        '| Receiver | Reconciled rows | Interior rows | Full-overlap selected events |',
        '|---|---:|---:|---:|']
    for side, s in r['sources'].items():
        lines.append(f'| {side} | {s["provenance"]["overlap_verified"]:,} | {s["interior_rows"]:,} | {s["provenance"]["selected_events"]:,} |')
    lines += ['', 'The same 128 endpoint rows are outside every lag/history population. Interior',
        'uncertain rows and no-event rows are retained. Every local raw setting, legacy',
        'click word and timestamp excursion reconciles. No remote outcome is used.', '',
        '## Cancellation: fixed local-history conditions', '',
        'Units are ppm of ALL interior rows. Bounds are weighted stratum contributions,',
        'not probabilities conditional on a detected photon. The table uses epsilon=0',
        'and the event-supported completion assumption. All other cases are in results.json.', '',
        '| Receiver | Local setting | Condition | Sum of stratum TVs | Hidden TV beyond pooled |',
        '|---|---:|---|---:|---:|']
    selected = [v for v in r['cancellation'] if v['period'] == 'full' and v['epsilon'] == 0 and v['envelope'] == 'event_supported']
    for v in selected:
        lines.append(f'| {v["receiver"]} | {v["local_setting"]} | {v["condition"]} | {interval(v["stratified_tv"])} | {interval(v["hidden_tv"])} |')
    ideal = [v for v in r['cancellation'] if v['epsilon'] == 0]
    positive = sum(v['hidden_tv'] is not None and v['hidden_tv'][0] > 0 for v in ideal)
    interactions = [(v, f) for v in ideal for f, c in v['interactions'].items() if c is not None and (c[0] > 0 or c[1] < 0)]
    lines += ['', f'Across all {len(ideal)} ideal-assignment condition/receiver/setting/period/envelope comparisons:',
        f'**{positive} positive hidden-TV lower bounds; {len(interactions)} signed feature interactions exclude zero.**',
        'Overlapping comparisons are covered jointly, not treated as independent replications.',
        'Failure to resolve cancellation is not proof that instantaneous effects vanish.', '',
        '### Observed condition occupancy', '',
        '| Receiver | Clock state 1, known rows | Recovery state 1, known rows | Recovery state unknown |',
        '|---|---:|---:|---:|']
    for side, s in r['sources'].items():
        lines.append(f'| {side} | {sum(s["group_membership"][2]):,} | {sum(s["group_membership"][4]):,} | {sum(s["group_unknown"][4]):,} |')
    lines += ['', 'These fixed states were not balanced or optimized after seeing their occupancy.',
        'Recovery means a detector record somewhere in the previous 64 rows; it is not a',
        'measured dead-time classification. Clock uses the completed i-2 to i-1 interval.', '',
        '## Temporal fingerprints', '',
        'Lag means O_i against X_(i+lag). Positive lags use later archived settings.',
        'All lag scores are descriptive 4*(k1-k0)/N; their missingness ranges are not',
        'confidence intervals. Unequal assignment frequencies can change these scores.', '',
        '| Receiver | Local setting | Lag | Any-event score (ppm) | Event-supported completion range (ppm) |',
        '|---|---:|---:|---:|---:|']
    for v in r['lag_scores']:
        if v['period'] == 'full' and v['feature'] == 'any' and v['envelope'] == 'event_supported':
            lines.append(f'| {v["receiver"]} | {v["local_setting"]} | {v["lag"]:+d} | {v["trusted_score"]*1e6:.2f} | {interval(v["completion_range"])} |')
    lines += ['', '### Separate future-setting freshness test', '',
        'This asks whether earlier local records predict a later remote bit under the',
        'stated per-history freshness bound. It does not test retrocausality. The stake',
        'mixture and all 64 tests are paid for; unknown rows receive worst-case factors.', '',
        '| Assumed eta | Envelope | Rejections / 64 | Largest log E lower bound |',
        '|---:|---|---:|---:|']
    for eta in ETAS:
        for envelope in ENVELOPES:
            vs = [v for v in r['future_tests_results'] if v['assumed_eta'] == eta and v['envelope'] == envelope]
            lines.append(f'| {eta:g} | {envelope} | {sum(v["reject_freshness"] for v in vs)} / {len(vs)} | {max(v["log_e_lower"] for v in vs):.3f} |')
    lines += ['', '### How much imbalance would the future test need?', '',
        'Fixed example: lag +1, any-event feature, eta=0, event-supported completion.',
        'The threshold holds known-event total and uncertain rows fixed. It is a',
        'sensitivity calculation, not an additional test or a fitted effect estimate.', '',
        '| Receiver | Local setting | Known events | Uncertain candidate rows | Observed absolute imbalance | Minimum rejecting imbalance |',
        '|---|---:|---:|---:|---:|---:|']
    for v in r['future_tests_results']:
        if v['lag'] == 1 and v['feature'] == 'any' and v['assumed_eta'] == 0 and v['envelope'] == 'event_supported':
            k0, k1 = v['known_event_counts']
            lines.append(f'| {v["receiver"]} | {v["local_setting"]} | {k0+k1} | {v["compatible_unknown_rows"]} | {abs(k1-k0)} | {v["minimum_rejecting_known_imbalance"]} |')
    lines += ['', '## Exact limits on interpretation', '',
        '- **Full-tag invisibility:** opposite offsets of 0.49 bins inside each rounding cell',
        '  give a 0.98-bin event-conditioned analog shift with identical recorded tags.',
        '  Under the adopted 78.125 ps/bin unit this illustrative shift is 76.5625 ps.',
        '  The event-conditioned latent TV is one; full-trial latent TV is the event probability.',
        '  This is a quantizer-model witness, not a measured apparatus delay or sensitivity limit.',
        '- **Coarse-summary invisibility:** moving between distinct recorded tags in the same',
        '  early/core/late category is invisible to this study even though the full tags differ.',
        '- **Cancellation:** the exact toy model has pooled TV zero but state-averaged TV',
        '  200 ppm. A pooled null therefore does not exclude structured effects.',
        '- **Lag ambiguity:** an ordinary one-row delayed response with persistent settings',
        '  produces current AND future associations. Equal observational tables can also',
        '  arise from direct influence or a latent common cause with different interventions.',
        '- **Selection artifact:** conditioning on O=X manufactures a perfect association',
        '  from independent variables. Our chosen states exclude current O and X.', '',
        'Exact fractions and four measured-frequency witness anchors are in results.json.',
        'No latent law was fitted to establish these equalities.', '',
        '## Decision', '',
        'This completes the declared three-way interrogation of existing records. Its',
        'deliverables are reproducible conditional bounds, separately justified temporal',
        'diagnostics, and constructive identifiability limits. The archive alone supplies',
        'no certified assignment/calibration premise that converts them to a new physics claim.',
        'See conclusions.md for the precise scope and which additional assumptions break',
        'each ambiguity; those are limits of the result, not unimplemented analysis items.', '',
        'See protocol.md, statistics.md, sources.md and validation.md. All fixed results,',
        'including halves and uncertainty sensitivities, are committed in results.json.', '']
    return '\n'.join(lines)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--directory', type=Path, default=HERE)
    p.add_argument('--check', action='store_true')
    a = p.parse_args()
    result = analyze(a.directory)
    outputs = {'results.json': artifact_json(result), 'report.md': report(result)}
    for name, content in outputs.items():
        path = a.directory/name
        if a.check:
            if path.read_text() != content:
                raise SystemExit(name+' differs')
        else:
            path.write_text(content)
    print('timing-structure analysis '+('verified' if a.check else 'written'))


if __name__ == '__main__':
    main()
