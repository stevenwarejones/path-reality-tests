#!/usr/bin/env python3
"""Run-preserving retrospective forecast of paired persistent readout excursions.

The target is an observed classification event during continuous measurement.
No IID-window inference or physical/QEC error probability is asserted.
"""
import argparse
import json
from pathlib import Path
import pickle

import numpy as np

from cqm_analysis import close
from cqm_sources import acquire

HERE = Path(__file__).resolve().parent
PERIOD = 0.7144833333333334
SAMPLE_PERIOD = 3e-6
PHASE_ORIGIN = -1.233 + 0.392 + 0.200124
WINDOW = 33
PERSISTENCE = 3
RUNS = 3100
BLOCK_SIZE = 50
TRAIN_BLOCKS = 31
SAMPLES = 524287  # Half-open dwell intervals; final boundary is not a new sample.
COLUMNS = ['all_windows', 'ground_A', 'ground_B', 'both_ground',
           'excursion_A_and_ground_A', 'excursion_B_and_ground_B',
           'double_excursion_and_both_ground', 'excursion_A_and_both_ground',
           'excursion_B_and_both_ground']


def native_scalar(dtype, raw):
    # Primitive native scalars reduce the otherwise very large pickle memo.
    # Never permit object, structured, or pointer-bearing scalar payloads.
    if not isinstance(dtype, np.dtype) or dtype.kind not in 'iuf' or dtype.itemsize > 8:
        raise ValueError('unsupported source scalar')
    return np._core.multiarray.scalar(dtype, raw).item()


class NumericUnpickler(pickle.Unpickler):
    def find_class(self, module, name):
        allowed = {('numpy.core.multiarray', '_reconstruct'): np._core.multiarray._reconstruct,
                   ('numpy', 'ndarray'): np.ndarray, ('numpy', 'dtype'): np.dtype,
                   ('numpy.core.multiarray', 'scalar'): native_scalar}
        if module == 'numpy._core.multiarray':
            module = 'numpy.core.multiarray'
        if (module, name) not in allowed:
            raise ValueError(f'non-numeric pickle global rejected: {module}.{name}')
        return allowed[module, name]


def dwell_states(value, samples=SAMPLES):
    """Validate intervals before numeric conversion; use deposited half-open edges."""
    if type(value) is not list or not value:
        raise ValueError('invalid dwell list')
    for row in value:
        if type(row) not in (tuple, list) or len(row) != 3 or any(type(x) is not int for x in row):
            raise ValueError('dwell entries must be integer triples')
        state, start, end = row
        if not 0 <= state <= 2 or not 0 <= start < end <= samples:
            raise ValueError('invalid dwell state or range')
    a = np.asarray(value, dtype=np.int64)
    if a[0, 1] != 0 or a[-1, 2] != samples or np.any(a[1:, 1] != a[:-1, 2]):
        raise ValueError('dwell intervals do not partition the acquisition')
    if np.any(a[1:, 0] == a[:-1, 0]):
        raise ValueError('adjacent dwell intervals have the same state')
    return np.repeat(a[:, 0].astype(np.int8), a[:, 2]-a[:, 1])


def window_events(states, window=WINDOW, persistence=PERSISTENCE):
    states = np.asarray(states)
    if states.ndim != 1 or states.dtype.kind not in 'iu':
        raise ValueError('classified states must be a one-dimensional integer array')
    if not 1 <= persistence <= window or len(states) <= window:
        raise ValueError('invalid window or persistence')
    if np.any((states < 0) | (states > 2)):
        raise ValueError('invalid classified state')
    n = (len(states)-1)//window
    starts = np.arange(n)*window
    future = states[1:1+n*window].reshape(n, window) != 0
    persistent = np.ones((n, window-persistence+1), dtype=bool)
    for j in range(persistence):
        persistent &= future[:, j:j+window-persistence+1]
    return starts, states[starts] == 0, np.any(persistent, axis=1)


def count_run(states, lag, bins):
    if len(states) != 2 or len(states[0]) != len(states[1]) or not np.isfinite(lag):
        raise ValueError('unmatched acquisition')
    start, ga, ea = window_events(states[0])
    other, gb, eb = window_events(states[1])
    if not np.array_equal(start, other):
        raise ValueError('window clocks differ')
    both = ga & gb
    phase = (start*SAMPLE_PERIOD + lag + PHASE_ORIGIN) % PERIOD
    indices = np.floor(phase/PERIOD*bins).astype(int)
    masks = [np.ones(len(start), dtype=bool), ga, gb, both, ga & ea, gb & eb,
             both & ea & eb, both & ea, both & eb]
    return np.column_stack([np.bincount(indices[m], minlength=bins) for m in masks])


def extract_blocks(cache, resolutions=(715, 1430, 3575)):
    base = cache/'mechanical/Zenodo/Fig4FigS6FigS7'
    with (base/'results_gef_2d').open('rb') as f:
        axis = NumericUnpickler(f).load()
        lags = NumericUnpickler(f).load()
        dwells = NumericUnpickler(f).load()
        if f.read():
            raise ValueError('trailing source object')
    if not isinstance(axis, np.ndarray) or axis.shape != (1000,) or axis.dtype != np.float64:
        raise ValueError('unexpected notebook time axis')
    if not np.allclose(axis, np.linspace(0, 524288*SAMPLE_PERIOD, 1000), rtol=0, atol=1e-12):
        raise ValueError('unexpected time-axis scaling')
    if type(lags) is not list or len(lags) != RUNS or not np.isfinite(lags).all():
        raise ValueError('invalid acquisition offsets')
    independent_lags = np.load(base/'lags.npy', allow_pickle=False)
    if independent_lags.shape != (RUNS,) or not np.array_equal(independent_lags, lags):
        raise ValueError('independent lag file disagrees')
    if type(dwells) is not list or len(dwells) != 2 or any(type(v) is not list or len(v) != RUNS for v in dwells):
        raise ValueError('unexpected paired acquisition schema')
    blocks = {bins: np.zeros((RUNS//BLOCK_SIZE, bins, 9), dtype=np.int64) for bins in resolutions}
    for run in range(RUNS):
        states = [dwell_states(dwells[q][run]) for q in range(2)]
        for bins in resolutions:
            blocks[bins][run//BLOCK_SIZE] += count_run(states, lags[run], bins)
        # Preserve acquisition boundaries, while releasing already processed objects.
        for q in range(2):
            dwells[q][run] = None
        if run % 500 == 0:
            print(f'Validated and counted paired acquisition {run}/{RUNS}', flush=True)
    return blocks


def fit_marginals(training, policy_bins=100):
    if training.ndim != 2 or training.shape[1] != 9 or np.any(training[:, 1:3] <= 0):
        raise ValueError('missing single-qubit calibration exposure')
    pa = training[:, 4]/training[:, 1]
    pb = training[:, 5]/training[:, 2]
    predicted = pa*pb
    n = len(training)
    if not 1 <= policy_bins <= n:
        raise ValueError('invalid policy length')
    values = []
    for i in range(n):
        selected = (i+np.arange(policy_bins)) % n
        # Both-ground exposure is an initial-state eligibility covariate.
        # No double-excursion counts or validation outcomes enter the policy.
        if training[selected, 3].sum() <= 0:
            raise ValueError('missing initial-state eligibility exposure')
        values.append(float(np.average(predicted[selected], weights=training[selected, 3])))
    chosen = int(np.argmin(values))
    return predicted, chosen


def score(counts, predicted, selected):
    if counts[selected, 3].sum() <= 0:
        raise ValueError('zero eligible paired exposure')
    exposure = int(counts[selected, 3].sum())
    observed = int(counts[selected, 6].sum())
    expected = float(counts[selected, 3] @ predicted[selected])
    return {'eligible_windows': exposure, 'observed_double_excursions': observed,
            'predicted_double_excursions': expected,
            'observed_per_million': observed/exposure*1e6,
            'predicted_per_million': expected/exposure*1e6,
            'observed_to_predicted_ratio': observed/expected if expected else None}


def intervention_clocks(cache):
    result = []
    for name in ['Fig2/774_t.npy', 'Fig2/782_t.npy', 'Fig6/833_t.npy', 'FigS9/321_t.npy']:
        a = np.load(cache/'mechanical/Zenodo'/name, allow_pickle=False)
        if a.shape != (32768,) or a.dtype != np.float64 or not np.allclose(a, np.arange(len(a))*.001, rtol=0, atol=1e-12):
            raise ValueError('unexpected intervention clock')
        result.append({'file': 'Zenodo/'+name, 'samples': len(a), 'sample_period_s': .001,
                       'span_s': float(a[-1]), 'can_directly_score_99us_target': False})
    return result


def analyze_blocks(blocks, clocks):
    counts = blocks[715]
    train = counts[:TRAIN_BLOCKS].sum(axis=0)
    test = counts[TRAIN_BLOCKS:].sum(axis=0)
    predicted, best = fit_marginals(train)
    chosen = (best+np.arange(100)) % 715
    all_bins = np.arange(715)
    validation_blocks = []
    for i, block in enumerate(counts[TRAIN_BLOCKS:]):
        validation_blocks.append({'acquisition_range_0based': [1550+50*i, 1599+50*i],
                                  'all_phases': score(block, predicted, all_bins),
                                  'fixed_quiet_policy': score(block, predicted, chosen)})
    sensitivity = []
    for bins, values in blocks.items():
        a = values[:TRAIN_BLOCKS].sum(0)
        b = values[TRAIN_BLOCKS:].sum(0)
        forecast, _ = fit_marginals(a, policy_bins=100*(bins//715))
        factor = bins//715
        selected = (best*factor + np.arange(100*factor)) % bins
        sensitivity.append({'phase_bins': bins, 'phase_bin_s': PERIOD/bins,
                            'all_phases': score(b, forecast, np.arange(bins)),
                            'same_fixed_policy': score(b, forecast, selected)})
    timing = []
    for shift in [-20, -10, -5, 0, 5, 10, 20]:
        timing.append({'phase_shift_s': shift*PERIOD/715,
                       **score(test, predicted, (chosen+shift) % 715)})
    drift = []
    for width in [1, 2, 31]:
        expectation = 0.
        for start in range(TRAIN_BLOCKS, len(counts), width):
            section = counts[start:min(start+width, len(counts))].sum(0)
            current, _ = fit_marginals(section)
            expectation += float(section[:, 3] @ current)
        drift.append({'calibration_acquisitions_per_block': width*50,
                      'expected_from_validation_marginals': expectation,
                      'observed_to_expected': int(test[:, 6].sum())/expectation,
                      'status': 'Post-validation diagnostic, not a held-out forecast or null test.'})
    blind_probability = float(train[:, 4].sum()/train[:, 1].sum() * train[:, 5].sum()/train[:, 2].sum())
    blind = np.full(715, blind_probability)
    joint_fitted = train[:, 6]/train[:, 3]
    comparators = {
        'phase_blind_independent': {'all_phases': score(test, blind, all_bins),
                                    'fixed_quiet_policy': score(test, blind, chosen),
                                    'scope': 'Algorithmic phase-information ablation, not full-source removal.'},
        'joint_outcome_fitted_benchmark': {'all_phases': score(test, joint_fitted, all_bins),
                                          'fixed_quiet_policy': score(test, joint_fitted, chosen),
                                          'scope': 'Uses training joint target counts; not a mechanism-derived prediction.'}}
    summary_all = score(test, predicted, all_bins)
    summary_quiet = score(test, predicted, chosen)
    return {'status': 'Retrospective readout-event forecast; full physical-mechanism gate remains unmet',
            'run_count': RUNS, 'samples_per_deposited_half_open_trace': SAMPLES,
            'sample_period_s': SAMPLE_PERIOD, 'window_s': WINDOW*SAMPLE_PERIOD,
            'persistence_consecutive_readouts': PERSISTENCE,
            'phase_period_s': PERIOD, 'phase_origin_s': PHASE_ORIGIN,
            'training_acquisitions_0based': [0,1549], 'validation_acquisitions_0based': [1550,3099],
            'uncertainty_status': 'Acquisition-block variation only; no IID-shot confidence or nominal future coverage.',
            'policy': {'start_bin': best, 'bins': 100, 'phase_start_s': best*PERIOD/715,
                       'duration_s': 100*PERIOD/715, 'cycle_duty_fraction': 100/715,
                       'validation_eligible_fraction': summary_quiet['eligible_windows']/summary_all['eligible_windows']},
            'training_all_phases': score(train, predicted, all_bins),
            'training_fixed_quiet_policy': score(train, predicted, chosen),
            'validation_all_phases': summary_all, 'validation_fixed_quiet_policy': summary_quiet,
            'validation_observed_rate_ratio_all_to_quiet': summary_all['observed_per_million']/summary_quiet['observed_per_million'],
            'validation_blocks': validation_blocks,
            'comparators': comparators,
            'post_validation_phase_resolution_check': sensitivity,
            'post_validation_clock_shift_check': timing,
            'post_validation_marginal_drift_check': drift,
            'intervention_clocks': clocks}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache', type=Path, required=True)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    acquire(args.cache, verify_only=True, mechanical_prediction=True)
    result = analyze_blocks(extract_blocks(args.cache), intervention_clocks(args.cache))
    output = HERE/'results/mechanical-prediction.json'
    if args.check:
        close(result, json.loads(output.read_text()))
        print('Paired acquisition counts, fixed forecasts and intervention clocks reproduced')
    else:
        output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
        print(output)
