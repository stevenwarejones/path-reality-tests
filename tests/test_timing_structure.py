"""Independent oracles and adversarial checks for the fixed archive interrogation."""
import importlib.util
from fractions import Fraction as F
from itertools import product
import json
from math import exp, expm1, log
from pathlib import Path
import sys
import tempfile
import unittest

import h5py
import numpy as np

HERE = Path(__file__).resolve().parents[1]/'studies/timing-structure'
sys.path.insert(0, str(HERE))
try:
    def load(name, file):
        spec = importlib.util.spec_from_file_location(name, HERE/file)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module
    extract = load('structure_extract_test', 'extract.py')
    inference = load('structure_inference_test', 'inference.py')
    witnesses = load('structure_witnesses_test', 'witnesses.py')
finally:
    sys.path.pop(0)


class StructureStatisticsTests(unittest.TestCase):
    def test_paid_families_and_count_inversion(self):
        self.assertEqual(inference.COUNT_LABELS, 160)
        self.assertEqual(inference.FUTURE_TESTS, 64)
        self.assertEqual(inference.COUNT_ALPHA+inference.FUTURE_ALPHA, .01)
        b = log(2*len(inference.LAMBDAS)*160*2/.005)
        lo, hi = inference.count_interval(300, 320, 1000000)
        for v in inference.LAMBDAS:
            self.assertLessEqual(v*300-expm1(v)*lo, b+1e-9)
            self.assertLessEqual(-v*320-expm1(-v)*hi, b+1e-9)

    def test_adaptive_freshness_mixture_expectation(self):
        # All paths, history-dependent propensity, and predictable event selection.
        for eta in (0., .1):
            expectation = 0.
            for bits in product((0, 1), repeat=9):
                probability = 1.
                counts = [0, 0]
                for i, bit in enumerate(bits):
                    q = .5+eta*(1 if sum(bits[:i]) % 2 else -1)
                    probability *= q if bit else 1-q
                    selected = i >= 1 and (bits[i-1] == 1 or i % 3 == 0)
                    if selected:
                        counts[bit] += 1
                expectation += probability*exp(inference.future_log_e(*counts, 0, eta))
            self.assertLessEqual(expectation, 1.+1e-12)

    def test_unknown_betting_is_pathwise_conservative(self):
        baseline = inference.future_log_e(2, 3, 4, .01)
        # Unknown rows can be no event, event with bit 0, or event with bit 1.
        for completion in product((0, 1, 2), repeat=4):
            actual = inference.future_log_e(2+completion.count(1), 3+completion.count(2), 0, .01)
            self.assertLessEqual(baseline, actual+1e-12)

    def test_sensitivity_and_missingness_do_not_increase_evidence(self):
        previous = inference.future_log_e(100, 200, 0, 0)
        for eta in (.0001, .001, .01, .4):
            current = inference.future_log_e(100, 200, 0, eta)
            self.assertLessEqual(current, previous)
            previous = current
        self.assertLess(inference.future_log_e(100, 200, 10, 0), inference.future_log_e(100, 200, 0, 0))

    def test_minimum_imbalance_is_exact_at_integer_boundary(self):
        for unknown in (0, 20, 200):
            n = 1000
            minimum = inference.minimum_imbalance(n, unknown, 0)
            self.assertIsNotNone(minimum)
            majority = (n+minimum)//2
            self.assertTrue(inference.future_result(n-majority, majority, unknown, 0)['reject_freshness'])
            self.assertFalse(inference.future_result(n-majority+1, majority-1, unknown, 0)['reject_freshness'])
        self.assertIsNone(inference.minimum_imbalance(0, 1, 0))

    def test_tv_and_hidden_tv_for_exact_reversing_laws(self):
        pooled = inference.tv([[0., 0.]]*4)
        # Weighted contribution of half the rows, event moves early->late.
        a = inference.tv([[0., 0.], [-.001, -.001], [0., 0.], [.001, .001]])
        b = inference.tv([[0., 0.], [.001, .001], [0., 0.], [-.001, -.001]])
        np.testing.assert_allclose(inference.cancellation(pooled, a, b)['hidden_tv'], [.002, .002])
        self.assertIsNone(inference.cancellation(None, a, b)['hidden_tv'])

    def test_assignment_bias_and_completions_widen_contrasts(self):
        initial = inference.contrast([(1000, 1000), (1200, 1200)], 10000000, 0)
        wide = inference.contrast([(950, 1050), (1150, 1250)], 10000000, .01)
        self.assertLessEqual(wide[0], initial[0])
        self.assertGreaterEqual(wide[1], initial[1])

    def test_descriptive_completion_enumeration(self):
        result = inference.descriptive_score(3, 5, 2, 4, 100)
        lo, hi = result['completion_range']
        for a in range(3):
            for b in range(5):
                score = 4*(5+b-3-a)/100
                self.assertLessEqual(lo, score)
                self.assertGreaterEqual(hi, score)

    def test_empty_and_invalid_inputs(self):
        self.assertAlmostEqual(inference.future_log_e(0, 0, 0, 0), 0.)
        self.assertEqual(inference.future_result(0, 0, 0, 0)['family_adjusted_p_upper'], 1.)
        for args in [(1, 0, 10), (0, 11, 10), (0, 0, 0)]:
            with self.assertRaises(ValueError):
                inference.count_interval(*args)
        with self.assertRaises(ValueError):
            inference.future_log_e(1, 1, -1, 0)
        with self.assertRaises(ValueError):
            inference.future_log_e(1, 1, 0, .5)


class StructureExtractionTests(unittest.TestCase):
    def test_history_excludes_current_and_uses_exact_window(self):
        values = np.zeros(130, dtype=int)
        values[[0, 64, 129]] = 1
        got = extract.prior_any(values)
        want = [values[i-64:i].any() for i in range(64, len(values))]
        np.testing.assert_array_equal(got, want)
        self.assertTrue(got[0])
        self.assertFalse(got[-1])  # record at current row 129 cannot enter its own history

    def test_unknown_remote_label_is_not_double_counted_in_future_test(self):
        k, e, u, future = extract.summarize(np.array([1]), np.array([0]), np.array([False]),
            np.array([2]), np.array([True]), np.array([0]))
        self.assertEqual(k.sum(), 0)
        self.assertEqual(u[0, :, 0, 1].sum(), 2)
        self.assertEqual(future[0, 0, 1], 1)

    def test_unknown_state_may_enter_both_groups_but_no_trusted_group(self):
        args = (np.array([1]), np.array([1]), np.array([False]), np.array([2]), np.array([True]), np.array([0]))
        for member in (False, True):
            k, e, u, f = extract.summarize(*args, np.array([member]), np.array([True]))
            self.assertEqual(e.sum(), 0)
            self.assertEqual(u[0, 0, 0, 1], 1)

    def test_full_aggregation_matches_row_oracle_and_chunk_boundaries(self):
        rng = np.random.default_rng(17)
        n = 273  # odd midpoint; every fixed lag fits the interior
        local = rng.integers(1, 3, n, dtype='u1')
        remote = rng.integers(1, 3, n, dtype='u1')
        local[86], remote[185] = 0, 3
        data = np.zeros(n, dtype=extract.CACHE_DTYPE)
        data['category'] = rng.integers(0, 4, n)
        data['detector'] = (data['category'] > 0) | (rng.random(n) < .1)
        data['clock'] = rng.integers(0, 2, n)
        bad = ~np.isin(local, (1, 2)) | ~np.isin(remote, (1, 2))
        bad[140:147] = True
        shapes = [(2, 2, 2, 4), (2, 2, 2), (2, 2, 2, 2), (2, 2, 2)]
        expected_g = [np.zeros((5,)+s, dtype='i8') for s in shapes]
        expected_l = [np.zeros((9,)+s, dtype='i8') for s in shapes]
        def add(targets, index, i, j, membership, state_unknown):
            hh = int(i >= n//2)
            uncertain = bool(bad[i] or bad[j] or state_unknown)
            if not uncertain and membership:
                x, y, c = int(remote[j])-1, int(local[i])-1, int(data['category'][i])
                targets[1][index, hh, x, y] += 1
                if c:
                    targets[0][index, hh, x, y, 0] += 1
                    targets[0][index, hh, x, y, c] += 1
            elif uncertain and (membership or state_unknown):
                for y in (0, 1):
                    if local[i] in (1, 2) and local[i] != y+1:
                        continue
                    targets[3][index, hh, y, 0] += 1
                    targets[3][index, hh, y, 1] += int(data['detector'][i])
                    for x in (0, 1):
                        if remote[j] in (1, 2) and remote[j] != x+1:
                            continue
                        targets[2][index, hh, x, y, 0] += 1
                        targets[2][index, hh, x, y, 1] += int(data['detector'][i])
        for i in range(64, n-64):
            clock, recovery = bool(data['clock'][i-1]), bool(data['detector'][i-64:i].any())
            for g in range(5):
                membership = True if g == 0 else (clock == (g == 2) if g < 3 else recovery == (g == 4))
                unknown = False if g == 0 else (bool(bad[i-2:i].any()) if g < 3 else bool(bad[i-64:i].any()))
                add(expected_g, g, i, i, membership, unknown)
            for j, lag in enumerate(extract.LAGS):
                add(expected_l, j, i, i+lag, True, False)
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            hdf = directory/'source.hdf5'
            with h5py.File(hdf, 'w') as h:
                h['alice/settings'], h['bob/settings'] = local, remote
                jumps = np.zeros((4, 2), dtype='i8')
                jumps[0] = [140, 145]
                jumps[3] = [-100000000, 100001600]
                h['alice/badSyncInfo'] = jumps
                h['bob/badSyncInfo'] = np.zeros((4, 0), dtype='i8')
            cache = directory/'alice-structure.npy'
            np.save(cache, data)
            (directory/'alice-structure.json').write_text(json.dumps(dict(rows=n,
                hdf_sha256=extract.decoder.sha256(hdf), cache_sha256=extract.decoder.sha256(cache))))
            baseline = None
            for chunk in (1, 17, 1000):
                got = extract.aggregate(hdf, directory, 'alice', chunk)
                for family, expected in [('groups', expected_g), ('lags', expected_l)]:
                    for key, values in zip(('counts', 'exposures', 'unknown', 'future_unknown'), expected):
                        np.testing.assert_array_equal(got[family][key], values)
                self.assertEqual(got['interior_rows'], n-128)
                self.assertEqual(sum(got['half_rows']), n-128)
                if baseline is not None:
                    self.assertEqual(got, baseline)
                baseline = got
            data['category'][100] ^= 1
            np.save(cache, data)
            with self.assertRaisesRegex(ValueError, 'hash mismatch'):
                extract.aggregate(hdf, directory, 'alice')


class StructureWitnessTests(unittest.TestCase):
    def test_quantization_exact_for_negative_and_positive_cells(self):
        for tag in range(-3, 4):
            for sign in (-1, 1):
                self.assertEqual((F(tag)+sign*F(49, 100)+F(1, 2)).__floor__(), tag)
        for p in (F(0), F(1, 1000), F(1)):
            w = witnesses.quantization(p)
            self.assertEqual(F(w['recorded_tv']), 0)
            self.assertEqual(F(w['latent_tv']), p)
        with self.assertRaises(ValueError):
            witnesses.quantization(offset=F(1, 2))

    def test_arbitrary_recorded_arm_laws_have_nonunique_lifts(self):
        w = witnesses.arbitrary_record_lifts()
        self.assertEqual(w['recorded_law_change_in_either_arm'], '0')
        self.assertEqual(F(w['unidentified_gap_difference_bins']), F(49, 50))
        self.assertNotEqual(*w['latent_event_mean_gaps'])

    def test_coarse_summary_has_distinct_fine_records(self):
        w = witnesses.coarse_graining()
        self.assertEqual(F(w['coarse_category_tv']), 0)
        self.assertGreater(F(w['fine_record_tv']), 0)

    def test_cancellation_lag_aliasing_and_collider(self):
        w = witnesses.all_witnesses()
        self.assertEqual(F(w['cancellation']['pooled_tv']), 0)
        self.assertGreater(F(w['cancellation']['hidden_tv']), 0)
        lag = w['temporal_aliasing']
        self.assertEqual(F(lag['observational_tv_between_causal_models']), 0)
        self.assertNotEqual(lag['direct_model_interventional_gap'], lag['common_cause_interventional_gap'])
        self.assertGreater(F(lag['lag_contrasts']['1']), 0)
        self.assertEqual(w['selection_collider']['selected_gap'], '1')


if __name__ == '__main__':
    unittest.main()
