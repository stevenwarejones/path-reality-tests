"""Actual-label inference: adaptive memory, missing rows and distribution bounds."""
import importlib.util
from itertools import product
from math import exp, expm1, log
from pathlib import Path
import sys
import unittest

import numpy as np

HERE = Path(__file__).resolve().parents[1] / 'studies/synthetic-timing-shifts'
sys.path.insert(0, str(HERE))
try:
    import actual
    import actual_extract
finally:
    sys.path.pop(0)


class ActualTimingTests(unittest.TestCase):
    def test_bernoulli_supermartingale_with_adaptive_memory(self):
        # Enumerate a strongly non-IID process: propensity depends on all past bits.
        for lam in (-2., -.08, .08, 2.):
            expectation = 0.
            crossing = 0.
            for bits in product((0, 1), repeat=8):
                probability, score, crossed = 1., 0., False
                for i, bit in enumerate(bits):
                    mu = .03 if sum(bits[:i]) % 2 else .91
                    probability *= mu if bit else 1-mu
                    score += lam*bit-expm1(lam)*mu
                    crossed |= score >= log(5.)
                expectation += probability*exp(score)
                crossing += probability*crossed
            self.assertLessEqual(expectation, 1.+1e-12)
            self.assertLessEqual(crossing, .2+1e-12)

    def test_count_inversion_matches_all_paid_boundaries(self):
        budget = log(2*len(actual.LAMBDAS)*actual.LABELS*actual.STARTS/actual.ALPHA)
        lo, hi = actual.count_interval(1000, 1030, 100000)
        for v in actual.LAMBDAS:
            self.assertLessEqual(v*1000-expm1(v)*lo, budget+1e-9)
            self.assertLessEqual(-v*1030-expm1(-v)*hi, budget+1e-9)

    def test_completion_and_assignment_uncertainty_only_widen(self):
        baseline = actual.effect_interval([(100, 100), (130, 130)], 1000000, 0)
        for arms, epsilon in [([(90, 120), (110, 150)], 0),
                              ([(100, 100), (130, 130)], .01)]:
            expanded = actual.effect_interval(arms, 1000000, epsilon)
            self.assertLessEqual(expanded[0], baseline[0])
            self.assertGreaterEqual(expanded[1], baseline[1])

    def test_biased_assignment_can_explain_apparent_difference(self):
        # Same potential response rate in both arms, strongly biased assignment.
        n = 100000000
        arms = [(150000, 150000), (350000, 350000)]
        ideal = actual.effect_interval(arms, n, 0)
        allowed = actual.effect_interval(arms, n, .1)
        self.assertGreater(ideal[0], 0)
        self.assertLessEqual(allowed[0], 0)
        self.assertGreaterEqual(allowed[1], 0)

    def test_large_fixed_effect_is_detectable_without_iid_test(self):
        interval = actual.effect_interval([(1000, 1000), (5000, 5000)], 100000000, 0)
        self.assertGreater(interval[0], 0)
        self.assertLess(interval[0], .00016)
        self.assertGreater(interval[1], .00016)

    def test_missing_rows_keep_noevents_and_compatible_contexts(self):
        local = np.array([1, 1, 2, 0, 2, 2])
        remote = np.array([1, 2, 1, 2, 0, 2])
        bad = np.array([False, False, True, False, False, False])
        category = np.array([0, 2, 1, 3, 0, 1])
        detector = np.array([False, True, True, True, False, True])
        half = np.array([0, 0, 0, 1, 1, 1])
        known, exposure, uncertainty, missing = actual_extract.summarize_rows(
            local, remote, bad, category, detector, half)
        self.assertEqual(missing, 3)
        self.assertEqual(exposure.sum(), 3)  # includes the trusted no-event row
        self.assertEqual(known[..., 0].sum(), 2)
        np.testing.assert_array_equal(known[..., 0], known[..., 1:].sum(axis=-1))
        self.assertEqual(uncertainty[..., 0].sum(), 5)  # ambiguity permits multiple contexts
        self.assertEqual(uncertainty[..., 1].sum(), 3)
        self.assertEqual(uncertainty[0, 0, 1, 0], 1)
        self.assertEqual(uncertainty[0, 1, 1, 0], 0)

    def test_tv_bounds_cover_known_probability_pairs(self):
        rng = np.random.default_rng(81)
        for _ in range(100):
            a, b = rng.dirichlet(np.ones(4), 2)
            delta = b-a
            features = [-delta[0], *delta[1:]]
            contrasts = [[max(-1, d-.03), min(1, d+.03)] for d in features]
            lo, hi = actual.tv_interval(contrasts)
            tv = .5*np.abs(delta).sum()
            self.assertLessEqual(lo, tv+1e-12)
            self.assertGreaterEqual(hi, tv-1e-12)
            exact = actual.tv_interval([[d, d] for d in features])
            np.testing.assert_allclose(exact, [tv, tv], atol=1e-12)

    def test_nuisance_budget_is_paid_once_and_invalid_inputs_fail(self):
        intervals = [[-.001, .001]]*4
        self.assertAlmostEqual(actual.tv_interval(intervals, .00005)[1], .00205)
        self.assertIsNone(actual.tv_interval([None]*4))
        for args in [(-1, 1, 10), (2, 1, 10), (0, 11, 10), (0, 0, 0)]:
            with self.assertRaises(ValueError):
                actual.count_interval(*args)
        with self.assertRaises(ValueError):
            actual.effect_interval([(0, 0)]*2, 10, .25)


if __name__ == '__main__':
    unittest.main()
