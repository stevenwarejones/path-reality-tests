"""Finite-sample null calibration and adversarial injection checks."""
import importlib.util
import itertools
import math
from pathlib import Path
import sys
import unittest

import numpy as np

PATH = Path(__file__).resolve().parents[1] / "studies/synthetic-timing-shifts/design.py"
SPEC = importlib.util.spec_from_file_location("synthetic_timing_design", PATH)
design = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = design
SPEC.loader.exec_module(design)


class SyntheticTimingTests(unittest.TestCase):
    def test_exact_binomial_null_is_superuniform(self):
        # Enumerate the finite null distribution using integer combinatorics,
        # independently of the scipy CDF used by the tested implementation.
        for n in range(1, 25):
            for alpha in (.001, .01, .05, .2, 1.):
                mass = sum(math.comb(n, k) for k in range(n+1)
                           if design.exact_score_pvalue(k, n-k) <= alpha) / 2**n
                self.assertLessEqual(mass, alpha + 1e-14)
        self.assertEqual(design.exact_score_pvalue(0, 0), 1)
        self.assertAlmostEqual(design.exact_score_pvalue(0, 8), 2/256)

    def test_null_calibration_after_learning_and_family_correction(self):
        # Enumerate independent fair labels for four training and eight test
        # clicks. Local timing records and regimes are fixed before assignments.
        # This checks the learned-feature/family composition, not just one CDF.
        cfg = design.Config(trials=1200, blocks=12, training_blocks=4, family_alpha=.2)
        trials = np.full((12, 2), 50, dtype=np.int64)
        rejected = 0
        for labels in itertools.product((0, 1), repeat=12):
            clicks = np.zeros((12, 2, len(design.TIMING_SIGN)), dtype=np.int64)
            for block, x in enumerate(labels):
                clicks[block, x, 6 if block % 3 else 11] = 1
            # The learner's denominators are fixed ancillary exposure counts.
            # The calibration proof requires only independence of test labels;
            # these are a conditional score fixture, not a generated trial table.
            pvalues = design.evaluate(cfg, trials, clicks)
            rejected += min(pvalues.values()) <= cfg.family_alpha / len(design.TESTS)
        self.assertLessEqual(rejected / 2**12, cfg.family_alpha)

    def test_no_click_conservation_and_tail_retention(self):
        cfg = design.Config(trials=10003)
        trials, clicks = design.simulate(cfg, design.Scenario("large_shift", shift_bins=100),
                                         np.random.default_rng(7))
        no_clicks = trials - clicks.sum(axis=2)
        self.assertEqual(int(trials.sum()), cfg.trials)
        self.assertTrue(np.all(no_clicks >= 0))
        self.assertEqual(int(no_clicks.sum() + clicks.sum()), cfg.trials)
        p = design.timing_probabilities(100, 2)
        self.assertEqual(float(p[-1]), 1)
        self.assertEqual(float(p.sum()), 1)

    def test_shift_changes_timing_but_not_click_probability(self):
        cfg = design.Config()
        for reverse in (False, True):
            case = design.Scenario("shift", shift_bins=1, reverse=reverse)
            for block in (0, 1):
                p0, t0 = design.cell_law(cfg, case, block, 0)
                p1, t1 = design.cell_law(cfg, case, block, 1)
                self.assertEqual(p0, p1)
                self.assertGreater(np.abs(t1-t0).sum(), .1)

    def test_reversing_distribution_cancels_before_regime_conditioning(self):
        cfg = design.Config()
        for case in (design.Scenario("timing", shift_bins=1, reverse=True),
                     design.Scenario("rate", rate_gap_fraction=.5, reverse=True)):
            laws = []
            for x in (0, 1):
                mixture = np.zeros(len(design.TIMING_SIGN)+1)
                for block in (0, 1):
                    p, times = design.cell_law(cfg, case, block, x)
                    mixture += np.r_[1-p, p*times] / 2
                laws.append(mixture)
            np.testing.assert_allclose(laws[0], laws[1], atol=1e-15, rtol=0)

    def test_frozen_detector_finds_reversal_and_exposes_transport_failure(self):
        cfg = design.Config(trials=640000)
        trials = np.full((32, 2), 10000, dtype=np.int64)
        clicks = np.zeros((32, 2, len(design.TIMING_SIGN)), dtype=np.int64)
        for block in range(32):
            for x in (0, 1):
                clicks[block, x, 10 if x == block % 2 else 7] = 20
        p = design.evaluate(cfg, trials, clicks)
        self.assertEqual(p["click_count"], 1)
        self.assertEqual(p["pooled_timing"], 1)
        self.assertLess(p["trained_timing"], 1e-20)
        self.assertLess(p["trained_histogram"], 1e-20)
        # Erase the regime contrast in training while leaving test data intact.
        for block in range(cfg.training_blocks):
            clicks[block] = clicks[0]
        failed = design.evaluate(cfg, trials, clicks)
        self.assertEqual(failed["trained_timing"], 1)
        self.assertEqual(failed["trained_histogram"], 1)

    def test_null_drift_has_no_assignment_dependence(self):
        cfg = design.Config()
        case = design.Scenario("drift", drift=True)
        for block in range(cfg.blocks):
            p0, t0 = design.cell_law(cfg, case, block, 0)
            p1, t1 = design.cell_law(cfg, case, block, 1)
            self.assertEqual(p0, p1)
            np.testing.assert_array_equal(t0, t1)
        self.assertNotEqual(design.cell_law(cfg, case, 0, 0)[0],
                            design.cell_law(cfg, case, 31, 0)[0])

    def test_zero_event_records_and_invalid_inputs(self):
        cfg = design.Config(trials=320)
        p = design.evaluate(cfg, np.full((32, 2), 5), np.zeros((32, 2, 18), dtype=int))
        self.assertEqual(p, dict.fromkeys(design.TESTS, 1.))
        for args in ({"trials": 0}, {"blocks": 31}, {"training_blocks": 32},
                     {"click_probability": 0}, {"timing_sigma_bins": float("nan")}):
            with self.assertRaises(ValueError):
                design.Config(**args)
        with self.assertRaises(ValueError):
            design.Scenario("bad", rate_gap_fraction=3)
        with self.assertRaises(ValueError):
            design.exact_score_pvalue(-1, 5)
        lo, hi = design.binomial_interval(0, 400)
        self.assertEqual(lo, 0)
        self.assertGreater(hi, 0)


if __name__ == "__main__":
    unittest.main()
