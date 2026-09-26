"""Coverage and power checks from direct binomial sums; all counts are synthetic."""
import importlib.util
import itertools
import math
from pathlib import Path
import re
import unittest

import numpy as np
from scipy.stats import binom

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('design', ROOT/'studies/phase-intervention-design/design.py')
design = importlib.util.module_from_spec(spec)
spec.loader.exec_module(design)


class BinomialTests(unittest.TestCase):
    def test_cp_marginal_coverage_and_monotone_endpoints(self):
        for n in (1, 7, 40):
            for family in (4, 5):
                intervals = [design.binomial_interval(k, n, .1, family) for k in range(n+1)]
                self.assertEqual(intervals[0][0], 0)
                self.assertEqual(intervals[-1][1], 1)
                self.assertTrue(all(a[0] <= b[0] and a[1] <= b[1]
                                    for a, b in zip(intervals, intervals[1:])))
                for p in np.linspace(0, 1, 101):
                    risk = sum(binom.pmf(k, n, p) for k, (lo, hi) in enumerate(intervals)
                               if p < lo or p > hi)
                    self.assertLessEqual(risk, .1/family + 1e-12)

    def test_unequal_counts_simultaneous_coverage_and_null_risk(self):
        # Enumerate independently: no radius from the implementation enters this sum.
        ns, ps, alpha = [4, 6, 5, 7], [.15, .6, .35, .4], .2
        for method in ('cp', 'hoeffding'):
            intervals = [[design.binomial_interval(k, n, alpha, 4, method)
                          for k in range(n+1)] for n in ns]
            marginal = [sum(binom.pmf(k, n, p) for k, (l, u) in enumerate(ci) if l <= p <= u)
                        for n, p, ci in zip(ns, ps, intervals)]
            self.assertGreaterEqual(math.prod(marginal), 1-alpha-1e-12)
            risk = 0
            for ks in itertools.product(*(range(n+1) for n in ns)):
                ci = [intervals[i][k] for i, k in enumerate(ks)]
                lower = max(ci[0][0]-ci[2][1], ci[2][0]-ci[0][1],
                            ci[1][0]-ci[3][1], ci[3][0]-ci[1][1])
                if lower > .2 + 1e-12:
                    risk += math.prod(binom.pmf(k, n, p) for k, n, p in zip(ks, ns, ps))
            self.assertLessEqual(risk, alpha + 1e-12)
        # Original unequal-count Hoeffding rule has the same decision as clipped intervals.
        for ks in ([0, 2, 4, 1], [4, 6, 0, 0], [1, 3, 2, 4]):
            self.assertEqual(design.classify(ks, ns, .2, alpha)['decision'],
                             design.interval_classify(ks, ns, .2, alpha, 'hoeffding')['decision'])

    def test_strict_interval_threshold(self):
        ks, ns = [90, 50, 10, 50], [100]*4
        for method in ('cp', 'hoeffding'):
            first = design.interval_classify(ks, ns, 0, method=method)
            threshold = max(first['lower_absolute_contrasts'])
            self.assertGreater(threshold, 0)
            self.assertEqual(design.interval_classify(ks, ns, threshold, method=method)['decision'],
                             'inconclusive')
            self.assertEqual(design.interval_classify(ks, ns, threshold-1e-8, method=method)['decision'],
                             'reject-calibrated-null')

    def test_power_envelope_tails_and_worst_case_including_occupation(self):
        ps, ns, occ = [.8, .5, .2, .5], [13000, 11000, 15000, 12000], (.5, 14000)
        for method in ('cp', 'hoeffding'):
            r = design.certified_power(ps, ns, .02, method=method, occupation=occ)
            self.assertTrue(r['certified'])
            for p, n, lo, hi in zip(ps+[occ[0]], ns+[occ[1]], r['count_lower'], r['count_upper']):
                self.assertLessEqual(binom.cdf(lo-1, n, p), .1/10 + 1e-12)
                self.assertLessEqual(binom.sf(hi, n, p), .1/10 + 1e-12)
            # Monotone endpoints imply the envelope extrema suffice, checked at all corners.
            for ks in itertools.product(*zip(r['count_lower'], r['count_upper'])):
                result = design.interval_classify(ks[:4], ns, .02, method=method,
                                                  occupation=(ks[4], occ[1]))
                self.assertEqual(result['decision'], 'reject-calibrated-null')

    def test_small_exact_power_exceeds_certificate(self):
        # A zero-probability opposite phase reduces exact power to a single binomial sum.
        for method in ('cp', 'hoeffding'):
            plan = design.interval_plan([.8, 0, 0, 0], .02, method=method)
            n = plan['trials_per_setting']
            exact = sum(binom.pmf(k, n, .8) for k in range(n+1)
                        if design.interval_classify([k, 0, 0, 0], [n]*4, .02, method=method)
                        ['decision'] == 'reject-calibrated-null')
            self.assertGreaterEqual(exact, plan['power_lower_bound'])
            self.assertGreater(plan['contrast_lower'], plan['null_upper'])

    def test_occupation_changes_the_null_and_uses_five_intervals(self):
        ks, ns = [8000, 5000, 2000, 5000], [10000]*4
        result = design.interval_classify(ks, ns, .02, occupation=(5000, 10000))
        self.assertEqual(result['occupation_interval'], design.binomial_interval(5000, 10000, .01, 5))
        self.assertAlmostEqual(result['null_upper_bound'], result['occupation_interval'][1]+.02)
        self.assertEqual(result['decision'], 'reject-calibrated-null')
        self.assertEqual(design.interval_classify(ks, ns, .02, occupation=(8000, 10000))['decision'],
                         'inconclusive')

    def test_invalid_intervals_and_designs(self):
        for args in [(0, 0, .01), (-1, 10, .01), (11, 10, .01), (True, 10, .01),
                     (1, 10, 0), (1, 10, .01, 0), (1, 10, .01, 4, 'unknown')]:
            with self.assertRaises(ValueError):
                design.binomial_interval(*args)
        with self.assertRaises(ValueError):
            design.interval_plan([.6, .5, .4, .5], .02, w=.5)

    def test_snapshot_tolerance_preserves_discrete_results_and_structure(self):
        expected = {'plan': [{'contrast': .6, 'tiny': 0., 'trials': 8002,
                              'certified': True, 'method': 'cp'}]}
        import copy
        perturbed = copy.deepcopy(expected)
        perturbed['plan'][0].update(contrast=.6+1e-12, tiny=5e-13)
        design.compare_snapshot(perturbed, expected)
        for key, value in [('trials', 8003), ('trials', 8002.),
                           ('certified', 1), ('method', 'hoeffding'),
                           ('contrast', .61), ('contrast', float('nan'))]:
            changed = copy.deepcopy(expected)
            changed['plan'][0][key] = value
            with self.assertRaisesRegex(ValueError, rf'root\.plan\[0\]\.{key}'):
                design.compare_snapshot(changed, expected)
        with self.assertRaisesRegex(ValueError, 'keys differ at root'):
            design.compare_snapshot({}, expected)
        with self.assertRaisesRegex(ValueError, 'length differs at root.plan'):
            design.compare_snapshot({'plan': []}, expected)

    def test_question_register_covered_by_dispositions(self):
        study = ROOT/'studies/wen-2026-propagator'
        registered = set(re.findall(r'^\| (Q\d+) \|', (study/'open-questions.md').read_text(), re.M))
        report = set()
        for cell in re.findall(r'^\| (Q[^|]+) \|', (study/'completion-report.md').read_text(), re.M):
            for first, last in re.findall(r'Q(\d+)(?:[–-]Q(\d+))?', cell):
                report.update(f'Q{i:02}' for i in range(int(first), int(last or first)+1))
        self.assertTrue(registered)
        self.assertEqual(registered, report)


if __name__ == '__main__':
    unittest.main()
