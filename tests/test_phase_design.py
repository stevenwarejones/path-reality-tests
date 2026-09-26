"""Prospective decision-rule checks; every fixture here is synthetic."""
import copy
import importlib.util
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
PHASE = ROOT / 'studies/phase-intervention-design'


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


design = load('phase_design', PHASE / 'design.py')
checks = load('phase_checks', PHASE / 'intervention_checks.py')
ident = load('identifiability', ROOT / 'studies/wen-2026-propagator/identifiability_checks.py')


class PhaseDesignTests(unittest.TestCase):
    def test_radius_and_strict_sample_plan(self):
        self.assertAlmostEqual(design.contrast_radius(10000, .01), .036564,
                               delta=.000001)
        plan = design.sample_plan(.1, .02)
        self.assertEqual(plan['trials_per_setting'], 6841)
        for d, b in ((.8, .02), (.1, .02), (.05, .02), (.03, .02)):
            n = design.sample_plan(d, b)['trials_per_setting']
            self.assertGreater(d, b + design.contrast_radius(n, .01) + design.contrast_radius(n, .1))
            self.assertLessEqual(d, b + design.contrast_radius(n-1, .01) + design.contrast_radius(n-1, .1))

    def test_invalid_designs_and_counts_fail(self):
        for d, b in ((.01, .02), (.02, .02), (1.1, .02), (.1, -.1), (math.nan, 0)):
            with self.assertRaises(ValueError):
                design.sample_plan(d, b)
        for n, alpha in ((0, .01), (True, .01), (1.2, .01), (10, 0), (10, 1), (10, math.nan)):
            with self.assertRaises(ValueError):
                design.contrast_radius(n, alpha)
        for counts, totals in (([1]*3, [10]*3), ([1]*4, [0]*4), ([11]*4, [10]*4),
                                ([1.5]*4, [10]*4), ([True]*4, [10]*4), ([-1]*4, [10]*4)):
            with self.assertRaises(ValueError):
                design.classify(counts, totals, .02)

    def test_real_and_imaginary_quadratures_and_unequal_trials(self):
        self.assertEqual(design.classify([8000, 4000, 0, 4000], [10000]*4, .02)['decision'],
                         'reject-calibrated-null')
        self.assertEqual(design.classify([4000, 0, 4000, 8000], [10000]*4, .02)['decision'],
                         'reject-calibrated-null')
        result = design.classify([50, 100, 150, 200], [100, 200, 300, 400], 0)
        self.assertEqual(result['decision'], 'inconclusive')
        self.assertEqual(result['estimates'], [.5]*4)
        self.assertGreater(result['per_setting_radius'][0], result['per_setting_radius'][3])

    def test_equality_at_threshold_is_inconclusive(self):
        baseline = design.classify([90, 50, 10, 50], [100]*4, 0)
        boundary = max(baseline['lower_absolute_contrasts'])
        self.assertEqual(design.classify([90, 50, 10, 50], [100]*4, boundary)['decision'],
                         'inconclusive')

    def test_exact_binomial_null_risk(self):
        # Independently enumerate the two pair statistics for a fixed iid null.
        n, alpha = 40, .1
        for p in (.05, .5, .95):
            pmf = [math.comb(n, k)*p**k*(1-p)**(n-k) for k in range(n+1)]
            pair_risk = sum(pmf[i]*pmf[j] for i in range(n+1) for j in range(n+1)
                            if abs(i/n-j/n) > math.sqrt(2*math.log(8/alpha)/n))
            four_setting_risk = 1-(1-pair_risk)**2
            self.assertLessEqual(four_setting_risk, alpha)

    def test_power_bound_against_independent_binomial_sum(self):
        plan = design.sample_plan(.8, .02)
        n = plan['trials_per_setting']
        threshold = plan['rejection_contrast_threshold']
        # Alternative p0=.8, ppi=0; this pair alone is sufficient for rejection.
        exact_pair_power = sum(math.comb(n, k)*.8**k*.2**(n-k)
                               for k in range(n+1) if k/n > threshold)
        self.assertGreaterEqual(exact_pair_power, plan['power_lower_bound'])

    def test_snapshot_check_rejects_changed_values_missing_keys_and_nan(self):
        expected = {'phase': [.1, .2], 'cases': 48, 'kind': 'synthetic'}
        for key, value in [('phase', [.1, .3]), ('phase', [.1, math.nan]), ('cases', True)]:
            changed = copy.deepcopy(expected)
            changed[key] = value
            with self.assertRaises(ValueError):
                checks.check_snapshot(expected, changed)
        with self.assertRaises(ValueError):
            checks.check_snapshot(expected, {'phase': [.1, .2]})

    def test_identifiability_counterexamples_are_exact_and_distinct(self):
        result = ident.analyze()
        self.assertEqual(result['shared_factor_pairing']['aligned_product_mean'], '101/100')
        self.assertEqual(result['shared_factor_pairing']['reversed_product_mean'], '99/100')
        self.assertEqual(result['phase_magnitude_pairing']['coherent_intensities'], ['8', '10'])

    def test_external_directory_checks_leave_snapshots_unchanged(self):
        scripts = [PHASE / 'design.py', PHASE / 'intervention_checks.py',
                   ROOT / 'studies/wen-2026-propagator/identifiability_checks.py']
        snapshots = [PHASE / 'results/design.json', PHASE / 'results/intervention.json',
                     ROOT / 'studies/wen-2026-propagator/results/identifiability.json']
        before = [p.read_bytes() for p in snapshots]
        with tempfile.TemporaryDirectory() as tmp:
            for script in scripts:
                result = subprocess.run([sys.executable, str(script), '--check'], cwd=tmp,
                                        capture_output=True, text=True)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(before, [p.read_bytes() for p in snapshots])


if __name__ == '__main__':
    unittest.main()
