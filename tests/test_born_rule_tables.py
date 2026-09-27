"""Scientific failure modes: mixture convention, physical certificates, rank and source parsing."""
import csv
import copy
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
import zipfile
import numpy as np

STUDY = Path(__file__).resolve().parents[1]/'studies/born-rule-identifiability'
sys.path.insert(0, str(STUDY))
from br_tables import (directions, h, model, noise_radius, centered_singular,
                       validate_fit, probabilities, spectral_boundary)
from br_table_sources import read_archive
from joint_tables import ideal_checks, check


class TablePhysicsTests(unittest.TestCase):
    def test_secondary_report_roundoff_does_not_change_rendering(self):
        from br_combination_report import render
        audit = json.loads((STUDY/'results/portfolio-audit.json').read_text())
        result = json.loads((STUDY/'results/combination-feasibility.json').read_text())
        changed = copy.deepcopy(result)
        for row in changed['qubit_equivalence']:
            row['max_exact_family_residual'] = 4.440892098500626e-16
        self.assertEqual(render(audit, result), render(audit, changed))
        reordered = copy.deepcopy(audit)
        reordered['ibm'] = dict(reversed(list(reordered['ibm'].items())))
        self.assertEqual(render(audit, result), render(reordered, result))
        changed['qubit_equivalence'][0]['max_exact_family_residual'] = 1e-5
        self.assertNotEqual(render(audit, result), render(audit, changed))

    def test_viviani_exact_signed_polynomial(self):
        r, m = directions('viviani')
        # Both signs, full declared parameter domain, fixed outcome/order.
        for theta in [-1, -.1, -.01, 0, .01, .1, 1]:
            p = (1+h(r@m.T, theta))/2
            self.assertAlmostEqual(np.linalg.det(np.c_[p, np.ones(5)]),
                                   -3*theta*(3*theta+16)/1024, places=13)
        self.assertAlmostEqual((1+r@m.T)[0, 0]/2, .75, places=13)

    def test_full_harmonic_rank_not_every_design(self):
        checks = ideal_checks()
        for t, values in checks['generic_singular_values'].items():
            self.assertEqual(np.count_nonzero(np.array(values) > 1e-10), 4 if float(t) == 0 else 11)
        # Repeating one direction cannot magically attain rank 11.
        x = np.ones((16, 16))*.3
        self.assertEqual(np.linalg.matrix_rank((1+h(x, .1))/2), 1)

    def test_antipodal_ensemble_is_not_nonlinear_mixed_probability(self):
        c, x, t = .4, .6, .2
        ensemble = ((1+c)/2*(1+h(x, t))+(1-c)/2*(1+h(-x, t)))/2
        self.assertAlmostEqual(ensemble, (1+c*h(x, t))/2)
        self.assertGreater(abs(ensemble-(1+h(c*x, t))/2), .001)

    def test_arbitrary_physical_qubit_has_centered_rank_three(self):
        rng = np.random.default_rng(64)
        r = rng.normal(size=(13, 3)); r /= np.maximum(1, np.linalg.norm(r, axis=1))[:, None]
        m = rng.normal(size=(7, 3)); m /= np.maximum(1, np.linalg.norm(m, axis=1))[:, None]
        p = (1+r@m.T)/2
        self.assertLess(np.linalg.norm(centered_singular(p)[3:]), 1e-14)

    def test_drift_mixture_can_increase_rank(self):
        rng = np.random.default_rng(924)
        tables = []
        for _ in range(2):
            r = rng.normal(size=(9, 3)); r /= np.linalg.norm(r, axis=1)[:, None]
            m = rng.normal(size=(4, 3)); m /= np.linalg.norm(m, axis=1)[:, None]
            tables.append((1+r@m.T)/2)
        self.assertTrue(all(centered_singular(p)[-1] < 1e-14 for p in tables))
        self.assertGreater(centered_singular(np.mean(tables, axis=0))[-1], .05)

    def test_subgaussian_radius_and_coverage_regression(self):
        shots = np.full((5, 4), 10000)
        radius = noise_radius(shots, .01)
        self.assertGreater(noise_radius(shots/4, .01), radius*1.99)
        rng = np.random.default_rng(905)
        draws = rng.binomial(shots, .5, size=(2000, 5, 4))/shots-.5
        # Conservative theorem, tested away from the nominal tail; not a proof of coverage.
        self.assertLessEqual(np.count_nonzero(np.linalg.norm(draws, axis=(1, 2)) > radius), 20)

    def test_published_count_and_physical_certificate_reproduction(self):
        audit = json.loads((STUDY/'results/table-audit.json').read_text())
        result = json.loads((STUDY/'results/joint-tables.json').read_text())
        check(result, audit)
        scan = audit['records']['scan']
        self.assertEqual(scan['missing_indices_within_recorded_range'], [10])
        self.assertEqual(scan['job_indices'][-1], 115)
        self.assertFalse(scan['job_list_complete'])
        p = probabilities(audit['records']['viviani']['counts_0_1_by_job'])
        self.assertAlmostEqual(np.linalg.det(np.c_[p, np.ones(5)]), -.0002982960473816028, places=14)
        for kind in ['scan', 'viviani']:
            row = result['records'][kind]
            fit = min(row['leakage'], key=lambda x: x['max_abs_residual'])
            p = np.array(row['mean_table'])
            q = validate_fit(fit, kind, p)
            self.assertLess(abs(q-p).max(), 2e-8)
            _, physical = model(fit['parameters'], kind, leakage=True)
            self.assertTrue(np.all(physical['leakage'] <= fit['leakage_cap']+1e-12))
            # Independent density/effect construction, not just model() re-evaluation.
            pauli = np.array([[[0, 1], [1, 0]], [[0, -1j], [1j, 0]], [[1, 0], [0, -1]]])
            rho, effect = [], []
            for c, r, lam in zip(physical['contrast'], physical['preparation_axes'], physical['leakage']):
                a = np.zeros((3, 3), complex)
                a[:2, :2] = (1-lam)*(np.eye(2)+c*np.einsum('k,kab->ab', r, pauli))/2
                a[2, 2] = lam
                self.assertGreaterEqual(np.linalg.eigvalsh(a).min(), -1e-12)
                self.assertAlmostEqual(np.trace(a).real, 1)
                rho.append(a)
            for lo, hi, m, z in zip(physical['readout_low'], physical['readout_high'], physical['measurement_axes'], physical['leak_response']):
                e = np.zeros((3, 3), complex)
                e[:2, :2] = ((lo+hi)*np.eye(2)+(hi-lo)*np.einsum('k,kab->ab', m, pauli))/2
                e[2, 2] = z
                self.assertGreaterEqual(np.linalg.eigvalsh(e).min(), -1e-12)
                self.assertLessEqual(np.linalg.eigvalsh(e).max(), 1+1e-12)
                effect.append(e)
            born = np.array([[np.trace(r@e).real for e in effect] for r in rho])
            np.testing.assert_allclose(born, p, atol=2e-8, rtol=0)
            corrupted = dict(fit, parameters=list(fit['parameters']))
            corrupted['parameters'][-1] = 1.5
            with self.assertRaises(ValueError):
                validate_fit(corrupted, kind, p)

    def test_parser_reversed_outcomes_and_reject_bad_counts(self):
        with tempfile.TemporaryDirectory() as td:
            path = Path(td)/'fixture.zip'
            def make(bad=False):
                stream = io.StringIO()
                writer = csv.writer(stream)
                writer.writerow(['', '1 0', '0 0', 'n', 'j'])
                for i in range(5):
                    for j in range(4):
                        writer.writerow([4*i+j, 3, 7.5 if bad else 7, i, j])
                with zipfile.ZipFile(path, 'w') as z:
                    z.writestr('fixture/wyniki_testy_0.csv', stream.getvalue())
            make()
            result = read_archive(path, 5, 1, 1, 10)
            self.assertEqual(result['counts_0_1_by_job'][0][0][0], [7, 3])
            make(True)
            with self.assertRaises(ValueError):
                read_archive(path, 5, 1, 1, 10)


if __name__ == '__main__':
    unittest.main()
