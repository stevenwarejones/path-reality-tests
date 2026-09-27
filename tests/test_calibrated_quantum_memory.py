"""Physical, inference and acquisition-failure tests; no network or source data."""
import importlib.util
import itertools
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

STUDY = Path(__file__).resolve().parents[1]/'studies/calibrated-quantum-memory'


def load(name):
    spec = importlib.util.spec_from_file_location('memory_'+name, STUDY/(name+'.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


A = load('audit')
C = load('certificate')


class CountAuditTests(unittest.TestCase):
    def rows(self):
        return {k: {'00': 40, '01': 10, '10': 20, '11': 30} for k in A.LABELS}

    def test_regrouping_is_lossless_involution(self):
        rng = np.random.default_rng(1)
        rows = {k: dict(zip(A.OUTCOMES, map(int, rng.multinomial(8000, [.2,.3,.1,.4]))))
                for k in sorted(A.LABELS)}
        archived = A.regroup(rows)
        self.assertEqual(A.regroup(archived), rows)
        self.assertEqual(sum(map(lambda c: sum(c.values()), archived.values())), 8000*324)
        self.assertGreater(len({sum(c.values()) for c in archived.values()}), 1)

    def test_schema_rejects_missing_outcome_and_row(self):
        rows = self.rows()
        del rows[next(iter(rows))]['11']
        with self.assertRaises(ValueError): A.regroup(rows)
        rows = self.rows()
        del rows[next(iter(rows))]
        with self.assertRaises(ValueError): A.regroup(rows)

    def test_schema_rejects_floats_negative_and_boolean_counts(self):
        for bad in (-1, .2, True):
            rows = self.rows()
            rows[next(iter(rows))]['00'] = bad
            with self.assertRaises(ValueError): A.regroup(rows)

    def test_past_diagnostic_does_not_reject_constant_marginal(self):
        result = A.past_bound(self.rows(), .05, 324)
        self.assertTrue(all(r['epsilon_lower'] == 0 for r in result))

    def test_analytic_distance_to_shared_marginal(self):
        rows = self.rows()
        rows['xp,x,xp,x'] = {'00': 900000, '01': 0, '10': 100000, '11': 0}
        rows['xp,x,xm,x'] = {'00': 100000, '01': 0, '10': 900000, '11': 0}
        bound = max(r['epsilon_lower'] for r in A.past_bound(rows,.05,324))
        self.assertGreater(bound, .397)
        self.assertLess(bound, .4)

    def test_hash_mismatch_fails_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory)/'test.json'; p.write_text('{}')
            with self.assertRaises(ValueError):
                A.verify_file(p, {'bytes':2, 'sha256':'0'*64})


class MemoryCertificateTests(unittest.TestCase):
    def test_kraus_complete_positivity_and_normalization(self):
        for lam in (-1/3,0.,.4,.9,1.):
            ks = C.depolarizing_kraus(lam)
            np.testing.assert_allclose(sum(k.conj().T@k for k in ks), C.I, atol=1e-14)
            j = C.choi(ks)
            self.assertGreaterEqual(np.linalg.eigvalsh(j).min(), -1e-14)
            np.testing.assert_allclose(np.trace(j.reshape(2,2,2,2),axis1=1,axis2=3), C.I,atol=1e-14)

    def test_physical_swap_and_bypass_are_equal_on_full_operator_basis(self):
        self.assertLess(C.physical_checks()['maximum_residual'],1e-14)

    def test_reset_probe_distinguishes_exact_pair(self):
        for lam in (.1,.6,1.):
            output=C.channel(C.STATES[0],C.depolarizing_kraus(lam))
            self.assertAlmostEqual(np.trace(C.STATES[0]@output).real,(1+lam)/2)
        self.assertAlmostEqual(C.calibration_probabilities(C.depolarizing_kraus(0))[0,0,0],.5)

    def test_single_outcome_bound_cannot_be_applied_to_hidden_flags(self):
        for rho in C.STATES:
            average = sum(p@rho@p.conj().T/4 for p in [C.I,*C.P])
            corrected = sum(p.conj().T@(p@rho@p.conj().T)@p/4 for p in [C.I,*C.P])
            np.testing.assert_allclose(average,C.I/2,atol=1e-14)
            np.testing.assert_allclose(corrected,rho,atol=1e-14)

    def test_continuum_of_classical_memory_models_obeys_eb_bound(self):
        # Arbitrary measurement direction and different conditional output direction.
        rng=np.random.default_rng(6)
        for _ in range(100):
            v=rng.normal(size=3);v/=np.linalg.norm(v)
            w=rng.normal(size=3);w/=np.linalg.norm(w)
            m=(C.I+sum(v[a]*C.P[a] for a in range(3)))/2
            out=(C.I+sum(w[a]*C.P[a] for a in range(3)))/2
            f=np.mean([np.trace(rho@(np.trace(m@rho)*out+
                       np.trace((C.I-m)@rho)*(C.I-out))).real for rho in C.STATES])
            self.assertAlmostEqual(f,.5+np.dot(v,w)/6)
            self.assertLessEqual(f,2/3+1e-14)

    def test_complete_source_removal_witnesses_and_joint_feasibility(self):
        data=json.loads((STUDY/'results/synthetic-certificate.json').read_text())
        for key,value in data['source_removal_witnesses'].items():
            if type(value) is bool: self.assertTrue(value,key)
        result=C.infer(data['calibration_successes'],data['dynamics_successes'])
        self.assertTrue(result['reject'])
        self.assertGreater(result['margin'],.04)
        self.assertAlmostEqual(result['margin'],data['certificate']['margin'])

    def test_unsupported_tables_fail_closed(self):
        for cal,dyn in [(np.zeros((3,3,1),dtype=int),np.zeros(6,dtype=int)),
                        (np.zeros((3,3,2)),np.zeros(6,dtype=int)),
                        (np.full((3,3,2),8001),np.zeros(6,dtype=int))]:
            with self.assertRaises(ValueError): C.infer(cal,dyn)

    def test_affine_inconsistent_calibration_is_not_memory_rejection(self):
        cal=np.full((3,3,2),4000,dtype=int)
        cal[0,0,:]=8000; cal[0,1,:]=0
        result=C.infer(cal,np.full(6,8000,dtype=int))
        self.assertEqual(result['status'],'calibration_affine_inconsistent')
        self.assertFalse(result['reject'])

    def test_measured_record_never_claims_memory_certificate(self):
        data=json.loads((STUDY/'results/nmn-audit.json').read_text())
        self.assertFalse(data['quantum_memory_certified'])
        self.assertEqual(data['simultaneous_rows'],3240)
        uq=[r for r in data['runs'] if r['source']=='NMN_lab_rslts.json'][0]
        self.assertAlmostEqual(uq['strongest_past_marginal_contrast']['epsilon_lower'],
                               .057844137838283474)


if __name__ == '__main__':
    unittest.main()
