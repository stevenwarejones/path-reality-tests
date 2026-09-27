"""Numerical diagnostic invariants, separate from the standard-library certificates."""
from pathlib import Path
import sys,unittest
try:
    import numpy as np
    import scipy
except ImportError:
    np=None

@unittest.skipIf(np is None,'Numerical diagnostic requires pinned scientific dependencies')
class DynamicsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        p=Path(__file__).resolve().parents[1]/'studies/collective-interference-identifiability'
        sys.path.insert(0,str(p))
        import dynamics, extended_dynamics, diagonal_dynamics, fresh_power
        cls.d=dynamics
        cls.e=extended_dynamics;cls.diag=diagonal_dynamics;cls.power=fresh_power
        sys.path.pop(0)

    def test_zero_time_and_no_tunnelling(self):
        q,D=self.d.prediction([0,0,0,0,.9],[0,2])
        for j in range(3):
            self.assertTrue(np.allclose(q[:,j,12*(4+j)+5],.9))
            self.assertTrue(np.allclose(q[:,j,-1],.1))
        self.assertTrue(np.allclose(D,0))

    def test_shared_fit_is_not_mislabelled_compatible(self):
        r=self.d.evaluate()
        self.assertEqual(r['singleton_shots'],4266)
        self.assertEqual(r['singleton_cell_interval_failure_count'],6)
        self.assertFalse(r['all_records_compatible'])
        self.assertLess(r['quadrature_61_vs_81_max_difference'],1e-10)


    def test_diagonal_zero_limit_matches_separable_channel(self):
        records=[dict(time_ms=t,input_y=j,y_start=-2,x_start=-2,width=5) for t in [.5,2.4] for j in [-1,0,1]]
        p=[120,110,5,-3,2,-4,1,2,-.1,.04,.96]
        q=self.e.prediction(p,records,order=9,half=6)
        q2=self.diag.prediction(p+[0],records,order=9,half=6)
        self.assertTrue(np.allclose(q,q2,atol=1e-12,rtol=1e-10))

    def test_fresh_test_size_and_monotonic_power(self):
        for q in [.001,.03,.2,.49]:
            self.assertLessEqual(self.power.power(80,100,2*q,q),2/320+1e-12)
        power=self.power.power
        self.assertGreater(power(200,200,.2,.03),power(200,200,.15,.03))
        self.assertGreater(power(200,200,.2,.03),power(200,200,.2,.05))
        self.assertGreater(power(200,200,.2,.03),power(200,200,.2,.03,2.1,.01))

    def test_label_gram_bounds_for_complex_product_labels(self):
        from itertools import permutations
        rng=np.random.default_rng(37)
        labels=np.eye(3,dtype=complex)+.03*(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))
        labels/=np.linalg.norm(labels,axis=1)[:,None]
        overlaps=labels.conj()@labels.T
        a=np.max(abs(overlaps-np.eye(3)))
        perms=list(permutations(range(3)))
        gram=np.array([[np.prod([overlaps[p[k],q[k]] for k in range(3)]) for q in perms] for p in perms])
        eig=np.linalg.eigvalsh(gram);r=3*a*a+2*a**3
        self.assertGreaterEqual(eig.min(),1-r-1e-12)
        self.assertLessEqual(eig.max(),1+r+1e-12)

    def test_july2_unique_records_and_fixed_point_cluster_envelope(self):
        import json
        records=json.loads((self.e.HERE/'results/july2-dynamics-records.json').read_text())
        self.assertEqual(len(records),4)
        self.assertEqual(sum(r['shots'] for r in records),1624)
        p=json.loads((self.e.HERE/'results/diagonal-dynamics-parameters.json').read_text())['July2']['parameters']
        for r in self.diag.july2_bunch(p):
            self.assertLessEqual(r['cluster_partition_upper_at_this_point'],2*r['distinguishable']+1e-12)
            self.assertLessEqual(r['boson'],6*r['distinguishable']+1e-12)
