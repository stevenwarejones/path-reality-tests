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
        import dynamics
        cls.d=dynamics
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
