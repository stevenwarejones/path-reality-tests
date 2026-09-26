"""Independent checks of boundary certificates, reconstruction and design scope."""
import importlib.util
from pathlib import Path
from fractions import Fraction as F
import sys
import unittest
import numpy as np
from scipy.optimize import linprog

ROOT=Path(__file__).resolve().parents[1]/'studies/path-compatibility'
sys.path.insert(0,str(ROOT))
import boundary as b


class CompatibilityTests(unittest.TestCase):
    def test_sharp_reference_thresholds(self):
        self.assertEqual(b.required_disturbance(F(16,25)),F(297,1225))
        self.assertEqual(b.required_disturbance(b.Q_MAX),F(1,50))
        self.assertEqual(b.required_disturbance(F(1)),F(1,50))
        self.assertIsNone(b.required_disturbance(b.G-F(1,100000)))
        self.assertEqual(b.CROSS,F(22223,25823))

    def test_exact_attainers_cover_branches_and_larger_disturbance(self):
        for q in [b.G,b.CROSS,b.Q_MAX,F(1),F(16,25)]+[F(i,1000) for i in range(540,1001)]:
            b.verify_model(b.attaining_model(q))
            b.verify_model(b.attaining_model(q,F(1)))

    def test_independent_kernel_linear_program(self):
        # 8 raw transition masses and d. No hard-coded attainer in this LP.
        def idx(l,ptr,j): return l*4+ptr*2+j
        mu=[float(b.F0),float(1-b.F0)]
        for q in [float(b.G),.64,float(b.CROSS),.875,.95,1.]:
            eq=[];rhs=[];ub=[];ur=[]
            for l in range(2):
                row=np.zeros(9);row[l*4:l*4+4]=1;eq.append(row);rhs.append(1)
                row=np.zeros(9);row[l*4:l*4+2]=1;ub.append(row);ur.append(q)
                row=np.zeros(9);row[idx(l,0,l)]=-1;row[idx(l,1,l)]=-1;row[8]=-1
                ub.append(row);ur.append(-1)
            for ptr in range(2):
                for j in range(2):
                    row=np.zeros(9)
                    for l in range(2): row[idx(l,ptr,j)]=mu[l]
                    eq.append(row);rhs.append(float(b.TARGET[ptr][j]))
            objective=np.zeros(9);objective[8]=1
            result=linprog(objective,A_ub=ub,b_ub=ur,A_eq=eq,b_eq=rhs,bounds=[(0,1)]*9,method='highs')
            self.assertTrue(result.success,result.message)
            expected=b.required_disturbance(F(str(q)))
            # Decimal rounding can place the lowest cap infinitesimally below G.
            if expected is None: expected=b.required_disturbance(b.G)
            self.assertAlmostEqual(result.fun,float(expected),places=9)

    def test_counterfeit_attainer_rejected(self):
        model=b.attaining_model(F(16,25));model['d']=F(13,320)
        with self.assertRaises(AssertionError): b.verify_model(model)

    def test_conditional_floor_and_reversed_sign(self):
        self.assertEqual(b.positive_floor(F(16,25),F(1,50),(b.F0,b.F0)),F(833,31250))
        self.assertEqual(b.positive_floor(1,1,(0,1)),0)
        self.assertGreater(b.positive_floor(F(16,25),F(1,50),(b.F0,b.F0)),b.B)
        self.assertEqual(b.positive_floor(F(16,25),F(1,50),(b.F0,b.F0),F(1,10)),0)


if __name__=='__main__':unittest.main()
