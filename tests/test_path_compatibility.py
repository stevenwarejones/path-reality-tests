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
import wen_reconstruction as w
import measurements as m


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

    def test_real_images_and_reference_slice_axes(self):
        x,y,frames,xp,idx,target,se,hashes=w.load_example()
        self.assertEqual(len(hashes),6)
        self.assertTrue(np.array_equal(x[idx],xp))
        self.assertEqual(frames['R'].shape,(321,101))
        self.assertEqual(target.shape,(17,));self.assertEqual(se.shape,(17,2))

    def test_synthetic_calibrated_reconstruction(self):
        psi=np.array([2.,3.]);z=np.array([.2+.4j,-.3+.1j]);C=psi*z
        base=10*np.ones(2)
        signals={'Plus':base+C.real/2,'Minus':base-C.real/2,
                 'R':base+C.imag/2,'L':base-C.imag/2,'PSI':psi**2}
        gains={ch:1.1 for ch in w.CHANNELS};backgrounds={ch:2. for ch in w.CHANNELS}
        recorded={ch:signals[ch]*gains[ch]+backgrounds[ch] for ch in w.CHANNELS}
        np.testing.assert_allclose(w.reconstruct(recorded,gains,backgrounds),z,atol=1e-14)

    def test_coordinate_enclosures_include_shared_calibrations(self):
        _,y,frames,_,idx,_,_,_=w.load_example();signals=w.reduce_frames(frames,y,idx,'all_y')
        box=w.conditional_box(signals,.05,1.)
        rng=np.random.default_rng(901)
        for _ in range(100):
            gains={ch:rng.uniform(.95,1.05) for ch in w.CHANNELS}
            background={ch:rng.uniform(0,1) for ch in w.CHANNELS}
            z=w.reconstruct(signals,gains,background)
            self.assertTrue(np.all(z.real>=box[:,0]-1e-12) and np.all(z.real<=box[:,1]+1e-12))
            self.assertTrue(np.all(z.imag>=box[:,2]-1e-12) and np.all(z.imag<=box[:,3]+1e-12))

    def test_zero_reference_is_reported_not_clipped(self):
        v={ch:np.ones(2) for ch in w.CHANNELS};v['PSI'][0]=0
        with self.assertRaises(ValueError):w.reconstruct(v)
        with self.assertRaises(ValueError):w.conditional_box(v)

    def test_full_pointer_probabilities_are_phase_gauge_invariant(self):
        p=1.3+.2j;k=.4-.7j;u=np.exp(.61j)
        # Every polarization analyzer depends on the same two amplitudes.
        for a,c in [(1,1),(1,-1),(1,1j),(1,-1j)]:
            self.assertAlmostEqual(abs(a*p+c*k)**2,abs(a*u*p+c*u*k)**2)
        self.assertNotAlmostEqual((u*k).imag,k.imag)

    def test_measurement_budget_and_loss_scoring(self):
        r=m.analyze();best=r['candidates'][1]
        self.assertEqual(best['eligible_main_trials'],22)
        self.assertLessEqual(best['type_I_and_II_upper_bound'],.01)
        self.assertIsNone(best['total_eligible_trials'])
        self.assertEqual(r['unrestricted_completion_result']['worst_case_separation'],0)
        # Erasures contribute zero instead of conditioning them away.
        self.assertEqual(F(4,5)*F(9,10)*F(9,10),F(best['signed_score_margin']))

    def test_analytic_sine_bounds_hold_at_family_extrema(self):
        self.assertGreaterEqual(np.sin(7*np.pi/40),.5)
        self.assertGreaterEqual(np.sin(11*np.pi/40),.5)
        self.assertGreaterEqual(np.sin(3*np.pi/8),.9)
        self.assertGreaterEqual(np.sin(21*np.pi/40),.9)
        self.assertLess(m.main_trial_budget(F(81,125))['eligible_main_trials'],
                        m.main_trial_budget(F(9,25))['eligible_main_trials'])

if __name__=='__main__':unittest.main()
