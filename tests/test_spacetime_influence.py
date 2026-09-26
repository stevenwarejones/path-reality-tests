"""Independent small-sample and countermodel checks for the spacetime study."""
import importlib.util
import itertools
import math
from pathlib import Path
import unittest

import numpy as np
from scipy.stats import binom

STUDY = Path(__file__).resolve().parents[1]/'studies/spacetime-causal-influence'
spec = importlib.util.spec_from_file_location('spacetime_design',STUDY/'design.py')
d = importlib.util.module_from_spec(spec)
spec.loader.exec_module(d)


class SpacetimeTests(unittest.TestCase):
    def test_cp_inversion_and_coverage(self):
        for n in (1,2,7,15):
            for p in (0,0.01,0.2,0.5,0.9,1):
                coverage = 0
                for k in range(n+1):
                    lo,hi = d.cp_interval(k,n,0.05)
                    if lo-1e-13 <= p <= hi+1e-13:
                        coverage += binom.pmf(k,n,p)
                    if k:
                        self.assertAlmostEqual(binom.sf(k-1,n,lo),0.025,places=10)
                    if k<n:
                        self.assertAlmostEqual(binom.cdf(k,n,hi),0.025,places=10)
                self.assertGreaterEqual(coverage,0.95-1e-12)

    def test_unequal_sample_difference_coverage(self):
        for p0,p1 in itertools.product((0,0.15,0.5,0.9,1),repeat=2):
            cover, false_reject = 0,0
            for k0 in range(5):
                for k1 in range(8):
                    mass=binom.pmf(k0,4,p0)*binom.pmf(k1,7,p1)
                    ci=d.iid_difference(k0,4,k1,7,0.1)
                    cover += mass*(ci['signed_lower']-1e-12 <= p1-p0 <= ci['signed_upper']+1e-12)
                    false_reject += mass*d.classify(ci,abs(p1-p0))['reject']
            self.assertGreaterEqual(cover,0.9-1e-12)
            self.assertLessEqual(false_reject,0.1+1e-12)

    def test_empty_arm_is_vacuous(self):
        ci=d.iid_difference(0,0,3,3)
        self.assertEqual(ci['absolute_lower'],0)
        self.assertEqual(ci['absolute_upper'],1)

    def test_absolute_interval_and_strict_boundary(self):
        self.assertEqual(d.absolute_interval(-0.2,0.1),(0,0.2))
        self.assertEqual(d.absolute_interval(0.2,0.3),(0.2,0.3))
        self.assertFalse(d.classify({'absolute_lower':0.1,'absolute_upper':0.3},0.1)['reject'])
        self.assertEqual(d.classify({'absolute_lower':0.0,'absolute_upper':0.3},0.1)['clean_absolute_upper'],0.4)

    def test_score_power_against_direct_binomial_sum(self):
        for n in (10,40,120):
            for gap,budget in ((0,0),(0.3,0.05),(0.8,0.1)):
                expected=sum(binom.pmf(k,n,(1+gap)/2) for k in range(n+1)
                             if d.classify(d.score_interval(k,n,0.2),budget)['reject'])
                self.assertAlmostEqual(d.exact_score_power(n,gap,budget,0.2),expected,places=12)

    def test_power_guarantees_and_unidentifiable_margin(self):
        for gap,budget in ((0.1,0.002),(0.03,0.01),(0.003,0.002)):
            n=d.sufficient_trials(gap,budget)
            self.assertIsInstance(n,int)
            self.assertGreater(gap,budget+d.score_radius(n,0.01)+math.sqrt(2*math.log(10)/n))
            self.assertGreaterEqual(d.exact_score_power(n,gap,budget),0.9)
        self.assertIsNone(d.sufficient_trials(0.01,0.01))
        self.assertIsNone(d.sufficient_trials(0.0,0.002))

    def test_memory_null_exact_enumeration(self):
        # Arbitrary deterministic receiver memory: B_i is previous X XOR previous B.
        # Enumerate all fresh fair setting sequences; no binomial assumption about B.
        n,alpha=12,0.25
        failures=0
        for xs in itertools.product((0,1),repeat=n):
            last_x,last_b,matches=0,0,0
            for x in xs:
                b=last_x ^ last_b
                matches += x==b
                last_x,last_b=x,b
            failures += d.classify(d.score_interval(matches,n,alpha),0)['reject']
        self.assertLessEqual(failures/2**n,alpha)

    def test_randomization_score_identity_and_bias(self):
        for p0,p1,pi in itertools.product((0,0.2,0.7,1),repeat=3):
            mean=(1-pi)*(1-2*p0)+pi*(2*p1-1)
            self.assertAlmostEqual(mean,(p1-p0)+(2*pi-1)*(p1+p0-1))
            self.assertLessEqual(abs(mean-(p1-p0)),2*abs(pi-0.5)+1e-12)

    def test_common_cause_and_selection(self):
        # Independent fair R,Y before selection. Conditional selected R=Y.
        joint={(r,y):0.25 for r,y in itertools.product((0,1),repeat=2)}
        for y in (0,1):
            self.assertEqual(joint[1,y]/sum(joint[r,y] for r in (0,1)),0.5)
            selected=sum(joint[r,y] for r in (0,1) if r==y)
            self.assertEqual((joint[1,y] if y==1 else 0)/selected,y)
        # X=R observationally; intervening on X leaves R fair in this model.
        self.assertEqual(joint[1,0]+joint[1,1],0.5)

    def test_coupling_and_contamination_budgets(self):
        for q,r0,r1,e0,e1 in itertools.product((0,0.25,1),repeat=5):
            p0=(1-e0)*q+e0*r0
            p1=(1-e1)*q+e1*r1
            self.assertLessEqual(abs(p1-p0),max(e0,e1)+1e-12)
        self.assertAlmostEqual(d.coupling_budget([0.001,0.002],[0.003,0.004]),0.01)
        # Arbitrary couplings can attain e0+e1: a common q=1/2, opposite flips.
        self.assertAlmostEqual(abs(0.7-0.3),0.2+0.2)

    def test_geometry_uncertainty_and_vacuum_boundary(self):
        g=d.geometry(30,0.5,0.5,0.1,-5,20,35,65)
        self.assertTrue(g['spacelike'])
        self.assertAlmostEqual(g['margin_m'],28.9-d.C_M_PER_NS*70)
        self.assertFalse(d.geometry(3,0.5,0.5,0.1,-5,20,35,65)['spacelike'])
        self.assertFalse(d.geometry(d.C_M_PER_NS*70,0,0,0,-5,20,35,65)['spacelike'])
        # An earlier setting leak destroys an otherwise good geometry.
        self.assertFalse(d.geometry(30,0.5,0.5,0.1,-100,20,35,65)['spacelike'])

    def test_coarse_binary_bound_not_full_tv(self):
        p=np.array([0.5,0.5,0]); q=np.array([0.5,0,0.5])
        self.assertEqual(abs(p[0]-q[0]),0)
        self.assertEqual(np.abs(p-q).sum()/2,0.5)

    def test_after_the_fact_filter_can_create_preresponse(self):
        # Centered smoothing imports future causal data into a past time bin.
        x=np.array([0.,0.,0.,1.,1.,1.,1.])
        self.assertEqual(x[2],0)
        self.assertGreater(np.convolve(x,np.ones(3)/3,mode='same')[2],0)

    def test_snapshot_exact_discrete_and_float_tolerance(self):
        d.compare({'p':0.3,'n':100},{'p':0.3+1e-13,'n':100})
        with self.assertRaisesRegex(AssertionError,r'\$\.n'):
            d.compare({'p':0.3,'n':100},{'p':0.3,'n':101})
        with self.assertRaisesRegex(AssertionError,r'\$\.p'):
            d.compare({'p':0.3},{'p':0.301})

    def test_input_validation(self):
        for k,n in ((True,10),(1,2.5),(-1,4),(5,4)):
            with self.assertRaises(ValueError): d.score_interval(k,n)
        for a in (0,1,float('nan')):
            with self.assertRaises(ValueError): d.score_radius(10,a)
        with self.assertRaises(ValueError): d.coupling_budget([-0.01],[])
        with self.assertRaises(ValueError): d.geometry(30,0,0,0,1,0,0,1)


if __name__ == '__main__':
    unittest.main()
