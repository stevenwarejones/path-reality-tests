"""Independent physical counterexamples and certificate corruption tests."""
import copy
from fractions import Fraction as F
import importlib.util
import itertools
import json
from pathlib import Path
import sys
import unittest

STUDY=Path(__file__).resolve().parents[1]/'studies/collective-interference-identifiability'
sys.path.insert(0,str(STUDY))
import certificate as collective_certificate
import witness as collective_witness
sys.path.pop(0)

class CollectiveTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        read=lambda n:json.loads((STUDY/'results'/n).read_text())
        cls.counts=read('counts.json');cls.cert=read('certificate.json');cls.w=read('physical-witness.json')

    def test_exact_certificate(self):
        r=collective_certificate.verify(self.cert,self.counts)
        self.assertEqual(r['leaf_count'],4)
        self.assertGreater(r['separation_margin'],.0003)

    def test_physical_witness(self):
        r=collective_witness.verify(self.w,self.counts,self.cert)
        self.assertGreater(r['contraction_slack_lower'],.04)
        self.assertTrue(r['joint_summary_feasible'])
        self.assertFalse(r['complete_many_body_record_feasibility_established'])

    def test_fake_extra_calibration_shots_rejected(self):
        c=copy.deepcopy(self.counts);c['single_shots']*=3
        with self.assertRaises(AssertionError):collective_certificate.verify(self.cert,c)

    def test_invalid_confidence_endpoint_rejected(self):
        c=copy.deepcopy(self.cert);c['bunch_lower']='1/2'
        with self.assertRaises(AssertionError):collective_certificate.verify(c,self.counts)

    def test_optimistic_upper_bound_rejected(self):
        c=copy.deepcopy(self.cert);c['polynomial_upper']='1/100'
        with self.assertRaises(AssertionError):collective_certificate.verify(c,self.counts)

    def test_no_missing_branch(self):
        c=copy.deepcopy(self.cert);del c['tree']['left']
        with self.assertRaises(KeyError):collective_certificate.verify(c,self.counts)

    def test_dual_corruption_rejected(self):
        c=copy.deepcopy(self.cert)
        def corrupt(n):
            if 'split' in n:corrupt(n['left']);corrupt(n['right'])
            else:n['lambda']={};n['mu']='0'
        corrupt(c['tree'])
        with self.assertRaises(AssertionError):collective_certificate.verify(c,self.counts)

    def test_exact_binomial_and_sqrt_arithmetic(self):
        import math
        for n in range(1,10):
            for k in range(n+1):
                p=F(2,7);num,den=collective_certificate.cdf_numerator(n,k,p)
                truth=sum(F(math.comb(n,j))*p**j*(1-p)**(n-j) for j in range(k+1))
                self.assertEqual(F(num,den),truth)
        for q in [F(0),F(1,931),F(3,7),F(100)]:
            l,u=collective_witness.sqrt_interval(q);self.assertLessEqual(l*l,q);self.assertGreaterEqual(u*u,q)

    def test_cluster_factor_attained_by_physical_contraction(self):
        # Three detected modes, T_ij=1/4. Operator norm=3/4<1, so loss dilation exists.
        # All six assignment amplitudes equal 1/64. H={e,(01)} has two members.
        perms=list(itertools.permutations(range(3)));a=F(1,64)
        def same_coset(s,t):
            return s==t or tuple(1-v if v<2 else v for v in s)==t
        p2=sum(a*a for s in perms for t in perms if same_coset(s,t))
        pd=6*a*a;pb=36*a*a
        self.assertEqual(p2,2*pd);self.assertEqual(pb,6*pd)
        self.assertLess(F(3,4),1)

    def test_common_channel_drift_breaks_product_of_averages(self):
        # For each hidden drift setting, three columns are distinct basis vectors
        # within its chosen row. Each map is an isometry; labels are orthogonal.
        channels=[[[int(o==3*r+j) for j in range(3)] for o in range(6)] for r in range(2)]
        for T in channels:
            for j in range(3):
                for k in range(3):
                    self.assertEqual(sum(T[o][j]*T[o][k] for o in range(6)),int(j==k))
        marginals=[[sum(F(T[3*r+x][j],2) for T in channels for x in range(3)) for r in range(2)] for j in range(3)]
        d=sum(marginals[0][r]*marginals[1][r]*marginals[2][r] for r in range(2))
        self.assertEqual(d,F(1,4));self.assertGreater(F(1),2*d)

if __name__=='__main__':unittest.main()
