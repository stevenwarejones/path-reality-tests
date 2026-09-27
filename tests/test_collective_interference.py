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
import cell_certificate as cells
import full_models as full
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

class ExpandedRegionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        read=lambda n:json.loads((STUDY/'results'/n).read_text())
        cls.counts=read('counts.json');cls.cert=read('cell-certificate.json')
        cls.w=read('physical-witness.json');cls.hist=read('parity-histogram.json')
        cls.region=read('full-region.json');cls.models={k:read('full-'+k+'.json') for k in ['boson','cluster']}

    def test_cell_certificate_and_statistical_tradeoff(self):
        r=cells.verify(self.cert,self.counts,self.w)
        self.assertEqual(r['cell_event_upper'],.0112)
        self.assertLess(r['row_only_margin_same_region'],0)
        self.assertGreater(r['margin'],.0012)

    def test_cell_dual_and_endpoint_corruption(self):
        for change in ['dual','endpoint']:
            c=copy.deepcopy(self.cert)
            if change=='dual':c['conditional_duals'][4]['upper']='0'
            else:c['conditional_intervals'][3][5]=['1','1']
            with self.assertRaises(AssertionError):cells.verify(c,self.counts,self.w)

    def test_conditional_collision_identity(self):
        a=[F(1,2),F(1,3),F(1,6)]+[F(0)]*9
        b=list(reversed(a));c=[F(1,12)]*12
        direct=sum(a[i]*b[j]*c[k] for i in range(12) for j in range(12) for k in range(12) if len({i,j,k})==3)
        formula=1-sum(a[x]*b[x]+a[x]*c[x]+b[x]*c[x] for x in range(12))+2*sum(a[x]*b[x]*c[x] for x in range(12))
        self.assertEqual(direct,formula)
        # Every true product satisfies all outer inequalities.
        v=a+b+c+[a[x]*b[x] for x in range(12)]+[a[x]*c[x] for x in range(12)]+[b[x]*c[x] for x in range(12)]+[a[x]*b[x]*c[x] for x in range(12)]
        rows,rhs,*_=cells.conditional_problem([[F(0)]*12 for _ in range(3)],[[F(1)]*12 for _ in range(3)])
        self.assertTrue(all(sum(co*v[j] for j,co in row.items())<=r for row,r in zip(rows,rhs)))

    def test_both_complete_region_models(self):
        for kind,m in self.models.items():
            r=full.verify(m,self.hist,self.region,self.counts,self.cert)
            self.assertTrue(r['full_parity_region_compatible'])
            self.assertGreater(r['margin'],.000625)
        # The higher-only model is not silently promoted into joint feasibility.
        with self.assertRaises(AssertionError):full.single_region(self.models['cluster'],self.cert)

    def test_nonphysical_and_fake_pattern_regions_rejected(self):
        m=copy.deepcopy(self.models['boson']);m['amplitude_denominator']=1
        with self.assertRaises(AssertionError):full.transport(m)
        r=copy.deepcopy(self.region);r['pattern_intervals']['1']=['0','1']
        # Zero lower is a valid outward endpoint: an invalid positive lower must fail.
        r['pattern_intervals']['1']=['1/2','1']
        with self.assertRaises(AssertionError):full.verify(self.models['boson'],self.hist,r,self.counts,self.cert)
        h=copy.deepcopy(self.hist);h['patterns'][0][1]+=1
        with self.assertRaises(AssertionError):full.histogram(h,self.counts)

    def test_loss_trace_against_explicit_isometry_enumeration(self):
        # 8 x 3 complex isometry: Fourier(4) columns, split by 3/5 and 4/5.
        # Two visible modes; six explicit loss modes. Independent block permanents.
        import math
        roots=[(1,0),(0,1),(-1,0),(0,-1)]
        U=[[(scale*roots[(r*j)%4][0],scale*roots[(r*j)%4][1]) for j in range(3)] for scale in [3,4] for r in range(4)]
        def amp(outs,cols):
            z=(0,0)
            for perm in itertools.permutations(cols):
                term=(1,0)
                for o,j in zip(outs,perm):term=full.mul(term,U[o][j])
                z=full.add(z,term)
            return z
        def pattern(outs):
            from collections import Counter
            return tuple(i for i,n in sorted(Counter(outs).items()) if i<2 and n%2)
        for kind in ['boson','cluster']:
            truth={p:F(0) for n in range(3) for p in itertools.combinations(range(2),n)}
            if kind=='boson':
                from collections import Counter
                for outs in itertools.combinations_with_replacement(range(8),3):
                    den=10**6*math.prod(math.factorial(n) for n in Counter(outs).values())
                    truth[pattern(outs)]+=F(full.norm(amp(outs,(0,1,2))),den)
            else:
                for pair in itertools.combinations_with_replacement(range(8),2):
                    prob=F(full.norm(amp(pair,(0,1))),10**4*(2 if pair[0]==pair[1] else 1))
                    for o in range(8):truth[pattern(pair+(o,))]+=prob*F(full.norm(U[o][2]),100)
            model=full.ParityModel(U[:2],10,kind)
            self.assertEqual(sum(truth.values()),1)
            for pat,p in truth.items():self.assertEqual(F(model.parity(pat),model.den),p)

if __name__=='__main__':unittest.main()
