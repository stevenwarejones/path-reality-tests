"""Adversarial checks of dimension, composition, physicality and artifact binding."""
import importlib.util,json,sys,unittest
from pathlib import Path
from fractions import Fraction as F
import numpy as np
D=Path(__file__).resolve().parents[1]/'studies/unmeasured-dynamics-gain'
sys.path.insert(0,str(D))
import dynamics_verify as dv
import resource_certificate as rc
from dynamics_common import baseline_models,BASE,region,GeneralEvaluator

class UnmeasuredDynamics(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rows=json.loads((BASE/'results/observations.json').read_text());cls.rows=cls.rows['rows'] if isinstance(cls.rows,dict) else cls.rows
        cls.R=region(cls.rows);cls.models=baseline_models()

    def test_all_certificates_offline(self):
        self.assertEqual(dv.run(),json.loads((D/'results/certificate.json').read_text()))

    def test_corrupted_retained_count_rejected(self):
        m=dict(self.models['joint']);m['probabilities']=m['probabilities'].copy();m['probabilities'][10]=0
        with self.assertRaisesRegex(ValueError,'retained constraint'):dv.audit(self.rows,m,self.R,['GST','RB'])

    def test_dropped_source_is_not_silently_retained(self):
        with self.assertRaisesRegex(ValueError,'retained constraint'):dv.audit(self.rows,self.models['RB_only'],self.R,['GST','RB'])

    def test_saved_probabilities_are_not_a_model(self):
        rows=[{'word':(0,0,2,2),'expression':'GiGiGyGy','length':4}]
        m={'x':self.models['joint']['x'],'probabilities':[.25]}
        with self.assertRaisesRegex(ValueError,'binding'):dv.binding(rows,m)

    def test_invalid_spam_rejected(self):
        m=dict(self.models['joint']);m['x']=m['x'].copy();m['x'][48]=1.2
        with self.assertRaises(ValueError):dv.audit(self.rows,m,self.R,['GST','RB'])

    def test_corrupted_external_source_rejected(self):
        import tempfile
        from source_audit import audit_extra
        with tempfile.TemporaryDirectory() as directory:
            (Path(directory)/'RB-clifford-counts.txt').write_text('corrupted')
            with self.assertRaisesRegex(ValueError,'checksum'):audit_extra(directory,[])

    def test_root_rounding_and_wrong_domain(self):
        for n in [0,1,2,3,4,10**20-1,10**20+1]:
            q=dv.ceilroot(n);self.assertGreaterEqual(q*q,n)
            if q:self.assertLess((q-1)**2,n)
        with self.assertRaisesRegex(ValueError,'domain'):dv.lower_bound(dv.GRID,dv.GRID,0)

    def test_qutrit_counterexample(self):
        V=np.roll(np.eye(3),1,axis=0);rho=np.diag([1,0,0]);E=np.diag([0,0,1]);g=[V@V,V]
        def p(w):
            s=rho.copy()
            for k in w:s=g[k]@s@g[k].T
            return np.trace(E@s)
        self.assertEqual([p(w) for w in [(),(0,0),(1,1),(0,0,1,1)]],[0,0,1,0])
        self.assertEqual(dv.lower_bound(0,0,dv.GRID),1)

    def test_trace_distance_inequality_random_qubits(self):
        rng=np.random.default_rng(927)
        for _ in range(300):
            r=rng.normal(size=3);r*=rng.random()/np.linalg.norm(r);s=rng.normal(size=3);s*=rng.random()/np.linalg.norm(s)
            a=(1+r[2])/2;b=(1+s[2])/2
            self.assertLessEqual(np.linalg.norm(r-s)/2,np.sqrt(a*(1-b))+np.sqrt(b*(1-a))+1e-14)

    def test_noncommuting_order(self):
        from sqd_general import unitary
        from sqd_models import PAULI
        from search import pack_channels
        # z then x versus x then z, with a y-sensitive readout rotated into xz plane.
        x=pack_channels([[unitary(np.array([0,0,np.pi/2]))],[unitary(np.array([0,np.pi/2,0]))],[np.eye(2)]],[1,0,1,np.pi/2])
        rows=[{'word':(0,1),'expression':'GiGx'},{'word':(1,0),'expression':'GxGi'}]
        a=GeneralEvaluator(rows)(x);b=dv.prior.independent_probabilities(rows,x)
        np.testing.assert_allclose(a,b,atol=1e-12);self.assertGreater(abs(a[0]-a[1]),.4)

    def test_resource_bad_rayleigh_and_length_rejected(self):
        with self.assertRaisesRegex(ValueError,'negative partial transpose'):rc.verify_resource([[1,0],[0,0],[0,0],[0,0]])
        with self.assertRaisesRegex(ValueError,'length'):rc.power(rc.eye(4),20001)
        with self.assertRaisesRegex(ValueError,'positive definiteness'):rc.positive_definite([[rc.iv.mpc(-1)]])

    def test_exact_similarity_word_cancellation(self):
        # Rational noncommuting matrices: algebraic equality, not an eigenvalue fit.
        A=[[F(1),F(0)],[F(1,7),F(2,3)]];B=[[F(1),F(0)],[F(-1,5),F(1,2)]]
        def mm(a,b):return [[sum(x*y for x,y in zip(row,col)) for col in zip(*b)] for row in a]
        S=[[F(1),F(0)],[F(0),F(1009,1000)]];T=[[F(1),F(0)],[F(0),F(1000,1009)]]
        self.assertEqual(mm(mm(mm(S,A),T),mm(mm(S,B),T)),mm(mm(S,mm(A,B)),T))

if __name__=='__main__':unittest.main()
