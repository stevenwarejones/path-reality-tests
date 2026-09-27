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

    def test_target_grid_rounding_preserves_full_error_budget(self):
        # 0.1 * 1e8 rounds to an integer although exact float(0.1)*1e8 does not.
        for z in [.1,.3,.7,np.nextafter(.5,0),np.nextafter(.5,1)]:
            low,high=dv.outward_target(z);center=F.from_float(float(z))
            self.assertLessEqual(low,center-dv.ERROR)
            self.assertGreaterEqual(high,center+dv.ERROR)

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


class TransferMethodBarrier(unittest.TestCase):
    def test_exact_complete_GST_constraints(self):
        import transfer_barrier as tb
        result=tb.run()
        self.assertFalse(result['physical_model'])
        self.assertEqual(result,json.loads((D/'results/transfer-barrier.json').read_text()))
        self.assertEqual(result['checked_two_sided_GST_constraints'],4663)
        self.assertGreater(result['minimum_GST_slack']['decimal'],.0028)

    def test_infinite_word_proof_guards(self):
        import transfer_barrier as tb
        tb.algebra_checks(F(21,100),F(1279,2500),F(1163,2500))
        with self.assertRaisesRegex(ValueError,'range argument'):tb.algebra_checks(F(19,100),F(1279,2500),F(1163,2500))
        with self.assertRaisesRegex(ValueError,'exceptional'):tb.algebra_checks(F(21,100),F(99,100),F(1,100))

    def test_transfer_assignment_exceptional_words(self):
        import transfer_barrier as tb
        rng=np.random.default_rng(210927)
        known={(0,0,2):{'k':22,'n':50},(0,0,2,2,2):{'k':26,'n':50},():{'k':0,'n':50}}
        words=[tuple(rng.integers(0,3,size=int(rng.integers(0,12)))) for _ in range(100)]
        words.extend([dv.TARGET,dv.TARGET[:-1],dv.TARGET+(2,),(),(0,0)])
        for endpoint in [0,1]:
            for u in words:
                for v in [(),dv.TARGET,dv.TARGET[:-1],(0,0)]:
                    for suffix in [(),(2,),(2,2),(0,2,2),dv.TARGET]:
                        p=float(tb.scalar_probability(u,known,endpoint));q=float(tb.scalar_probability(v,known,endpoint))
                        diff=abs(float(tb.scalar_probability(u+suffix,known,endpoint)-tb.scalar_probability(v+suffix,known,endpoint)))
                        self.assertLessEqual(diff,np.sqrt(p*(1-q))+np.sqrt(q*(1-p))+1e-14)

    def test_gram_containment_on_physical_gate_set(self):
        from gram_probe import contraction_matrix
        from dynamics_common import operations
        gates,initial,effect,_=operations(baseline_models()['joint']['x'])
        rng=np.random.default_rng(927);states=rng.normal(size=(3,8));states/=2*np.maximum(1,np.linalg.norm(states,axis=0))
        for g in range(3):
            T=gates[0,g,1:4,1:4];t=gates[0,g,1:4,0];after=T@states+t[:,None]
            vectors=np.c_[states,after,t];Q=vectors.T@vectors;M=contraction_matrix(Q,[(j,j+8) for j in range(8)],16)
            np.testing.assert_allclose(M,states.T@(np.eye(3)-T.T@T)@states,atol=1e-15)
            self.assertGreaterEqual(np.linalg.eigvalsh(M).min(),-1e-14)

    def test_nonCP_map_passes_contraction_relaxation(self):
        # Transposition is a ball isometry but its Choi matrix is the swap, with eigenvalue -1.
        T=np.diag([1.,-1.,1.]);np.testing.assert_array_equal(T.T@T,np.eye(3))
        swap=np.array([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]],float)
        self.assertEqual(np.linalg.eigvalsh(swap).min(),-1)

if __name__=='__main__':unittest.main()
