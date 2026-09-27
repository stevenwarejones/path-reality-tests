"""Independent physics, assumption-violation, inference and source-parser tests."""
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np
from scipy.stats import binom
from scipy.spatial.transform import Rotation

STUDY=Path(__file__).resolve().parents[1]/'studies/shared-quantum-dynamics'
sys.path.insert(0,str(STUDY))
from sqd_sources import parse_word,read_counts,external_cache,split
from sqd_models import (parameters,operations,kraus,certificate,Evaluator,
                        density_probability,reference_probability,dormant_flag_probability,IDEAL)
from sqd_equivalence import (fiber_interval,fiber_point,invariants,formula_probability,
                            leakage_population,monitored_loss,clock_probability,synthetic_certificate)
from sqd_inference import deviance,simultaneous_intervals


class SharedDynamics(unittest.TestCase):
    def test_parser_order_repetition_and_rejection(self):
        self.assertEqual(parse_word('Gx(GyGi)^3Gx'),(1,2,0,2,0,2,0,1))
        self.assertEqual(parse_word('{}'),())
        for bad in ['Gz','Gx;print(1)','(Gx)^0','(Gx)^20001','((Gx)^2)^2','Gx garbage','']:
            with self.assertRaises(ValueError):parse_word(bad)

    def test_denominator_and_schema_fail_closed(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'counts.txt';header='## Columns = plus count, count total\n'
            p.write_text(header+'Gx 3.0 285 0 0\n')
            self.assertEqual(read_counts(p,'RB')[0]['n'],285)
            for row in ['Gx 286 285 0 0','Gx -1 50 0 0','Gx 1.2 50 0 0','Gx 1 0 0 0','Gx nan 50 0 0','Gx 1 50 1 0']:
                p.write_text(header+row+'\n')
                with self.assertRaises(ValueError):read_counts(p,'RB')
            p.write_text('plus total\nGx 1 50 0 0\n')
            with self.assertRaises(ValueError):read_counts(p,'RB')
        with self.assertRaises(ValueError):external_cache(STUDY/'cache')

    def test_all_models_native_against_complex_kraus(self):
        rng=np.random.default_rng(771);words=[(),(1,),(1,1),(1,2,0)]+[tuple(rng.integers(0,3,25)) for _ in range(5)]
        E=Evaluator([{'word':w} for w in words])
        for model in ['qubit','leakage','alternating','quasistatic']:
            x,_,_=parameters(model);x[:9]=rng.normal(0,.4,9)
            if model=='leakage':x[15:]=[2.,3.,.8]
            p=E(x,model)
            np.testing.assert_allclose(p,[density_probability(w,x,model) for w in words],atol=1e-11,rtol=0)

    def test_cptp_entangled_input_and_povm(self):
        for model in ['qubit','leakage','alternating','quasistatic']:
            x,_,_=parameters(model)
            if model=='leakage':x[15:]=[4,7,.9]
            for b in range(2):
                for g in range(3):
                    Ks=kraus(x,model,b,g);d=Ks[0].shape[0]
                    bell=np.eye(d).reshape(-1)/np.sqrt(d);rho=np.outer(bell,bell)
                    result=sum(np.kron(K,np.eye(d))@rho@np.kron(K.conj().T,np.eye(d)) for K in Ks)
                    self.assertAlmostEqual(np.trace(result).real,1,places=12)
                    self.assertGreaterEqual(np.linalg.eigvalsh(result).min(),-1e-12)
            c=certificate(x,model);self.assertLess(c['trace_preservation_max_error'],1e-12)
            self.assertGreaterEqual(c['choi_min_eigenvalue'],-1e-12)
        x,_,_=parameters('qubit');x[10]=-1
        with self.assertRaises(ValueError):operations(x,'qubit')

    def test_nominal_outcome_and_noncommuting_order(self):
        x,_,_=parameters('qubit');x[9:12]=0;x[12:15]=[0,1,1]
        self.assertAlmostEqual(density_probability((),x,'qubit'),0)
        self.assertAlmostEqual(density_probability((1,1),x,'qubit'),1)
        a=density_probability((1,2,1,2),x,'qubit');b=density_probability((1,1,2,2),x,'qubit')
        self.assertGreater(abs(a-b),.4)

    def test_exact_fiber_formula_and_two_state_flag(self):
        rng=np.random.default_rng(710);x,_,_=parameters('leakage');x[:9]=rng.normal(0,.3,9);x[9:12]=[4,5,6];x[15:]=[1,2,.85]
        lo,hi=fiber_interval(x)
        for lam in np.linspace(lo,hi,7):
            y=fiber_point(x,lam)
            for w in [(),(1,),(1,2,0),(2,0,1,1)*12]:
                a=density_probability(w,x,'leakage');b=density_probability(w,y,'leakage')
                self.assertAlmostEqual(a,b,places=11)
                self.assertAlmostEqual(a,formula_probability(w,y),places=11)
                self.assertAlmostEqual(a,dormant_flag_probability(w,y),places=11)
        for lam in [lo-1e-6,hi+1e-6]:
            with self.assertRaises(ValueError):fiber_point(x,lam)

    def test_sharp_fiber_bound_violations(self):
        x,_,_=parameters('leakage');x[9:12]=[4,5,6];x[15:]=[1,2,.85];v=invariants(x);lo,hi=fiber_interval(x,False)
        # Just below lower bound requires leakage response above one.
        self.assertGreater(v['a']-v['B']*v['s']/(lo*.9),1)
        # Upper endpoint is return=0 for this example; above requires negative return.
        self.assertLess(v['s']-hi*1.1,0)
        y=fiber_point(x,lo);self.assertAlmostEqual(y[17],1)

    def test_population_monitor_breaks_equivalence(self):
        c=synthetic_certificate();self.assertLess(c['terminal_probability_max_error'],1e-10)
        self.assertGreater(abs(c['one_gate_monitor_a']-c['one_gate_monitor_b']),1e-4)
        for name in ['leakage_parameters_a','leakage_parameters_b']:
            x=c[name];s=invariants(x)['s']
            for N in [1,7,1000]:self.assertAlmostEqual(monitored_loss(leakage_population(N,x),N,s),x[15]*1e-4,places=13)

    def test_clock_and_memory_separated_only_by_declared_wait(self):
        x,_,_=parameters('alternating');x[15]=3
        for w in [(1,1),(2,0,1),(1,2,0)*4]:
            self.assertAlmostEqual(reference_probability(w,x,'alternating'),clock_probability(w,x),places=13)
        self.assertGreater(abs(clock_probability((1,1),x)-clock_probability((1,1),x,0)),.0008)

    def test_coherent_leakage_violates_flag_assumption(self):
        # Identical block diagonals, opposite computational/leakage coherence.
        plus=np.array([1,0,1],complex)/np.sqrt(2);minus=np.array([1,0,-1],complex)/np.sqrt(2)
        H=np.array([[1,0,1],[0,np.sqrt(2),0],[1,0,-1]],complex)/np.sqrt(2)
        rho=np.outer(plus,plus.conj());sig=np.outer(minus,minus.conj())
        np.testing.assert_allclose(np.diag(rho),np.diag(sig))
        self.assertGreater(abs((H@rho@H.conj().T)[0,0]-(H@sig@H.conj().T)[0,0]),.99)
        # Thus discarding coherences is essential, not a general qutrit theorem.

    def test_gate_dependent_loss_breaks_uniform_closed_form(self):
        x,_,_=parameters('leakage');x[15:]=[2,3,.9];M,r,E,_=operations(x,'leakage');word=(1,2)*20
        modified=M.copy();new=x.copy();new[15]=8;N,_,_,_=operations(new,'leakage');modified[:,2]=N[:,2]
        state=r.copy()
        for g in word:state=modified[0,g]@state
        self.assertGreater(abs(E@state-formula_probability(word,x)),1e-4)

    def test_observable_gauge_not_coordinate_difference(self):
        x,_,_=parameters('qubit');x[:9]=np.arange(9)*.03
        rz=Rotation.from_rotvec([0,0,.001]).as_matrix();y=x.copy()
        y[:9]=((IDEAL+x[:9].reshape(3,3)*.001)@rz.T-IDEAL).reshape(-1)/.001
        for w in [(1,2,0,2,1),(1,1,2)*5]:
            self.assertAlmostEqual(density_probability(w,x,'qubit'),density_probability(w,y,'qubit'),places=12)
        self.assertGreater(np.max(abs(x-y)),.1)

    def test_exact_interval_coverage_and_dependence_violation(self):
        # Directly enumerate a binomial law rather than mirror the beta formula.
        for p in [.05,.3,.8]:
            coverage=0.
            for k in range(21):
                lo,hi=simultaneous_intervals(np.array([k]),np.array([20]))
                if lo[0]<=p<=hi[0]:coverage+=binom.pmf(k,20,p)
            self.assertGreaterEqual(coverage,.95-1e-12)
        # Perfectly dependent Bernoulli shots all equal Z~Bernoulli(.5).
        for k in [0,10000]:
            lo,hi=simultaneous_intervals(np.array([k]),np.array([10000]))
            self.assertFalse(lo[0]<=.5<=hi[0])

    def test_saturation_and_frozen_split(self):
        self.assertLess(np.max(abs(deviance(np.array([0,3,10]),np.array([10,10,10]),np.array([0,.3,1])))),1e-10)
        rows=[{'family':'GST','length':1024},{'family':'GST','length':1025},{'family':'RB','length':1000},{'family':'RB','length':1001}]
        self.assertEqual(split(rows).tolist(),[True,False,True,False])

    def test_offline_results_reproduce(self):
        # Explicit import avoids collisions with other study modules called report.
        spec=importlib.util.spec_from_file_location('sqd_report',STUDY/'report.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
        text=m.make_report();self.assertEqual(text,(STUDY/'results/report.md').read_text())


class GeneralDynamicsEscalation(unittest.TestCase):
    def test_general_channels_include_nonunital_and_independent_propagation(self):
        from sqd_general import pack,channel_kraus,GeneralEvaluator,density_probability as dp,certificate as cert
        from sqd_models import PAULI
        from sqd_transfer import amplitude_damping
        rng=np.random.default_rng(1931)
        # Construct process factors from independent amplitude-damping Kraus maps.
        # Identity nominal gate isolates the nonunital part.
        kk=amplitude_damping(.31)
        coeff=np.array([[np.trace(P@K)/2 for P in PAULI] for K in kk]).T
        L=np.linalg.cholesky(coeff@coeff.conj().T+1e-15*np.eye(4))
        reconstructed=channel_kraus(pack(L),0)
        rho=np.diag([.2,.8])
        actual=sum(K@rho@K.conj().T for K in reconstructed)
        expected=sum(K@rho@K.conj().T for K in kk)
        np.testing.assert_allclose(actual,expected,atol=2e-14)
        x=np.r_[np.tile(pack(L),3),.83,.12,.94,.73]
        words=[(0,0,0),(),(0,),(0,0),(0,1),(0,0),(1,2,0,0,1),tuple(rng.integers(0,3,80))]
        pp=GeneralEvaluator([{'word':w} for w in words])(x)
        np.testing.assert_allclose(pp,[dp(w,x) for w in words],atol=1e-12)
        self.assertLess(cert(x)['trace_preservation_max_error'],1e-12)
        self.assertGreaterEqual(cert(x)['choi_min_eigenvalue'],-1e-12)
        with self.assertRaises(ValueError):channel_kraus(np.zeros(16),0)

    def test_repeat_bound_random_arbitrary_cptp_and_spam(self):
        from sqd_repeat import repeat_upper
        rng=np.random.default_rng(312)
        for _ in range(300):
            A=rng.normal(size=(4,2,2))+1j*rng.normal(size=(4,2,2))
            vals,V=np.linalg.eigh(sum(K.conj().T@K for K in A));inv=(V*vals**-.5)@V.conj().T
            kk=[K@inv for K in A]
            B=rng.normal(size=(2,2))+1j*rng.normal(size=(2,2));rho=B@B.conj().T;rho/=np.trace(rho)
            Q,_=np.linalg.qr(rng.normal(size=(2,2))+1j*rng.normal(size=(2,2)))
            E=(Q*rng.uniform(size=2))@Q.conj().T;p=[]
            for j in range(3):
                p.append(float(np.trace(E@rho).real));rho=sum(K@rho@K.conj().T for K in kk)
            self.assertLessEqual(p[2],repeat_upper(p[0],p[1])+1e-12)
        # Near the informative boundary: pure qubit coherent rotations.
        for angle in np.linspace(0,.4,30):
            self.assertLessEqual(np.sin(2*angle)**2,repeat_upper(0,np.sin(angle)**2)+1e-12)

    def test_repeat_bound_assumption_violations(self):
        from sqd_repeat import repeat_upper
        # Three-level cyclic permutation, binary effect |2><2|.
        U=np.roll(np.eye(3),1,axis=0);rho=np.diag([1.,0,0]);p=[]
        for _ in range(3):p.append(float(rho[2,2]));rho=U@rho@U.T
        self.assertEqual(p,[0.,0.,1.]);self.assertGreater(p[2],repeat_upper(p[0],p[1]))
        # A two-state clock: first call changes flag; second applies X.
        flag=0;bit=0;p=[0.]
        for _ in range(2):
            bit^=flag;flag^=1;p.append(float(bit))
        self.assertEqual(p,[0.,0.,1.])
        with self.assertRaises(ValueError):repeat_upper(-.1,.2)

    def test_exact_outward_binomial_and_root_certificates(self):
        import math
        from sqd_repeat import exact_upper,rational_repeat,GRID,repeat_upper
        for n in [12,50]:
            for k in range(n+1):
                a=exact_upper(k,n,120)
                if k<n:
                    num=sum(math.comb(n,j)*a**j*(GRID-a)**(n-j) for j in range(k+1))
                    self.assertLessEqual(120*num,GRID**n)
                self.assertGreaterEqual(rational_repeat(a,a)/GRID+1e-15,repeat_upper(a/GRID,a/GRID))
        # Exact binomial coverage, including boundaries, with refitted nuisance absent.
        for p in [.001,.03,.2,.5,.9]:
            miss=sum(binom.pmf(k,12,p) for k in range(13) if exact_upper(k,12,120)/GRID<p)
            self.assertLessEqual(miss,1/120+1e-12)

    def test_repeat_certificate_loses_coverage_under_unrestricted_shot_dependence(self):
        from sqd_repeat import exact_upper,rational_repeat,GRID
        # Ordinary qubit with E=I/2 gives p0=p1=p2=1/2. One shared fair
        # latent bit copied to every shot yields all-zero rows with probability
        # 1/2, on which an incorrectly assumed binomial certificate excludes p2.
        bound=rational_repeat(exact_upper(0,50,120),exact_upper(0,285,120))/GRID
        self.assertLess(bound,.5)

    def test_nonunital_gate_dependent_similarity_and_return_constraint(self):
        from transfer_audit import setup
        from sqd_transfer import similarity,apply,transform,amplitude_damping
        old,new,locked,state,E,E2,kk,kp,tp,lower=setup();T=similarity(.6)
        np.testing.assert_allclose(new,[T@G@np.linalg.inv(T) for G in old],atol=1e-14)
        self.assertGreater(np.max(abs(new-locked)),1e-6)
        rng=np.random.default_rng(615)
        for _ in range(15):
            word=rng.integers(0,3,25);a=state.copy();b=state.copy();rho=np.diag([1.,0,0]).astype(complex);sigma=rho.copy()
            for g in word:a=old[g]@a;b=new[g]@b;rho=apply(kk[g],rho);sigma=apply(kp[g],sigma)
            self.assertAlmostEqual(E@a,E2@b,places=12)
            self.assertAlmostEqual(np.trace(np.diag([.02,.98,.4])@rho).real,E@a,places=12)
            self.assertAlmostEqual(np.trace(np.diag([.02,.98,1/3])@sigma).real,E2@b,places=12)
            self.assertAlmostEqual(sigma[2,2].real,.6*rho[2,2].real,places=12)
        # Generic polarized return also has a physical neighborhood of kappa=1.
        from sqd_transfer import matrix,TAU
        original_return=np.diag([.6,.4]);ks=amplitude_damping(.02)
        ks2,l2,e2,t2=transform(ks,.02,.1,.95,original_return)
        tt=similarity(.95)
        np.testing.assert_allclose(matrix(ks2,l2,e2,t2),tt@matrix(ks,.02,.1,original_return)@np.linalg.inv(tt),atol=1e-14)
        # Large nonunital translation makes the formal similarity nonphysical.
        with self.assertRaises(ValueError):transform(amplitude_damping(.8),.001,.001,.5)
        with self.assertRaises(ValueError):transform(amplitude_damping(.1),.001,.001,1.1)

    def test_nonunital_two_state_memory_realizes_leakage(self):
        from sqd_transfer import amplitude_damping,qutrit_kraus,memory_probability,apply
        from sqd_general import unitary
        from sqd_models import IDEAL
        channels=[[K@unitary(IDEAL[g]) for K in amplitude_damping(.03+.01*g)] for g in range(3)]
        loss=[.02,.04,.01];ret=[.1,.2,.3];returns=[np.diag([.6,.4]),np.diag([.3,.7]),np.diag([.5,.5])]
        qq=[qutrit_kraus(channels[g],loss[g],ret[g],returns[g]) for g in range(3)]
        rho=np.diag([.7,.3]);E=np.array([[.2,.1],[.1,.8]]);z=.3
        word=(1,2,0,1,1,2)*7;state=np.zeros((3,3),complex);state[:2,:2]=rho
        effect=np.zeros((3,3),complex);effect[:2,:2]=E;effect[2,2]=z
        for g in word:state=apply(qq[g],state)
        self.assertAlmostEqual(memory_probability(word,channels,loss,ret,returns,rho,E,z),np.trace(effect@state).real,places=12)

    def test_exact_leakage_lift_and_sharp_cp_boundary(self):
        from sqd_transfer import lift_qubit,amplitude_damping,reset_kraus,TAU,apply,qutrit_kraus
        ks=[np.sqrt(.9)*K for K in amplitude_damping(.1)]+[np.sqrt(.1)*K for K in reset_kraus(TAU)]
        J=sum(np.outer(K.reshape(-1,order='F'),K.reshape(-1,order='F').conj()) for K in ks);budget=2*np.linalg.eigvalsh(J).min()
        base,tau=lift_qubit(ks,budget*.8,.2);qq=qutrit_kraus(base,budget*.8,.2,tau)
        rho=np.array([[.7,.1j],[-.1j,.3]]);state=np.zeros((3,3),complex);state[:2,:2]=rho
        E=np.array([[.2,.1],[.1,.8]]);effect=np.zeros((3,3),complex);effect[:2,:2]=E;effect[2,2]=np.trace(E@TAU)
        for _ in range(40):
            rho=apply(ks,rho);state=apply(qq,state)
            np.testing.assert_allclose(state[:2,:2]+state[2,2]*TAU,rho,atol=1e-12)
            self.assertAlmostEqual(np.trace(E@rho).real,np.trace(effect@state).real,places=12)
        self.assertGreater(state[2,2].real,0)
        with self.assertRaises(ValueError):lift_qubit(ks,budget*1.01,.2)
        with self.assertRaises(ValueError):lift_qubit(ks,budget*.8,.001)
        # A unitary is on the boundary: no positive replacement-loss budget.
        with self.assertRaises(ValueError):lift_qubit([np.eye(2)],.001,.2)

    def test_escalation_offline_report(self):
        import escalation_report
        self.assertEqual(escalation_report.make_report(),(STUDY/'results/escalation.md').read_text())


if __name__=='__main__':unittest.main()
