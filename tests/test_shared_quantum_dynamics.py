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


if __name__=='__main__':unittest.main()
