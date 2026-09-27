"""Finite-sample calibration, exact normalization identity, and stability failures."""
import sys,json,math,unittest
from pathlib import Path
from fractions import Fraction as F
import numpy as np
from mpmath import iv
D=Path(__file__).resolve().parents[1]/'studies/detector-normalized-dynamics'
sys.path.insert(0,str(D))
import adequacy as ad
import normalized as nm
import stability as st
from dynamics_common import baseline_models,operations
from sqd_models import PAULI

class DetectorNormalizedDynamics(unittest.TestCase):
    def test_finite_mgf_matches_independent_binomial_enumeration(self):
        for n in [1,4,8]:
            for p in [.1,.5,.83]:
                t=F(1,16);value=ad.log_mgf(n,iv.mpf(p),t)
                total=0.
                for k in range(n+1):
                    dev=0.
                    if k:dev+=2*k*math.log(k/(n*p))
                    if n-k:dev+=2*(n-k)*math.log((n-k)/(n*(1-p)))
                    total+=math.comb(n,k)*p**k*(1-p)**(n-k)*math.exp(float(t)*dev)
                self.assertAlmostEqual(float(value.mid),math.log(total),places=12)

    def test_chernoff_bound_dominates_exact_small_sample_tail(self):
        n=8;p=.2;values=[];masses=[]
        for k in range(n+1):
            values.append(float(ad.deviance(k,n,iv.mpf(p)).mid));masses.append(math.comb(n,k)*p**k*(1-p)**(n-k))
        for k in range(n+1):
            out=ad.evaluate([{'k':k,'n':n}],[p],F(1,8));exact=sum(m for m,d in zip(masses,values) if d>=values[k]-1e-12)
            self.assertGreaterEqual(out['tail_upper']['decimal']+1e-12,exact)

    def test_precision_is_restored_for_inherited_certificates(self):
        old=iv.dps
        ad.evaluate([{'k':0,'n':2}],[.5],F(1,8))
        self.assertEqual(iv.dps,old)

    def test_perfect_within_row_dependence_breaks_binomial_calibration(self):
        # K=0 or n, each with probability 1/2, while marginal shot probability is 1/2.
        bound=ad.evaluate([{'k':0,'n':20}],[.5],F(1,4))['tail_upper']['decimal']
        self.assertLess(bound,.025)  # Actual dependent-law tail for these two extremes is one.

    def test_degenerate_probability_and_invalid_tilt_are_rejected(self):
        with self.assertRaises(ValueError):ad.enclosure(0)
        with self.assertRaises(ValueError):ad.log_mgf(10,iv.mpf(.5),F(1,2))
        with self.assertRaises(ValueError):ad.deviance(11,10,iv.mpf(.5))

    def test_normalization_identity_for_nonunital_channel(self):
        gamma=.3;K=[np.diag([1,np.sqrt(1-gamma)]),np.array([[0,np.sqrt(gamma)],[0,0]])]
        E=np.array([[.15,.1j],[-.1j,.85]]);a=.17;b=.6
        dual=lambda A:sum(k.conj().T@A@k for k in K)
        diameter=lambda A:np.ptp(np.linalg.eigvalsh(A))
        np.testing.assert_allclose(dual(a*np.eye(2)+b*E),a*np.eye(2)+b*dual(E),atol=1e-15)
        self.assertAlmostEqual(diameter(dual(E))/diameter(E),diameter(dual(a*np.eye(2)+b*E))/diameter(a*np.eye(2)+b*E),places=13)

    def test_normalized_snapshot_and_path_binding(self):
        self.assertEqual(nm.HERE,D)
        self.assertEqual(nm.run(),json.loads((D/'results/normalized.json').read_text()))

    def test_changed_gate_invalidates_affine_witness_identity(self):
        x=baseline_models()['joint']['x'];y=x.copy();y[2]+=.01
        with self.assertRaises(ValueError):nm.affine_relation(x,y)

    def test_response_range_tightness_has_physical_kraus_realization(self):
        a=.0452128544;b=.7681878105;K=[]
        for i,p in enumerate([a,b]):
            for j,q in enumerate([1-p,p]):
                k=np.zeros((2,2));k[j,i]=np.sqrt(q);K.append(k)
        np.testing.assert_allclose(sum(k.T@k for k in K),np.eye(2),atol=1e-15)
        E=np.diag([0.,1.]);dual=sum(k.T@E@k for k in K)
        np.testing.assert_allclose(dual,np.diag([a,b]),atol=1e-15)
        self.assertAlmostEqual(np.ptp(np.linalg.eigvalsh(dual)),.7229749561,places=14)

    def test_unitary_errors_do_not_change_normalized_contrast(self):
        E=np.diag([.1,.8]);U=(np.eye(2)-1j*PAULI[1])/np.sqrt(2)
        self.assertAlmostEqual(np.ptp(np.linalg.eigvalsh(U.conj().T@E@U))/np.ptp(np.linalg.eigvalsh(E)),1,places=14)
        rho=np.diag([1.,0.]);self.assertGreater(abs(np.trace(E@rho)-np.trace(E@U@rho@U.conj().T)),.3)

    def test_overlap_stability_countermodel_and_snapshot(self):
        out=st.run();self.assertEqual(out,json.loads((D/'results/stability.json').read_text()))
        self.assertEqual((out['matched_words'],out['matched_rows']),(8,23))
        E0=np.diag([.01,0]);E1=np.diag([.01,1]);rho=np.diag([1.,0.])
        self.assertEqual(np.trace(E0@rho),np.trace(E1@rho));self.assertEqual(np.linalg.norm(E1-E0,2),1)

if __name__=='__main__':unittest.main()
