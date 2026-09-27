"""Independent matrix, probability, approximation and statistical checks."""
import importlib.util
import itertools
import json
from pathlib import Path
import sys
import unittest
import numpy as np
from scipy.linalg import expm
from scipy.stats import binom

ROOT=Path(__file__).resolve().parents[1]/'studies'/'continuum-finite-models'
spec=importlib.util.spec_from_file_location('continuum_model',ROOT/'model.py')
m=importlib.util.module_from_spec(spec); sys.modules[spec.name]=m; spec.loader.exec_module(m)


class ContinuumFinite(unittest.TestCase):
    def test_cyclic_matrix_all_small_sizes(self):
        ring=m.Ring(length=5.3,mass=2.7,hbar=.43)
        for n in range(1,18):
            f=m.dft(n); h=ring.kinetic(n)
            np.testing.assert_allclose(f.conj().T@f,np.eye(n),atol=4e-14)
            np.testing.assert_allclose(h@f,f*ring.lattice(np.arange(n),n),atol=1e-12)
            np.testing.assert_allclose(np.linalg.eigvalsh(h),np.sort(ring.lattice(np.arange(n),n)),atol=1e-12)
        self.assertEqual(ring.kinetic(1)[0,0],0)
        self.assertLess(ring.kinetic(2)[0,1],0)

    def test_independent_matrix_evolution_and_spectral_counterexample(self):
        rng=np.random.default_rng(9041); ring=m.Ring(); n=13
        modes=m.check_modes(np.arange(-3,4),n); f=m.dft(n)
        for _ in range(8):
            c=rng.normal(size=7)+1j*rng.normal(size=7); c/=np.linalg.norm(c)
            full=np.zeros(n,complex); full[modes%n]=c
            q,_=np.linalg.qr(rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)))
            t=float(rng.uniform(-3,3)); evolved=expm(-1j*t*ring.kinetic(n))@(f@full)
            expected=f[:,modes%n]@m.evolve(c,ring.lattice(modes,n),t)
            np.testing.assert_allclose(evolved,expected,atol=2e-13)
            labels=np.where(np.arange(n)<=n//2,np.arange(n),np.arange(n)-n)
            hs=f@np.diag(ring.continuum(labels))@f.conj().T
            spectral=expm(-1j*t*hs)@(f@full)
            continuum=f[:,modes%n]@m.evolve(c,ring.continuum(modes),t)
            # Arbitrary basis measurement and independently generated failure rates.
            efficiencies=rng.uniform(0,1,n)
            def detect(state):
                raw=np.abs(q.conj().T@state)**2
                return np.r_[efficiencies*raw,np.sum((1-efficiencies)*raw)]
            np.testing.assert_allclose(detect(spectral),detect(continuum),atol=5e-13)
            self.assertAlmostEqual(detect(spectral).sum(),1)
        # Exact spectral couplings are not nearest-neighbor.
        self.assertGreater(abs(hs[0,2]),.01)

    def test_cosine_remainder_and_second_order_rate(self):
        ring=m.Ring(); previous=None
        for n in (16,32,64,128,256):
            j=np.arange(-3,4); diff=ring.continuum(j)-ring.lattice(j,n)
            upper=(ring.length/n)**2*ring.k(j)**4/24
            self.assertTrue(np.all(diff >= -1e-14))
            self.assertTrue(np.all(diff <= upper+1e-13))
            if previous is not None: self.assertTrue(3.7 < previous/diff[-1] < 4.1)
            previous=diff[-1]
        for x in np.linspace(-20,20,4001):
            gap=x*x/2-(1-np.cos(x))
            self.assertGreaterEqual(gap,-1e-14)
            self.assertLessEqual(gap,x**4/24+1e-13)

    def test_tail_and_probability_bound(self):
        ring=m.Ring(); rng=np.random.default_rng(182)
        for _ in range(30):
            c=rng.normal(size=11)+1j*rng.normal(size=11); c/=np.linalg.norm(c)
            modes=np.arange(-5,6); keep=np.abs(modes)<=2
            projected,tail,norm_bound=m.project_normalize(c,keep)
            self.assertAlmostEqual(np.linalg.norm(c-projected),norm_bound,places=12)
            self.assertAlmostEqual(np.sqrt(1-abs(np.vdot(c,projected))**2),np.sqrt(tail),places=12)
            t=.4; n=128
            a=m.evolve(c,ring.continuum(modes),t)
            b=m.evolve(projected,ring.lattice(modes,n),t)
            q,_=np.linalg.qr(rng.normal(size=(11,11))+1j*rng.normal(size=(11,11)))
            self.assertLessEqual(m.tv(abs(q@a)**2,abs(q@b)**2),m.approximation(ring,n,2,t,tail)['tv_bound']+1e-13)
        out,tail,_=m.project_normalize([0,1],[True,False])
        self.assertIsNone(out); self.assertEqual(tail,1)

    def test_readout_from_born_rule(self):
        for theta in np.linspace(-7,7,19):
            for eta,v in itertools.product((0,.3,1),(0,.6,1)):
                rho=np.array([[1,v*np.exp(1j*theta)],[v*np.exp(-1j*theta),1]])/2
                plus=eta*np.ones((2,2))/2
                minus=eta*np.array([[1,-1],[-1,1]])/2
                failure=(1-eta)*np.eye(2)
                p=np.array([np.trace(rho@e).real for e in (plus,minus,failure)])
                np.testing.assert_allclose(p,m.readout(theta,eta,v),atol=1e-15)

    def test_blind_spots_and_calibration_degeneracies(self):
        ring=m.Ring(); n=16
        self.assertEqual(m.phase_gap(ring,n,-2,2,9),0)
        self.assertEqual(m.phase_gap(ring,n,0,2,0),0)
        delta=m.phase_gap(ring,n,0,2,1)
        wrapped=delta*2*np.pi/abs(delta)
        np.testing.assert_allclose(m.readout(wrapped),m.readout(0),atol=1e-14)
        for eta,v in ((0,1),(.8,0)):
            np.testing.assert_array_equal(m.readout(.91,eta,v),m.readout(0,eta,v))
        self.assertAlmostEqual(m.tv(m.readout(.7),m.readout(.7)),0) # offset absorbs gap
        scale=float(ring.continuum(1)/ring.lattice(1,n))
        self.assertAlmostEqual(scale*ring.lattice(1,n),ring.continuum(1))
        self.assertNotAlmostEqual(scale*ring.lattice(2,n),ring.continuum(2))

    def test_aliasing_rejected_not_refitted(self):
        for modes,n in (([-4,4],8),([0,1.5],8),([0,0],8),([],8)):
            with self.assertRaises(ValueError): m.check_modes(modes,n)
        for n in (0,True,3.5):
            with self.assertRaises(ValueError): m.check_sites(n)

    def test_continuous_nuisance_enclosure(self):
        rng=np.random.default_rng(2026)
        for phase,r in ((0,.1),(np.pi,.3),(2*np.pi,.2),(17,9)):
            lo,hi=m.outcome_box(phase,r,.7,.1,.42,.12,.03)
            for _ in range(100):
                eta=rng.uniform(.6,.8); w=rng.uniform(.30,.54)
                p=m.readout(rng.uniform(phase-r,phase+r),eta,w/eta)
                self.assertTrue(np.all(p>=lo-1e-15)); self.assertTrue(np.all(p<=hi+1e-15))
        self.assertEqual(m.cosine_interval(-.1,.1)[1],1)
        self.assertEqual(m.cosine_interval(3,4)[0],-1)
        self.assertEqual(m.cosine_interval(0,7),(-1,1))

    def test_iid_joint_tv_and_test_error(self):
        p=np.array([.4,.35,.25]); q=np.array([.5,.4,.1])
        pn=qn=np.array([1.])
        for n in range(1,6):
            pn=np.kron(pn,p); qn=np.kron(qn,q)
            distance=m.tv(pn,qn)
            self.assertLessEqual(distance,n*m.tv(p,q)+1e-14)
            reject=qn>pn; alpha=pn[reject].sum(); beta=qn[~reject].sum()
            self.assertAlmostEqual(alpha+beta,1-distance)

    def test_strict_budget_and_exact_binomial_power(self):
        for gap in (.01,.1,.4,.9):
            n=m.required_trials(gap,72)
            self.assertLess(m.radius(n,72,.025)+m.radius(n,72,.1),gap)
            self.assertGreaterEqual(m.radius(n-1,72,.025)+m.radius(n-1,72,.1),gap)
        n=m.required_trials(.2,3); threshold=.2+m.radius(n,3,.025)
        cutoff=int(np.floor(n*threshold))
        self.assertLessEqual(binom.sf(cutoff,n,.2),.025)
        self.assertGreaterEqual(binom.sf(cutoff,n,.4),.9)
        self.assertIsNone(m.required_trials(0))

    def test_register_snapshot_and_provenance(self):
        result=json.loads((ROOT/'results/sensitivity.json').read_text())
        self.assertEqual(len(result['rows']),660)
        self.assertEqual(result['data_inputs'],[])
        for r in result['rows']:
            if r['efficiency']==0 or r['visibility']==0:
                self.assertEqual(r['certified_coordinate_gap'],0)
            self.assertIsNone(r['total_source_attempts'])
            if r['trials_per_setting'] is not None:
                self.assertEqual(r['total_eligible_trials'],24*r['trials_per_setting']+3*r['calibration_trials_per_stratum'])
        register=json.loads((ROOT/'questions.json').read_text())
        self.assertEqual({r['id'] for r in register},{f'CF{i:02d}' for i in range(1,13)})
        for row in register:
            self.assertTrue((ROOT/row['evidence']).exists())
        report=(ROOT/'README.md').read_text()
        for row in register: self.assertIn(row['id'],report)
        for path in ROOT.rglob('*'):
            if path.is_file() and '__pycache__' not in path.parts:
                self.assertIn(path.suffix,{'.md','.py','.json','.svg','.txt'})


if __name__=='__main__': unittest.main()
