"""Independent operator and joint-law checks for the formal connections.

Finite matrices test concrete instances; these tests do not certify infinite
Hilbert-space theorems or the floating interval study's external assumptions.
"""
import importlib.util
from pathlib import Path
import sys
import unittest
import numpy as np
from scipy.linalg import expm

ROOT = Path(__file__).resolve().parents[1] / 'studies' / 'continuum-finite-models'
spec = importlib.util.spec_from_file_location('continuum_connection_model', ROOT / 'model.py')
m = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = m
spec.loader.exec_module(m)


def complete_operators(embedding, reference):
    """Actual Kraus matrices, including the entire orthogonal complement."""
    size = embedding.shape[0]
    row = np.array([1., np.exp(1j * reference)]) / np.sqrt(2)
    plus = np.outer([1., 0.], row)
    minus = np.outer([1., 0.], row * [1., -1.])
    return [embedding @ plus @ embedding.conj().T,
            embedding @ minus @ embedding.conj().T,
            np.eye(size) - embedding @ embedding.conj().T]


def noisy_law(operators, state, eta, visibility):
    raw = np.array([np.linalg.norm(k @ state)**2 for k in operators])
    transition = np.array([[eta*(1+visibility)/2, eta*(1-visibility)/2, 0],
                           [eta*(1-visibility)/2, eta*(1+visibility)/2, 0],
                           [1-eta, 1-eta, 1]])
    return transition @ raw


class ContinuumConnections(unittest.TestCase):
    def test_complete_transport_on_and_off_accessible_subspace(self):
        rng = np.random.default_rng(91)
        for n in (4, 7, 13):
            embedding = m.dft(n)[:, [0, 1]]
            operators = complete_operators(embedding, .37)
            np.testing.assert_allclose(sum(k.conj().T @ k for k in operators),
                                       np.eye(n), atol=2e-14)
            # Completeness must hold outside the two prepared modes too.
            for _ in range(5):
                state = rng.normal(size=n) + 1j*rng.normal(size=n)
                state /= np.linalg.norm(state)
                law = noisy_law(operators, state, .73, .82)
                self.assertTrue(np.all(law >= -1e-14))
                self.assertAlmostEqual(float(law.sum()), 1.)
            outside = m.dft(n)[:, 2]
            np.testing.assert_allclose(noisy_law(operators, outside, 1., 1.), [0, 0, 1], atol=2e-14)

    def test_physical_schrodinger_derivative_arbitrary_site_states(self):
        rng = np.random.default_rng(430)
        ring = m.Ring(length=5.7, mass=2.3, hbar=.61)
        for n in (1, 2, 4, 9):
            h = ring.kinetic(n)
            u = rng.normal(size=n) + 1j*rng.normal(size=n)
            u /= np.linalg.norm(u)
            t = .23
            propagator = lambda time: expm(-1j*time*h/ring.hbar) @ u
            delta = 1e-5
            # Centered finite differences of a matrix exponential, independent
            # of the implementation's diagonal-phase evolution.
            derivative = (propagator(t+delta)-propagator(t-delta))/(2*delta)
            np.testing.assert_allclose(1j*ring.hbar*derivative, h@propagator(t), atol=3e-8, rtol=3e-8)
            f = m.dft(n)
            spectral = f @ m.evolve(f.conj().T@u, ring.lattice(np.arange(n), n), t, ring.hbar)
            np.testing.assert_allclose(spectral, propagator(t), atol=3e-13)
            np.testing.assert_allclose(propagator(0), u, atol=1e-14)

    def test_noisy_n4_witness_and_all_blind_time_bins(self):
        ring = m.Ring()
        n = 4
        embedding = m.dft(n)[:, [0, 1]]
        u = embedding @ (np.ones(2)/np.sqrt(2))
        continuum_h = embedding @ np.diag(ring.continuum([0, 1])) @ embedding.conj().T
        lattice_h = ring.kinetic(n)
        gap = .5-4/np.pi**2
        t = np.pi/gap
        ops = complete_operators(embedding, t/2)
        for eta, v in ((1., 1.), (.8, .9), (0., .9), (.8, 0.)):
            p = noisy_law(ops, expm(-1j*t*continuum_h)@u, eta, v)
            q = noisy_law(ops, expm(-1j*t*lattice_h)@u, eta, v)
            np.testing.assert_allclose(p, [eta*(1+v)/2, eta*(1-v)/2, 1-eta], atol=1e-13)
            np.testing.assert_allclose(q, [eta*(1-v)/2, eta*(1+v)/2, 1-eta], atol=1e-13)
            self.assertAlmostEqual(np.abs(p-q).sum()/2, eta*v)
            if eta*v == 0:
                np.testing.assert_allclose(p, q, atol=1e-13)
        for k in (-2, -1, 0, 1, 2):
            time = 2*np.pi*k/gap
            for reference in (.19, 1.4):
                ops = complete_operators(embedding, reference)
                p = noisy_law(ops, expm(-1j*time*continuum_h)@u, .7, .6)
                q = noisy_law(ops, expm(-1j*time*lattice_h)@u, .7, .6)
                np.testing.assert_allclose(p, q, atol=2e-13)

    def test_shared_scale_exact_one_mode_and_two_mode_obstruction(self):
        ring = m.Ring(length=8., mass=2., hbar=.4)
        for n in (5, 8, 32, 128):
            ec = ring.continuum(np.array([1, 2]))
            ea = ring.lattice(np.array([1, 2]), n)
            scale = ec[0]/ea[0]
            self.assertGreater(scale, 0)
            self.assertAlmostEqual(scale*ea[0], ec[0])
            self.assertAlmostEqual(ec[1]/ec[0], 4.)
            self.assertAlmostEqual(ea[1]/ea[0], 4*np.cos(np.pi/n)**2)
            self.assertGreater(ec[1]-scale*ea[1], 0)
            for time in (0., .7, 4.):
                np.testing.assert_allclose(np.exp(-1j*time*ec[0]/ring.hbar),
                                           np.exp(-1j*time*scale*ea[0]/ring.hbar), atol=1e-14)

    def test_heterogeneous_products_zero_and_unequal_counts(self):
        # Different outcome alphabets; coordinates within a draw are mutually
        # exclusive, never independent Bernoulli observations.
        p = [np.array([.6, .4]), np.array([.2, .5, .3]), np.array([1.])]
        q = [np.array([.55, .45]), np.array([.22, .48, .3]), np.array([1.])]
        for counts in ((0, 0, 0), (1, 0, 0), (1, 3, 2), (2, 1, 0)):
            jp = jq = np.array([1.])
            budget = 0.
            for a, b, count in zip(p, q, counts):
                budget += count*np.abs(a-b).sum()/2
                for _ in range(count):
                    jp, jq = np.kron(jp, a), np.kron(jq, b)
            tv = np.abs(jp-jq).sum()/2
            self.assertLessEqual(tv, min(1., budget)+1e-14)
            # Best deterministic rejection event attains the TV error bound.
            rejection = (jq > jp).astype(float)
            alpha_beta = jp@rejection + 1-jq@rejection
            self.assertAlmostEqual(alpha_beta, 1-tv)
            randomized = np.linspace(0., 1., len(jp))
            self.assertGreaterEqual(jp@randomized + 1-jq@randomized, 1-tv-1e-14)


if __name__ == '__main__':
    unittest.main()
