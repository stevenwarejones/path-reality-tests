"""Adversarial tests of physical normalization, source readers and certificate scope."""
import importlib.util
import io
from pathlib import Path
import sys
import unittest

import numpy as np
from scipy import constants as co
from scipy.integrate import quad
from scipy.special import erf

STUDY = Path(__file__).resolve().parents[1]/"studies/collapse-compatibility"
# Isolate common module names from the other studies' tests.
def local_module(name):
    path = STUDY/(name+".py")
    spec = importlib.util.spec_from_file_location("collapse_"+name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module

r = local_module("response")
s = local_module("spectral_bounds")
i = local_module("inference")
data = local_module("io_data")
_previous = {name: sys.modules.get(name) for name in ["response", "io_data"]}
try:
    sys.modules["response"] = r
    sys.modules["io_data"] = data
    combo = local_module("combinations")
    baseline = local_module("baselines")
finally:
    for name, value in _previous.items():
        if value is None:
            sys.modules.pop(name, None)
        else:
            sys.modules[name] = value


class CollapseResponses(unittest.TestCase):
    def test_cube_force_torque_limit_and_nonnegative(self):
        force, torque = combo.cube_response(1e-7)
        self.assertGreater(force, 0)
        self.assertAlmostEqual(torque/force/(.046**2/6), 1, places=4)
        with self.assertRaises(ValueError): combo.cube_response(.1)

    def test_radiation_neutrality_high_energy_and_colored_suppression(self):
        self.assertAlmostEqual(float(combo.radiation_charge_factor(0.)), 0., places=10)
        self.assertAlmostEqual(float(combo.radiation_charge_factor(1e9))/(54**2+54), 1., places=6)
        energies = np.linspace(1, 30, 2001)
        white = combo.radiation_response(energies, np.ones_like(energies), 1e-7, 0)
        colored = combo.radiation_response(energies, np.ones_like(energies), 1e-7, 1e-12)
        self.assertLess(colored/white, 1e-12)

    def test_count_level_white_collapse_injection(self):
        rng = np.random.default_rng(1729)
        x = np.linspace(0, 10*r.D, 410)
        coefficient = r.triangular_white_coefficient(170000., 158., 1e-7)
        for rate in [0., 1e-8]:
            expected = .4*np.exp(-rate*coefficient)
            counts = rng.poisson(1000*(1+expected*np.cos(2*np.pi*x/r.D+.3)))
            fit = baseline.harmonic_fit(x, counts)
            self.assertLess(abs(fit["visibility"]-expected), 5*fit["visibility_se_gaussian"])

    def test_psd_level_gamma_injection(self):
        rng = np.random.default_rng(31415)
        frequency = np.linspace(3515., 3550., 1468)
        true_b = 2.6
        predicted = baseline.psd_model(frequency, [1., true_b, .1, 3532.75, 1.], 1e6)
        measured = rng.gamma(80, predicted/80)
        fit = baseline.fit_psd(np.array([frequency, measured, measured/np.sqrt(80)]).T)
        self.assertLess(abs(fit["B_phi0_squared_per_hz"]/(true_b*1e-18)-1), .1)
        self.assertEqual(len(fit["omitted_zero_based_indices"]), 6)

    def test_benchmark_escape_preserves_finite_covariance(self):
        g = r.sphere_static_coefficient(1e-7, 2200., 1e-7, 1e-7)
        last = None
        for tau in [1e4, 1e6, 1e8]:
            lam = 10/(g*r.ou_time(.01, tau))
            self.assertAlmostEqual(lam*g*r.ou_time(.01, tau), 10.)
            response = lam*r.ou_spectrum(2*np.pi*.003, tau)
            if last is not None: self.assertLess(response, last/90)
            last = response
        self.assertAlmostEqual(lam/(2*tau)/(10/(g*.01**2)), 1., places=8)

    def test_white_trajectory_independent_time_integral(self):
        mass, velocity, rc = 170000., 158., 1e-7
        time = r.L/velocity
        sep = co.h*time/(mass*co.m_u*r.D)
        independent = 2*mass**2*quad(lambda t: -np.expm1(-(sep*t/time)**2/(4*rc**2)), 0, time, epsabs=1e-14)[0]
        self.assertAlmostEqual(r.triangular_white_coefficient(mass, velocity, rc)/independent, 1, places=10)

    def test_small_sphere_point_limit(self):
        radius, density, separation, rc = 1e-11, 2200., 1e-7, 1e-7
        mass = 4*np.pi/3*radius**3*density
        point = (mass/co.m_u)**2*(-np.expm1(-separation**2/(4*rc**2)))
        self.assertAlmostEqual(r.sphere_static_coefficient(radius, density, separation, rc)/point, 1, places=7)

    def test_ou_time_white_zero_short_and_independent_integral(self):
        self.assertEqual(r.ou_time(.01, 0), .01)
        self.assertEqual(r.ou_time(0, 1), 0)
        for tau in [1e-3, .01, 1., 1e8]:
            expected = quad(lambda u: (.01-u)*np.exp(-u/tau)/tau, 0, .01, epsabs=1e-25)[0]
            self.assertAlmostEqual(r.ou_time(.01, tau)/expected, 1, places=10)

    def test_moving_response_is_not_static(self):
        v, rc, tau = 158., 1e-7, 1e-6
        independent = quad(lambda u: np.exp(-2*rc*u/(v*tau)-u*u)*2*rc/(v*tau), 0, 10)[0]
        self.assertAlmostEqual(r.moving_local_ou_factor(v, rc, tau)/independent, 1, places=10)
        self.assertLess(r.moving_local_ou_factor(v, rc, tau), .002)
        self.assertGreater(r.ou_time(.01, tau)/.01, .99)

    def test_prescribed_static_paths_frequency_kernel(self):
        nodes, weights = np.polynomial.legendre.leggauss(60)
        time, rc, d = .01, 1e-7, 2e-7
        t, w = (nodes+1)*time/2, weights*time/2
        a = np.zeros((len(t), 3))
        b = a+np.array([d, 0, 0])
        for omega in [0, 50, 1000]:
            expected = (1-np.exp(-d*d/(4*rc*rc)))*time*time*np.sinc(omega*time/(2*np.pi))**2
            self.assertAlmostEqual(r.prescribed_path_kernel(t, w, a, b, omega, rc)/expected, 1, places=10)
        self.assertEqual(r.prescribed_path_kernel(t, w, a, a, 100, rc), 0.)

    def test_layer_force_independent_fourier_integral(self):
        rc, thickness = 1e-7, 370e-9
        edges = np.arange(48)*thickness
        density = np.array([7170 if k%2==0 else 2200 for k in range(47)])
        def integrand(z):
            q = z/rc
            numerator = np.sum(density*(np.exp(1j*q*edges[1:])-np.exp(1j*q*edges[:-1])))
            return np.exp(-z*z)*abs(numerator)**2/rc
        iz = 2*quad(integrand, 0, 10, limit=2000, epsabs=.01, epsrel=1e-9)[0]
        def transverse(length):
            # Independent real-space autocorrelation integral.
            return 2*np.sqrt(np.pi)/rc*quad(lambda u: (length-u)*np.exp(-u*u/(4*rc*rc)), 0, min(length, 20*rc), epsabs=1e-24)[0]
        independent = co.hbar**2*rc**3/(np.pi**1.5*co.m_u**2)*iz*transverse(113e-6)*transverse(82e-6)
        self.assertAlmostEqual(r.layered_force_coefficient(rc)/independent, 1, places=7)

    def test_zero_grating_and_density_scaling(self):
        s0, s1 = r.optical_coefficients(np.array([170000.]), np.array([158.]), [0, 0, 0])
        np.testing.assert_allclose(s0, 1)
        np.testing.assert_allclose(s1, 0)
        g = r.sphere_static_coefficient(1e-7, 2200, 1e-7, 1e-7)
        self.assertAlmostEqual(r.sphere_static_coefficient(1e-7, 4400, 1e-7, 1e-7)/g, 4.)


class CollapseCertificates(unittest.TestCase):
    def example(self):
        return dict(knots=[0,1,2], kernels=[[1,1,0],[0,1,1]], macro=[1,2,1],
                    tail_slopes=[0,0], macro_tail_slope=0, weights=[1,1], bounds=[2,3])

    def test_exact_continuum_and_lp(self):
        args = self.example()
        self.assertEqual(s.affine_certificate(**args)["upper_bound"], "5")
        result = s.finite_bound(args["kernels"], args["bounds"], args["macro"])
        self.assertAlmostEqual(result["value"], 5)

    def test_unmeasured_band_escapes(self):
        result = s.finite_bound([[1,0,0],[0,1,0]], [2,3], [1,1,1])
        self.assertEqual(result["status"], "unbounded")
        self.assertEqual(result["escape_column"], 2)

    def test_negative_weights_and_failed_domination_rejected(self):
        for weights in [[-1,3], [1,0]]:
            args = self.example(); args["weights"] = weights
            with self.assertRaises(ValueError): s.affine_certificate(**args)

    def test_tail_omission_and_tail_violation_rejected(self):
        for slope in [None, 1]:
            args = self.example(); args["macro_tail_slope"] = slope
            with self.assertRaises(ValueError): s.affine_certificate(**args)

    def test_between_grid_counterexample_rejected(self):
        # The endpoint grid [0,2] passes, but the actual declared function peaks at 1.
        args = self.example(); args["macro"] = [1,3,1]
        with self.assertRaises(ValueError): s.affine_certificate(**args)

    def test_float_and_nonfinite_inputs_rejected(self):
        args = self.example(); args["weights"] = [1.,1]
        with self.assertRaises(ValueError): s.affine_certificate(**args)
        with self.assertRaises(ValueError): s.finite_bound([[np.nan]], [1], [1])

    def test_combining_need_not_help(self):
        alone = s.finite_bound([[1,1]], [2], [1,1])
        joint = s.finite_bound([[1,1],[1,1]], [2,100], [1,1])
        self.assertEqual(alone["value"], joint["value"])

    def test_gaussian_calibration_and_injection_recovery(self):
        result = i.calibration_demo(20260926, 5000, .05)
        self.assertEqual(result["cases"][0]["false_exclusions_of_true_response"], 0)
        injected = result["cases"][1]
        self.assertLess(injected["false_exclusions_of_true_response"]/5000, .065)
        np.testing.assert_allclose(injected["recovered_spectrum_from_mean"], [.7,.4], atol=.01)

    def test_parser_rejects_coordinate_and_value_corruption(self):
        for content in [b"1 2\n0 3\n", b"0 1\n1 nan\n", b"0 1 2\n1 2 3\n"]:
            with self.assertRaises(ValueError): data.table(content, 2)


if __name__ == "__main__":
    unittest.main()
