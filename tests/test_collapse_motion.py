"""Independent physical checks for motion, finite observations and shared-frame scope."""
import importlib.util
from pathlib import Path
import sys
import unittest
import numpy as np
from scipy import constants as co
from scipy.integrate import quad
from scipy.signal.windows import blackman

STUDY=Path(__file__).resolve().parents[1]/'studies/collapse-compatibility'
def load(name):
    spec=importlib.util.spec_from_file_location('motion_test_'+name,STUDY/(name+'.py'))
    m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m);return m
r=load('response')
previous=sys.modules.get('response');sys.modules['response']=r
try:m=load('moving_response')
finally:
    if previous is None:sys.modules.pop('response',None)
    else:sys.modules['response']=previous


class MovingPhysics(unittest.TestCase):
    def test_moving_sphere_static_matches_independent_fourier(self):
        for radius in [1e-7,3e-7]:
            g=m.sphere_overlap(0,radius,2200,1e-7)-m.sphere_overlap(1e-7,radius,2200,1e-7)
            self.assertAlmostEqual(g/r.sphere_static_coefficient(radius,2200,1e-7,1e-7),1,places=10)
            self.assertAlmostEqual(m.moving_sphere_exponent_per_a(.01,np.inf,0,radius,2200,1e-7,1e-7)/(g*.01**2/2),1,places=12)

    def test_sphere_force_closed_form_matches_radial_integral(self):
        radius,rc=25e-7,1e-7;a=radius/rc;mass=4*np.pi*radius**3*7430/3
        from scipy.special import spherical_jn
        integral=quad(lambda q:q**4*np.exp(-q*q)*(3*spherical_jn(1,a*q)/(a*q))**2,0,10,epsabs=1e-14,limit=500)[0]
        independent=4/(3*np.sqrt(np.pi))*(co.hbar*mass/r.M0/rc)**2*integral
        self.assertAlmostEqual(m.sphere_force_q0(radius,7430,rc)/independent,1,places=8)

    def test_finite_tau_quadrature_vs_fourier_not_narrow_grid(self):
        length,rc,v,omega,tau=5.,1.,2.,.7,3.
        i0=m.top_hat_covariance(length,0,rc)
        direct=quad(lambda k:np.exp(-k*k*rc*rc)*length**2*np.sinc(k*length/(2*np.pi))**2/(1+tau*tau*(omega-v*k)**2),-12,12,points=[omega/v],epsabs=1e-10,limit=400)[0]/i0
        self.assertAlmostEqual(m.lateral_ou_factor(omega,tau,v,length,rc)/direct,1,places=7)

    def test_distributional_limit_and_zero_speed_are_different(self):
        omega,v,length,rc=.018,30000.,.046,1e-7
        asym=m.lateral_frozen_factor(omega,v,length,rc)
        self.assertLess(abs(m.lateral_ou_factor(omega,1.,v,length,rc)/asym-1),length/v)
        self.assertLess(1e12*m.lateral_ou_factor(omega,1e12,0,length,rc),asym*.01)
        with self.assertRaises(ValueError):m.lateral_frozen_factor(omega,0,length,rc)

    def test_moving_benchmark_not_static_normalization(self):
        stationary=m.moving_sphere_exponent_per_a(.01,1e6,0,1e-7,2200,1e-7,1e-7)
        moving=m.moving_sphere_exponent_per_a(.01,1e6,30000,1e-7,2200,1e-7,1e-7)
        self.assertGreater(stationary/moving,1e8)
        frozen=m.moving_sphere_exponent_per_a(.01,np.inf,30000,1e-7,2200,1e-7,1e-7)
        self.assertAlmostEqual(moving/frozen,1,places=10)

    def test_benchmark_orientation_matters(self):
        perpendicular=m.moving_sphere_exponent_per_a(.01,np.inf,30000,1e-7,2200,1e-7,1e-7,0.)
        parallel=m.moving_sphere_exponent_per_a(.01,np.inf,30000,1e-7,2200,1e-7,1e-7,1.)
        self.assertGreater(perpendicular/parallel,1e6)

    def test_blackman_dc_formula_matches_direct_samples(self):
        for sym in [False,True]:
            n,dt,f=1001,.1,.137
            w=blackman(n,sym=sym);direct=2*dt*abs(w@np.exp(-2j*np.pi*f*dt*np.arange(n)))**2/(w@w)
            self.assertAlmostEqual(m.blackman_dc_gain(n,dt,f,sym)/direct,1,places=7)

    def test_window_dc_and_detrend_injection(self):
        b=m.periodogram_weights(256,1.,.013,detrend=False)
        bd=m.periodogram_weights(256,1.,.013,detrend=True)
        q0,a=3.,.2
        rng=np.random.default_rng(20260927)
        amplitudes=rng.normal(0,np.sqrt(a*q0/2),size=20000)
        measured=np.mean(abs(amplitudes*np.sum(b))**2)
        self.assertLess(abs(measured/(a*m.frozen_periodogram_per_a(q0,b))-1),.03)
        self.assertLess(m.frozen_periodogram_per_a(q0,bd),1e-23*m.frozen_periodogram_per_a(q0,b))

    def test_delayed_pair_can_cancel_even_at_large_equal_time_separation(self):
        # A translated frozen field reaches the second body later. Low omega
        # cancels; half a temporal cycle is anticorrelated, not independent.
        self.assertEqual(m.relative_pair_factor(0),0)
        self.assertAlmostEqual(float(m.relative_pair_factor(np.pi)),1)
        self.assertLess(float(m.relative_pair_factor(.003*2*np.pi*.376/30000)),1e-12)

    def test_uniform_geometry_bound_dominates_resolved_axis(self):
        length,rc=.046,1e-7
        for v in [.001,1.,30000.]:
            # A separable torque-axis covariance with arbitrary amplitude.
            q0=1.;tv=1e20
            upper=m.moving_psd_upper_per_a(q0,tv,np.sqrt(3)*length,rc,v)
            for omega in [.0,.018,1.]:
                exact=2*q0*m.lateral_frozen_factor(omega,v,length,rc)
                self.assertGreaterEqual(upper,exact)

    def test_prescribed_two_time_kernel_checks_local_map(self):
        # Temporarily expose response only to the function's local import.
        prev=sys.modules.get('response');sys.modules['response']=r
        try:exact=m.triangle_path_exponent_per_a(171283.8,158,30000,1e-7,1e6)
        finally:
            if prev is None:sys.modules.pop('response',None)
            else:sys.modules['response']=prev
        local=np.sqrt(np.pi)*1e-7/30000*r.triangular_white_coefficient(171283.8,158,1e-7)
        self.assertAlmostEqual(exact/local,1,places=6)

    def test_bad_bound_and_zero_speed_fail_closed(self):
        with self.assertRaises(ValueError):m.covariance_length_bound(1.,1.,1.,1.)
        with self.assertRaises(ValueError):m.moving_psd_upper_per_a(1.,1e20,1.,1.,0.)


if __name__=='__main__':unittest.main()
