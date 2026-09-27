"""Independent controls for the feasibility bounds, including failure directions."""
import json
from pathlib import Path
import sys
import unittest
import numpy as np
from numpy.polynomial.hermite import hermgauss
from scipy.linalg import expm
from scipy.integrate import quad

STUDY=Path(__file__).resolve().parents[1]/'studies/moving-noise-boundary'
sys.path.insert(0,str(STUDY))
from noise_kernels import (quantum_tv_bound,conditional_tv_bound,harmonic_visibility_bound,
    sphere_spatial_difference,sphere_difference_lower,band_held_lower,held_kernel,
    line_periodogram_gain,separated_band_gain_upper,blackman_norm,response_covariance,
    scalar_phase_exponent)


class MovingNoiseBoundary(unittest.TestCase):
    def test_noncommuting_exact_gaussian_quantum_evolution(self):
        # Independent Gaussian-Hermite integration of actual matrix exponentials.
        x=np.array([[0,1],[1,0]],complex);z=np.diag([1.,-1.]);rho=np.array([[1,0],[0,0]],complex)
        nodes,weights=hermgauss(48)
        for sigma,T in [(.03,.4),(.2,.7),(.6,.8)]:
            h0=1.7*x+.4*z;u0=expm(-1j*T*h0);expected=u0@rho@u0.conj().T
            noisy=np.zeros((2,2),complex)
            for n,w in zip(nodes,weights):
                u=expm(-1j*T*(h0+np.sqrt(2)*sigma*n*z))
                noisy+=w/np.sqrt(np.pi)*(u@rho@u.conj().T)
            distance=np.linalg.norm(noisy-expected,ord='nuc')/2
            self.assertLessEqual(distance,quantum_tv_bound(sigma*sigma,1,T)*(1+1e-12))
            self.assertGreater(distance,0)

    def test_commuting_gaussian_normalization_and_saturation(self):
        for variance in [1e-9,.03,.5,100]:
            exact=.5*(-np.expm1(-2*variance))
            self.assertLessEqual(exact,quantum_tv_bound(variance,1,1))
        self.assertEqual(quantum_tv_bound(1e20,1,1),1.)
        self.assertEqual(quantum_tv_bound(0,1,1,.02),.02)
        with self.assertRaises(ValueError): quantum_tv_bound(-1,1,1)

    def test_postselection_is_not_unconditional_distance(self):
        p=np.array([.001,.009,.99]);q=np.array([.002,.009,.989])
        eps=np.abs(p-q).sum()/2
        conditional=np.abs(p[:2]/p[:2].sum()-q[:2]/q[:2].sum()).sum()/2
        self.assertGreater(conditional,eps*50)
        self.assertLessEqual(conditional,conditional_tv_bound(eps,.01))

    def test_signed_harmonic_fit_bound_with_nonuniform_positions(self):
        rng=np.random.default_rng(432);x=np.linspace(0,2*np.pi,41)+rng.normal(0,.03,41)
        design=np.column_stack([np.ones(41),np.cos(x),np.sin(x)])
        inverse=np.linalg.pinv(design);row_l1=np.abs(inverse).sum(axis=1)
        mean=.08;visibility=.35;beta=np.array([mean,-mean*visibility,0]);eps=1e-4
        bound=harmonic_visibility_bound(eps,mean,visibility,row_l1)
        for _ in range(80):
            perturbed=inverse@(design@beta+rng.choice([-eps,eps],41))
            changed=np.hypot(*perturbed[1:])/perturbed[0]
            self.assertLessEqual(abs(changed-visibility),bound)

    def test_spatial_certificate_against_independent_prior_fourier_integral(self):
        sys.path.insert(0,str(STUDY.parent/'collapse-compatibility'))
        from moving_response import moving_sphere_exponent_per_a
        for rc in [.5e-7,1e-7,2e-7]:
            exact=sphere_spatial_difference(1e-7,2200,1e-7,rc)
            independent=2*moving_sphere_exponent_per_a(.01,np.inf,0,1e-7,2200,1e-7,rc)/.01**2
            self.assertAlmostEqual(exact/independent,1.,places=8)
            self.assertLess(sphere_difference_lower(1e-7,2200,1e-7,rc),exact)
        self.assertGreaterEqual(float(held_kernel(1,.01)),band_held_lower(1,.01))

    def test_discrete_line_kernel_against_direct_double_covariance(self):
        n=93;dt=.02;t=np.arange(n)*dt;f=8.321;nu=.137123
        covariance=np.cos(2*np.pi*nu*(t[:,None]-t[None,:]))
        for symmetric in [False,True]:
            divisor=n-1 if symmetric else n
            w=.42-.5*np.cos(2*np.pi*np.arange(n)/divisor)+.08*np.cos(4*np.pi*np.arange(n)/divisor)
            b=w*np.exp(-2j*np.pi*f*t)
            direct=2*dt*(b.conj()@covariance@b).real/(w@w)
            self.assertAlmostEqual(direct/float(line_periodogram_gain(f,nu,n,dt,symmetric)),1,places=7)
            self.assertAlmostEqual(w@w,blackman_norm(n,symmetric),places=11)

    def test_uniform_bound_covers_atoms_between_bins_and_narrow_spectra(self):
        n=8192;dt=1e-4;fmin=3520.;fmax=3540.
        bound=separated_band_gain_upper(fmin,fmax,1,n,dt)
        for nu in [.1,.137123456,.5,1.,1.-1e-12]:
            for symmetric in [False,True]:
                values=line_periodogram_gain(np.linspace(fmin,fmax,131),nu,n,dt,symmetric)
                self.assertLessEqual(float(values.max()),bound)
        # A grid excluding an exact line center badly misses a narrowing response.
        widths=np.array([1e-1,1e-4,1e-8]);offset=.012345
        offgrid=widths/(offset**2+widths**2)
        self.assertLess(offgrid[-1],offgrid[0]*1e-5)
        self.assertAlmostEqual(quad(lambda x:1/(1+x*x),-np.inf,np.inf)[0],np.pi,places=10)

    def test_orthogonal_detrending_bound_and_DC_removal(self):
        n=91;dt=.03;t=np.arange(n)*dt
        design=np.column_stack([np.ones(n),t]);P=np.eye(n)-design@np.linalg.pinv(design)
        w=np.blackman(n);b=P@(w*np.exp(-2j*np.pi*.73*t))
        for nu in [0,.1,.89123]:
            c=np.cos(2*np.pi*nu*(t[:,None]-t[None,:]))
            val=2*dt*(b.conj()@c@b).real/(w@w)
            self.assertLessEqual(val,2*n*dt+1e-12)
            if nu==0:self.assertLess(abs(val),1e-14)

    def test_doppler_sign_acceleration_and_delayed_correlations(self):
        t=np.linspace(0,2,41);v=3.;k=np.array([[2.,0,0]]);omega=np.array([5.]);weights=[.7]
        c=response_covariance(t,np.column_stack([v*t,t*0,t*0]),k,omega,weights)
        expected=.7*np.cos((5-2*v)*(t[:,None]-t[None,:]))
        np.testing.assert_allclose(c,expected,atol=2e-15)
        accelerated=np.column_stack([v*t+.3*t*t,t*0,t*0])
        ca=response_covariance(t,accelerated,k,omega,weights)
        self.assertGreater(np.max(abs(ca-c)),.1)
        self.assertGreaterEqual(np.linalg.eigvalsh(ca).min(),-1e-13)
        # Two samples at different times and places can be exactly correlated.
        delayed=response_covariance([0,2],[[0,0,0],[6,0,0]],[[2,0,0]],[6],[1])
        self.assertAlmostEqual(delayed[0,1],1.)

    def test_finite_time_frozen_and_fast_motion_limits(self):
        t=np.linspace(0,1,2001);q=np.ones(len(t))/2000;q[[0,-1]]/=2
        k=np.array([[1.,0,0]]);omega=np.array([0.]);weights=np.array([1.]);d=np.array([.4,0,0])
        static=np.zeros((len(t),3))
        zero=scalar_phase_exponent(t,static,static+d,k,omega,weights,q)
        self.assertAlmostEqual(zero,1-np.cos(.4),places=12)
        moving=static.copy();moving[:,0]=2*np.pi*t
        fast=scalar_phase_exponent(t,moving,moving+d,k,omega,weights,q)
        self.assertLess(fast,1e-25)
        white=quad(lambda f:2*float(held_kernel(f,1)),0,200,limit=1000)[0]
        self.assertLess(abs(white-1),2/(np.pi**2*200))

    def test_correlated_background_invalidates_component_bound(self):
        # field = 100 Z; background = -100 Z + .1 X. PSD matrix remains positive.
        covariance=np.array([[10000.,-10000.],[-10000.,10000.01]])
        self.assertGreater(np.linalg.eigvalsh(covariance).min(),0.)
        total=np.ones(2)@covariance@np.ones(2)
        self.assertAlmostEqual(total,.01,places=10)
        self.assertGreater(covariance[0,0]/total,999999)

    def test_claim_gates_and_removal_do_not_enable_exclusion(self):
        r=json.loads((STUDY/'results/boundary.json').read_text())
        self.assertFalse(r['claim_flags']['empirical_joint_survivor'])
        self.assertFalse(r['claim_flags']['exclusion'])
        self.assertFalse(r['claim_flags']['genuine_identification_gain'])
        self.assertEqual(r['source_removal']['sensor_addition_gain_in_selected_band'],1.)
        self.assertEqual(r['sodium']['ensemble_bound_for_prior_SSB_variance'],1.)
        self.assertLess(r['sensor']['fraction_of_minimum_measured_psd_upper'],1.)
        self.assertLess(r['sodium']['maximum_visibility_change_over_SE_upper'],1.)

class ColdAtomKernelTests(unittest.TestCase):
    def test_free_flight_atoms_and_limits(self):
        from cal_response import free_flight_integral
        for T in [.017,.157,1.]:
            for f in [0.,1e-10,.137123456789,1.,7.81,1000.]:
                direct=quad(lambda t:(T-t)*np.cos(2*np.pi*f*t),0,T,epsabs=1e-12,limit=2000)[0]
                direct-=1j*quad(lambda t:(T-t)*np.sin(2*np.pi*f*t),0,T,epsabs=1e-12,limit=2000)[0]
                self.assertAlmostEqual(abs(free_flight_integral(f,T)-direct),0.,places=11)
            self.assertAlmostEqual(float(abs(free_flight_integral(0.,T))**2),T**4/4)
        # Different free-flight times are genuinely different temporal kernels.
        ratios=[abs(free_flight_integral(f,.157)/free_flight_integral(f,.017))**2 for f in [0.,40.]]
        self.assertGreater(max(ratios)/min(ratios),10.)

    def test_relative_force_and_doppler_injection(self):
        from cal_response import cloud_line_kernel,free_flight_integral,HBAR_OVER_M0
        k=np.array([1e6,0,0]);r=np.array([[0,0,0],[np.pi/1e6,0,0]])
        weights=np.array([.3,.7]);v=np.array([2e-5,0,0]);T=.157;f=1.7
        result=cloud_line_kernel(k,f,v,r,weights,T,[1,0,0])
        self.assertAlmostEqual(result['observed_frequency_hz'],f-20/(2*np.pi))
        self.assertAlmostEqual(result['width']/result['single_particle'],.84)
        self.assertAlmostEqual((result['width']+result['centroid'])/result['single_particle'],1.)
        common=cloud_line_kernel(k,f,v,np.zeros((2,3)),weights,T,[1,0,0])
        self.assertEqual(common['width'],0.)
        # Exact Gaussian two-quadrature injection, including centroid refitting.
        from numpy.polynomial.hermite import hermgauss
        nodes,w=hermgauss(4);ensemble=0.
        amplitude=1j*HBAR_OVER_M0*k[0]*free_flight_integral(result['observed_frequency_hz'],T)*np.exp(1j*(r@k))
        for i,x in enumerate(nodes):
            for j,y in enumerate(nodes):
                displacement=np.sqrt(2)*(x*amplitude.real+y*amplitude.imag)
                displacement+=123. # Ordinary common acceleration/pointing nuisance.
                centered=displacement-np.dot(weights,displacement)
                ensemble+=w[i]*w[j]/np.pi*np.dot(weights,centered**2)
        self.assertAlmostEqual(ensemble/result['width'],1.,places=9)

if __name__=='__main__':unittest.main()
