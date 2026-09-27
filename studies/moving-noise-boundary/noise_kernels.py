"""Finite-observation responses and conservative, non-Markov quantum bounds.

Frequencies are Hz; measures carry variance in s^-2 for unit mass m0.
See theory.md for quantifiers and normalization. No fitting or source access here.
"""
import numpy as np
from scipy.integrate import quad


def quantum_tv_bound(variance, mass_u, duration, tail=0.):
    """Half diamond-norm bound on the ensemble channel, including a bad tail.

    Valid for independent input, Gaussian Fourier field and deterministic controls.
    The min(1, ...) is essential for heavy tails / large perturbations.
    """
    if variance < 0 or mass_u < 0 or duration < 0 or not 0 <= tail <= 1:
        raise ValueError('Nonnegative physical inputs and tail in [0,1] required')
    x = 2 * variance * mass_u**2 * duration**2
    core = 1. if x >= np.log(3.) else .5 * np.expm1(x)
    return float(tail + (1-tail)*core)


def conditional_tv_bound(epsilon, survival_lower):
    if epsilon < 0 or not 0 < survival_lower <= 1:
        raise ValueError('Invalid distance or survival bound')
    return 1. if epsilon >= survival_lower else min(1., 2*epsilon/(survival_lower-epsilon))


def harmonic_visibility_bound(epsilon, mean_probability, visibility, row_l1):
    """OLS visibility perturbation from a uniform unconditional bin-probability bound.

    row_l1 are the absolute row sums of the pseudoinverse [1,cos,sin].
    This is a deterministic sensitivity bound, not a confidence interval.
    """
    d0 = epsilon * row_l1[0]
    da = epsilon * np.hypot(*row_l1[1:])
    if mean_probability <= d0:
        return float('inf')
    return float((da+abs(visibility)*d0)/(mean_probability-d0))


def sphere_spatial_difference(radius, density, separation, rc, m0=1.66053906660e-27):
    """G = g(0)-g(d), dimensionless; positive Gaussian overlap quadrature.

    For covariance sigma^2 exp(-|x-y|^2/(4 rc^2)), D=G sigma^2 T^2 at DC.
    """
    def g(d):
        def integrand(u):
            r = radius*u
            overlap = np.pi*(4*radius+r)*(2*radius-r)**2/12
            if d == 0:
                angle = 2*np.exp(-r*r/(4*rc*rc))
            elif r == 0:
                angle = 2*np.exp(-d*d/(4*rc*rc))
            else:
                angle = 2*rc*rc/(r*d)*(np.exp(-(r-d)**2/(4*rc*rc))-np.exp(-(r+d)**2/(4*rc*rc)))
            return 2*np.pi*r*r*overlap*angle*radius
        return density**2/m0**2*quad(integrand,0,2,epsabs=1e-60,epsrel=2e-11)[0]
    return g(0)-g(separation)


def held_kernel(frequency_hz, duration):
    """Time factor for two held arms, normalized to one at zero frequency."""
    return np.sinc(np.asarray(frequency_hz)*duration)**2


def band_held_lower(fmax, duration):
    """Analytic continuum lower bound via sin x >= x-x^3/6 (x <= 1)."""
    x = np.pi*fmax*duration
    if not 0 <= x <= 1:
        raise ValueError('This explicit lower bound requires pi*fmax*T <= 1')
    return float((1-x*x/6)**2)


def geometric_sum(cycles, n):
    x = (np.asarray(cycles)+.5) % 1 - .5
    return n*np.sinc(n*x)/np.sinc(x)*np.exp(1j*np.pi*(n-1)*x)


def blackman_transform(frequency, n, dt, symmetric=False):
    divisor = n-1 if symmetric else n
    coefficients = {0:.42, 1:-.25, -1:-.25, 2:.04, -2:.04}
    return sum(c*geometric_sum(np.asarray(frequency)*dt+j/divisor,n)
               for j,c in coefficients.items())


def blackman_norm(n, symmetric=False):
    if n < 8:
        raise ValueError('Use n >= 8')
    return .3046*(n-1 if symmetric else n)


def line_periodogram_gain(observed_hz, line_hz, n, dt, symmetric=False):
    """Expected one-sided PSD / variance for real stationary line, in seconds.

    No grid approximation: includes both sidebands and all sample aliases.
    """
    a = blackman_transform(np.asarray(observed_hz)-line_hz,n,dt,symmetric)
    b = blackman_transform(np.asarray(observed_hz)+line_hz,n,dt,symmetric)
    return dt*(abs(a)**2+abs(b)**2)/blackman_norm(n,symmetric)


def separated_band_gain_upper(fmin, fmax, bandmax, n, dt):
    """Uniform bound for both Blackman conventions and every atom in [0,bandmax].

    Deliberately loose triangle bound, with no cancellation or grid assumption.
    Require every shifted sideband below Nyquist, so sin is monotone.
    """
    shift = 2/((n-1)*dt)
    delta = fmin-bandmax-shift
    if delta <= 0 or (fmax+bandmax+shift)*dt >= .5:
        raise ValueError('Separated sub-Nyquist bands required')
    return float(2*dt/(.3046*(n-1)*np.sin(np.pi*delta*dt)**2))


def mechanical_gain(f, f0, q):
    """|chi(f)/chi(0)|^2; dimensionless."""
    return f0**4/((f0*f0-np.asarray(f)**2)**2+(np.asarray(f)*f0/q)**2)


def response_covariance(times, positions, k_vectors, omega, weights):
    """Finite-time covariance of a real field, including arbitrary trajectories.

    C_ij=sum weight cos(k.(x_i-x_j)-omega*(t_i-t_j)). Each row represents
    the symmetric +/- pair of a nonnegative spectral measure. Not a PSD slice.
    """
    t,x,k,w,a = map(np.asarray,(times,positions,k_vectors,omega,weights))
    if np.any(a < 0) or not np.isfinite(a).all():
        raise ValueError('Nonnegative finite spectral weights required')
    phase = x@k.T-t[:,None]*w
    return (np.cos(phase)*a)@np.cos(phase).T+(np.sin(phase)*a)@np.sin(phase).T


def scalar_phase_exponent(times, path_a, path_b, k_vectors, omega, weights, quadrature):
    """Exact Gaussian prescribed-path exponent for discrete quadrature, not optics."""
    t = np.asarray(times)
    a = np.exp(1j*(np.asarray(path_a)@np.asarray(k_vectors).T-t[:,None]*omega))
    b = np.exp(1j*(np.asarray(path_b)@np.asarray(k_vectors).T-t[:,None]*omega))
    z = np.asarray(quadrature)@(a-b)
    return float(.5*np.dot(weights,abs(z)**2))


def sphere_difference_lower(radius, density, separation, rc, m0=1.66053906660e-27):
    """Continuum lower certificate by integrating only |k| <= 1/radius.

    mu(k)/M >= .9 and 1-cos(k.d) >= 11(k.d)^2/24 on this ball.
    No numerical spatial quadrature enters the certificate.
    """
    if min(radius,density,separation,rc) <= 0 or separation > radius:
        raise ValueError('Positive parameters and separation <= radius required')
    mass = 4*np.pi*radius**3*density/3
    return ((mass/m0)**2 * rc**3/np.pi**1.5 * np.exp(-(rc/radius)**2)
            * .9**2 * 11/24 * separation**2 * 4*np.pi/(15*radius**5))
