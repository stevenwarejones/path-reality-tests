"""Free particle on a circle: finite alternatives and conservative design bounds.

SI or consistent dimensionless units; exact spectral time evolution (no timestep).
Every readout has plus, minus, and failure outcomes. Synthetic predictions only.
"""
from dataclasses import dataclass
from math import ceil, cos, floor, isfinite, log, pi, sin, sqrt
import numpy as np


@dataclass(frozen=True)
class Ring:
    length: float = 2 * pi
    mass: float = 1.0
    hbar: float = 1.0

    def __post_init__(self):
        if not all(isfinite(x) and x > 0 for x in (self.length, self.mass, self.hbar)):
            raise ValueError('positive finite physical scales required')

    def k(self, j):
        return 2 * pi * np.asarray(j) / self.length

    def continuum(self, j):
        return self.hbar**2 * self.k(j)**2 / (2 * self.mass)

    def lattice(self, j, sites):
        check_sites(sites)
        a = self.length / sites
        # Stable form of 1-cos(x), including very fine grids.
        return 2 * self.hbar**2 / (self.mass * a*a) * np.sin(self.k(j)*a/2)**2

    def kinetic(self, sites):
        check_sites(sites)
        a = self.length / sites
        shift = np.roll(np.eye(sites), 1, axis=0)
        # Add both directed shifts even when N=2 and they coincide.
        return self.hbar**2/(2*self.mass*a*a) * (2*np.eye(sites)-shift-shift.T)


def check_sites(sites):
    if isinstance(sites, bool) or not isinstance(sites, (int, np.integer)) or sites < 1:
        raise ValueError('sites must be a positive integer')


def check_modes(modes, sites):
    check_sites(sites)
    modes = np.asarray(modes)
    if modes.ndim != 1 or len(modes) == 0 or not np.all(np.isfinite(modes)):
        raise ValueError('nonempty one dimensional finite mode list required')
    if not np.all(modes == modes.astype(int)) or len(set(modes)) != len(modes):
        raise ValueError('distinct integer modes required')
    if 2 * np.max(np.abs(modes)) >= sites:
        raise ValueError('strict no-aliasing band N > 2 J required')
    return modes.astype(int)


def dft(sites):
    check_sites(sites)
    return np.exp(2j*pi*np.outer(np.arange(sites), np.arange(sites))/sites)/sqrt(sites)


def evolve(coefficients, energies, time, hbar=1.0):
    c, e = np.asarray(coefficients, complex), np.asarray(energies, float)
    if c.ndim != 1 or c.shape != e.shape or not np.all(np.isfinite(c)) or not np.all(np.isfinite(e)):
        raise ValueError('finite equal length vectors required')
    if not isfinite(time) or not isfinite(hbar) or hbar <= 0:
        raise ValueError('finite time and positive hbar required')
    if not np.isclose(np.vdot(c, c).real, 1, atol=1e-12, rtol=0):
        raise ValueError('normalized state required')
    return c * np.exp(-1j * time * e / hbar)


def project_normalize(coefficients, keep):
    c = np.asarray(coefficients, complex)
    if not np.isclose(np.vdot(c,c).real, 1, atol=1e-12, rtol=0):
        raise ValueError('normalized state required')
    keep = np.asarray(keep, bool)
    if c.shape != keep.shape:
        raise ValueError('mask shape mismatch')
    retained = float(np.vdot(c[keep],c[keep]).real)
    tail = float(np.vdot(c[~keep],c[~keep]).real)
    if retained == 0:
        return None, tail, 2.0
    out = np.where(keep, c/sqrt(retained), 0)
    # Rationalized form preserves a nonzero error when tail is below machine epsilon.
    # The tolerance above admits last-bit normalization drift; clamp only at endpoints.
    tau=min(1.0,max(0.0,tail))
    return out, tail, sqrt(2*tau/(1+sqrt(1-tau)))


def calibration_phase_bias(phase_radius):
    """One calibration interval's worst-case contrast bias, not center drift."""
    if not isfinite(phase_radius) or phase_radius < 0:
        raise ValueError('finite nonnegative calibration phase radius required')
    return min(2., phase_radius**2/2)


def approximation(ring, sites, cutoff, time, tail=0.0):
    """Separate the derived pure-state trace bound from the proved Lean norm bound.

    tv_bound uses sqrt(tail); formal_tv_bound uses the conservative 2*sqrt(tail).
    Both share the globally proved 1/24 retained-band coefficient.
    """
    check_modes(np.arange(-cutoff,cutoff+1), sites)
    if not 0 <= tail <= 1:
        raise ValueError('tail probability outside [0,1]')
    a, k = ring.length/sites, abs(float(ring.k(cutoff)))
    epsilon = abs(time)*ring.hbar*a*a*k**4/(24*ring.mass)
    # Pure-state trace distance for normalized projection is sqrt(tail).
    return {'amplitude_band_bound': epsilon, 'tv_bound': min(1.0,sqrt(tail)+epsilon),
            'formal_tv_bound': min(1.0,2*sqrt(tail)+epsilon),
            'tv_bound_basis': 'pure-state trace-distance derivation; not the formal tail theorem'}


def readout(phase, efficiency=1., visibility=1.):
    if not all(isfinite(x) for x in (phase,efficiency,visibility)):
        raise ValueError('finite parameters required')
    if not (0 <= efficiency <= 1 and 0 <= visibility <= 1):
        raise ValueError('efficiency and visibility in [0,1] required')
    contrast = visibility*cos(phase)
    return np.array([efficiency*(1+contrast)/2,efficiency*(1-contrast)/2,1-efficiency])


def phase_gap(ring, sites, j0, j1, time):
    check_modes([j0,j1],sites)
    delta = ring.lattice([j0,j1],sites)-ring.continuum([j0,j1])
    return float(time*(delta[1]-delta[0])/ring.hbar)


def tv(p,q):
    p,q=np.asarray(p,float),np.asarray(q,float)
    if p.shape != q.shape or np.any(p < 0) or np.any(q < 0):
        raise ValueError('equal shape nonnegative distributions required')
    if not np.isclose(p.sum(),1,atol=1e-12,rtol=0) or not np.isclose(q.sum(),1,atol=1e-12,rtol=0):
        raise ValueError('normalized distributions required')
    return float(np.abs(p-q).sum()/2)


def cosine_interval(lo,hi):
    if not isfinite(lo) or not isfinite(hi) or lo > hi:
        raise ValueError('finite ordered phase interval required')
    if hi-lo >= 2*pi:
        return (-1.,1.)
    values=[cos(lo),cos(hi)]
    if ceil(lo/(2*pi)) <= floor(hi/(2*pi)):
        values.append(1.)
    if ceil((lo-pi)/(2*pi)) <= floor((hi-pi)/(2*pi)):
        values.append(-1.)
    return min(values),max(values)


def outcome_box(phase, phase_radius, eta, eta_radius, contrast, contrast_radius, contamination=0.):
    """Analytic enclosure of a continuous nuisance box, not a grid optimizer.

    contrast is eta*visibility. Rectangular relaxation ignores eta/contrast
    correlations, enlarging the null. contamination is complete-outcome TV.
    """
    if min(phase_radius,eta_radius,contrast_radius,contamination) < 0:
        raise ValueError('nonnegative uncertainty radii required')
    if not 0 <= contrast <= eta <= 1 or contamination > 1:
        raise ValueError('0 <= contrast <= efficiency <= 1 required')
    el,eh=max(0.,eta-eta_radius),min(1.,eta+eta_radius)
    wl,wh=max(0.,contrast-contrast_radius),min(1.,contrast+contrast_radius)
    cl,ch=cosine_interval(phase-phase_radius,phase+phase_radius)
    products=[w*c for w in (wl,wh) for c in (cl,ch)]
    ql,qh=min(products),max(products)
    lower=np.array([(el+ql)/2,(el-qh)/2,1-eh])-contamination
    upper=np.array([(eh+qh)/2,(eh-ql)/2,1-el])+contamination
    return np.maximum(0,lower),np.minimum(1,upper)


def box_gap(left,right):
    """Certified coordinate gap. Zero alone does not establish joint overlap."""
    al,ah=left; bl,bh=right
    return float(max(0.,np.max(al-bh),np.max(bl-ah)))


def radius(n,coordinates,error):
    if isinstance(n,bool) or not isinstance(n,int) or n < 1 or coordinates < 1 or not 0 < error < 1:
        raise ValueError('positive integer count and valid multiplicity/error required')
    return sqrt(log(2*coordinates/error)/(2*n))


def required_trials(gap,coordinates=3,alpha=.025,beta=.1):
    """Strict sufficient per-stratum n for simultaneous fixed-count Hoeffding.

    Reject if any empirical coordinate is farther than r_alpha from its null
    interval. Guaranteed power >=1-beta when r_alpha+r_beta < gap.
    """
    if not isfinite(gap) or gap < 0 or gap > 1 or coordinates < 1:
        raise ValueError('invalid gap or multiplicity')
    if not 0 < alpha < 1 or not 0 < beta < 1:
        raise ValueError('errors must lie in (0,1)')
    if gap == 0:
        return None
    threshold=(sqrt(log(2*coordinates/alpha))+sqrt(log(2*coordinates/beta)))**2/(2*gap**2)
    n=max(1,floor(threshold)+1)
    # Correct last-bit rounding against the actual strict criterion.
    while radius(n,coordinates,alpha)+radius(n,coordinates,beta) >= gap:
        n+=1
    while n > 1 and radius(n-1,coordinates,alpha)+radius(n-1,coordinates,beta) < gap:
        n-=1
    return n


def calibration_trials(probability_radius,error=.025):
    if not 0 < probability_radius < 1 or not 0 < error < 1:
        raise ValueError('calibration radius/error must lie in (0,1)')
    # Three fixed strata: total detection, plus at phase zero and phase pi.
    n=ceil(log(6/error)/(2*probability_radius**2))
    return n,3*n
