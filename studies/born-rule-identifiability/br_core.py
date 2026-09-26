"""Conditional optical models; no universal Born-rule test is enabled."""
import math
import numpy as np

SIGNS = np.array([-1, 1, 1, -1, 1, -1, -1, 1])
BITS = ((np.arange(8)[:, None] >> np.arange(3)) & 1).astype(float)


def statistics(y):
    y = np.asarray(y, float)
    if y.shape[-1] != 8 or not np.isfinite(y).all():
        raise ValueError('Require finite eight-setting intensities')
    p = y-y[..., :1]
    singles = p[..., [1, 2, 4]]
    if np.any(singles <= 0):
        raise ValueError('Peres requires positive background-subtracted singles')
    pair = np.stack([p[..., 3]-p[..., 1]-p[..., 2], p[..., 5]-p[..., 1]-p[..., 4],
                     p[..., 6]-p[..., 2]-p[..., 4]], axis=-1)
    den = np.abs(pair).sum(axis=-1)
    if np.any(den <= 1e-12):
        raise ValueError('Unstable Sorkin denominator')
    c = pair/(2*np.sqrt(singles[..., [0, 0, 1]]*singles[..., [1, 2, 2]]))
    return {'peres': (c*c).sum(axis=-1)-2*np.prod(c, axis=-1),
            'epsilon_V': y @ SIGNS, 'denominator_V': den, 'correlations_AB_AC_BC': c}


def validate_gram(g):
    g = np.asarray(g)
    if g.shape != (3, 3) or not np.isfinite(g).all() or not np.allclose(g, g.conj().T):
        raise ValueError('Require finite Hermitian 3 by 3 matrix')
    if np.linalg.eigvalsh(g).min() < -1e-12:
        raise ValueError('Not positive semidefinite')


def gram_from_lower(y):
    p = np.asarray(y, float)-y[0]
    g = np.diag(p[[1, 2, 4]])
    for i, j, s in [(0, 1, 3), (0, 2, 5), (1, 2, 6)]:
        g[i, j] = g[j, i] = (p[s]-g[i, i]-g[j, j])/2
    validate_gram(g)
    return g


def imaginary_limit(g):
    validate_gram(g)
    if np.iscomplexobj(g) or np.linalg.eigvalsh(g).min() <= 0:
        raise ValueError('Require real positive definite reference')
    return math.sqrt(float(np.linalg.det(g)/g[2, 2]))


def complete_gram(g, u):
    if not math.isfinite(u) or abs(u) > imaginary_limit(g)*(1+1e-12):
        raise ValueError('Imaginary coherence outside PSD interval')
    out = np.array(g, complex)
    out[0, 1] += 1j*u
    out[1, 0] -= 1j*u
    validate_gram(out)
    return out


def table(g, delta=0., theta=0., phase=0., background=0., leakage=0j, quadratic=0.):
    """Theta is additive all-open intensity in V, NOT a probability law.

    Context phase acts on A in the all-open setting only. Fixed coherent
    leakage and quadratic detector response are separate synthetic imperfections.
    """
    validate_gram(g)
    if not np.isfinite([delta, theta, phase, background, leakage, quadratic]).all():
        raise ValueError('Nonfinite model parameter')
    v = BITS.astype(complex)+(1-BITS)*leakage
    v[7, 0] *= np.exp(1j*(delta+phase))
    q = np.einsum('si,ij,sj->s', v.conj(), g, v).real
    q[7] += theta
    if q.min() < -1e-12:
        raise ValueError('Negative optical intensity')
    if np.min(1+2*quadratic*q) <= 0:
        raise ValueError('Nonmonotone detector response')
    return background+q+quadratic*q*q


def phase_shift(g, u, delta):
    return 2*float(g[0, 1]+g[0, 2])*(math.cos(delta)-1)+2*u*math.sin(delta)


def phase_envelope(g, radius):
    """Sharp image over |delta|<=radius, |u|<=U, for radius<=pi/2."""
    if not 0 <= radius <= math.pi/2:
        raise ValueError('Phase radius outside certified domain')
    u = imaginary_limit(g)
    r = float(g[0, 1]+g[0, 2])
    values = []
    for v in [-u, u]:
        stationary = math.atan2(v, r)
        angles = [0., radius]+[stationary+k*math.pi for k in [-1, 0, 1]
                               if 0 <= stationary+k*math.pi <= radius]
        values += [phase_shift(g, v, d) for d in angles]
    return min(values), max(values)


def cancellation_witness(g, delta, theta=0., error=0., **kwargs):
    y = table(g, delta, theta, **kwargs)
    yp = table(g, delta, theta, phase=math.pi+error, **kwargs)
    return (y[7]+yp[7])/2-y[1]-y[6]+y[0]


def phase_error_bound(g, radius):
    validate_gram(g)
    if not 0 <= radius <= math.pi:
        raise ValueError('Invalid phase error radius')
    bc = float((g[1, 1]+g[2, 2]+2*g[1, 2].real).real)
    return 2*math.sqrt(max(0., float(g[0, 0].real)*bc))*math.sin(radius/2)


def leakage_response_bound(g, amplitude_bound, quadratic_bound, max_power=1.):
    """Triangle envelope for W; not a calibration inferred from archive means."""
    validate_gram(g)
    if not 0 <= amplitude_bound <= 1 or quadratic_bound < 0 or max_power <= 0:
        raise ValueError('Invalid calibration')
    a = np.abs(g)*max_power
    ell = amplitude_bound
    leak = 4*ell*(a[0, 1]+a[0, 2])+ell**2*(a[1:, 1:].sum()+a[0, 0]+a.sum())
    return float(leak+4*quadratic_bound*a.sum()**2)


def budget(gap, voltage_range=1., alpha=.01, power=.9):
    if not math.isfinite(gap) or voltage_range <= 0 or not 0 < alpha < 1 or not 0 < power < 1:
        raise ValueError('Invalid design')
    if gap <= 0:
        return None
    n = math.ceil(3.5*voltage_range**2*(math.sqrt(math.log(2/alpha))+
                  math.sqrt(math.log(1/(1-power))))**2/(2*gap**2))
    return {'per_setting': n, 'total': 5*n}


def hac_covariance(values, ids, lag):
    """Descriptive Bartlett covariance of sample mean, retaining original gaps."""
    x = np.asarray(values, float)
    if x.ndim == 1:
        x = x[:, None]
    ids = np.asarray(ids)
    if len(x) < 2 or len(ids) != len(x) or np.any(np.diff(ids) <= 0) or lag < 0:
        raise ValueError('Invalid ordered cycle series')
    z = x-x.mean(0)
    out = z.T @ z
    byid = {int(c): i for i, c in enumerate(ids)}
    for k in range(1, lag+1):
        pairs = [(i, byid[int(c)+k]) for i, c in enumerate(ids) if int(c)+k in byid]
        if pairs:
            a, b = np.array(pairs).T
            cross = z[a].T @ z[b]
            out += (1-k/(lag+1))*(cross+cross.T)
    return out/(len(z)*(len(z)-1))
