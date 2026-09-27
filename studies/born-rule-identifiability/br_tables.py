"""Physical table models, constructive upper bounds and conservative spectral lower bounds."""
import numpy as np
from functools import lru_cache
from scipy.optimize import least_squares
from scipy.special import xlogy


def rotation(angle):
    z = np.diag(np.exp(np.array([-1, 1])*1j*angle/2))
    return z.conj().T @ np.array([[1, -1j], [-1j, 1]]) @ z / np.sqrt(2)


def bloch(v):
    return np.array([2*(v[0].conj()*v[1]).real, 2*(v[0].conj()*v[1]).imag,
                     abs(v[0])**2-abs(v[1])**2])


@lru_cache(maxsize=2)
def directions(kind):
    """Independently reconstruct S_theta=Z_theta^dagger SX Z_theta in execution order."""
    if kind == 'viviani':
        a = np.zeros(5)
        b = np.pi*np.array([.25, -.25, .75, -.75, 0])
        t, f = np.zeros(4), -b[:4]
    elif kind == 'scan':
        eta, p = 1.23095942, np.arange(5)*2*np.pi/5
        a = np.r_[0, eta-np.pi, eta+np.pi+2*np.pi/3, eta+np.pi-2*np.pi/3, p]
        b = np.r_[0, 0, 2*np.pi/3, -2*np.pi/3, p+np.pi/2]
        t = np.r_[np.pi, np.pi/2, np.pi/2+2*np.pi/3, np.pi/2-2*np.pi/3]
        f = np.r_[0, np.pi, np.pi+2*np.pi/3, np.pi-2*np.pi/3]
    else:
        raise ValueError('Unknown geometry')
    zero = np.array([1, 0])
    r = np.array([bloch(rotation(y) @ rotation(x) @ zero) for x, y in zip(a, b)])
    m = np.array([bloch((rotation(x) @ rotation(y)).conj().T @ zero) for x, y in zip(t, f)])
    return r, m


def h(x, theta):
    if abs(theta) > 1:
        raise ValueError('Outside normalized monotone theory domain')
    return x + theta*x*(1-x*x)/2


def frame(v):
    # Deterministic tangent basis, including poles (choose least-aligned coordinate).
    e = np.eye(3)[np.argmin(abs(v), axis=1)]
    a = np.cross(v, e)
    a /= np.linalg.norm(a, axis=1)[:, None]
    return np.stack([a, np.cross(v, a)], axis=-1)


def moved(v, coordinates):
    w = v + np.einsum('nki,ni->nk', frame(v), coordinates)
    return w / np.linalg.norm(w, axis=1)[:, None]


def model(parameters, kind, theta=0, leakage=False):
    r, m = directions(kind)
    n = len(r)
    v = np.asarray(parameters)
    expected = 4*n+20 if leakage else 3*n+16
    if v.shape != (expected,):
        raise ValueError('Parameter dimension')
    rr, mm = moved(r, v[:2*n].reshape(n, 2)), moved(m, v[2*n:2*n+8].reshape(4, 2))
    c = v[2*n+8:3*n+8]
    offset = 4*n+8 if leakage else 3*n+8
    lo, hi = v[offset:offset+4], v[offset+4:offset+8]
    # Mixture of antipodal PURE preparations, weights (1+c)/2 and (1-c)/2.
    # Never f_theta of a mixed-state Born probability.
    q = lo + (hi-lo)*(1+c[:, None]*h(rr @ mm.T, theta))/2
    lam = v[3*n+8:4*n+8] if leakage else np.zeros(n)
    leak_response = v[-4:] if leakage else np.zeros(4)
    p = (1-lam[:, None])*q+lam[:, None]*leak_response
    return p, {'preparation_axes': rr, 'measurement_axes': mm, 'contrast': c,
               'readout_low': lo, 'readout_high': hi, 'leakage': lam,
               'leak_response': leak_response, 'qubit_table': q}


def fit_table(p, kind, theta=0, tangent_box=.15, leakage_cap=None, starts=2):
    """Construct feasible models. Nonconvex fits are UPPER bounds, never global exclusions."""
    p = np.asarray(p)
    n = len(p)
    leakage = leakage_cap is not None
    low = np.r_[np.full(2*n+8, -tangent_box), np.zeros(n)]
    high = np.r_[np.full(2*n+8, tangent_box), np.ones(n)]
    if leakage:
        low, high = np.r_[low, np.zeros(n)], np.r_[high, np.full(n, leakage_cap)]
    low, high = np.r_[low, np.zeros(12 if leakage else 8)], np.r_[high, np.ones(12 if leakage else 8)]
    best = None
    for seed in range(starts):
        rng = np.random.default_rng(730+seed)
        v = np.r_[rng.uniform(-min(.01, tangent_box/2), min(.01, tangent_box/2), 2*n+8), np.full(n, .99)]
        if leakage:
            v = np.r_[v, rng.uniform(.1, .9, n)*leakage_cap]
        v = np.r_[v, np.full(4, .01), np.full(4, .99)]
        if leakage:
            v = np.r_[v, rng.uniform(.1, .9, 4)]
        result = least_squares(lambda v: (model(v, kind, theta, leakage)[0]-p).ravel(),
                               v, bounds=(low, high), max_nfev=1200,
                               gtol=1e-12, ftol=1e-12, xtol=1e-12)
        if best is None or np.linalg.norm(result.fun) < np.linalg.norm(best.fun):
            best = result
    prediction, physical = model(best.x, kind, theta, leakage)
    r, m = directions(kind)
    angles = np.arccos(np.clip(np.r_[np.sum(r*physical['preparation_axes'], axis=1),
                                     np.sum(m*physical['measurement_axes'], axis=1)], -1, 1))
    return {'theta': theta, 'tangent_box': tangent_box, 'leakage_cap': leakage_cap,
            'parameters': best.x.tolist(), 'prediction': prediction.tolist(),
            'max_abs_residual': float(abs(prediction-p).max()),
            'frobenius_residual': float(np.linalg.norm(prediction-p)),
            'max_axis_angle_rad': float(angles.max()),
            'max_leakage': float(physical['leakage'].max()),
            'optimizer_success': bool(best.success), 'evaluations': int(best.nfev),
            'status': 'constructive numerical upper bound; global optimality not certified'}


def validate_fit(fit, kind, p, tolerance=2e-8):
    v = np.array(fit['parameters'])
    n = len(p)
    if np.any(abs(v[:2*n+8]) > fit['tangent_box']+tolerance):
        raise ValueError('Axes outside declared box')
    if np.any(v[2*n+8:] < -tolerance) or np.any(v[2*n+8:] > 1+tolerance):
        raise ValueError('Nonphysical probability/contrast')
    prediction, physical = model(v, kind, fit['theta'], fit['leakage_cap'] is not None)
    if fit['leakage_cap'] is not None and physical['leakage'].max() > fit['leakage_cap']+tolerance:
        raise ValueError('Leakage exceeds cap')
    if abs(prediction-np.array(fit['prediction'])).max() > tolerance:
        raise ValueError('Prediction is not the stated physical model')
    if abs(np.linalg.norm(prediction-p)-fit['frobenius_residual']) > tolerance:
        raise ValueError('Residual certificate mismatch')
    if abs(abs(prediction-p).max()-fit['max_abs_residual']) > tolerance:
        raise ValueError('Supremum residual certificate mismatch')
    return prediction


def centered_singular(p):
    return np.linalg.svd(p-np.mean(p, axis=0), compute_uv=False)


def noise_radius(shots, alpha):
    """Independent bounded Bernoulli trials, possibly unequal means; joint sub-Gaussian norm bound."""
    v = 1/(4*np.asarray(shots, dtype=float))
    t = np.log(1/alpha)
    return float(np.sqrt(v.sum()+2*np.sqrt(np.square(v).sum()*t)+2*v.max()*t))


def spectral_boundary(counts, alpha=.025):
    a = np.asarray(counts)
    summed = a.sum(axis=0)
    shots = summed.sum(axis=2)
    p = summed[:, :, 0]/shots
    singular = centered_singular(p)
    tail = float(np.linalg.norm(singular[3:]))
    noise = noise_radius(shots, alpha)
    jobs = len(a)
    # Allows arbitrary within-job and cross-setting dependence. Independent jobs only.
    job_radius = float(np.sqrt(p.size*np.log(2*p.size/alpha)/(2*jobs)))
    byjob = a[:, :, :, 0]/a.sum(axis=3)
    return {'centered_singular_values': singular.tolist(),
            'mean_table_distance_lower': tail/np.sqrt(p.size),
            'iid_shot_noise_frobenius_radius': noise,
            'iid_shot_distance_lower': max(0, tail-noise)/np.sqrt(p.size),
            'independent_job_noise_frobenius_radius': job_radius,
            'independent_job_distance_lower': max(0, tail-job_radius)/np.sqrt(p.size),
            'alpha': alpha, 'job_se_max': float((byjob.std(axis=0, ddof=1)/np.sqrt(jobs)).max()),
            'within_job_tail_norms': [float(np.linalg.norm(centered_singular(q)[3:])) for q in byjob],
            'job_probability_range_max': float(np.ptp(byjob, axis=0).max())}


def deviance(counts, prediction):
    a = np.asarray(counts).sum(axis=0)
    n = a.sum(axis=2)
    p = np.clip(prediction, 1e-12, 1-1e-12)
    return float(2*np.sum(xlogy(a[:, :, 0], a[:, :, 0]/(n*p)) +
                         xlogy(a[:, :, 1], a[:, :, 1]/(n*(1-p)))))


def probabilities(counts):
    s = np.asarray(counts).sum(axis=0)
    return s[:, :, 0]/s.sum(axis=2)


def contextual_response(p, q):
    """Minimum reset-to-0/1 probability per cell for a specified base Q (not minimized over Q)."""
    return np.where(p >= q, (p-q)/np.maximum(1-q, 1e-15), (q-p)/np.maximum(q, 1e-15))
