"""Diagnostic free-cloud kernels, not a calibrated CAL inference.

Frequencies are cyclic Hz; k is rad/m. The covariance measure is the SAME
field variance measure as theory.md, with Fourier phase k.y - 2*pi*f*t.
The cloud paths here have common constant velocity and fixed offsets. A
lensed/interacting expanding CAL cloud requires its actual Green function.
"""
import math
import numpy as np

HBAR_OVER_M0 = 1.054571817e-34 / 1.66053906660e-27


def free_flight_integral(frequency_hz, duration_s):
    """Integral_0^T (T-t) exp(-2*pi*i*f*t) dt, including off-grid atoms."""
    if duration_s < 0:
        raise ValueError('Negative observation duration')
    x = 2 * np.pi * np.asarray(frequency_hz, dtype=float) * duration_s
    small = np.abs(x) < 0.1
    safe = np.where(small, 1., x)
    direct = (-np.expm1(-1j * safe) - 1j * safe) / safe**2
    series = sum((-1j * x)**n / math.factorial(n + 2) for n in range(12))
    return duration_s**2 * np.where(small, series, direct)


def cloud_line_kernel(k_vector, field_frequency_hz, velocity, offsets, weights,
                      duration_s, measurement_axis):
    """Width variance m² per field variance s^-2, to prescribed-path order.

Subtracting the cloud centroid is essential. Displacement and centroid
kernels are also returned to expose the common-force blind spot. No atoms
are declared independent merely because their equal-time separation is large.
"""
    k = np.asarray(k_vector, dtype=float)
    v = np.asarray(velocity, dtype=float)
    r = np.asarray(offsets, dtype=float)
    w = np.asarray(weights, dtype=float)
    axis = np.asarray(measurement_axis, dtype=float)
    if k.shape != (3,) or v.shape != (3,) or axis.shape != (3,):
        raise ValueError('Vectors must have three coordinates')
    if r.shape != (len(w), 3) or np.any(w < 0) or not np.isclose(w.sum(), 1.):
        raise ValueError('A normalized nonnegative cloud is required')
    if not np.isclose(np.linalg.norm(axis), 1.):
        raise ValueError('Measurement axis must be a unit vector')
    observed = field_frequency_hz - np.dot(k, v)/(2*np.pi)
    phase = np.exp(1j * (r @ k))
    form = np.dot(w, phase)
    # The centered sum is stable and nonnegative even for almost uniform force.
    relative = float(np.dot(w, np.abs(phase-form)**2))
    single = float(HBAR_OVER_M0**2 * np.dot(k, axis)**2 *
                   abs(free_flight_integral(observed, duration_s))**2)
    return {'width': single*relative, 'centroid': single*abs(form)**2,
            'single_particle': single, 'observed_frequency_hz': float(observed)}
