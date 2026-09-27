"""Observation identities; no universal cross-chip coefficients are assumed."""
import math


def odd_probability(intensity):
    if not math.isfinite(intensity) or intensity < 0:
        raise ValueError("integrated tunneling intensity must be finite and nonnegative")
    return -math.expm1(-2 * intensity) / 2


def any_probability(intensity):
    if not math.isfinite(intensity) or intensity < 0:
        raise ValueError("integrated tunneling intensity must be finite and nonnegative")
    return -math.expm1(-intensity)


def any_bounds(contrast, max_intensity=None):
    """Sharp mixed-Poisson bounds; unbounded-support upper end is a supremum."""
    if not math.isfinite(contrast) or not 0 <= contrast <= 1:
        raise ValueError("contrast must be in [0,1]")
    low = contrast / (1 + math.sqrt(1 - contrast))
    if max_intensity is None:
        return low, contrast
    if not math.isfinite(max_intensity) or max_intensity < 0:
        raise ValueError("invalid intensity cap")
    if contrast == 1 or contrast > -math.expm1(-2 * max_intensity) + 1e-14:
        raise ValueError("intensity cap incompatible with contrast")
    return low, contrast / (1 + math.exp(-max_intensity))


def two_point_witness(contrast, intensity):
    """Mix Lambda=0 and Lambda=intensity, reproducing the given parity contrast."""
    if not 0 < contrast < 1 or intensity <= 0:
        raise ValueError("interior contrast and positive intensity required")
    weight = contrast / (2 * odd_probability(intensity))
    if weight > 1 + 1e-14:
        raise ValueError("intensity too small")
    weight = min(weight, 1.0)  # Remove a possible roundoff-only negative weight.
    return {"intensities": [0.0, intensity], "weights": [1 - weight, weight]}


def corrected_contrast(observed_odd, background_odd):
    """Independent XOR background correction, before readout errors."""
    if not 0 <= observed_odd <= 0.5 or not 0 <= background_odd < 0.5:
        raise ValueError("invalid probability")
    result = (2 * observed_odd - 2 * background_odd) / (1 - 2 * background_odd)
    if result < 0:
        raise ValueError("negative signal contrast")
    return result
