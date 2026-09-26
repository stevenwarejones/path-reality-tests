#!/usr/bin/env python3
"""Synthetic checks of intervention and window identities; no observed data input."""
import argparse
import json
import hashlib
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
PHASES = np.arange(4) * np.pi / 2


def probability(effect, state):
    value = np.trace(effect @ state)
    assert abs(value.imag) < 1e-12
    assert -1e-12 <= value.real <= 1 + 1e-12
    return float(value.real)


def unitary(rng, dimension):
    matrix = rng.normal(size=(dimension, dimension)) + 1j * rng.normal(size=(dimension, dimension))
    return np.linalg.qr(matrix)[0]


def phase_check(rho, effect, projection):
    complement = np.eye(len(rho)) - projection
    dephased = projection @ rho @ projection + complement @ rho @ complement
    baseline = probability(effect, dephased)
    cross = np.trace(effect @ projection @ rho @ complement)
    direct = []
    residuals = []
    for theta in PHASES:
        phase = complement + np.exp(1j * theta) * projection
        observed = probability(effect, phase @ rho @ phase.conj().T)
        predicted = baseline + 2 * (np.exp(1j * theta) * cross).real
        residuals.append(abs(observed - predicted))
        np.testing.assert_allclose(observed, predicted, atol=1e-12, rtol=0)
        np.testing.assert_allclose(probability(effect, phase @ dephased @ phase.conj().T),
                                   baseline, atol=1e-12, rtol=0)
        direct.append(observed)
    np.testing.assert_allclose([direct[0] - direct[2], direct[1] - direct[3]],
                               [4 * cross.real, -4 * cross.imag], atol=1e-12, rtol=0)
    return direct, baseline, max(residuals)


def analyze():
    rng = np.random.default_rng(3212026)
    maximum_error = 0.0
    for dimension in (2, 3, 7):
        for _ in range(16):
            raw = rng.normal(size=(dimension, dimension)) + 1j * rng.normal(size=(dimension, dimension))
            rho = raw @ raw.conj().T
            rho /= np.trace(rho)
            basis = unitary(rng, dimension)
            effect = basis @ np.diag(rng.uniform(0, 1, dimension)) @ basis.conj().T
            # Noncommuting mixed states and effects catch sign/conjugation mistakes.
            projection = np.diag([1] * (dimension // 2) + [0] * (dimension - dimension // 2))
            _, _, error = phase_check(rho, effect, projection)
            maximum_error = max(maximum_error, error)

    h = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)
    initial = np.array([1, 0], dtype=complex)
    region = np.diag([0, 1])
    psi = h @ initial
    rho = np.outer(psi, psi.conj())
    eta = 0.8
    effect_zero = h.conj().T @ np.diag([eta, 0]) @ h
    pzero, null, _ = phase_check(rho, effect_zero, region)
    effect_one = h.conj().T @ np.diag([0, eta]) @ h
    pone, _, _ = phase_check(rho, effect_one, region)
    failure = np.full(4, 1 - eta)
    np.testing.assert_allclose(np.array(pzero) + np.array(pone) + failure, 1, atol=1e-12)

    # A coherent state can have zero contrast for a detector insensitive to phase.
    insensitive, _, _ = phase_check(rho, eta * np.eye(2), region)
    np.testing.assert_allclose(insensitive, eta, atol=1e-12)

    # Final containment does not certify harmless intermediate truncation.
    keep_zero = np.diag([1, 0])
    free = h @ h @ initial
    projected = h @ keep_zero @ h @ initial
    survival = float(np.vdot(projected, projected).real)
    np.testing.assert_allclose(abs(free) ** 2, [1, 0], atol=1e-12)
    np.testing.assert_allclose(abs(projected) ** 2, [.25, .25], atol=1e-12)
    np.testing.assert_allclose(abs(projected[0]) ** 2 / survival, .5, atol=1e-12)

    # Independent random checks of the multi-plane telescoping error bound.
    max_ratio = 0.0
    for _ in range(32):
        dimension = 5
        free_state = unitary(rng, dimension)[:, 0]
        projected_state = free_state.copy()
        error_bound = 0.0
        for step in range(4):
            propagation = unitary(rng, dimension)
            free_state = propagation @ free_state
            projected_state = propagation @ projected_state
            if step < 3:
                projection = np.diag([1, 1, 1, 0, 0])
                error_bound += np.linalg.norm(free_state - projection @ free_state)
                projected_state = projection @ projected_state
        error = np.linalg.norm(free_state - projected_state)
        assert error <= error_bound + 1e-12
        max_ratio = max(max_ratio, float(error / error_bound))

    alpha, n, target_width = .01, 10000, .02
    contrast_radius = math.sqrt(2 * math.log(8 / alpha) / n)
    report = {
        "kind": "synthetic algebra/precision diagnostics; no empirical fit or power claim",
        "source": {"filename": "intervention_checks.py",
                   "sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
        "seed": 3212026,
        "mixed_state_cases": 48,
        "phase_identity_max_error": maximum_error,
        "lossy_balanced_example": {"efficiency": eta, "selected_bin": pzero,
                                  "other_bin": pone, "failure": failure.tolist(),
                                  "dephased_selected_bin": null},
        "coherent_zero_contrast_counterexample": insensitive,
        "intermediate_window_counterexample": {
            "free_selected_bin": float(abs(free[0]) ** 2),
            "projected_selected_bin_unconditioned": float(abs(projected[0]) ** 2),
            "projected_survival": survival,
            "projected_selected_bin_conditioned": float(abs(projected[0]) ** 2 / survival)},
        "window_bound_random_cases": 32,
        "largest_norm_error_over_bound": max_ratio,
        "precision_design": {"family_error": alpha, "trials_per_setting": n,
                             "contrast_radius": contrast_radius,
                             "target_contrast_radius": target_width,
                             "required_trials_per_setting": math.ceil(2 * math.log(8 / alpha) / target_width ** 2)},
        "limitations": ["Numerical checks are not formal proofs.",
                        "No calibration/nuisance bound has been measured.",
                        "Hoeffding calculation assumes independent fixed-probability Bernoulli trials."]
    }
    return report


def check_snapshot(expected, actual, path='result'):
    """Match structure exactly and numerical values within roundoff tolerance."""
    if isinstance(expected, dict):
        if not isinstance(actual, dict) or expected.keys() != actual.keys():
            raise ValueError(f'{path}: keys differ')
        for key in expected:
            check_snapshot(expected[key], actual[key], f'{path}.{key}')
    elif isinstance(expected, list):
        if not isinstance(actual, list) or len(expected) != len(actual):
            raise ValueError(f'{path}: list lengths differ')
        for i, (left, right) in enumerate(zip(expected, actual)):
            check_snapshot(left, right, f'{path}[{i}]')
    elif isinstance(expected, float):
        if isinstance(actual, bool) or not isinstance(actual, (int, float)) or not math.isfinite(actual):
            raise ValueError(f'{path}: invalid number')
        if not math.isclose(expected, actual, rel_tol=1e-10, abs_tol=1e-12):
            raise ValueError(f'{path}: number differs')
    elif type(expected) is not type(actual) or expected != actual:
        raise ValueError(f'{path}: value differs')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true', help='Validate against reviewed results without writing')
    args = parser.parse_args()
    report = analyze()
    if args.check:
        check_snapshot(report, json.loads((ROOT / 'results/intervention.json').read_text()))
        print('Intervention algebra, full loss outcomes and window counterexamples passed')
    else:
        (ROOT / "results" / "intervention.json").write_text(json.dumps(report, indent=2) + "\n")
        print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
