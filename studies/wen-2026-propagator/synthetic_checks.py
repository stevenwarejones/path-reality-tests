#!/usr/bin/env python3
"""Independent synthetic calculations, NOT a fit or replication of observed data."""
import argparse
import json
import math
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent


def smape(a, b):
    denominator = np.abs(a) + np.abs(b)
    # Both zero contributes zero; no silent division of a nonzero difference by zero.
    return float(np.mean(np.divide(2 * np.abs(a - b), denominator,
                                  out=np.zeros_like(denominator, dtype=float), where=denominator != 0)))


def enumerate_amplitudes(kernels, start):
    m, n, n2 = kernels.shape
    assert n == n2
    coordinates = np.array(np.unravel_index(np.arange(n ** m), (n,) * m)).T
    previous = np.column_stack((np.full(n ** m, start), coordinates[:, :-1]))
    amp = np.prod(kernels[np.arange(m), coordinates, previous], axis=1)
    return coordinates, previous, amp


def check_composition():
    # Independent implementations: enumerate routes versus multiply transfer matrices.
    rng = np.random.default_rng(63017)
    k = rng.normal(size=(3, 3, 3)) + 1j * rng.normal(size=(3, 3, 3))
    paths, _, amp = enumerate_amplitudes(k, 1)
    route_sum = np.array([amp[paths[:, -1] == i].sum() for i in range(3)])
    transfer = (k[2] @ k[1] @ k[0])[:, 1]
    np.testing.assert_allclose(route_sum, transfer, rtol=1e-13, atol=1e-13)
    return float(np.max(np.abs(route_sum - transfer)))


def compare_snapshot(actual, expected, path="root"):
    """Compare declared numerical outputs with tolerances, never image bytes."""
    if isinstance(expected, dict):
        if not isinstance(actual, dict) or actual.keys() != expected.keys():
            raise ValueError(f"Snapshot keys differ at {path}")
        for key in expected:
            compare_snapshot(actual[key], expected[key], f"{path}.{key}")
    elif isinstance(expected, list):
        if not isinstance(actual, list) or len(actual) != len(expected):
            raise ValueError(f"Snapshot length differs at {path}")
        for index, (a, e) in enumerate(zip(actual, expected)):
            compare_snapshot(a, e, f"{path}[{index}]")
    elif isinstance(expected, float):
        if not isinstance(actual, (int, float)) or not math.isfinite(actual) or not math.isclose(
                actual, expected, rel_tol=1e-10, abs_tol=1e-12):
            raise ValueError(f"Snapshot number differs at {path}: {actual} vs {expected}")
    elif type(actual) is not type(expected) or actual != expected:
        raise ValueError(f"Snapshot value differs at {path}: {actual} vs {expected}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="Check the committed numerical snapshot")
    parser.add_argument("--output-dir", type=Path, help="Write regenerated outputs to this directory")
    args = parser.parse_args()
    n, m, wavelength, dz = 17, 5, 795e-9, 15e-3
    paths, previous, _ = enumerate_amplitudes(np.ones((m, n, n)), 8)
    differences = paths - previous
    action_index = np.sum(differences ** 2, axis=1)
    endpoint = paths[:, -1] - 8
    report = {"kind": "synthetic-only; no experimental workbooks loaded",
              "paths": n ** m, "paths_per_endpoint": n ** (m - 1),
              "full_kernel_entries": m * n ** 2,
              "used_kernel_entries_fixed_start": n + (m - 1) * n ** 2,
              "composition_max_abs_error": check_composition(), "coordinates": []}
    curves = []
    for dx in (5.73e-6, 48.72e-6 / 8):
        phase = np.pi * dx ** 2 / (wavelength * dz) * action_index
        relative = phase - np.pi * dx ** 2 / (wavelength * dz) * endpoint ** 2 / m
        assert np.min(relative) > -1e-12
        amp = np.exp(1j * phase)
        coherent = np.array([abs(amp[endpoint == i].sum()) ** 2 for i in range(-8, 9)])
        coherent /= coherent.sum()  # declared shape normalization, NOT absolute detection rate
        flat = np.ones(n) / n
        report["coordinates"].append({"dx_m": dx,
            "max_transverse_slope": float(16 * dx / dz),
            "max_relative_action_over_pi_hbar": float(relative.max() / np.pi),
            "fraction_relative_action_at_least_2pi": float(np.mean(relative >= 2 * np.pi)),
            "normalized_shape": coherent.tolist(),
            "tv_coherent_vs_uniform": float(np.abs(coherent - flat).sum() / 2)})
        curves.append(coherent)

    alpha_sd, beta_sd = 0.0518, 0.0459 * np.pi
    noise_results = []
    for seed in range(8):
        rng = np.random.default_rng(seed)
        alpha = rng.normal(0, alpha_sd, (m, n, n))
        beta = rng.normal(0, beta_sd, (m, n, n))
        # Every occurrence of a given kernel entry shares the SAME noisy estimate.
        magnitude = np.prod(1 + alpha[np.arange(m), paths, previous], axis=1)
        phase_error = np.sum(beta[np.arange(m), paths, previous], axis=1)
        probability = magnitude ** 2
        noise_results.append({"seed": seed, "path_probability_smape": smape(probability, probability.mean()),
                              "path_phase_error_sd_over_pi": float(np.std(phase_error) / np.pi)})
    report["independent_kernel_noise"] = noise_results
    report["iid_segment_phase_sd_over_pi"] = float(np.sqrt(m) * beta_sd / np.pi)
    # Exact counterexample: a common multiplicative calibration factor affects every path equally.
    report["common_mode_counterexample"] = {
        "probability_multiplier_for_alpha_005": 1.05 ** (2 * m),
        "within_run_probability_smape": 0,
        "explanation": "A shared gain error changes absolute scale while leaving within-run path uniformity exact."}
    report["coordinate_convention_shape_tv"] = float(np.abs(curves[0] - curves[1]).sum() / 2)
    if args.check:
        compare_snapshot(report, json.loads((ROOT / "results" / "synthetic.json").read_text()))
        print("Synthetic numerical snapshot and independent composition check passed")
        if args.output_dir is None:
            return
    out = args.output_dir if args.output_dir is not None else ROOT / "results"
    out.mkdir(parents=True, exist_ok=True)
    (out / "synthetic.json").write_text(json.dumps(report, indent=2) + "\n")
    fig, ax = plt.subplots(figsize=(8, 4.6), constrained_layout=True)
    for curve, label in zip(curves, ("5.73 μm centers", "6.09 μm centers")):
        ax.plot(range(-8, 9), curve, marker="o", label=label)
    ax.axhline(1 / n, color="gray", linestyle="--", label="Incoherent, normalized")
    ax.set(xlabel="Final coordinate index", ylabel="Probability after shape normalization",
           title="Synthetic finite-grid sum — not measured data")
    ax.legend()
    fig.savefig(out / "synthetic-shapes.png", dpi=180)
    plt.close(fig)
    print(json.dumps({"coordinates": report["coordinates"],
                      "noise": noise_results, "composition_error": report["composition_max_abs_error"]}, indent=2))


if __name__ == "__main__":
    main()
