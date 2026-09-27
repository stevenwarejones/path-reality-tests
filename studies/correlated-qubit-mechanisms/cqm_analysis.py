#!/usr/bin/env python3
"""Reproduce source summaries and deterministic feasibility diagnostics."""
import argparse
import csv
from fractions import Fraction
import io
import itertools
import json
from pathlib import Path
import pickle
import zipfile

import numpy as np
from scipy.optimize import linprog

from cqm_sources import acquire
from cqm_response import any_bounds, any_probability, odd_probability, two_point_witness

HERE = Path(__file__).resolve().parent


def rows(z, name):
    return list(csv.reader(io.StringIO(z.read(name).decode())))[1:]


def affine_certificate(data, envelope=4):
    """Exact annihilating contrasts, plus independently solved minimax LP.

    Every reported rational certificate annihilates 1 and x exactly. It is
    a global lower bound over ALL affine functions, not an optimizer failure.
    The LP supplies an attainable upper bound with nonnegative slope/intercept.
    """
    rational = [[Fraction(v) for v in row] for row in data]
    best_k = (Fraction(0), None, None)
    best_delta = (Fraction(0), None, None)
    for size in (2, 3):
        for idx in itertools.combinations(range(len(data)), size):
            r = [rational[i] for i in idx]
            if size == 2:
                if r[0][0] != r[1][0]:
                    continue
                w = [Fraction(1), Fraction(-1)]
            else:
                x = [a[0] for a in r]
                w = [x[1] - x[2], x[2] - x[0], x[0] - x[1]]
                if not any(w):
                    continue
            assert sum(w) == 0 and sum(wi * ri[0] for wi, ri in zip(w, r)) == 0
            residual = abs(sum(wi * ri[1] for wi, ri in zip(w, r)))
            error = sum(abs(wi) * ri[2] for wi, ri in zip(w, r))
            k = residual / error
            delta = (residual - envelope * error) / sum(map(abs, w))
            if k > best_k[0]:
                best_k = k, idx, w
            if delta > best_delta[0]:
                best_delta = delta, idx, w
    a = np.asarray(data, dtype=float)
    x, y, error = a.T
    if np.any(error <= 0) or not np.isfinite(a).all():
        raise ValueError("invalid rate table")
    design = np.column_stack([np.ones(len(x)), x])
    fit = np.linalg.lstsq(design / error[:, None], y / error, rcond=None)[0]
    lp_k = linprog([0, 0, 1], A_ub=np.vstack([
        np.column_stack([design, -error]), np.column_stack([-design, -error])]),
        b_ub=np.r_[y, -y], bounds=[(0, None)] * 3, method="highs")
    lp_delta = linprog([0, 0, 1], A_ub=np.vstack([
        np.column_stack([design, -np.ones(len(x))]),
        np.column_stack([-design, -np.ones(len(x))])]),
        b_ub=np.r_[y + envelope * error, -y + envelope * error],
        bounds=[(0, None)] * 3, method="highs")
    if not lp_k.success or not lp_delta.success:
        raise ValueError("minimax computation failed")
    # A directly checked feasible point gives an upper bound, irrespective
    # of the solver's claim of optimality. Round it outward by 1e-9.
    k_upper = max(0.0, float(np.max(abs(design @ lp_k.x[:2] - y) / error))) + 1e-9
    delta_upper = max(0.0, float(np.max(abs(design @ lp_delta.x[:2] - y) - envelope * error))) + 1e-9
    if k_upper < float(best_k[0]) or delta_upper < float(best_delta[0]):
        raise ValueError("primal/certificate inconsistency")

    def certificate(best):
        value, idx, weights = best
        return {"lower_bound": float(value), "exact_lower_bound": str(value),
                "row_indices_0based": list(idx) if idx else [],
                "exact_weights": [str(w) for w in weights] if weights else []}

    return {"rows": len(data), "dose_range_m_minus2": [float(x.min()), float(x.max())],
            "weighted_affine_intercept_s_minus1": float(fit[0]),
            "weighted_affine_slope_m2_s_minus1": float(fit[1]),
            "required_error_multiplier": certificate(best_k),
            "feasible_multiplier_upper": k_upper,
            "additive_discrepancy_s_minus1": certificate(best_delta),
            "feasible_discrepancy_upper_s_minus1": delta_upper,
            "feasible_discrepancy_affine_parameters": lp_delta.x[:2].tolist(),
            "fits_four_reported_errors": delta_upper < 1e-8}


def gamma_analysis(path):
    curves = {}
    footprint = []
    with zipfile.ZipFile(path) as z:
        for chip in ("nonCu", "1umCu"):
            q = 2 if chip == "nonCu" else 5
            curves[chip + "_charge"] = affine_certificate(rows(z, f"Figure02/Figure2b_{chip}_Q{q}.csv"))
            table = rows(z, f"Figure03/Figure3b_{chip}.csv")
            for i, q in enumerate([3, 4, 6] if chip == "nonCu" else range(1, 7)):
                selected = [[r[0], r[1 + 2*i], r[2 + 2*i]] for r in table]
                curves[chip + f"_parity_Q{q}"] = affine_certificate(selected)
            for i, r in enumerate(rows(z, f"Figure05/Figure5c_{chip}.csv")):
                distance_native, contrast, error = map(float, r)
                # Errors are preserved only as a diagnostic envelope, not CI.
                lo = max(0, contrast - 4 * error)
                hi = min(1, contrast + 4 * error)
                if not 0 <= contrast <= 1 or lo > hi:
                    raise ValueError("invalid footprint contrast")
                footprint.append({"chip": chip, "row": i,
                    "distance_native_header_says_mm": distance_native,
                    "distance_mm_interpreted_as_native_um": distance_native / 1000,
                    "unit_resolution": "CSV header says mm but values are consistent with micrometers: 2020 maps to 2.020 mm on the stated 8 mm chip. Inference from paper geometry, not an author-verified correction; unused by probability bounds.",
                    "reported_poisoning_contrast": contrast, "reported_error": error,
                    "any_tunnel_bounds_at_point_estimate": list(any_bounds(contrast)),
                    "any_tunnel_bounds_four_error_envelope": [any_bounds(lo)[0], hi]})
    return {"curves": curves, "footprint": footprint,
            "total_rate_points": sum(v["rows"] for v in curves.values()),
            "interpretation": "Exploratory deterministic consistency test. Four reported errors is NOT a calibrated confidence region."}


class PrimitiveUnpickler(pickle.Unpickler):
    def find_class(self, module, name):
        raise ValueError(f"non-primitive pickle object rejected: {module}.{name}")


def mutual_information(counts):
    counts = np.asarray(counts, dtype=float).reshape(2, 2)
    if np.any(counts < 0) or counts.sum() <= 0:
        raise ValueError("invalid contingency counts")
    p = counts / counts.sum()
    independent = p.sum(axis=0)[None, :] * p.sum(axis=1)[:, None]
    selected = p > 0
    return float(np.sum(p[selected] * np.log2(p[selected] / independent[selected])))


def mechanical_analysis(cache):
    p = cache / "mechanical/Zenodo/Fig3aceFig5Fig7FigS11/data/statistics_wrapped_period"
    with p.open("rb") as f:
        a = np.asarray(PrimitiveUnpickler(f).load(), dtype=np.int64)
        if f.read():
            raise ValueError("unexpected extra object")
    if a.shape != (714, 4) or np.any(a < 0):
        raise ValueError("unexpected phase-count schema")
    count = a.sum(axis=0)
    weights = a.sum(axis=1) / a.sum()
    phase = a.reshape(-1, 2, 2) / a.sum(axis=1)[:, None, None]
    px = phase[:, 1, :].sum(axis=1)
    py = phase[:, :, 1].sum(axis=1)
    total_cov = count[3]/count.sum() - ((count[2]+count[3])/count.sum())*((count[1]+count[3])/count.sum())
    within_cov = float(weights @ (phase[:, 1, 1] - px*py))
    phase_cov = float(weights @ (px*py) - (weights @ px)*(weights @ py))
    return {"phase_bins": 714, "phase_bin_s": 0.001, "nominal_period_s_from_notebook": 0.7144833333333334,
            "total_counts": int(a.sum()), "contingency_counts": count.tolist(),
            "pooled_MI_bits_plugin": mutual_information(count),
            "weighted_phase_conditional_MI_bits_plugin": float(sum(w*mutual_information(c) for w, c in zip(weights, a))),
            "total_covariance": float(total_cov), "within_phase_covariance": within_cov,
            "between_phase_covariance": phase_cov,
            "between_phase_covariance_fraction": phase_cov/total_cov,
            "covariance_identity_residual": float(total_cov - within_cov - phase_cov),
            "interpretation": "Descriptive reconstruction from phase-aggregated counts; no independent-trial CI or microscopic source identification."}


def iaia_analysis(path):
    with zipfile.ZipFile(path) as z:
        a = np.asarray(rows(z, "Figure4b.dat"), dtype=float)
        b = np.asarray(rows(z, "Figure2a.dat"), dtype=float)
    if a.shape != (7, 16) or b.shape[1] != 8:
        raise ValueError("unexpected Iaia figure schema")
    return {"observed_nonCu_to_Cu_rate_ratios": {
        "single_qubit": (a[:3, 1]/a[:3, 9]).tolist(),
        "twofold": (a[3:6, 1]/a[3:6, 9]).tolist(),
        "threefold": float(a[6, 1]/a[6, 9])},
        "injection_Figure2a_rows": len(b),
        "nonCu_peak_delay_us": float(b[np.argmax(b[:, 5]), 4]),
        "interpretation": "Ratios of published observed rates, not uncertainty intervals or independent event-level reconstructions."}


def compute(cache):
    contrast = 0.2
    fixed = -np.log1p(-contrast)/2
    witness = two_point_witness(contrast, 5.0)
    return {"status": "FEASIBILITY GATE: no certified empirical combination gain; mission not achieved",
        "gamma": gamma_analysis(cache / "gamma.zip"),
        "mechanical": mechanical_analysis(cache), "iaia": iaia_analysis(cache / "iaia.zip"),
        "synthetic_observation_witness": {"contrast": contrast,
            "homogeneous_intensity": float(fixed), "homogeneous_any_probability": any_probability(float(fixed)),
            "heterogeneous_model": witness,
            "heterogeneous_any_probability": sum(w*any_probability(l) for w, l in zip(witness["weights"], witness["intensities"])),
            "heterogeneous_contrast": sum(w*2*odd_probability(l) for w, l in zip(witness["weights"], witness["intensities"])),
            "scope": "Synthetic scalar parity-observation witnesses, NOT models verified against complete source datasets."}}


def close(actual, expected, path="root"):
    if isinstance(expected, dict):
        if set(actual) != set(expected):
            raise AssertionError(path + ": keys differ")
        for k in expected:
            close(actual[k], expected[k], path + "." + k)
    elif isinstance(expected, list):
        if len(actual) != len(expected):
            raise AssertionError(path + ": lengths differ")
        for i, (a, e) in enumerate(zip(actual, expected)):
            close(a, e, path + f"[{i}]")
    elif isinstance(expected, float):
        if not np.isclose(actual, expected, rtol=2e-8, atol=2e-10):
            raise AssertionError(f"{path}: {actual} != {expected}")
    elif actual != expected:
        raise AssertionError(path + ": values differ")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", required=True, type=Path)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    acquire(args.cache, verify_only=True)
    result = compute(args.cache)
    target = HERE / "results/analysis.json"
    if args.check:
        close(result, json.loads(target.read_text()))
        print("All measured summaries, rational certificates and observation bounds reproduced")
    else:
        target.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
        print(target)
