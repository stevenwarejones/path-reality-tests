"""Deterministic report and compact scientific plots from reviewed derived artifacts."""
from pathlib import Path
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent


def generate(results,output):
    output = Path(output); output.mkdir(parents=True,exist_ok=True)
    b = results["baselines"]["data"]
    c = results["compatibility"]["data"]
    selected = c["selected"]["rows"]
    first = [r for r in selected if r["benchmark"]=="100 nm silica sphere" and r["rc_m"]==1e-7]
    white = first[0]; slow = next(r for r in first if r["tau_s"]==1)
    f = b["force"]; n = b["sodium"]; m = b["macroscopicity"]
    rad = c["runner_up"]; rw = rad["rows"][0]
    lines = [
        "# Collapse compatibility: audited baselines and conditional spectral feasibility",
        "",
        "No new physical exclusion is established. This report retains the original stationary comparison. The [moving-field extension](moving-report.md) now recomputes a joint ensemble-dephasing witness using one frame, a moving benchmark, trajectory metadata, assembly bounds, finite windows and the sodium mass/velocity mixture. Its empirical and approximation limits are explicit.",
        "",
        "## Observed-source reproduction",
        "",
        f"- Sodium: {n['scan_count']} scans, {n['count_bins']} count bins; interpolated mean mass {n['mean_mass_u']:.1f} u. Independent linear versus author nonlinear visibility fits differ by at most {n['max_visibility_fit_discrepancy']:.3g}. The RMS fixed-optics visibility residual is {n['rms_quantum_visibility_residual']:.5f}.",
        f"- Deposited macroscopicity posterior convention: tau_e 5% quantile {m['tau5_s']:.6g} s, log10(tau_e)={m['log10_tau5']:.6f}. This fixes empirical phase, source rate and calibration; it is not a new confidence limit.",
        f"- Sensor workbook: reproduced force estimate {f['literal_workbook']['force_estimate_N2_per_hz']:.6g} N²/Hz and quoted 95% upper {f['literal_workbook']['upper95_N2_per_hz']:.6g} N²/Hz.",
        f"- Independent B-error-only temperature regression differs from the published orthogonal-fit intercept by {f['B0_discrepancy_in_published_se']:.3f} published standard errors. Missing T/Q covariance prevents exact orthogonal-fit replay.",
        f"- All 14 spectra are fitted with explicit six-bin masks and assumed Qprime. Largest absolute B discrepancy is {max(abs(r['relative_B_discrepancy']) for r in f['spectra']):.2%}. These are conditional reproductions, not exact replicas.",
        "- Workbook audit: M26 sums error contributions; N23 uses stiffness uncertainty in a derivative that calls for stiffness. The independent-quadrature recalculation is a sensitivity check, not a replacement confidence bound. The paper supplement's reported uncertainty also differs from its main text and workbook.",
        f"- XENONnT 2022 subset: {rad['audit']['events']} events; no-background-subtraction Poisson signal-count upper {rad['audit']['signal_counts_upper95_no_background_subtraction']:.3f}. The newer 1–140 keV collapse analysis is not reproduced by this 1–30 keV release.",
        "",
        "## Selected combination: kHz layer sensor + mHz space envelope",
        "",
        "At rc=100 nm, for a static 100 nm silica sphere separated by 100 nm, 10 ms duration and required coherence factor exp(-10):",
        "",
        "| OU correlation time | Required rate /s | Sensor alone upper /s | Space alone upper /s | Joint upper /s | Compatible |",
        "|---|---:|---:|---:|---:|---|",
    ]
    for row in first:
        lines.append(f"| {row['tau_s']:.3g} s | {row['required_lambda']:.3g} | {row['force_only_lambda_upper']:.3g} | {row['space_only_lambda_upper']:.3g} | {row['joint_lambda_upper']:.3g} | {row['benchmark_feasible']} |")
    lines += [
        "",
        f"At white noise the sensor dominates. At tau=1 s the space envelope improves on the sensor alone by {slow['force_only_lambda_upper']/slow['space_only_lambda_upper']:.3g}, but still permits the declared benchmark. Removing either dataset recovers the corresponding single-source column. This is complementary physics, not increased event statistics.",
        "",
        "All selected benchmark/grid combinations remain compatible, including half/double envelope stress factors. These stress factors are not calibrated uncertainty coverage. Layer-thickness sensitivity is retained in compatibility.json.",
        "",
        "The selected values are deterministic comparisons against quoted envelopes. They have no asserted joint 95% coverage, use an isolated sensor load, and assume stationary-body spectral factorization. The common noise rest frame is not fixed by the records.",
        "",
        "## Runner-up: sensor + radiation",
        "",
        f"At white noise and rc=100 nm the approximate 1–30 keV radiation recast gives rate uppers {rw['radiation_only_upper_alpha_range'][0]:.3g}–{rw['radiation_only_upper_alpha_range'][1]:.3g} /s across the published atomic-distance range. The detector efficiency and charge cancellation are included; resolution/migration is not. This is a feasibility result, not a certified benchmark exclusion.",
        "",
        "The runner-up JSON retains radiation alone, sensor alone and their minimum at every correlation time. Radiation loses leverage rapidly when the spectrum suppresses photon frequencies; the white bound is never copied into a colored analysis. Efficiency and atomic-distance sensitivity are explicit.",
        "",
        "## Conditional limitation and scope",
        "",
        "Choosing lambda(tau)=C/[G H(T,tau)] keeps stationary macroscopic suppression at C. Every finite positive-frequency response tends to zero as tau grows; the covariance amplitude lambda/(2tau) stays finite. Saved witnesses meet the benchmark exactly. The exact affine certificate checker separately rejects missing tails and failed between-grid domination.",
        "",
        "This is not a demonstrated universal-model escape. Moving through spatially correlated noise shifts the sampled frequencies; the local sodium diagnostic differs by nearly three orders of magnitude from a stationary holding-time approximation at a 1 microsecond memory time. Actual instrument leakage may also measure nominally unobserved frequencies.",
        "",
        "## Novelty, validation and next input",
        "",
        "The broad strategy, mechanical spectral robustness, and force/torque noise ratios have substantial prior art. No breakthrough or novel physical theorem is claimed. This PR contributes reproducible audits, executable candidate comparison, and explicit reasons the ambitious exclusion is not yet enabled.",
        "",
        "Validation includes independent time/Fourier geometry calculations, the inspected author-module oracle, synthetic count/PSD injection tests, simultaneous-bound calibration, deliberate spectral gaps and invalid-tail certificates. Offline checks recompute conditional outputs and regenerate this report; full-source checks additionally reread hash-pinned originals.",
        "",
        "The moving-field extension supplies a conditional response calculation. The remaining high-value inputs are a calibrated approximation-error envelope for sodium and a reconstructed LPF torque/attitude likelihood. A slowly held massive superposition at several durations would directly address the remaining zero-frequency response.",
        "",
        "See [derivation](../derivation.md), [candidate comparison](../combinations.md), [source dictionary](../data-dictionary.md), and [literature audit](../literature.md).",
        "",
        "![Observed baselines and conditional comparisons](diagnostics.svg)",
        ""
    ]
    (output/"report.md").write_text("\n".join(lines))
    plt.rcParams.update({"svg.hashsalt":"collapse-compatibility-v1","svg.fonttype":"none","font.size":9})
    fig,axes = plt.subplots(2,2,figsize=(11,7),layout="constrained")
    scans = n["scans"]
    x = [r["g2_power_mw"] for r in scans]
    axes[0,0].scatter(x,[r["visibility"] for r in scans],s=9,label="fitted measurements")
    axes[0,0].scatter(x,[r["quantum_visibility"] for r in scans],s=8,marker="x",label="fixed optics")
    axes[0,0].set(xlabel="G2 power (mW)",ylabel="visibility",title="Sodium: measured baseline")
    axes[0,0].legend(fontsize=8)
    inp = f["temperature_fit_inputs"]
    axes[0,1].errorbar([r["T_over_Q_K"]*1e9 for r in inp],[r["B"]*1e18 for r in inp],
                       yerr=[r["se_B"]*1e18 for r in inp],fmt="o",ms=3)
    xx = np.linspace(0,1500,100); fit = f["weighted_y_only_fit"]
    axes[0,1].plot(xx,(fit["B0"]+fit["B1"]*xx*1e-9)*1e18,label="T >=100 mK fit")
    axes[0,1].set(xlabel="T/Q (nK)",ylabel="B (1e-18 Phi0²/Hz)",title="Sensor: thermal baseline")
    axes[0,1].legend(fontsize=8)
    nonzero = first[1:]
    for key,label in [("force_only_lambda_upper","sensor"),("space_only_lambda_upper","space"),("required_lambda","benchmark")]:
        axes[1,0].loglog([r["tau_s"] for r in nonzero],[r[key] for r in nonzero],label=label)
    axes[1,0].set(xlabel="OU correlation time (s)",ylabel="rate (1/s)",title="Conditional stationary envelope model")
    axes[1,0].legend(fontsize=8)
    rr = [r for r in rad["rows"] if r["rc_m"]==1e-7 and r["tau_s"]>0]
    axes[1,1].loglog([r["tau_s"] for r in rr],[r["radiation_only_upper_alpha_range"][1] for r in rr],label="radiation, approximate response")
    axes[1,1].loglog([r["tau_s"] for r in rr],[r["force_only_upper"] for r in rr],label="sensor")
    axes[1,1].set(xlabel="OU correlation time (s)",ylabel="rate upper (1/s)",title="Runner-up: colored-noise sensitivity")
    axes[1,1].legend(fontsize=8)
    fig.savefig(output/"diagnostics.svg",metadata={"Date":None})
    svg = output/"diagnostics.svg"
    svg.write_text("\n".join(line.rstrip() for line in svg.read_text().splitlines())+"\n")
    plt.close(fig)
