"""Deterministic review report; source originals remain external."""
def render(out):
    rep=out['reproduction.json']; ident=out['identifiability.json']
    lines=['# Optical identifiability: limited result','',
        'The archived mean tables admit an explicit ordinary partially coherent optical model with a shutter-context phase shift, indistinguishable from an additive all-open intensity diagnostic. This is a conditional mean-table construction, not a full temporal likelihood fit or a Born-rule violation. No new fundamental theorem or priority is claimed.','',
        '## Reproduction','',
        'Three Zenodo 5497283 originals are pinned by size, MD5 and SHA-256 in [manifest](../manifest.json). The analysis retains whole cycles and covariance. Intensities are detector volts; independent voltage-to-power calibration is absent.','',
        '| Label | Cycles | Reproduced F | Published F | Reproduced κ | Published κ |',
        '|---|---:|---:|---:|---:|---:|']
    for t,d in rep['datasets'].items():
        lines.append(f"| {t}°C | {d['retained_cycles']} | {d['mean_F']:.8f} | {d['published_summary']['F']:.4f} | {d['kappa_ratio_of_means']:.8f} | {d['published_summary']['kappa']:.5f} |")
    lines+=['','κ uses the ratio of means (paper B.3) and reproduces the reported precision. The ~0.0001 F differences remain unresolved without author processing code; no correction is tuned to force agreement. Published exclusions are one-based; zero-based and unfiltered sensitivities remain in the JSON. The 23°C temperature array is absent. The 30°C housing array has one NaN. Contrast data use four-setting cycles and have an unresolved temperature scale.','',
        'Bartlett HAC lags 0,5,10,20,40 preserve original cycle gaps. These descriptive estimates do not have asserted finite-sample coverage. Drift slopes are reported without detrending. No measured discovery p-value is calculated.','',
        '## Exact mean-table construction','',
        'Six lower-setting intensities fix real H. A PSD completion G has imaginary AB entry u with |u|≤sqrt(det(H)/H_CC). A phase δ acting on A in the all-open context changes the intensity by 2(H_AB+H_AC)(cosδ−1)+2u sinδ. This occupies exactly the same cell as θ·tr(H). The full complex quadratic form independently verifies each pair.','',
        '| Label | Ordinary δ (rad) | Alternative θ | Maximum mismatch (V) |','|---|---:|---:|---:|']
    for t,d in ident['datasets'].items():
        lines.append(f"| {t}°C | {d['ordinary_phase_rad']:.6f} | {d['alternative_theta']:.7f} | {d['max_table_residual_V']:.2g} |")
    lines+=['','These are allowed nuisance models, not an inferred mechanism or the paper’s particular crosstalk model. Partial coherence and context phase would need independent bounds. The diagnostic is not a normalized alternative quantum theory.','',
        '## Separate versus joint','',
        'Individual regions are in [identifiability.json](identifiability.json). Intersections below share only diagnostic θ; each temperature keeps independent coherence and phase nuisances. Real H is fixed to the estimated mean. These are conditional feasibility regions, not confidence intervals.','',
        '| Phase radius (rad) | Joint θ endpoints | Nonempty | Contains zero |','|---:|---|---|---|']
    for r in ident['joint']:
        lo,hi=r['theta_interval'];lines.append(f"| {r['phase_radius_rad']:.2f} | [{lo:.7f}, {hi:.7f}] | {r['nonempty']} | {r['contains_zero']} |")
    lines+=['','Reversed endpoints indicate an empty intersection. Opposite residual signs reject the restricted zero-phase shared diagnostic, not quantum mechanics. Ordinary θ=0 survives the independently varying .03-radian regions. No AION/NIST/Wen likelihood is multiplied into this optical result.','',
        '## Additional observable and adversarial tests','',
        'Add an all-open exposure advancing A by π. W=(I_ABC(0)+I_ABC(π))/2−I_A−I_BC+I_empty cancels the unknown A phase and coherence in the declared model, while retaining a phase-independent intensity defect. One extra setting suffices for the constructed pair; this is not universal Born-rule identifiability.','',
        'Phase error e contributes at most 2 sqrt(G_AA q_BC) sin(e/2), and an independently certified per-reading bias r adds 4r. Prospective budgets use five independent bounded setting means, alpha .01 and power at least .9; worst-case separation is signal minus twice the systematic bound. These are conditional design costs, not achieved apparatus sensitivity.','',
        '| Raw-intensity simulation | Rejections / 4000 | Monte Carlo 95% interval |','|---|---:|---|']
    for r in out['simulations.json']['rows']:
        lo,hi=r['binomial_95_interval'];lines.append(f"| {r['scenario']} | {r['rejected']} | [{lo:.4f}, {hi:.4f}] |")
    lines+=['','Ignoring leakage/nonlinearity inflates false positives. An analytic synthetic envelope prevents these exclusions but can erase power for the target signal. Gaussian known-noise simulations include common block power/phase drift; arbitrary within-setting drift or dependence remains outside their scope.','',
        '## Scope, novelty and next measurement','',
        'Quadratic interference, coherence tomography, phase cycling and apparatus-induced higher-order signals are known. The contribution is the reproducible source audit, archive-scale explicit pairs, sharp conditional region and costed extra setting. See [derivation](../derivation.md), [theory](../theory.md), [literature](../literature.md) and the expanded [combination search](../combinations.md).','',
        'The most useful next optical measurement is an independently calibrated antipodal all-open exposure with detector, leakage and drift bounds. No substantive new formal physics theorem is established, so a Lean companion is not opened.','']
    return '\n'.join(lines)
