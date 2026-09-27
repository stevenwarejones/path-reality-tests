"""Generate the measured combination report from checked outputs."""
import numpy as np


def render(result, data):
    joint = result['joint']
    lines = ['# Measured scan–Viviani joint analysis', '',
        '**Finding:** two distinct measured preparation–measurement designs constrain different',
        'patterns of the proposed probability deformation, but both full mean tables admit an',
        'explicit ordinary qutrit explanation with small preparation-dependent leakage.',
        'This is a bounded identifiability result, not a measured Born-rule violation.', '',
        'The numerical deformation profiles are feasible fits, not globally certified exclusions.',
        'The independent-shot lower bounds are conditional; allowing arbitrary within-job',
        'dependence removes their statistical separation. See [assumptions and derivations](../joint-tables.md).', '',
        '## Distinct acquisitions', '',
        '| Acquisition | Measured table | Jobs | Shots per setting per job | Recorded metadata jobs / job-list rows |',
        '|---|---:|---:|---:|---:|']
    for kind in ['scan', 'viviani']:
        row = data['records'][kind]
        n = len(row['counts_0_1_by_job'][0])
        lines.append(f"| {kind} | {n} × 4 | {row['actual_jobs']} | {row['shots_per_circuit']*row['repetitions_per_setting_job']:,} | {row['metadata'].get('jobs')} / {row['job_list_rows']} |")
    lines += ['', 'Nairobi scan: December 2022, [7470893](https://zenodo.org/records/7470893).',
              'Nairobi Viviani: October 2023, [21775462 v7](https://zenodo.org/records/21775462).',
              'Counts/settings are read from every circuit row. Metadata job totals and job lists',
              'are incomplete; no missing jobs are imputed. Apparatus parameters are independent.', '',
              '## Physical compatibility boundary', '',
              'Metric: maximum absolute probability change needed to reach any stable ordinary',
              'qubit table. Lower bounds use centered rank ≤3; upper bounds use explicit physical states/effects.', '',
              '| Acquisition | Exact-mean lower bound | 95% joint conditional lower bound | Physical qubit residual upper | Smallest tested leakage cap fitting every mean (<2e-8) |',
              '|---|---:|---:|---:|---:|']
    for kind in ['scan', 'viviani']:
        r = result['records'][kind]
        b = r['boundary']
        caps = [f['leakage_cap'] for f in r['leakage'] if f['max_abs_residual'] < 2e-8]
        lines.append(f"| {kind} | {b['mean_table_distance_lower']:.8g} | {b['iid_shot_distance_lower']:.8g} | {r['ordinary']['max_abs_residual']:.8g} | {min(caps) if caps else 'none'} |")
    lines += ['', 'Leakage caps are feasible upper constructions, **not** estimated leakage or proven minima.',
              'The quantum model is block diagonal: a qubit plus one leakage level; its leakage population',
              'depends on preparation, while each measurement has one leakage response independent of preparation.',
              'All states/effects, full predictions and numerical residuals are saved in `joint-tables.json`.',
              'This also supplies an ordinary persistent classical-flag explanation. Constant leakage alone',
              'would not supply the extra rank. No calibration record here bounds this preparation dependence.', '',
              f"The joint exact-mean lower bound is {joint['ordinary_mean_table_distance_lower']:.8g};",
              f"the simultaneous independent-shot lower bound is {joint['ordinary_iid_shot_distance_lower_95pct']:.8g}.",
              f"The independent-job-only lower bound is {joint['ordinary_independent_job_distance_lower_95pct']:.8g}.", '',
              '## Separate, joint and removal profiles', '',
              'Pure-direction errors, antipodal preparation mixtures and asymmetric measurement readout',
              'are fitted independently in each acquisition. Tangent coordinates are bounded by ±0.15;',
              'this is a declared sensitivity family, not an externally certified control calibration.', '',
              '| Analysis | Best sampled theta | Interpretation |', '|---|---:|---|',
              f"| Scan alone / remove Viviani | {joint['remove_viviani_theta']:.6g} | Extra scan settings constrain deformation shape |",
              f"| Viviani alone / remove scan | {joint['remove_scan_theta']:.6g} | Independent preparation–measurement consistency |",
              f"| Joint | {joint['selected_grid_theta']:.6g} | One common theta; independent apparatus |",
              '| Add/remove Lima rotation family | unchanged identification | Exact angle-warp ambiguity survives |', '',
              'These are grid minimizers of equally weighted acquisition mean-square residuals, **not**',
              'theta estimates with confidence intervals. They neither prove the alternatives incompatible',
              'nor establish unique identification. General preparation ensembles are outside this restricted fit.', '',
              '| theta | Scan maximum residual | Viviani maximum residual | Joint mean-square objective |',
              '|---:|---:|---:|---:|']
    for i, theta in enumerate(joint['grid_theta']):
        a = result['records']['scan']['deformation_profile'][i]
        b = result['records']['viviani']['deformation_profile'][i]
        lines.append(f"| {theta:g} | {a['max_abs_residual']:.6g} | {b['max_abs_residual']:.6g} | {joint['equal_acquisition_mse_profile'][i]:.6g} |")
    lines += ['', '## Frozen held-out job-block check', '',
              'The first two thirds of archive-numbered jobs train the model; the final third is held out.',
              'This retrospective split follows public-data/paper inspection and is not a prospective blind test.',
              f"The train-only joint grid selects theta={joint['train_selected_grid_theta']:g}. Predictions are then frozen.", '',
              '| Acquisition | Held-out max error, common deformation | Held-out max error, ordinary qutrit |',
              '|---|---:|---:|']
    for kind in ['scan', 'viviani']:
        lines.append(f"| {kind} | {joint['holdout'][kind]['max_abs_error']:.6g} | {result['records'][kind]['holdout_train_leakage_max_error']:.6g} |")
    lines += ['', 'Both train-fitted models miss the held-out job block: a mean-table explanation is not a',
              'validated stationary model of the complete acquisition. Deviances are descriptive;',
              'no Wilks or chi-square calibration is used.', '', '## Drift and ordinary adversaries', '']
    for kind in ['scan', 'viviani']:
        fits = result['records'][kind]['job_leakage_constructions']
        residual = max(f['max_abs_residual'] for f in fits)
        passed = sum(f['max_abs_residual'] < 2e-8 for f in fits)
        lines += [f"For {kind}, job-specific ordinary qubit-plus-leakage constructions with cap .005",
                  f"fit {passed}/{len(fits)} individual job tables to <2e-8; maximum residual is {residual:.5g}."]
    lines += ['These constructions allow independently varying instruments between jobs. They are',
              'compatibility certificates, not a fitted time-evolution law or successful held-out predictions.',
              'Complete acquisition timing is missing, so the job-block split is not certified chronological.', '']
    drift = result['aria_drift_control']
    lines += [f"The separate measured Aria control has pooled determinant {drift['pooled_determinant']:.7g}",
              f"and mean job determinant {np.mean(drift['job_determinants']):.7g}; maximum job probability range",
              f"is {drift['boundary']['job_probability_range_max']:.5g}. This is a published drift issue, not a new anomaly.",
              'Within-job singular values, circuit-order residual diagnostics and job variation are retained.',
              'Circuit order is recorded; shot chronology and complete timestamps are unavailable.', '',
              'Count-level simulations use actual job/cell shot totals and physically fitted contrasts.',
              'A rejection below rejects the **stable qubit model class**, not ordinary quantum theory.',
              'Drift, leakage and memory intentionally violate that class and expose false Born-rule attribution.', '',
              '| Acquisition | Generating explanation | Stable-qubit rejections / simulations | Known-nuisance theta projection mean |',
              '|---|---|---:|---:|']
    for kind, cases in result['simulations'].items():
        for name, r in cases.items():
            lines.append(f"| {kind} | {name} | {r['stable_qubit_class_rejections']}/{r['replications']} | {r['theta_projection_mean_known_nuisance']:.6g} |")
    lines += ['', 'Projection recovery assumes the generating apparatus is known; it is not recovery after',
              'profiling all nuisances. Models with arbitrary within-job dependence are not excluded.', '',
              '## Added value and next measurement', '',
              'The implemented combination goes beyond the quarter-turn family by using two different',
              'measured preparation geometries. Removal changes the restricted deformation profile.',
              'Its defensible negative result is the joint physical qubit-plus-one-level construction',
              'and explicit probability-distance bracket, with dependence sensitivity—not novelty of rank witnesses.',
              'The most valuable missing observation is preparation-resolved leakage discrimination,',
              'interleaved with the same randomized tables and time-tagged repeats. It must constrain',
              'both leakage populations and leakage-level readout; another uncalibrated determinant is insufficient.', '',
              'Reproduce: `python studies/born-rule-identifiability/joint_tables.py --cache /tmp/born-rule-sources --refit`.',
              'Check saved constructions and raw-source reduction: use `--check` with the same cache.', '']
    return '\n'.join(lines)
