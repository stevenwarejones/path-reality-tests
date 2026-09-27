# Measured scan–Viviani joint analysis

**Finding:** two distinct measured preparation–measurement designs constrain different
patterns of the proposed probability deformation, but both full mean tables admit an
explicit ordinary qutrit explanation with small preparation-dependent leakage.
This is a bounded identifiability result, not a measured Born-rule violation.

The numerical deformation profiles are feasible fits, not globally certified exclusions.
The independent-shot lower bounds are conditional; allowing arbitrary within-job
dependence removes their statistical separation. See [assumptions and derivations](../joint-tables.md).

## Distinct acquisitions

| Acquisition | Measured table | Jobs | Shots per setting per job | Recorded metadata jobs / job-list rows |
|---|---:|---:|---:|---:|
| scan | 9 × 4 | 115 | 800,000 | 20 / 20 |
| viviani | 5 × 4 | 9 | 1,500,000 | 1 / 1 |

Nairobi scan: December 2022, [7470893](https://zenodo.org/records/7470893).
Nairobi Viviani: October 2023, [21775462 v7](https://zenodo.org/records/21775462).
Counts/settings are read from every circuit row. Metadata job totals and job lists
are incomplete; no missing jobs are imputed. Apparatus parameters are independent.

## Physical compatibility boundary

Metric: maximum absolute probability change needed to reach any stable ordinary
qubit table. Lower bounds use centered rank ≤3; upper bounds use explicit physical states/effects.

| Acquisition | Exact-mean lower bound | 95% joint conditional lower bound | Physical qubit residual upper | Smallest tested leakage cap fitting every mean (<2e-8) |
|---|---:|---:|---:|---:|
| scan | 0.00010661328 | 3.5803653e-05 | 0.00015246344 | 0.001 |
| viviani | 0.00032222679 | 0.00011911093 | 0.00038385007 | 0.002 |

Leakage caps are feasible upper constructions, **not** estimated leakage or proven minima.
The quantum model is block diagonal: a qubit plus one leakage level; its leakage population
depends on preparation, while each measurement has one leakage response independent of preparation.
All states/effects, full predictions and numerical residuals are saved in `joint-tables.json`.
This also supplies an ordinary persistent classical-flag explanation. Constant leakage alone
would not supply the extra rank. No calibration record here bounds this preparation dependence.

The joint exact-mean lower bound is 0.00032222679;
the simultaneous independent-shot lower bound is 0.00011911093.
The independent-job-only lower bound is 0.

## Separate, joint and removal profiles

Pure-direction errors, antipodal preparation mixtures and asymmetric measurement readout
are fitted independently in each acquisition. Tangent coordinates are bounded by ±0.15;
this is a declared sensitivity family, not an externally certified control calibration.

| Analysis | Best sampled theta | Interpretation |
|---|---:|---|
| Scan alone / remove Viviani | 0 | Extra scan settings constrain deformation shape |
| Viviani alone / remove scan | 0.01 | Independent preparation–measurement consistency |
| Joint | 0.004 | One common theta; independent apparatus |
| Add/remove Lima rotation family | unchanged identification | Exact angle-warp ambiguity survives |

These are grid minimizers of equally weighted acquisition mean-square residuals, **not**
theta estimates with confidence intervals. They neither prove the alternatives incompatible
nor establish unique identification. General preparation ensembles are outside this restricted fit.

| theta | Scan maximum residual | Viviani maximum residual | Joint mean-square objective |
|---:|---:|---:|---:|
| -0.03 | 0.0020711 | 0.00137873 | 2.11988e-06 |
| -0.02 | 0.0013401 | 0.00103764 | 1.09583e-06 |
| -0.01 | 0.000622056 | 0.000709789 | 4.32642e-07 |
| -0.005 | 0.000333247 | 0.000546445 | 2.30843e-07 |
| 0 | 0.000152463 | 0.00038385 | 1.15196e-07 |
| 0.001 | 0.000174484 | 0.000336074 | 9.48679e-08 |
| 0.002 | 0.000228213 | 0.000288309 | 7.97631e-08 |
| 0.003 | 0.000300078 | 0.000240556 | 6.98378e-08 |
| 0.004 | 0.0003617 | 0.000192814 | 6.49865e-08 |
| 0.005 | 0.000431319 | 0.000145084 | 6.52416e-08 |
| 0.006 | 0.000500885 | 9.73648e-05 | 7.06161e-08 |
| 0.008 | 0.00063986 | 1.96019e-06 | 9.67153e-08 |
| 0.01 | 0.000778622 | 2.81641e-12 | 1.37096e-07 |
| 0.015 | 0.00112457 | 0.000107182 | 2.79121e-07 |
| 0.02 | 0.0014691 | 0.00027217 | 5.03621e-07 |
| 0.03 | 0.00215452 | 0.000596566 | 1.20314e-06 |

## Frozen held-out job-block check

The first two thirds of archive-numbered jobs train the model; the final third is held out.
This retrospective split follows public-data/paper inspection and is not a prospective blind test.
The train-only joint grid selects theta=0.004. Predictions are then frozen.

| Acquisition | Held-out max error, common deformation | Held-out max error, ordinary qutrit |
|---|---:|---:|
| scan | 0.017255 | 0.0172477 |
| viviani | 0.00482012 | 0.00499311 |

Both train-fitted models miss the held-out job block: a mean-table explanation is not a
validated stationary model of the complete acquisition. Deviances are descriptive;
no Wilks or chi-square calibration is used.

## Drift and ordinary adversaries

For scan, job-specific ordinary qubit-plus-leakage constructions with cap .005
fit 113/115 individual job tables to <2e-8; maximum residual is 0.0013266.
For viviani, job-specific ordinary qubit-plus-leakage constructions with cap .005
fit 9/9 individual job tables to <2e-8; maximum residual is 8.0576e-11.
These constructions allow independently varying instruments between jobs. They are
compatibility certificates, not a fitted time-evolution law or successful held-out predictions.
Complete acquisition timing is missing, so the job-block split is not certified chronological.

The separate measured Aria control has pooled determinant 0.0009953881
and mean job determinant 0.000887462; maximum job probability range
is 0.6284. This is a published drift issue, not a new anomaly.
Within-job singular values, circuit-order residual diagnostics and job variation are retained.
Circuit order is recorded; shot chronology and complete timestamps are unavailable.

Count-level simulations use actual job/cell shot totals and physically fitted contrasts.
A rejection below rejects the **stable qubit model class**, not ordinary quantum theory.
Drift, leakage and memory intentionally violate that class and expose false Born-rule attribution.

| Acquisition | Generating explanation | Stable-qubit rejections / simulations | Known-nuisance theta projection mean |
|---|---|---:|---:|
| scan | stable_qubit | 0/200 | 1.24265e-07 |
| scan | deformation_0.01 | 200/200 | 0.0100089 |
| scan | ordinary_qutrit | 200/200 | -0.000465087 |
| scan | ordinary_correlated_drift | 200/200 | -0.0615433 |
| scan | ordinary_preparation_memory | 200/200 | -0.00742847 |
| viviani | stable_qubit | 0/200 | -9.77252e-06 |
| viviani | deformation_0.01 | 200/200 | 0.00998519 |
| viviani | ordinary_qutrit | 200/200 | 0.00284758 |
| viviani | ordinary_correlated_drift | 200/200 | -0.0193347 |
| viviani | ordinary_preparation_memory | 200/200 | -0.00601398 |

Projection recovery assumes the generating apparatus is known; it is not recovery after
profiling all nuisances. Models with arbitrary within-job dependence are not excluded.

## Added value and next measurement

The implemented combination goes beyond the quarter-turn family by using two different
measured preparation geometries. Removal changes the restricted deformation profile.
Its defensible negative result is the joint physical qubit-plus-one-level construction
and explicit probability-distance bracket, with dependence sensitivity—not novelty of rank witnesses.
The most valuable missing observation is preparation-resolved leakage discrimination,
interleaved with the same randomized tables and time-tagged repeats. It must constrain
both leakage populations and leakage-level readout; another uncalibrated determinant is insufficient.

Reproduce: `python studies/born-rule-identifiability/joint_tables.py --cache /tmp/born-rule-sources --refit`.
Check saved constructions and raw-source reduction: use `--check` with the same cache.
