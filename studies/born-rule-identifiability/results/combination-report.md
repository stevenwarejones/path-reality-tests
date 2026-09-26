# Expanded dataset-combination findings

Main result: a normalized single-qubit probability deformation is exactly indistinguishable from one ordinary control-axis warp throughout the audited integer gate-depth circuit family. Combining depths cannot identify it. A complementary-preparation witness can separate the pair ideally, but contrast/calibration uncertainty removes that prospective gain. No measured Born-rule exclusion is obtained. The waveguide runner-up supplies a separate exact optical ambiguity and phase-cycling design.

## Public source audit

| Archive | Evidence | Backend | Jobs | Circuit rows | Largest first-harmonic residual |
|---|---|---|---:|---:|---:|
| lima1-5 | observed | ibmq_lima | 290 | 29000 | 0.0016146 |
| lima1-5s | simulation | qasm_simulator | 290 | 29000 | 0.0002546 |
| lagos_bench | observed | ibm_lagos | 63 | 18396 | not mapped |

Hardware residuals are descriptive, with job-level errors in [portfolio-audit.json](portfolio-audit.json). They are not a theta fit or discovery significance. The simulator is not independent hardware evidence. Missing benchmark gate labels are not inferred from order.

| Photon archive | Outcome matrix | Available event settings | Max event/P discrepancy | Max POVM row-sum error |
|---|---|---:|---:|---:|
| 8-bin TMD | [14, 9] | 14 | 0.0000000 | 0.0025641 |
| PNR SNSPD | [122, 5] | 26 | 0.0069035 | 0.0000000 |

Zero outcomes are retained. The supplied detector POVMs are derived under assumed coherent-state statistics; they are not independent Born-free calibration. PNR classified events and fitted probability products differ and cannot be pooled as independent replications.

## Exact single-versus-joint result

For fθ(p)=p+θp(1−p)(2p−1), hθ(x)=2fθ((1+x)/2)−1, choose cos(gθ(φ))=hθ(cosφ). Since sin(nπ/2) is 0 or ±1, fθ(p_n(φ))=p_n(gθ(φ)) for every integer n. Independent complex circuit multiplication verifies the convention. Both datasets individually and jointly remain ambiguous with the same one-parameter warp. This is an exact conditional observation-map theorem, not an assertion that all real residuals follow this model.

| θ | Largest numerical equality error over checked depths | Ideal complementary squared radius |
|---:|---:|---:|
| -0.1 | 2.2e-16 | 0.9506250 |
| 0.01 | 2.2e-16 | 1.0050063 |
| 0.1 | 2.2e-16 | 1.0506250 |

## Raw-count feasibility and removal

Each of four prospective preparations uses 17,400,000 shots, matching the hardware pooled per-angle depth-1 scale. These new preparations were **not** measured in the public archive. Simultaneous exact binomial intervals use family alpha 0.01.

| Scenario | Rejections / 1000 | Rejections with radius 1.01 | MC 95% interval |
|---|---:|---:|---|
| ordinary | 0 | 0 | [0.0000, 0.0037] |
| probability_deformation | 1000 | 0 | [0.9963, 1.0000] |
| ordinary_common_axis_warp | 0 | 0 | [0.0000, 0.0037] |
| deformation_99_percent_contrast | 0 | 0 | [0.0000, 0.0037] |
| deformation_hardware_reference_contrast | 0 | 0 | [0.0000, 0.0037] |
| negative_deformation | 0 | 0 | [0.0000, 0.0037] |

Removing the complementary preparation restores exact ambiguity at any shot count. Relaxing the preparation/measurement consistency removes the unit-disk certificate; radius 1.01 is a sensitivity scenario, not a measured bound. Loss of contrast and negative theta can remain inside the ordinary disk, so the witness is not a general identifiability result.

The photon-calibration negative control separately verifies (P T)R_d=P(T R_d) with independent responses for two detectors; removing either detector does not remove this ordinary stochastic gauge. Treating T as a Born-rule modification would be a category error.

## What is and is not new

Known ingredients include quantum control gauges, phase cycling, PSD coherence and concentration bounds. This PR provides an explicit all-depth counterexample, public byte-level audit, complementary-preparation cost/sensitivity, and an optical runner-up. Physics novelty is not established. The expanded [candidate inventory and bridges](../combinations.md) gives four proposals and precise reasons the available combinations do not yet support a probability-rule exclusion.

Highest-value next input: independently certified complementary ±X/±Y preparations, with common channel/readout and contrast bounds. The optical alternative is a calibrated antipodal all-open phase exposure. No author contact, new cloud experiment, or raw redistribution was performed.
