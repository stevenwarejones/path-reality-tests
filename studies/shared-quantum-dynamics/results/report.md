# Shared quantum dynamics: measured combination and resource obstruction

**Question:** Can one ordinary set of operations predict both acquisition families, and what physical resources do their records identify?

**Result:** The interleaved GST/RB combination removes explicit ordinary-model ambiguities that either family leaves open. It does not identify leakage versus a bounded dormant classical flag. Within the declared incoherent leakage family, an exact all-word transformation trades loss, return, readout response and depolarization without changing any terminal probability. The sharp conditional compatibility interval is derived analytically.

**Prediction limitation:** None of the saved joint restricted-model fits passes every long-GST simultaneous cell interval. RB predictions are much better. This is not a complete temporal explanation, a rejection of all stationary qubit channels, or a necessary nonzero leakage/memory cost.

The endpoint is a restricted operational obstruction with explicit measured witnesses. Publication-level novelty is unestablished; the original paper already studied alternating errors. See [proof and model scope](../theory.md), [source audit](../sources.md), and [prior art](../prior-art.md).

## Claim status

| Claim | Status |
|---|---|
| Counts, words, denominators, preserved row order | Measured; hash-pinned and independently parsed |
| Same acquisition: GST/RB interleaving | Primary-paper assertion; timestamps absent from condensed files |
| Shared-gate physical fits and source-removal examples | Conditional numerical compatibility examples; restricted models, local optimizers |
| Exact leakage fiber and dormant-flag equality | Analytically certified within stated channel class; independently checked Kraus realizations |
| Fitted leakage populations are actual device populations | Unavailable; no operational subspace monitor |
| Globally necessary leakage or finite memory | Not established; compatible qubit training examples remain |
| Correlated drift and injected-error power | Simulated ordinary controls at actual exposures; not experimental discovery |
| Broader theory violation, novel-priority claim | Not claimed |

## Matched measured geometry

| Family | Rows | Shots | Train rows | Validation rows |
|---|---:|---:|---:|---:|
| GST | 4657 | 232850 | 3175 | 1482 |
| RB | 2004 | 563280 | 610 | 1394 |

Training is GST primitive length ≤1,024 and RB length ≤1,000. Longer sequences are held out for frozen prediction. No random split of duplicate circuits is used. Both files use Gi/Gx/Gy, one preparation, and binary terminal readout. Gates, preparation and readout are shared across families within the declared stationarity assumption. This does not connect Nairobi epochs.

The independent Gate A example is Gy⁴: GST records 0/50 plus counts, RB 3/285 in its first matching row. Local weighted Jacobian diagnostics at the same generic 18-parameter leakage point are:

| Design | Thresholded rank | Parameter count |
|---|---:|---:|
| GST | 15 | 18 |
| RB | 12 | 18 |
| joint | 15 | 18 |

The relative cutoff is 1e-6 of the largest singular value. Units and finite differences matter. These are diagnostics, not global rank certificates. The exact all-word fiber is the structural result.

## Frozen source-removal predictions

| Model | Fit source | GST validation deviance | RB validation deviance | GST / RB simultaneous-cell violations |
|---|---|---:|---:|---:|
| alternating | GST | 2672.44 | 1502.22 | 9 / 0 |
| alternating | RB | 24326.58 | 1522.85 | 372 / 0 |
| alternating | joint | 2734.80 | 1501.54 | 13 / 0 |
| leakage | GST | 2674.26 | 1513.51 | 9 / 0 |
| leakage | RB | 26944.83 | 1608.93 | 433 / 0 |
| leakage | joint | 2759.39 | 1498.26 | 9 / 0 |
| quasistatic | GST | 2672.72 | 1502.31 | 9 / 0 |
| quasistatic | RB | 24794.15 | 1520.96 | 380 / 0 |
| quasistatic | joint | 2738.56 | 1500.10 | 13 / 0 |
| qubit | GST | 2672.78 | 1502.22 | 9 / 0 |
| qubit | RB | 27894.83 | 1520.36 | 374 / 0 |
| qubit | joint | 2737.06 | 1500.08 | 13 / 0 |

Deviances are descriptive, not Wilks-calibrated global exclusions. Each displayed validation-family cell check is separately 95% simultaneous under iid Bernoulli shots per circuit (no independence between circuits is needed for its union bound). It excludes that fixed prediction vector when violated, not the full model class. Training/validation conditional independence would be needed to treat the fitted predictor as fixed for a prospective validation test. These data are retrospective and lack chronology.

## Identification gain beyond a shifted optimizer

An exploratory bank of 30 fully physical shared-operation models contains the twelve primary fits and fixed alternating-amplitude profiles. Nuisances are refitted using training counts only. A single 95% joint count region uses 3,785 exact binomial cell intervals and eight Hoeffding length-bin averages, splitting alpha equally. Source removal drops constraints from that same region; it never changes the physical model.

| Available source constraints | Bank models retained | Maximum attainable long-GST probability spread | Maximum attainable long-RB probability spread |
|---|---:|---:|---:|
| GST | 19 | 0.460122 | 0.730486 |
| RB | 18 | 0.909256 | 0.072290 |
| joint | 12 | 0.051962 | 0.063498 |

Adding RB removes 7 bank examples permitted by GST; adding GST removes 6 permitted by RB. These are explicit physical witnesses of complementary constraints, not just two separate best fits.

**Scope:** These spreads describe an attainable subset of the full prediction region. They are lower bounds on its diameter, not outer confidence bounds, and their reduction does not prove that all physical predictions lie in the displayed joint spread. The sampling bounds can exclude the individual saved prediction vectors even when the bank was chosen adaptively, because the count region is simultaneous. They do not make the selected bank exhaustive.

## Sharp conditional leakage fibers

For each fitted leakage observable vector, the following interval is analytically necessary and sufficient along its declared fixed-coordinate fiber. The parameter box is imposed for numerical sensitivity, not measured hardware calibration. Distinct points keep all measured-word probabilities equal and also agree on every unmeasured word.

| Fit source | Loss-rate interval per primitive | Population after 1,000 gates, example A / B | Max probability difference over all 6,661 rows |
|---|---|---|---:|
| GST | [4.2681137e-05, 5.4520221e-05] | 0.0268800 / 0.0334987 | < 1e-9 |
| RB | [5.9760313e-05, 6.2210162e-05] | 0.0580580 / 0.0601957 | < 1e-9 |
| joint | [2.2513796e-06, 1.6497358e-05] | 0.0029394 / 0.0156557 | < 1e-9 |

These are **not confidence intervals on actual leakage**. They condition on one local fitted observable vector; other fits and model classes remain possible. The joint example is a physical training-compatible model whose long-GST predictions fail. Its equivalence persists whether or not that particular vector fits the hardware. A separate two-state classical dormant flag reproduces the same channel recursion and all terminal probabilities. No arbitrary circuit-label memory is used.

A calibrated leakage projector after one primitive measures the loss rate directly and breaks this fiber. A known-contrast monitor also works with propagated calibration uncertainty. More terminal binary words, another device’s leakage data, or an uncalibrated extra detector output do not supply that information.

## Ordinary adversaries and calibration failure

All controls generate counts at the 3,785 actual training-row exposures. Every qubit fit refits all qubit nuisance coordinates; generating-class fits also refit nuisances. Each mechanism has 20 simulations. A 39-replicate null calibration at the measured joint qubit fit is deliberately compared with a separate-pilot, nuisance-adapted 39-replicate calibration. One start is used for these exploratory simulations; local failure is never an exclusion certificate.

| Injected ordinary mechanism | Fixed-reference exceedances / 20 | Pilot-adapted exceedances / 20 |
|---|---:|---:|
| alternating_memory | 20 | 20 |
| imperfect_reset | 11 | 1 |
| leakage | 20 | 20 |
| quasistatic_drift | 20 | 20 |
| readout_bias | 10 | 0 |

The fixed-reference threshold is not portable across SPAM nuisance values: it produces false alarms even for imperfect reset and readout that are inside the fitted qubit class. The pilot-adapted column is a conditional empirical sensitivity estimate, not a uniformly calibrated test. With only 20 trials, 0/20 or 20/20 has a two-sided exact 95% interval approximately [0,.168] or [.832,1]; intermediate rates have similarly substantial uncertainty. Generating-class fit scores and convergence flags are retained, so residuals are not classified as new physics.

A further 20 simulations share a coherent sign across hypothetical blocks of 64 file rows. Their qubit deviance range is 13132.5–20846.3. These blocks are an adversary, not reconstructed experimental jobs. Marginal drift-fit scores are retained; their independent-shot likelihood is not a correct joint likelihood for these correlated blocks.

The weakest structural direction is exact: along the leakage fiber the count distributions are identical under the same sampling law, so the power of any terminal-count test is at most its size at every exposure. A synthetic pair with different loss rates and a calibrated monitor explicitly verifies recovery of the loss rate from the added observation. No known-nuisance-only recovery is offered as evidence of terminal identification.

## Dependence and the stopping boundary

There are 8 words measured in both files. A simultaneous shared-word context-cost certificate gives lower bound 0.0 under its iid-within-pool assumption. It does not force context dependence, leakage, or memory.

No timestamps, independent job blocks, or reset logs are deposited in these selected records. Huge total exposure does not establish independence. With arbitrary dependence within the acquisition, no nontrivial confidence statement is claimed; the certified statistical lower bound is zero. The deterministic channel equalities and the conditional physical interval survive because they do not use sampling independence.

**What is proved impossible:** distinguishing the declared leakage-fiber points or their two-state dormant-flag encodings using any sequence of the specified gates, reset, and binary terminal readout alone. **What is merely not found:** an adequate frozen long-sequence fit in the four restricted numerical families, a globally necessary resource cost, and a matched leakage monitor for this acquisition. General CPTP fits, coherent leakage, richer bounded memory, and real acquisition drift are not ruled out.

A second exact ambiguity equates the alternating bit with an external two-phase slot clock under an assumed one-slot-per-primitive schedule. The records do not establish that schedule. A calibrated identity wait advancing the clock but not the gate-count flag separates the constructed pair; an unknown idle would not.

## Reproduction limits

Offline checks independently verify Kraus/Choi physics, synthetic and fitted equivalence witnesses, every saved count score, source-removal membership, and this report. They use committed reduced counts and saved predictions. The separate full-source route re-downloads hash-pinned bytes, reconstructs words/denominators, recomputes all primary and bank predictions, verifies independent density-matrix propagation, and regenerates source geometry. Numerical refits and simulation calibrations are explicit additional commands and are not promised bitwise unique.
