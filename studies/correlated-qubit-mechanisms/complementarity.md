# Does mechanical phase add information beyond qubit history?

**A strong qubit-history baseline removes most of the apparent need for the
mechanical channel.** On a fixed subset of the same validation acquisitions,
qubit history predicts 2,912 paired events against 2,927 observed. Adding the
deposited mechanical phase improves binary log loss by 0.39%, with mixed block
results. This does not establish indispensable mechanical information, a
practical mitigation gain, or a universal identifiability result.

A separate finite-window generative comparison favors a shared fluctuating
latent environment over the particular additional joint-event component tested
here. A common observed 1 ms flag is also available in existing off/on/shock
records. Thus incompatible sampling of the original 99 µs target is **not** a
reason to conclude that every useful intervention comparison needs new data.

The scientific target remains: predict an operationally meaningful outcome
requiring combined information, or construct empirically adequate competing
physical models explaining why that outcome remains unidentified. Neither has
been completed. Radiation and injection measurements still do not constrain
these forecasts jointly. The copper loss-channel construction remains
conditional on a feasible starting model, not a fit to all measured records.

## A stronger, matched retrospective comparison

The [protocol](complementarity-protocol.json) was recorded before computing these
new outputs. This is a retrospective extension after seeing the earlier results,
not preregistration or newly blind validation. We retain acquisitions 0–1549 for
training and 1550–3099 for validation; no hyperparameters are selected using the
validation outcomes.

All predictors use exactly the same eligible windows: start at sample 250,000
(0.75 s), then every 330 samples (0.99 ms), with both qubits currently classified
G. The target remains a persistent excursion in **both** qubits during the next
33 readouts (99 µs), requiring three consecutive non-G classifications per qubit.
This fixed sparse grid permits a richer model comparison with 1,253,978 training
and 1,258,687 validation windows. Its counts are not replacements for the previous
all-window counts of 48,528 predicted and 56,438 observed.

Each qubit-history feature uses only classifications at or before the window
start, within the current acquisition. Features include non-G fractions over
approximately 0.1, 1, 10 and 100 ms; fifteen preceding approximately 50 ms segments;
and capped time since the last non-G classification. These allow the history
model to learn periodic activity over more than one mechanical cycle without
receiving mechanical offsets or an explicit mechanical-period feature. All
models receive elapsed time in the acquisition. The phase models additionally
receive the sine and cosine of the deposited phase.

The three nonconstant predictors use the same histogram gradient-boosting
classifier, tree budget, regularization and training outcomes. They all fit the
**paired** target, unlike the earlier product-of-marginals rule. This is a
predictive benchmark comparison, not a mechanism-derived forecast. Model
settings, feature definitions and extraction digests are in the
[output](results/complementarity.json).

| Predictor | Predicted count | Binary log loss | Brier score |
|---|---:|---:|---:|
| Constant training joint frequency | 3,006.25 | 0.01642468 | 0.002320035 |
| Qubit history only | 2,911.68 | 0.00993937 | 0.002212225 |
| Deposited phase and elapsed time | 3,096.56 | 0.00995860 | 0.002207522 |
| Qubit history plus deposited phase | 2,812.06 | **0.00990014** | 0.002209225 |
| Observed | **2,927** | — | — |

Lower scores are better. The combined model improves log loss relative to
history alone by 0.00003922 per eligible window, or **0.39%**. It improves on phase
alone by 0.59%, but the phase model has the better Brier score. Combined log loss
is better than history alone in 21 of 31 contiguous 50-acquisition blocks, and
better than phase alone in 17 of 31. These are paired descriptive comparisons;
block independence, nominal confidence and future generalization are not
established. No family-wide optimality or exclusion is claimed for the fixed
learner and finite history representation.

For each predictor, a selection threshold is the 15th percentile of its training
scores. We apply that fixed threshold to validation, including tied scores.
Coverage therefore need not equal 15%, and differs between predictors.

| Predictor | Selected validation windows | Coverage | Observed paired events | Predicted events |
|---|---:|---:|---:|---:|
| History only | 236,513 | 18.79% | 0 | 1.58 |
| Phase | 188,585 | 14.98% | 2 | 0.29 |
| Combined | 224,044 | 17.80% | 1 | 1.43 |

These small counts and unequal coverages do not rank prospective controllers.
Zero observed events does not imply zero risk. They do show that the original
quiet-selection result cannot establish that a mechanical channel is necessary:
a substantive qubit-history alternative also finds a quiet observed subset.
The constant model assigns one score to all windows and retains all under the
tied-threshold rule; it supplies no selective policy.

## Causal availability of mechanical phase

[Kono et al.](https://arxiv.org/html/2305.02591v1), sections II.3–II.4, describe
aligning recorded vibration traces to a reference by maximizing cross-correlation,
and applying analogous alignment to continuous-monitoring traces. The inspected
dwell notebook loads completed offsets; it does not reconstruct their estimator
or certify that each offset used only samples preceding our forecast window.
The separately stored `lags.npy` agrees entry-for-entry, which checks association
but not causality or phase uncertainty.

Consequently, history features pass a strict past-only test within the released
classification stream; this does not certify the original classifier's latency.
The phase comparison uses **deposited retrospective information of unverified causal
availability**. It is not an equal-online-information competition. It tests the
value of those deposited phase labels beyond observed past qubit history. The
small gain cannot establish practical predictive complementarity, even with
this advantage to phase-assisted models.

The control archive has a separate triggered acceleration record (`Fig6/832.npy`),
but it is not the missing per-acquisition accelerometer stream for the 3,100 dwell
records. We do not substitute it. A causal phase implementation would need a
past-only estimator, known reference construction, latency, and timing-error
validation on the matched acquisitions.

## Competing finite-window generative models

We also predict structure beyond the total paired-event count. Each local future
window has seven possible marks: no persistent excursion, or one of three bins
for the first qualifying onset (readouts 1–11, 12–22, 23–31), crossed with whether
at least 12 of the 33 future readouts are non-G. The latter is **36 µs of observed
non-G occupancy**, not an uninterrupted dwell duration or a microscopic lifetime.
A full dwell-duration and whole-trace fit remains open.

A conditioning cell is one of 36 mechanical-phase bins crossed with four bins
of pooled preceding-10-ms non-G activity. Thus environmental variation is not
represented only by a fixed periodic mean. For each cell, let a and b be the
local mark probabilities. The models generate the paired marks as follows:

| Model | Probability law | Free parameters per cell |
|---|---|---:|
| Independent local marks | `a(i)b(j)` | 12 |
| Shared environmental state | `sum_z w(z) a_z(i)b_z(j)`, two latent states | 25 |
| Added joint-event component | `(1-lambda)a(i)b(j) + lambda d(i) 1[i=j>0]` | 18 |
| Unrestricted joint benchmark | Arbitrary normalized 7 × 7 table | 48 |

The shared state can produce paired excess through simultaneous changes in local
propensities, without a direct joint-event component. The alternative adds a
common event with identical **coarse marks**. That restriction is a particular
coupling hypothesis, not all possible common disturbances. Neither model
represents the latent environment's temporal evolution across windows, readout
backaction, or microscopic energy transport.

The shared state can also absorb unresolved deterministic phase variation within
the approximately 19.85 ms conditioning bins. Its success does not distinguish
such variation from a stochastic bath. These coarser model cells are separate
from the earlier fine-phase-bin sensitivity check.

All four use the same training counts and the same 100 pseudo-observations per
cell, distributed according to the global training joint table with Jeffreys
0.5 per category. The mixture models use 300 EM iterations; the shared-environment
fit chooses among three specified initializations by its **training** objective.
EM is not a global optimizer. The unrestricted table is a stronger-capacity
benchmark with identical retained smoothing, not a mechanistic explanation.

| Window model | Predicted paired count | Same onset-bin pairs | Both ≥36 µs non-G occupancy | 49-category log loss |
|---|---:|---:|---:|---:|
| Independent local | 1,783.17 | 616.36 | 350.94 | 0.15122402 |
| Shared environment | 2,871.56 | 999.03 | 552.88 | **0.15103561** |
| Added joint component | 1,982.41 | 857.48 | 459.93 | 0.15124188 |
| Joint-table benchmark | 2,930.98 | 1,006.70 | 581.78 | 0.15144348 |
| Observed | **2,927** | **988** | **508** | — |

The shared-environment fit improves validation log loss relative to the tested
joint component in 29/31 blocks and relative to the unrestricted table in 31/31.
Its onset-bin and count predictions are close, while its long-occupancy count is
8.8% high. A more flexible table predicts the overall count more closely but
has worse held-out mark log loss. This is why fitting a coincidence total alone
would be inadequate.

The result favors this shared-state predictor among the specified fits. It does
not statistically reject all joint-event mechanisms, establish a calibrated
aggregate physical fit, or furnish two models compatible with every released
constraint. In particular, the latent-state explanation is not proof of a
specific fluctuating physical environment, and correlated classification noise
remains possible. The earlier 16.3% discrepancy alone does not identify a cause.

## A common observable already present in intervention records

The [control audit](results/intervention-observable.json) verifies more than clock
agreement. Fig. 2 and Fig. 6 contain **byte-identical G/E/F calibration arrays**
617/618/619, and both notebooks project readout using PCA fitted to the same
combined G/E calibration. Fig. 6 divides by two recorded positive normalization
constants; the audit reverses those constants before using the common definition
`u=(r-mu_G)/(mu_E-mu_G)` and flagging `u>0.5`. We read the source helper as text;
we do not execute the authors' notebooks or helper code.

Two compatible *observed* quantities can therefore be counted at the native
1 ms clocks: simultaneous flags, and simultaneous flags at the next readout
conditional on both current readouts being unflagged. Adjacent-readout counts
never cross acquisition boundaries. The numbers below are descriptive, not an
intervention prediction fitted without the intervention outcome.

| Record | Acquisitions × samples | Simultaneous flags | Eligible unflagged starts | Paired flags 1 ms later |
|---|---:|---:|---:|---:|
| PT off, 774 | 1 × 32,768 | 0 | 32,623 | 0 |
| PT on, 782 | 1 × 32,768 | 44 | 31,749 | 35 |
| Controlled shock, 833 | 50 × 32,768 | 32 | 1,626,568 | 29 |

Shared reference calibration and sampling do not establish equal physical
observation response. The prepared-G reference itself flags approximately 1.32%
and 1.52% in the two channels, mixing preparation, dynamics and separation effects;
these cannot be treated as independently known false-positive rates. The saved
Fig. 2 records contain projected quadratures rather than raw IQ, and the original
projection orientation is taken as deposited. Readout duty, preparation, state-
dependent backaction, intervention timing and run drift must enter any physical
transfer model. One off record and one on record do not establish independent
randomized intervention trials. The controlled-shock acceleration file supplies
its own trigger, not a calibration bridge to every other acquisition.

**Decision:** retain this coarser native-clock observable as the next intervention
target. Do not downsample the continuous record and declare equivalence, and do
not yet conclude that new measurements are necessary. Before fitting transfer,
recover the preparation/readout settings and shock alignment sufficiently to
predict this common observed target under a withheld intervention. If those
cannot be constrained, demonstrate that specific observation ambiguity with
empirically adequate models.


## Subsequent intervention test

The [transfer analysis](transfer.md) now makes the native-1-ms intervention the
primary target. It finds same-number forcing/qubit records for shorter PT traces,
while auditing the shock reference as a separate run. A dynamic rate model and
an exactly PT-equivalent saturating extension both fail transfer. Their explicit
time and duration predictions do not make them adequate physical countermodels.
The finite-window model ranking above is therefore not promoted to a physical
mechanism claim, and further optimization of the 0.39% gain is not the priority.
