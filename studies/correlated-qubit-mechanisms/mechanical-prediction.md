# A retrospective forecast of short paired readout excursions

**The phase-conditioned point forecast improves substantially over a phase-blind
one, but does not establish a physical mechanism or a transferable intervention.**
Training on individual-qubit event frequencies predicts 48,528 paired events in
a separate half of the mechanical acquisitions; 56,438 occur. A fixed quiet-phase
policy selected without fitting the paired-event outcome contains 14 events in
3,424,600 eligible validation windows. These are observed classifications during
continuous measurement, not microscopic tunneling events or logical QEC errors.

This is a concrete prediction test after the gamma causal-attribution route was
stopped. It advances beyond pooled covariance, but the high-bar research goal
remains unmet: there is no globally certified joint physical model, complete
source-only witnesses, or full-record pair of competing physical mechanisms.

## Three different targets

| Target | What this analysis can say |
|---|---|
| Response among a defined observed population | Counts persistent readout excursions in windows starting with both qubits classified G; makes an ordered retrospective forecast |
| Response to a common incident-impact population | Not identified; the mechanical records have no impact-by-impact ionizing/non-ionizing classification or microscopic observation calibration |
| Causal effect of copper rather than other device differences | Not addressed by this apparatus; same-cooldown gamma chips do not provide a randomized copper intervention |

The gamma loss-channel construction is stated in
[matched-gamma.md](matched-gamma.md#causal-attribution-decision). No new
selection-sensitivity curve is substituted for an empirical bound on acceptance.

## Exact source and chronology

The source is [Kono et al.'s version-v0 archive](https://zenodo.org/records/11034817).
The acquisition script downloads only selected ZIP members and verifies their
pinned lengths and SHA-256 values. The full 3.1 GB ZIP remains unverified.

| Record | Contents and role |
|---|---|
| `Zenodo/Fig4FigS6FigS7/results_gef_2d` | Three consecutive pickle objects: 1000-point plotting axis, 3100 phase offsets, and two lists of 3100 dwell sequences. Each sequence contains `(state,start,end)` intervals, with states G/E/F = 0/1/2. This is processed classified data, not raw IQ. |
| `Zenodo/Fig4FigS6FigS7/lags.npy` | Independently stored 3100 offsets; every value equals the corresponding offset in the dwell file |
| `Zenodo/Fig4FigS6FigS7/Fig_Dwelltime_analysis.ipynb` | Supplies the 3 µs sample period, 0.7144833333 s mechanical period and phase-coordinate constants; read as text, not executed |
| Fig. 2/6/S9 notebooks and `774_t.npy`, `782_t.npy`, `833_t.npy`, `321_t.npy` | Establish the 1 ms readout clocks for the pulse-tube-off/on and controlled-shock candidates. These cannot directly score the 99 µs target. |

Every deposited dwell sequence partitions `[0,524287)` without gaps or overlap.
We use those half-open boundaries and omit the terminal boundary as a new sample.
The acquisition index is preserved throughout. Absolute wall timestamps and a
validated independent-block length are not supplied by these files. An array
index is not promoted into a wall-clock timestamp.

These are **different acquisitions** from the old 714-bin residual-population
summary. We do not join their event timestamps. The new core is the paired
continuous-readout record plus its synchronized, accelerometer-derived offsets;
these are dependent measurement channels from the same acquisitions, not
independent experiments. The original paper's analog readout, phase dependence
and correlation results remain prior art.

We use the deposited phase offsets and verify their one-to-one correspondence
with the dwell records. Their extraction from raw accelerometer traces is not
independently reconstructed here; agreement between the two stored copies does
not calibrate phase-estimation error.

## Target and fixed prediction rule

The [protocol](mechanical-protocol.json) was recorded locally before computing
this target. It was not externally preregistered; source figures and related
outcomes had already been inspected. The validation is retrospective.

1. Use acquisitions 0–1549 for calibration and 1550–3099 for validation. Preserve
   contiguous groups of 50 acquisitions for variation diagnostics.
2. Start disjoint 33-interval windows on the deposited 3 µs clock. A window is
   eligible when both qubits are classified G at its first readout.
3. For each qubit, count an excursion when at least three consecutive readouts
   among the following 33 are non-G. An E-to-F change does not interrupt the
   run. This persistence rule suppresses isolated classification spikes; it is
   **not** a calibrated correction for all readout errors.
4. Convert each window's start to mechanical phase using the deposited offset
   and the notebook constants. Divide the period into 715 bins, about 0.9993 ms
   each. This phase-bin width differs from the 99 µs target window.
5. Fit each qubit's marginal excursion frequency conditional on its own G start
   in each training phase bin. Do not fit the paired-excursion counts.

The tested observation model assumes independence of the two local readout
processes conditional on mechanical phase, with stable phase response between
calibration and validation. Local temporal memory is not fitted as an exponential
dwell law. If `p_A(j)` and `p_B(j)` are the two fitted conditional frequencies,
the predicted paired probability is `p_A(j) p_B(j)`. Conditional independence
of complete local paths would justify this product after conditioning on both
G starts; unresolved common phase/amplitude variation can invalidate it.

Select one contiguous 100-bin quiet interval minimizing the predicted product,
weighted by training both-G exposure. Initial-state eligibility is used;
training paired future outcomes and all validation outcomes are excluded from
policy selection. The selected interval is **50.963–150.891 ms** in the deposited
phase convention, with duration 99.928 ms, or 13.986% of the cycle. Validation
predictions condition on its measured initial-state eligibility exposure.

## Validation results and comparators

| Predictor | All phases: predicted count | Fixed quiet policy: predicted count |
|---|---:|---:|
| Independent, phase blind | 5,377.84 | 766.26 |
| Independent, phase conditioned; primary rule | 48,527.96 | 7.11 |
| Phase-conditioned benchmark fitted to training **paired outcomes** | 58,148.09 | 13.05 |
| Observed validation count | **56,438** | **14** |
| Eligible validation windows | 24,034,855 | 3,424,600 |

The primary all-phase point forecast underpredicts by 16.3% relative to its
prediction. Phase information corrects much of the phase-blind model's error,
but does not remove the residual. A benchmark that fits the target itself does
better and is included to avoid overstating the mechanism-derived forecast.
The phase-blind comparison is an **algorithmic ablation**, not proof that every
analysis of the full qubit-only source is weaker: those records may themselves
contain information about the periodic drive.

The measured validation rates are 2348.17 per million eligible windows over all
phases and 4.088 per million in the fixed quiet interval, a **574.4-fold observed
rate ratio**, retaining 14.248% of eligible windows. This is retrospective phase
selection in an unchanged measurement stream. It is not an experimentally
implemented intervention, a demonstrated online controller, or an intrinsic
error-suppression factor.

## Uncertainty and adversarial checks

There is no IID-shot confidence interval or nominal future-coverage claim.
Windows share boundary readouts, sequences contain memory, and consecutive
acquisitions can share slow drift. Choosing 50-acquisition blocks preserves
some dependence but does not prove that those blocks are independent.
[The output](results/mechanical-prediction.json) reports every validation block.

- The all-phase observed/predicted ratio ranges from 0.767 to 1.608 across the
  31 validation blocks. Thus the pooled 16.3% discrepancy is not a certified
  rejection of the full phase-driven physical model family.
- Quiet-policy validation block rates range from 0 to 35.69 per million;
  the training range is 0 to 26.72. The training block range therefore does
  **not** cover every validation block. Only nine validation blocks have a
  nonzero quiet-policy count. Small counts do not justify a precise future rate.
- A post-validation fivefold refinement to 3575 phase bins predicts 48,624.38
  all-phase events, still 16.1% below the observed count relative to prediction.
  The same fixed quiet interval still contains 14 events. Coarse binning alone
  does not explain the point-forecast mismatch.
- Refitting only the marginal frequencies within 50-acquisition validation
  blocks predicts 48,207.83 paired events. This post-validation diagnostic still
  leaves a 17.1% observed/predicted excess. It is not a held-out test or proof
  against every form of run variation.
- Shifting the fixed policy by approximately ±20 ms gives between 11 and 22
  observed events among the seven tested offsets. This is an exploratory clock
  sensitivity check, not a measured bound on online timing error.

The common excess could involve unresolved mechanical variation, common physical
bursts, correlated readout disturbances, selection by the initial state, or a
combination. We do not assign it to QPs, TLSs or radiation. The pulse-tube-on
record alone cannot validate the counterfactual with the pulse tube off.

## Transfer decision

The deposited off/on and controlled-shock clock vectors have 32,768 samples
separated by **1 ms**, rather than 3 µs. Their readout duty cycle and backaction
also differ. Fitting a free observation conversion between these protocols would
remove the very predictive restriction being tested. We therefore do not label
this a successful withheld-intervention prediction.

Likewise, the deposited offsets align records retrospectively. Their availability
does not establish that an equivalent phase estimate can be generated causally
before each proposed quiet window. An online phase estimator would need its own
validation. The measured readout-classification target cannot yet be translated
into free-evolution or QEC errors.

The next useful control is now specific: matched 3 µs acquisition and state-readout
calibration during pulse-tube-on, pulse-tube-off and controlled-shock trials,
with a causal phase estimate logged before the target window. Rare-event readout
controls must resolve a rate of a few events per million eligible windows if the
quiet-policy target is retained. The approximately 7 GHz microwave monitor from
the original work is not, by itself, a calibration of rare joint classification
errors or all noise near the qubit transition frequencies.

This investigation stops short of a full physical exclusion. The surviving
mechanisms have not been fitted to every retained source constraint, so no
source-removal witness or general impossibility theorem is claimed. The concrete
result is a reproducible retrospective forecast, its residual, and an identified
protocol mismatch that blocks the proposed intervention transfer.
