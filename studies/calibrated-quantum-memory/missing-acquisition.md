# Exact acquisition needed to test the declared cut

**Scope:** adapt the q0 February 21 storage acquisition to a middle
measure/reprepare intervention, or recover equivalent already-acquired records.
The inspected deposits do not verify such records. This is a requested-data
specification, not a request sent to authors and not a claim of impossibility.
For NMN, historical pulse/job and regrouping records must additionally resolve
the conditional semantics in [the audit](acquisition-semantics.md).

## Required settings and sufficient raw record schema

Freeze one physical setting y: measure Z at the cut and prepare a fixed output
state sigma_b for each outcome b in {0,1}; no averaging away of b. As a concrete
reference one may choose sigma_0=|0><0| and sigma_1=|1><1| (ideal dephasing),
or prepare |0><0| in both branches with documented conditional correction.
Either flagged reference is EB. Use exactly the same pulse/control bundle in
calibration and dynamics. This choice specifies a new protocol; it does not
assert that existing terminal QPT performed it.

| Family | Required input/output settings | Retained counts | Why needed |
|---|---|---|---|
| Isolated intermediate instrument calibration | Six inputs X±,Y±,Z±; three terminal Pauli axes x,y,z: **18 settings per epoch**, for the one fixed y | For each setting retain n(b,c) for all four b,c pairs and total N, including zeros; 72 cells | Determines each trace-nonincreasing output map, including postmeasurement state and input coherences. Existing g/e IQ calibration supplies neither. |
| Two-interval dynamics with same middle instrument | Same six initial states; same three terminal axes: **18 settings per delay pair**, with fixed intervention and recovery | All four b,c counts, N; six matching-axis success counts are derived from these same shots | Provides F6 and complete retained terminal/flag tables for joint feasibility and source-only fits. Original four-input terminal QPT has no intervening flag. |
| Independent SPAM and classifier validation | Probe-state and terminal-effect anchors covering x/y/z and all declared levels, or a justified external error bound | Raw counts, independent classifier validation IQ, frozen threshold/training split | Separates state preparation, axis rotation and assignment errors. Z assignment alone is not trusted x/y tomography. |
| Leakage and reset/acceptance records | All distinguishable noncomputational outputs and each available reset-history flag | Counts for every joint record, including failures and rejected trials | A hidden level or available flag can carry information past the nominal break. Conditioning only on accepted trials is not a TP protocol. |

If r contains more than one binary flag, replace the four-cell tables by
`|R| × 2` joint outcomes and increase M accordingly. If reliable qutrit
readout is needed, include leakage outcomes and change the model dimension;
**do not** coerce them into c=0/1 or discard them. Measurements proposed only
for diagnosis must be distinguished from those entering the exclusion.

A minimally sufficient count file has one row per outcome:

```text
source_id, epoch_id, block_id, device_id, system_id, environment_cut_id,
compiled_operation_hash, physical_setting_y, preparation_id, final_axis,
delay_before_ns, delay_after_ns, instrument_record_r, final_outcome_c,
count, attempted_shots, retained_shots, failed_shots, run_start, run_end
```

A separate settings file maps each ID to physical pulses/angles, register
wires, bit endianness, classifier version, classical control logic, and which
registers persist, are discarded, or drive later operations. Declare logical
output labels separately from physical circuit IDs. Supply native job IDs,
backend/software versions and immutable pulse/configuration hashes. Keep
unrounded integers; normalization and mitigation are analysis products.
At the binary baseline, all four counts sum to the same attempted N if no
extra failure outcome exists. If trials can fail, record that outcome and
model it or justify a prespecified missingness rule before analysis.

## Timing and transfer controls

1. Randomize/interleave all preparation and final-axis settings within short
   blocks; save ordering, per-block counts and start/end times. This tests
   context-dependent drift instead of assuming away NMN's past-marginal issue.
2. Bracket each dynamics block with **the same** instrument tomography,
   using unchanged classifier/pulse version and documented scheduling.
   Save actual runtime configuration, not only a post-run snapshot.
   Adjacent node IDs or a 16-second gap establish no numerical drift bound.
3. Freeze discriminator training before the held-out calibration shots. Save
   training and validation identifiers, thresholds, rotations and unit
   conversions. Alternatively preregister a valid selection-aware confidence
   construction; do not use in-sample confusion frequencies as fixed-rule trials.
4. During isolated instrument calibration, justify decoupling of the external
   environment. During dynamics, test/calibrate the effect of restoring its
   coupling and pulse timing. A quantitative transfer allowance d must follow
   from measurements plus a declared physical regularity model, or a justified
   hardware bound. Finite bracketing samples alone cannot bound arbitrary
   unobserved drift or instrument–environment coupling.
5. Restrict control access to the preparation label. State information should
   enter through q0 only; no later lookup of the classical preparation ID.
   Log feed-forward inputs/outputs and latency. Hardware attenuation/isolation
   or explicit extra-system models must support any leakage/coupling bound.

No universal shot count guarantees a positive result. As a **planning value**,
8000 attempted shots for each of 18 calibration and 18 dynamics settings gives
288000 trials per matched epoch/delay pair, plus SPAM/transfer controls.
This is not an empirical power estimate. For the binary baseline, intervals
for 72 calibration cells + 72 dynamics cells + 6 F6 success events give
M=150 and radius sqrt(log(2M/alpha)/(2N)) at per-family common N. Allocate
additional error budget to any estimated SPAM/transfer constraints before
looking at a margin. For input k+, success means c=0 along k; for k− it means c=1, summed over
all r. Success events are functions of retained cells; the
union bound remains valid despite their dependence. Do not reallocate freed
alpha after deleting a source.

## What becomes testable, and the completion gate

These tables make the [flagged inequality](flagged-instrument.md) evaluable:
compute a confidence-region upper bound epsilon_U on the **full flagged**
instrument distance to an EB reference; compute a certified dynamics-required
lower bound max(0,F_L−2/3); compare with calibration plus justified transfer
and SPAM allowances. A linear functional stronger than F6 may be considered,
but must be proved valid for the declared arbitrary-classical-memory null.
No failed nonconvex search constitutes a lower-bound certificate.

Before an empirical claim, deliver all of:

* A joint physical quantum model satisfying every retained probability
  constraint, CP/TP and causal constraints, including calibration and flags.
* A physical dynamics-only model fitting **all 18 joint tables** with allowed
  imperfect intervention and only classical environmental memory.
* A physical calibration-only model fitting **all 18 instrument tables**
  without the claimed environmental resource.
* Globally valid confidence coverage and a positive certified exclusion
  margin robust to the declared transfer/SPAM allowances.

For a negative identifiability result, instead construct two models fitting
**all complete observed records** with distinct memory interpretations and
state which added measurement distinguishes them. Merely showing that g/e
assignment probabilities do not characterize postmeasurement states is a
scoped observability gap, not such a complete empirical construction.

Current status: measured readout/QPT matching and schema reconstruction are
complete for the pinned subset; the intervention tables, transfer bounds,
joint feasibility and source-only witnesses required above remain unverified.
