# Feasibility on the actual local timing background

This follow-up measures recovery of injected effects on the pinned NIST Bob
record. It uses **fresh artificial assignments** and does not compare the actual
remote settings or outcomes. The [report](empirical-report.md) and
[results](empirical-results.json) are distinct from the original Gaussian study.

## Question and scope

Can frozen features recover timing perturbations against the empirical peak,
tails, count sparsity and chronological variation? The null experiment randomly
assigns independent fair bits to a *fixed* local background. Actual local device
memory is retained in that background, but the artificial assignments cannot have
caused it. Conditional null calibration is therefore valid without claiming that
the real experimental assignments satisfy the same assumptions.

This is conditional algorithmic sensitivity. It does not validate the real
random number generators, independence of experimental training and test data,
physical acquisition completeness, spacelike timing, or absence of ordinary
leakage. It supplies no p-value for real remote-setting dependence.

## Sources and extraction

The [existing manifest](../nist-bell-causal-audit/manifest.json) pins the HDF5 and
Bob ZIP. The extractor verifies both full hashes before use, discovers the local
overlap offset, and checks **every archived Bob setting and legacy click word**
against the raw decoder. It verifies the complete expected local timestamp-jump
list. No actual Alice setting or outcome is accessed by the analysis.

The five-value decoder API is unchanged by default. An optional `include_events`
argument exposes detector row indices, pulse numbers and phases computed by the
same period-correction code. This avoids a second timing implementation.

The local record is defined before artificial assignments: detector events in
baseline pulse bits 4, 5 and 6, with **no phase-radius cut**. All such events are
retained in an external NPZ cache. All valid rows without such events remain
no-event trials. Raw prefixes, suffixes and the final open interval are recorded.

Local timestamp-excursion support (including the following guard interval) and
ambiguous local settings are excluded *before* artificial assignment; counts are
reported. For Bob this removes 2,236 of 107,109,596 overlap rows and two selected
events. This differs deliberately from the causal audit's worst-case completion:
the target here is recovery on a declared retained population, not a bound for
the complete physical experiment. Unknown lost physical trials remain unknown.

The retained background has 31,532 selected events, each in a different trial,
and 107,075,828 no-event trials. The extractor supports multiclicks and keeps them
grouped. The score uses the first selected event in original record order, so
multiple photons could never count as independent assignment trials. That feature
is only one aspect of the local record, not an exhaustive full-record test.

## Frozen population and injections

Rows are divided into 32 chronological blocks by `floor(32*row/overlap_rows)`.
The first eight train; the remaining twenty-four test. Features are learned
separately for local setting × block parity (four predefined strata). The
alternating parity regime is deliberately supplied to match the injected pattern;
we have not discovered a comparable regime in real remote-setting data.

Each repetition draws one independent fair artificial bit per event-bearing
trial. Assignment totals for the remaining no-event trials are sampled exactly
by binomial draws. No-event trials stay in all exposure denominators. These
compressed assignments are equivalent to independent fair bits over every
retained row. All injections in a repetition share those bits for paired
comparisons; repetitions are independent artificial randomizations of the **same
physical background**, not new physical experiments.

The injection displaces copied detector timestamps by a common offset per trial.
Sync timestamps and the original pulse/trial membership stay fixed. Every
selected event survives; no new event is selected from outside the baseline pulse
population. Bins have both overflow tails, so no injected event disappears at a
histogram edge. This frozen-population model does not reproduce how a physical
shift might change a hardware gate, reassignment at a pulse boundary, detector
response or photon production. Click-count contrasts remain unchanged by design.

Three displacement models separate an important digitization ambiguity:

| Model | Perturbation | Interpretation |
|---|---|---|
| Continuous | Arm means move by -delta/2 and +delta/2 in decoded phase | Numerical sensitivity; fractional shifts of stored values can exploit discrete mass near thresholds |
| Digital | Symmetric whole-tag offsets; arm difference 2 or 4 bins | Exact integer edits to copied time tags, with no latent sub-bin model |
| Redigitized | Uniform latent u in [-1/2,1/2); integer increment floor(u + displacement + 1/2) | An explicit hypothetical rounding model, not a measured within-bin distribution |
| Centered latent | Original true times at tag-bin centers; round after the shift | An adversarial counterexample: small shifts can leave every recorded tag unchanged |

A tag bin is 78.125 ps. Redigitization uses the same latent positions across
injections within a repetition; they are independent of the artificial labels.
Sub-bin recovery depends on this uncalibrated latent distribution. Neither a
high continuous nor a high redigitized recovery rate establishes hardware timing
resolution or a picosecond physical exclusion limit.
The centered-latent model is not asserted to describe this apparatus; it shows
why the measured integer tags alone do not identify sub-bin response to a
physical perturbation. It supplies the opposite extreme to a uniform latent law.

Fixed and block-parity-reversing versions are included. A null local-clock
control moves both arms equally by ±1 bin according to block parity. A transfer
failure injects the same direction into both training regimes but opposite
directions at test. Real regime weights need not balance exactly; pooled
cancellation or total transfer failure is therefore not guaranteed.

## Statistics and reproducibility

The four scores are the same types as the Gaussian study: click count, pooled
early/late timing relative to the configured peak, trained timing signs and
trained histogram signs. Training alone determines the learned weights. All
tests use the same final three quarters and exact two-sided fair-sign p-values.
Bonferroni allocates 0.01/4 to each test. Under the artificial-label null, fixed
records and frozen features make all nonzero assignment scores independent fair
signs. This proof applies to these artificial bits, not to the original bits.

Each scenario has 400 randomizations. Clopper-Pearson intervals are pointwise
Monte Carlo intervals for recovery conditional on the one retained background
(and the specified latent model where relevant). They do not account for
between-run variability or for changing an injection model after seeing results.
This entire study is exploratory, not preregistered. Null non-rejection does not
establish physical no-signaling or certify the real assignment process.

```sh
python -m pip install -r studies/synthetic-timing-shifts/requirements.txt
python studies/nist-bell-causal-audit/fetch_data.py \
  --directory /tmp/nist-inputs --keys hdf5 bob
python studies/synthetic-timing-shifts/empirical_extract.py \
  --hdf5 /tmp/nist-inputs/hdf5.hdf5 --archive /tmp/nist-inputs/bob.zip \
  --output /tmp/nist-output/bob-timing.npz
python studies/synthetic-timing-shifts/empirical.py \
  --background /tmp/nist-output/bob-timing.npz --check
```

The full-source CI regenerates the external cache directly from both pinned
inputs and reproduces the committed results and report. Offline tests cover
population conservation, no-click exposures, grouped multiclicks, integer
rounding, frozen features and compatibility of the optional decoder output.
The source ZIP/HDF5 and event cache are never committed. Only derived aggregate
descriptions and injection-recovery summaries are published.

## Decision gate

The committed run supports continued method development. Whole-tag reversing
injections are recovered in 400/400 repetitions by the trained tests while the
pooled early/late score misses them. Under uniform-latent redigitization, trained
timing recovers the reversing 0.25-bin arm difference in 398/400 repetitions
(99.5%); trained histogram recovers 400/400. At 0.125 bin those rates fall to
205/400 (51.25%) and 182/400 (45.5%). Thus realistic rounding substantially reduces
the apparently effortless sensitivity of direct floating-point shifts.

The centered-latent version of the 0.25-bin shift produces **exactly the null
records**, and the training-to-test transfer-failure case is also undetected in
this run. Together these give concrete limits: recoverability depends on the
unobserved timing response and on a usable persistent feature. The observed
null any-test rate is 3/400; the local-clock control is 0/400. Neither establishes
physical no-signaling or precise experimental error calibration.

Proceed to real-label inference only after an explicit assignment/dependence
model and its calibration have been justified. First extend recovery checks to
another receiver/run and features selected without the deliberately supplied
reversal schedule. Acquisition failures, pulse-boundary selection and time-tag
rounding need their own sensitivity envelopes. This study can justify further
method development; it cannot yet justify a new physical exclusion claim.
