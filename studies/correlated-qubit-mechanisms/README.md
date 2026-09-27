# Correlated qubit mechanisms: research checkpoint

**Status: the consequential research goal is not achieved.** This study audits
whether radiation, controlled phonon injection and mechanical interventions can
jointly predict an intervention without fitting its outcome. It reproduces
selected measured summaries, certifies a narrow mathematical incompatibility
in a proposed rate model, and derives sharp observation bounds. It does **not**
establish empirical combination gain, a new mechanism exclusion, a transferable
poisoning law, or a non-ionizing background lower bound.

The gate decision is to stop the proposed cross-apparatus joint fit at present.
The unresolved requirements are a calibrated transfer from injection to the
radiation response, an event/readout likelihood with run-level uncertainty, and
a response measurement that distinguishes tunneling intensity from parity
contrast. These are requirements for this approach, not an impossibility result
for the underlying physics or all available methods.

## Primary target: intervention transfer

The [forcing-to-response test](transfer.md) now uses short PT on/off acquisitions
with matched acceleration and trigger records. It predicts native-1-ms shock
response, onset lags and flagged-run duration without fitting shock outcomes.
The simple transfer fails: 32 paired shock flags are observed, versus about
648,620 under a linear response and 1,073 under a saturating extension. Both
extensions are exactly identical on all retained PT input trajectories, but
neither is an adequate physical countermodel. After a separately labeled pre-trigger baseline diagnostic, 74.86% of shock
proxy values lie outside the PT calibration range, and its acceleration reference is a separate
run from the 50 qubit acquisitions.

This is a scoped negative transfer test. It supersedes further optimization of
the small history/phase benchmark gain as the primary research direction. A
transferable observation model and concrete combination payoff remain open.

## Qubit-history and intervention follow-up

The [matched information comparison](complementarity.md) tests qubit-history-only,
phase-only and combined predictors on the same retrospective windows. History
alone predicts 2,912 paired events against 2,927 observed. Adding deposited phase
improves log loss by only 0.39%, with mixed block performance and unverified causal
phase availability. Mechanical information is not established as indispensable.

A shared latent-environment window model predicts held-out count and coarse onset
structure more accurately than the tested additional joint-event component.
These are observed mark models, not full physical compatibility witnesses. The
control audit also finds byte-identical calibration references and a common
native-1-ms quadrature flag across off/on/shock records. The earlier 99 µs sampling
mismatch does not rule out useful coarser intervention work.

The next goal is a meaningful outcome requiring the combined information, or
empirically adequate competing physical models demonstrating nonidentification.
This checkpoint does not claim that goal is complete.

## Mechanical prediction follow-up

The [run-preserving mechanical analysis](mechanical-prediction.md) now uses
3,100 paired dwell records and their synchronized phase offsets. A rule fitted
to individual-qubit training frequencies predicts 48,528 joint persistent
excursions in separate validation acquisitions; 56,438 occur. Its fixed quiet
interval contains 14 events in 3.42 million eligible windows. The apparent
574-fold rate reduction is retrospective selection of **observed readout events**,
not a demonstrated physical mitigation or an online scheduling result.

The phase-only point forecast remains imperfect, and the off/on and controlled-
shock records use 1 ms sampling rather than the required 3 µs. We do not bridge
that observation change with fitted free parameters. The report distinguishes
selected readout response, common-impact physics and causal copper attribution;
none is silently substituted for another. The scientific completion gate remains
open. See the recorded [protocol](mechanical-protocol.json) and
[results](results/mechanical-prediction.json).

## Matched gamma revision

The [matched gamma follow-up](matched-gamma.md) recovers all 12 footprint count
denominators from the paper, corrects the direct use of net-parity bounds for
an any-switch readout, and tests the proposed multiple-window route. The
released windows do not sample a calibrated common latent population. Under
explicit observation assumptions, four distant pairs retain a contrast gap
through a deterministic error envelope; the strongest needs combined charge/
mask retention above **70.2438%** to imply a common-population difference. That
retention is not established. This is a quantified calibration requirement,
not the requested full-record mechanism separation.

## What is actually reproduced

- The complete Iaia and Larson figure archives were downloaded and every member
  hashed. The 3.1 GB Kono archive's 343-entry directory and 11 selected members
  were downloaded by exact byte ranges. Its full archive hash was **not** checked.
- Iaia's observed two-qubit coincidence suppression ratios are 81.25, 75 and
  114; the three-qubit ratio is 200. These reproduce figure values, not a new
  event-level analysis or a confidence statement.
- Kono's phase-aggregated data contain 12,697,600 observations. Pooled plug-in
  mutual information is 0.0027806635 bits. Weighted phase-conditional plug-in
  mutual information is 0.0000880256 bits. Between-phase modulation accounts
  for 84.1751% of the **covariance**, a different statistic. This is consistent
  with the original interpretation; it does not identify the microscopic cause
  or demonstrate conditional independence.
- The proposed affine rate model fails an envelope of four supplied error bars
  for every one of the nine Larson parity curves. Exact rational contrasts
  certify required error multipliers between 30.58 and 194.87. The two charge
  curves require only 1.42 and 1.92. These multipliers are **not sigma values**:
  the deposited error columns are not a verified full uncertainty model.
- In a mixed-Poisson tunneling model, a background-corrected parity contrast
  of 0.20 permits an any-tunneling probability between 0.10557 and 0.20 in
  closure. A scalar observation does not select a unique point in that interval.
  This is a model-conditional observation bound, not a new interpretation
  certified against all of the radiation experiment's data.

The precise gate, alternatives, and calibration requirements are in
[gate.md](gate.md). Definitions and proofs are in [model.md](model.md).
The [literature comparison](literature.md) identifies the closest prior work
and separates fitted agreement from held-out prediction. Exact results and
witness coefficients are in [analysis.json](results/analysis.json).

## Reproduce

From the repository root, with Python 3.12:

```bash
python -m pip install -r studies/correlated-qubit-mechanisms/requirements.txt
python studies/correlated-qubit-mechanisms/cqm_sources.py --cache /tmp/correlated-qubit-data
python studies/correlated-qubit-mechanisms/cqm_analysis.py --cache /tmp/correlated-qubit-data --check
python -m unittest discover -s tests -p 'test_correlated_qubit_mechanisms.py' -v
python scripts/check_repository.py
```

Routine unit tests need no network or external data. The source reproduction
route downloads approximately 18 MB, using member ranges instead of the entire
mechanical archive. Full-archive verification is separately available, but was
not run in this investigation:

```bash
python studies/correlated-qubit-mechanisms/cqm_sources.py --cache /tmp/correlated-qubit-data --full-mechanical
```

Use `--verify-only` to check existing files without downloading. Without
`--check`, the analysis regenerates its result file. The check compares exact
certificate strings exactly, and floating-point results with explicit tolerances.
The source acquisition verifies SHA-256, lengths, ZIP membership and member
hashes. It does not run source notebooks; the one analyzed pickle contains only
primitive lists and integers, and loading arbitrary classes is rejected.

[manifest.json](manifest.json) pins records, licenses, original filenames,
lengths, hashes and acquisition status. The three measured archives are CC BY
4.0 and retain their authors' rights. Source files remain external to git.
The [mechanical inventory](mechanical-inventory.json) is an archive directory,
not a claim to have inspected all 343 payloads.

The reported measurement reconstructions are deterministic. Twenty scientific and
integrity tests check Poisson sums, sharp extremizers, heterogeneous distributions,
independent-background correction, global affine certificates, phase mixing,
corruption rejection, pre-conversion integer/range/exposure validation, rational
upper residuals, background-path enumeration, selection extremizers, dwell
partition validation, a direct counting oracle, and protection against fitting
the joint validation target. No full G4CMP, microscopic transport fit, calibrated
coverage simulation, held-out intervention test or QEC simulation was performed.

## Reproduce the mechanical prediction

This optional route acquires ten additional selected members (about 68 MB
compressed), including the 393 MB uncompressed dwell record. It does not fetch
the full 3.1 GB archive. A numeric-only restricted unpickler and primitive scalar
conversion avoid importing arbitrary pickle classes and reduce memory use.
The source hashes are verified before deserialization. Allow roughly 2 GB RAM.

```bash
python studies/correlated-qubit-mechanisms/cqm_sources.py --cache /tmp/correlated-qubit-data --mechanical-prediction
python studies/correlated-qubit-mechanisms/cqm_mechanical.py --cache /tmp/correlated-qubit-data --check
python -m unittest discover -s tests -p 'test_cqm_mechanical.py' -v
```

The original `cqm_analysis.py --check` still reproduces the earlier gate without
requiring the additional dwell data. The new route processes every acquisition,
retains its block assignment, and verifies the separate phase-offset and
intervention-clock files. Refinement and benchmark diagnostics are explicitly
labeled as post-validation analyses; no nominal confidence level is attached.

## Reproduce the information and control comparisons

Install the pinned study requirements, including scikit-learn 1.8.0. The information
comparison uses the same dwell sources as above; feature extraction starts afresh
on every invocation and writes a reconstructible scratch cache outside the repo.
Allow roughly 3 GB RAM for extraction and fitting. Set thread counts as below for
the reference reproduction.

```bash
python -m pip install -r studies/correlated-qubit-mechanisms/requirements.txt
OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 python studies/correlated-qubit-mechanisms/cqm_complementarity.py --cache /tmp/correlated-qubit-data --features /tmp/cqm-history-features.npz --check
python studies/correlated-qubit-mechanisms/cqm_sources.py --cache /tmp/correlated-qubit-data --mechanical-prediction --mechanical-interventions
python studies/correlated-qubit-mechanisms/cqm_interventions.py --cache /tmp/correlated-qubit-data --check
python -m unittest discover -s tests -p 'test_cqm_complementarity.py' -v
```

The control route adds 11 selected payloads (32.4 MB compressed). The new tests
check future-data isolation, marks against a direct oracle, valid normalized
mixture laws, and acquisition-safe intervention counting. Neither forecast nor
control script attaches IID uncertainty to the recorded windows.

## Reproduce the intervention-transfer test

Six additional selected members provide the short matched PT records (about
1.23 MB compressed). The transfer route also verifies the preceding prediction
and intervention sources. No authors' notebook or helper is executed.

```bash
python studies/correlated-qubit-mechanisms/cqm_sources.py --cache /tmp/correlated-qubit-data --mechanical-transfer
OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 python studies/correlated-qubit-mechanisms/cqm_transfer.py --cache /tmp/correlated-qubit-data --check
python -m unittest discover -s tests -p 'test_cqm_transfer.py' -v
```

The tests compare transition matrices to an independent matrix exponential,
check causal forcing access, enumerate complete paths to verify duration
expectations, and prevent cross-acquisition onsets. The model's flag rates are
not reported as independently calibrated intrinsic qubit rates.
