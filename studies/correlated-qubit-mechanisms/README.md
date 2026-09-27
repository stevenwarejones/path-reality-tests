# Correlated qubit mechanisms: feasibility gate

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

The reported measurement reconstructions are deterministic. Eight scientific and
integrity tests check Poisson sums, sharp extremizers, heterogeneous distributions,
independent-background correction, global affine certificates, phase mixing,
corruption rejection, pre-conversion integer/range/exposure validation, rational
upper residuals, background-path enumeration and selection extremizers. No full G4CMP, microscopic transport fit, calibrated
coverage simulation, held-out intervention test or QEC simulation was performed.
