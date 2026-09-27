# Shared quantum dynamics across measured sequence families

**Can one physically consistent ordinary model predict different experimental
sequences, and do the records identify a necessary cost in leakage or memory?**

This study combines **4,657 GST sequences and 2,004 randomized-benchmarking records**
from the same interleaved 30 March 2015 trapped-ion acquisition. It preserves
[PR #11's table-level baseline](../born-rule-identifiability/README.md) and does not
share Nairobi apparatus parameters across epochs.

The [general-CPTP escalation](results/escalation.md) adds arbitrary stationary
qubit channels, five optimizer-independent outer prediction bounds, and an exact
leakage equivalence with gate-dependent transfer, nonunital dynamics and unknown
return-state polarization. An exact lift now attaches zero-leakage and positive-leakage
realizations to the **same compatible all-acquisition qubit prediction vector**,
with a sharp CP budget for this construction. Its [proof and limits](escalation-theory.md) separate
certified bounds from local fitting failures. Strict improvement of the complete
GST-only prediction region and a globally necessary resource cost remain unresolved.

The original result has two parts. The measured combination rules out explicit physical
models that either family alone permits. But in a declared shared-gate leakage
family, a continuous change in loss, return, leakage readout and depolarization
leaves **every** binary terminal sequence probability unchanged. An exact formula
gives the sharp conditional physical interval, and a two-state dormant classical
flag gives the same probabilities. This is a bounded, operational obstruction;
its resolution is a calibrated leakage-sensitive observation, not more of the
same binary terminal sequences.

**Limit:** the restricted and general-channel local training fits still miss long GST observations; budget-limited optimization is not an exclusion certificate. An all-length general qubit fit lies inside the declared conservative joint
count region, so that region does not require positive leakage or memory. No
validated complete device explanation, globally necessary nonzero leakage, failure of
ordinary quantum mechanics, or publication-level novelty is established.
The original experimental paper already studied an alternating-memory explanation.

Start with the [generated measured report](results/report.md),
[models and proof](theory.md), [source audit](sources.md), and
[claim-to-prior-art table](prior-art.md).

| Claim | Status |
|---|---|
| Shared-operation words and counts from two measured families | Measured; immutable source hashes |
| Complementary source-removal examples | Conditional physical witnesses; finite bank is not an exhaustive prediction region |
| All-word leakage equivalence and conditional compatibility interval | Analytic within the stated incoherent, uniform-transfer model |
| Two-state memory realization | Explicit CPTP construction; not arbitrary circuit-label memory |
| Actual device leakage population or minimum memory | Unavailable; no matched calibrated monitor |
| Adversary sensitivity / nuisance calibration | Simulated at actual exposures; retrospective, local numerical refits |
| General stationary qubit CPTP baseline | Implemented; an all-length model is compatible with the declared joint count region; held-out training fits still fail |
| Five doubled-word outer prediction bounds | Analytic full-class bound; exact outward binomial-tail arithmetic |
| Gate-dependent, nonunital leakage equivalence | Analytic with explicit positivity conditions and unknown return-state polarization |
| Adequate prediction of every held-out family | Not achieved by saved training fits; no global exclusion |

## Reproduce offline

Python 3.12 and a C++11 compiler (`g++`) are used. The compiler only builds the
small reviewed local propagation kernels into a temporary directory; no downloaded
code or stored pickle is executed. NumPy density-matrix calculations independently
check its physics.

```sh
python -m pip install -r studies/shared-quantum-dynamics/requirements.txt
python studies/shared-quantum-dynamics/report.py --check
python studies/shared-quantum-dynamics/escalation_report.py --check
python studies/shared-quantum-dynamics/general_compatibility.py --check
python -m unittest discover -s tests -p 'test_shared_quantum_dynamics.py' -v
```

Offline checks use reduced counts, word hashes and saved predictions. They verify
all scores, physical constructions, source-removal memberships and report text.
They do not reconstruct full primitive words without the originals.

## Full external-source reproduction

Use a cache outside the checkout. The downloader verifies SHA-256 before parsing
and fails closed if a source changes. Original archives, notebooks and scripts
are never committed or executed.

```sh
python studies/shared-quantum-dynamics/sqd_sources.py --cache /tmp/shared-dynamics-sources --download
python studies/shared-quantum-dynamics/run.py --cache /tmp/shared-dynamics-sources --check
python studies/shared-quantum-dynamics/certificates.py --cache /tmp/shared-dynamics-sources --check
python studies/shared-quantum-dynamics/general_fit.py --cache /tmp/shared-dynamics-sources --check
python studies/shared-quantum-dynamics/general_geometry.py --cache /tmp/shared-dynamics-sources --check
python studies/shared-quantum-dynamics/lift_certificate.py --cache /tmp/shared-dynamics-sources --check
python studies/shared-quantum-dynamics/polish_general.py --cache /tmp/shared-dynamics-sources --check
python studies/shared-quantum-dynamics/sqd_repeat.py --cache /tmp/shared-dynamics-sources --check
python studies/shared-quantum-dynamics/transfer_audit.py --cache /tmp/shared-dynamics-sources --check
python studies/shared-quantum-dynamics/observable_audit.py --cache /tmp/shared-dynamics-sources --check
python studies/shared-quantum-dynamics/candidate_audit.py --cache /tmp/shared-dynamics-candidates --download --check
python studies/shared-quantum-dynamics/report.py --check
python studies/shared-quantum-dynamics/escalation_report.py --check
python studies/shared-quantum-dynamics/general_compatibility.py --check
```

This checks selected source words/counts, all twelve primary physical predictions,
all thirty bank predictions, independent Kraus calculations, geometry diagnostics,
matched-source hashes, and runner-up field audits. The manual
`shared-dynamics-reproduce.yml` workflow runs this route independently of ordinary CI.

## Numerical refits and simulated controls

The [protocol](protocol.json) fixed the primary models, length splits and statistics
before fitting. Papers and some source rows had already been inspected; this is
retrospective. The witness-bank profiles and pilot-adapted calibration are explicitly
exploratory extensions, not independent confirmation or new primary selections.
No validation counts are used for nuisance fits or bank selection.

```sh
OPENBLAS_NUM_THREADS=1 python studies/shared-quantum-dynamics/run.py --cache /tmp/shared-dynamics-sources --refit
OPENBLAS_NUM_THREADS=1 python studies/shared-quantum-dynamics/removal.py --cache /tmp/shared-dynamics-sources
python studies/shared-quantum-dynamics/certificates.py --cache /tmp/shared-dynamics-sources
OPENBLAS_NUM_THREADS=1 python studies/shared-quantum-dynamics/controls.py --cache /tmp/shared-dynamics-sources
OPENBLAS_NUM_THREADS=1 python studies/shared-quantum-dynamics/dependence.py --cache /tmp/shared-dynamics-sources
OPENBLAS_NUM_THREADS=1 python studies/shared-quantum-dynamics/calibrate_controls.py --cache /tmp/shared-dynamics-sources
python studies/shared-quantum-dynamics/format_artifacts.py
python studies/shared-quantum-dynamics/report.py
```

Refits can yield different equivalent coordinates and need not be byte-identical.
Primary fits use two starts; exploratory control/bootstrap fits use one start.
Their failures are never exclusion certificates. Deviance exceedances are empirical
sensitivity diagnostics; nuisance-adapted calibration is not uniform finite-sample
coverage. Simulations retain complete convergence flags and nuisance parameters
where used for generating-class recovery.

The exact missing observation and conditions that would refute the interpretation
are stated in [the theory](theory.md). No acquisition was run, no hardware access
was bought, and no source author was contacted. No Lean companion is opened for
these recurrence identities.

## Escalation refits

`general_fit.py --cache CACHE` runs missing GST/RB/joint fits; `--all-lengths`
adds the separately labeled in-sample diagnostic. Saved fits are reused, so move
`results/general-fits.json` aside to refit from scratch. `polish_general.py --cache
CACHE` runs the exploratory training-only continuation, with absolute finite
differences and three redundant factor scales fixed. It does not use validation
outcomes. All initial fits use the frozen two-start budget; continuation is an
explicitly recorded numerical follow-up. `escalation_report.py` regenerates its
report and verifies exact bound arithmetic offline. The full-source route checks
all saved general predictions and all 6,661 generalized-equivalence controls.

The general-channel kernel reuses shared word prefixes for speed; tests compare
it with the original full propagation and independent density matrices. The
initial 120-evaluation continuation is preserved separately from the extended
1,200-evaluation budget, which terminated after 997 evaluations. Neither run
was selected for validation performance.
