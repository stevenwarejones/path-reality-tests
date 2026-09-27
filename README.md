# Path reality tests

This repository tests when claims that quantum paths or trajectories are “real”
can be checked by experiment. Studies connect explicit models to observable
predictions and document the data, assumptions and reproducible calculations
needed to compare them. Studies provide independent, reproducible checks.

| Study | Question | Status |
|---|---|---|
| [Fine timing information](studies/fine-timing-structure/README.md) | What do pulse identity and finer timing add beyond four categories? | Conditional archive bounds and paired recovery; no measured departure |
| [Continuum and finite models](studies/continuum-finite-models/README.md) | Which specified finite dynamics differ from a continuum ring? | Synthetic conditional design; formal completion and apparatus calibration pending |
| [Wen 2026 propagator](studies/wen-2026-propagator/README.md) | What can reconstructed propagators and endpoint statistics tell us about path descriptions? | Public-data audit available; full reconstruction needs additional records |
| [Path contextuality](studies/path-contextuality/README.md) | Can a finite weak path probe reject explicit noncontextual representations? | Prospective conditional design; exact instrument and calibrated-trial budget |
| [Phase-intervention design](studies/phase-intervention-design/README.md) | Which calibrated dephased models would a four-phase experiment exclude? | Proposed; synthetic checks and conditional power |
| [NIST Bell-record causal audit](studies/nist-bell-causal-audit/README.md) | What conditional influence bounds do public trial records support? | Full raw reconstruction, measured marginals and explicit calibration limits |
| [Born-rule identifiability](studies/born-rule-identifiability/README.md) | Can additional datasets separate probability changes from control and detector errors? | Measured scan–Viviani combination, ordinary qutrit compatibility certificates and conditional distance bounds; no identified probability-rule violation |
| [Synthetic timing shifts](studies/synthetic-timing-shifts/README.md) | Could sparse records reveal timing shifts or effects that reverse sign? | Synthetic and real-background injections; artificial assignments and held-out scores |
| [Spacetime causal influence](studies/spacetime-causal-influence/README.md) | Can a remote intervention change a complete local record outside its future light cone? | Conditional protocol, formal companion and synthetic sensitivity study |
| [Collapse compatibility](studies/collapse-compatibility/README.md) | What can interference, mechanical, space and radiation data jointly determine about collapse noise? | Reproduced sources; two conditional feasibility studies and explicit inference limits |

The Wen study now includes reproduced Figure 3/4 summaries, separate endpoint
and theory comparisons, sensitivity checks, and an explicit list of missing
inputs. Its [method comparison](studies/wen-2026-propagator/method-comparison.md)
explains which numerical differences affect reproducibility and which tested
choices preserve the broad comparison. The [completion report](studies/wen-2026-propagator/completion-report.md)
connects every audit requirement to evidence and identifies the remaining inputs.

The companion [ontology-separation](https://github.com/stevenwarejones/ontology-separation)
repository contains the formal Lean framework. Its
[path-interference case study](https://github.com/stevenwarejones/ontology-separation/blob/ab1eeb81119ff3fdcb46a771d9bcdf5e39ea3d29/docs/PATH_INTERFERENCE_CASE_STUDY.md)
explains the distinction between equivalent descriptions and observable model
separation.

The Wen study includes the pinned CC0 Dryad originals for reproducible offline
checks, with the official source, version and SHA-256 manifest retained. Paper
PDFs and complete workbook extracts are not committed. The NIST audit keeps its
source inputs external and hash-pinned; committed tables are derived analysis
artifacts. New source datasets or event extracts require explicit authorization
and source/license review before inclusion.
Repository code is distributed under the [Apache 2.0 license](LICENSE); external
data and papers retain their stated licenses. Changes are reviewed through pull
requests before merging.

See [sharp compatibility boundaries](studies/path-compatibility/README.md) for exact attainers and the two-witness decision extension.

See [timing structure and identifiability](studies/timing-structure/README.md) for fixed cancellation and temporal diagnostics on existing NIST records, with exact invisibility witnesses.
