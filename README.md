# Path reality tests

This repository tests when claims that quantum paths or trajectories are “real”
can be checked by experiment. Studies connect explicit models to observable
predictions and document the data, assumptions and reproducible calculations
needed to compare them. Studies provide independent, reproducible checks.

| Study | Question | Status |
|---|---|---|
| [Wen 2026 propagator](studies/wen-2026-propagator/README.md) | What can reconstructed propagators and endpoint statistics tell us about path descriptions? | Public-data audit available; full reconstruction needs additional records |
| [Path contextuality](studies/path-contextuality/README.md) | Can a finite weak path probe reject explicit noncontextual representations? | Prospective conditional design; exact instrument and calibrated-trial budget |
| [Phase-intervention design](studies/phase-intervention-design/README.md) | Which calibrated dephased models would a four-phase experiment exclude? | Proposed; synthetic checks and conditional power |
| [Spacetime causal influence](studies/spacetime-causal-influence/README.md) | Can a remote intervention change a complete local record outside its future light cone? | Conditional protocol, formal companion and synthetic sensitivity study |

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
PDFs and complete workbook extracts are not committed. Other datasets require
explicit source/license review before inclusion.
Repository code is distributed under the [Apache 2.0 license](LICENSE); external
data and papers retain their stated licenses. Changes are reviewed through pull
requests before merging.
