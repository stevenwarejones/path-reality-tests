# Phase-intervention design

The baseline null is the quantum P/Q-dephased class. Rejecting it in a heralded
Mach–Zehnder interferometer is the expected confirmation of arm coherence. This
study provides a calibrated, loss-complete template, and a substantive conditional
test of **definite-region models with local phase action and calibrated occupation**:
their complete-outcome contrast is at most the occupied-arm probability w.

**Status: analytical design and synthetic checks only.** The contribution is
formal and statistical re-analysis of established interference physics. No new
experimental mechanism, calibrated apparatus or achieved power is claimed.

| Material | Purpose |
|---|---|
| [Protocol](follow-up-protocol.md) | Four randomized phase settings, complete trial outcomes, null class and rejection rule |
| [Physical correspondence](physical-bridge.md) | Preparation, intervention, detector and finite-window assumptions |
| [Literature comparison](literature-comparison.md) | Closest theoretical and experimental precedents |
| [Open questions](open-questions.md) | Explicit feasibility, calibration and stronger-witness obligations |
| [Prospective sample planning](design.py) | Hoeffding and exact binomial intervals, occupation calibration and conditional power |
| [Validation](results/validation.md) | Numerical checks, counterexamples and limits |

The [protocol](follow-up-protocol.md#local-phase-definite-region-null) states the
rejected conjunction and its surviving alternatives. Neither test distinguishes
path sums from equivalent transfer matrices, excludes all trajectories or tests
faster-than-light or backward-time influence. A zero fringe need not mean dephasing.

## Reproduce the synthetic checks

With Python 3.12, from this study folder:

```sh
python -m pip install -r requirements.txt
python intervention_checks.py --check
python intervention_checks.py
python design.py --check
```

The `--check` commands compare against reviewed numerical snapshots without
writing them. `intervention_checks.py` without that flag regenerates
`results/intervention.json`; `design.py` regenerates `results/design.json`.
Paths resolve relative to the scripts, so
invocation from another working directory also works. No source data or workbook
inputs are required; numerical outputs come from the identified scripts.

The [pinned Lean case study](https://github.com/stevenwarejones/ontology-separation/blob/ab1eeb81119ff3fdcb46a771d9bcdf5e39ea3d29/docs/PATH_INTERFERENCE_CASE_STUDY.md)
certifies an ideal finite special case. It does not certify this experiment's
generalized instrument model or statistical analysis. `provenance.json` identifies
the scientific sources. The [general phase and local-region certificate](https://github.com/stevenwarejones/ontology-separation/blob/ffbbd0d3c0eb2066ea9b32b0dfd18e9c2d50bd4d/docs/PHASE_INTERVENTION.md)
proves the occupation bound and sharp lossy-family boundary; verification status
is recorded in [validation](results/validation.md).

## Relationship to other studies

The [Wen 2026 study](../wen-2026-propagator/README.md) audits published data and
their inference chain. This study proposes a new acquisition protocol; it neither
requires that audit to find a discrepancy nor inherits an empirical validation
from it. Each study has its own review scope, evidence and open questions.
