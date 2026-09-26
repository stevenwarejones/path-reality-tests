# Phase-intervention design

Can a controlled phase change in a selected spatial region distinguish coherent
propagation from a specified model that has lost coherence between that region
and its complement? This study supplies a proposed experiment, explicit model
predictions, calibration obligations, a finite-sample decision rule and a
conservative conditional power calculation.

**Status: proposed design, with synthetic checks only.** No experimental data have
been collected. Apparatus calibration, achieved statistical power and experimental
novelty remain unestablished. Phase interference is established physics. The intended
contribution is a reproducible test of clearly stated model classes, including
loss and nuisance bounds.

| Material | Purpose |
|---|---|
| [Protocol](follow-up-protocol.md) | Four randomized phase settings, complete trial outcomes, null class and rejection rule |
| [Physical correspondence](physical-bridge.md) | Preparation, intervention, detector and finite-window assumptions |
| [Literature comparison](literature-comparison.md) | Prior work, actual review depth and staged searches |
| [Open questions](open-questions.md) | Explicit feasibility, calibration and stronger-witness obligations |
| [Prospective sample planning](design.py) | Precision, power and a count-based decision rule |
| [Validation](results/validation.md) | Numerical checks, counterexamples and limits |

The null is the declared P/Q-dephased class with setting-independent preparation,
downstream response and loss, allowing a separately justified nuisance budget.
A contrast can reject that class. It cannot distinguish equivalent path-sum and
transfer descriptions, exclude every definite-trajectory theory, establish
literal occupation of every path, or demonstrate faster-than-light or
backward-time influence. A zero contrast can also occur for a coherent state.

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
the scientific sources. The newer formal phase-intervention development is linked in [validation](results/validation.md).

## Relationship to other studies

The [Wen 2026 study](../wen-2026-propagator/README.md) audits published data and
their inference chain. This study proposes a new acquisition protocol; it neither
requires that audit to find a discrepancy nor inherits an empirical validation
from it. Each study has its own review scope, evidence and open questions.
