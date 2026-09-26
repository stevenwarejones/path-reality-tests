# Requirement-to-evidence report

This is a conditional protocol/resource study and formal verification of known
causal reasoning. It supplies neither an executed spacelike test nor a measured
bound. The [research report](README.md) states the rejected conjunction and
surviving model classes. The [question register](open-questions.md) records all
ST01–ST12 dispositions, including externally blocked apparatus and source inputs.

Formal source links below point to the companion review branch. Both repositories
build independently; no unmerged cross-repository dependency is installed.

- [Finite model and countermodels](https://github.com/stevenwarejones/ontology-separation/blob/study/spacetime-influence/OntologySeparation/Experiments/SpacetimeInfluence.lean)
- [Algebra, calibration, intervals and geometry](https://github.com/stevenwarejones/ontology-separation/blob/study/spacetime-influence/OntologySeparation/Experiments/SpacetimeInfluenceBounds.lean)
- [Lean examples/tests](https://github.com/stevenwarejones/ontology-separation/blob/study/spacetime-influence/Tests/SpacetimeInfluence.lean)
- [Generated theorem report](https://github.com/stevenwarejones/ontology-separation/blob/study/spacetime-influence/examples/spacetime-influence.html)
- [Complete axiom snapshot](https://github.com/stevenwarejones/ontology-separation/blob/study/spacetime-influence/docs/AXIOM_AUDIT.txt)

Every theorem name below is in `OntologySeparation.SpacetimeInfluence`.

| Obligation | Evidence | Status / limit |
|---|---|---|
| Compare at least two protocols and select one | [protocol.md](protocol.md), comparison table and earlier-record design | Done; optical modulation selected for in-cone positive control. Earlier-record alternative has a precise disposition. |
| Distinguish interventions, correlations and complete outcomes | [physical-bridge.md](physical-bridge.md); `sharedBit_correlated`, `postselection_counterexample`, `dependentPreparation_gap` | Done; observational conditioning requires independence/consistency to estimate do-effects. |
| Derive operational null, normalize and show nonemptiness | `Model.joint`, `Model.observed`, `Model.marginal_formula`, `Model.no_influence`, `sharedBit` | Done for a finite sufficient causal class. |
| Normalized tunable alternative and separator | `alternative`, `alternative_gap`, `alternative_excluded`; g=0, g=1 and strict-positive tests | Done; excludes all finite models in the defined class for exact g>0. Not a field-theory alternative. |
| Access equivalence | `receiver_access_equivalent` | Done for the shared-bit/independent-sender construction and receiver-only access; universal hidden-path equivalence outside scope. |
| Physical-to-operational bridge | `local_operation_fixes_effect`; detector primary sources in [evidence table](evidence-table.md) | Conditional algebra done. Real-device localization/commutation and response calibration externally blocked, ST04/ST06. No continuum proof. |
| Premise countermodels | [physical bridge](physical-bridge.md), formal examples and [numerical tests](../../tests/test_spacetime_influence.py) | Done for preparation, receiver channel, selection, supports, storage, predictable settings and constant-effect interpretation. |
| Geometry/timing and backward record order | `spacelike_of_budget`, `earlier_record`; [config.json](config.json) and geometry tests | Done for supplied support intervals; actual support evidence absent. Vacuum c used. |
| Nuisance in witness units | `coupling_gap`, `three_event_bound`, `calibrated_gap`, `contamination_gap`; coupling-flip counterexample | Derived by coupling, union and triangle bounds; numerical eₓ values require calibration. |
| Binary CI, absolute upper bound, random counts | [statistics.md](statistics.md); `binary_influence`, `difference_interval`, `absolute_interval_upper`, `strict_interval_exclusion` | Done under stated assumptions; exponential/binomial coverage is derived in prose, not Lean. |
| Memory, losses, conditioning, multiplicity and stopping | Fixed-horizon score/bias derivation; `fair_score_mean`, `biased_score_mean`, `rejection_risk`; protocol all-trial rules | Done. No optional-stopping theorem; binary-TV interpretation of memory interval needs a constant effect. |
| Power/resources over gaps and nuisances | [design.py](design.py), [config.json](config.json), [results.json](results.json) | Done, explicitly synthetic/prospective, with independent checks and strict rounding. Achieved sensitivity blocked, ST07. |
| Literature stages, primary versions, adverse evidence and closest-work comparison | [search log](literature-search-log.md), [evidence table](evidence-table.md) | Done through 2026-09-26 UTC with access gaps disclosed in ST11. Closest binary methods comparison: Albanese 2026. No priority claim. |
| Experimental protocol, equipment and controls | [protocol.md](protocol.md) | Done as a proposal. Calibration records and lab specifications absent; no invented prices or availability. |
| Distinguish proof, computation and observation | [README evidence/status table](README.md#evidence-and-contribution) | Done. No new source datasets/PDFs or measurements included. |
| Full final-head verification and generated equality | Repository CI plus the checks below | Readiness depends on the actual PR head's green checks; links and exact tested revisions are recorded in the PR descriptions. |

## Reproducible verification

The empirical workflow runs the complete test suite, input allowlist, local-link
checker, all existing study checks and the new `design.py --check`. The new tests
independently verify CP tail inversion/small-count coverage, strict binomial
rejection probabilities, memory-null assignments, coupling budgets, timing
failures, selection and coarse-outcome limitations. Existing data and outputs
must remain byte-identical.

The formal workflow runs `sh scripts/check.sh`, regenerates the axiom snapshot
and theorem HTML, and requires generated-file equality. All new trust roots are
registered; only `propext`, `Classical.choice` and `Quot.sound` are admitted.
A targeted build alone is insufficient to establish full repository readiness.

Exact base revisions are recorded in the question register. Final-head CI run
links belong to the PR descriptions because inserting a commit's own identifier
into its contents would change that identifier. Apparatus premises remain
untested regardless of software CI status.
