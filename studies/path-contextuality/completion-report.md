# Requirement-to-evidence report

The mathematical design is a finite specialization of a known contextuality
theorem. No apparatus premises have been experimentally certified. The formal
companion's review state and CI are tracked in
[PR #87](https://github.com/stevenwarejones/ontology-separation/pull/87).

| Obligation | Evidence | Disposition |
|---|---|---|
| Prior-art and apparatus comparison before selecting a witness | `literature.md`, `evidence.json`, search log; KLP Theorem 3 selected | Done at recorded retrieval scope; peripheral supplement/version gaps retained |
| Explicit stochastic ontology and operational/ontic distinction | `derivation.md`; formal `Model`, `ResponseCap`, `Disturbance`, `cap_of_measurement_equivalence` | Implemented; formal CI pending |
| Bound with a nonempty null | `Model.bound`, `full_bound`, `null_nonempty`; independent row derivation | Implemented; formal CI pending |
| Normalized violating quantum instrument | Kraus completeness, full channel/effect identities, pulled-back POVM and `quantum_realizes_table`; independent matrix/Fraction tests | Numerical/exact arithmetic checks passed; formal CI pending |
| Adversarial disturbance/context/selection models | `exact_countermodels`, `PathContextualityCountermodels.lean`, six constructions in derivation | Exact Python checks passed; formal CI pending |
| No unjustified finite reduction | Universal proof over every finite cardinality; no LP/enumeration or polytope claim | Done; measurable extension is separately scoped and not Lean-certified |
| Complete outcomes and normalization convention | 9-outcome lossy table, bypass versus joint distinction, selection counterexample | Done for specified noise channel |
| Robustness with imperfect equivalences | `Model.robust_bound`, `ceiling_mono`, `Model.full_bound_of_upper`; explicit near-model or relaxed representation premise | Conditional result; operational closeness does not imply ontic closeness |
| Accessible calibration family | 44 synthetic contexts, exact spanning-coordinate determinant, protocol mapping | Ideal construction done; actual control availability/characterization externally blocked |
| Finite statistics and calibration costs | CP rule, Hoeffding certificate, fixed/random context counts, seeds, all-trial budgets and unit tests | Done conditional on IID, valid representations and characterized controls; not an end-to-end apparatus budget |
| Optimization and sample comparison | Exact rational 1,215-point grid, complete outputs | Grid minimum only; global optimization neither needed nor claimed |
| Motivating path/velocity/time bridge | Interface non-identifiability by freely appended history labels; primary observable comparison | Projector claim scoped; dynamic bridge outside this experiment |
| Distinction from #85 | Different rejected conjunction; PH05/PH07/PH08 mapped to PC01–PC03 | Done; original study untouched |
| No new-data or priority claims | Synthetic labels, primary citations, explicit closest prior theorem | Done |
| Final-head verification | PR #87 and companion CI | In progress; no review-ready formal claim yet |

## Claims table

| Category | Conclusion |
|---|---|
| Excluded in the ideal mathematical example | Finite stochastic ontic models matching joint and bypass probabilities while satisfying the negative-response cap and identity-plus-disturbance representation; hence their definite-occupancy subclass |
| Still compatible | Models dropping measurement/transformation noncontextuality, identity representation or common preparation; contextual/invasive trajectory models; generic Bohmian-style interpretations |
| Unknown empirically | Whether any proposed device meets the representation/calibration premises, trial stability, losses, mode completeness and required power |
| Depends on another physical bridge | Photon velocity, relativistic worldline claims, integrated excitation time and backward chronological motion |

Every PC question has a disposition in `open-questions.md`. A missing apparatus
calibration is not called a completed experiment. Source/supplement access gaps
cannot affect the exact instrument or conditional statistical calculations,
because no experimental parameters from those sources enter them.
