# General CPTP escalation: measured bounds and a broader obstruction

Can shared ordinary dynamics predict the acquisition, and what costs do the combined records identify? The answer remains split: we now have general-channel local fits, an optimizer-independent outer certificate, and a broader exact leakage obstruction. An all-length stationary qubit construction is compatible with a declared conservative joint count region. We do **not** have a successful held-out stationary-qubit predictor or a global exclusion.

## General CPTP fits

All three gates are arbitrary stationary qubit CPTP maps, with arbitrary physical preparation/readout up to a common unitary gauge. The original GST≤1024 / RB≤1000 training split is retained. The all-length fit is an explicitly in-sample diagnostic. Optimizer success is a termination flag, not a global optimality certificate.

| Fit source | Training/all-fit deviance | Evaluations | Terminated | Optimality diagnostic | Long GST violations | Long RB violations |
|---|---:|---:|---|---:|---:|---:|
| GST | 3234.204809 | 180 | False | 393 | 14 | 0 |
| RB | 662.758882 | 180 | False | 39.06 | 305 | 0 |
| all_joint | 7288.334963 | 180 | False | 1083 | 0 | 0 |
| joint | 3917.780551 | 180 | False | 60.77 | 10 | 0 |
| joint_continuation_120 | 3914.207563 | 120 | False | 18.7 | 19 | 0 |
| joint_continuation | 3914.136583 | 997 | True | 2.553 | 19 | 0 |

The observable Jacobian retains gauge redundancies: 52 factor coordinates describe at most 31 generic observable degrees of freedom. Relative 10⁻⁶ diagnostic ranks are GST 31, RB 23, joint 30. Each threshold scales with its own largest singular value, so these numbers are not a monotone dimension test or global identification theorem. In particular, no missing structural degree of freedom is demonstrated by this local rank comparison.

**Constructive full-acquisition compatibility:** the all-length model has zero cell violations in a simultaneous region allocating α=0.025 to all 6,661 cells and α=0.025 to ten family/length means. Minimum cell slack is 0.00513052; minimum slack across all constraints is 0.00250084, compared with a 2×10⁻⁸ numerical verification tolerance. Cell endpoints use exact outward binomial-tail arithmetic. This explicit stationary qubit model has no leakage level and no nontrivial classical memory. **No strictly positive leakage or memory cost follows from this region.** This is an in-sample, numerically verified physical compatibility example, not held-out validation, a globally optimal likelihood fit, or adequacy under every stronger statistic. The region/mean audit is an explicitly retrospective diagnostic; it does not replace the frozen predictive test.

Violations use the earlier simultaneous Clopper–Pearson intervals separately within each evaluation family. For all_joint these are in-sample diagnostics, **not held-out predictions**. Budget-limited or poorly stationary local solutions cannot settle whether a better general model exists. Even disappearance of cell violations would not prove a full sampling-model goodness of fit.

## Full-class outer prediction certificate

For a whole word channel Φ_w, common reset and binary effect, every stationary qubit CPTP model obeys p(ww) ≤ min(1,(2√p(w)+√p(empty))²). See the independent proof in [escalation-theory.md](../escalation-theory.md). Six exact outward binomial bounds, with Bonferroni tail 1/120, give simultaneous 95% coverage under within-row binomial sampling. There is no intermediate reset in ww. These doubled sequences are predictions, not newly acquired observations.

| RB row (zero-based) | Word length | Plus / shots | Doubled length | Joint outer upper bound | Either source removed: this relaxation |
|---|---:|---:|---:|---:|---:|
| 2 | 8 | 1/285 | 16 | 0.3728347088 | 1 |
| 40 | 38 | 4/285 | 76 | 0.5005422008 | 1 |
| 82 | 156 | 5/285 | 312 | 0.5364587006 | 1 |
| 116 | 319 | 7/285 | 638 | 0.6031255268 | 1 |
| 472 | 976 | 9/285 | 1952 | 0.6649465203 | 1 |

GST empty circuit: 0/50; outward upper bound 0.0913086895. Exact integer binomial tails and integer square-root rounding certify the printed bounds at denominator 10¹⁰.

**What removal establishes:** both sources are required by this particular certificate. **What it does not establish:** that either full single-source feasible set reaches one, or that these bounds improve on every implicit constraint of GST alone. The earlier 30-model attainable spread is retained only as historical exploratory evidence; it is not promoted to a full-class region. With unrestricted within-row dependence the reported statistical upper bound is one.

## Why an observable-basis relaxation was insufficient

A nominally complete four-fiducial GST matrix has nonzero point-estimate singular values, but the entrywise simultaneous count box contains a singular matrix. Two saved rational matrices inside that box have opposite **exact integer determinant signs**, so their segment crosses determinant zero. Repeated circuit entries are shared. Uniform inversion of this box is therefore impossible. This certifies failure of this relaxation, not impossibility of stronger CPTP-aware bounds: the singular witness need not have a common physical realization or satisfy other records.

## Gate-dependent, nonunital leakage obstruction

The exact similarity now permits gate-dependent loss/return, arbitrary computational CPTP channels, anisotropy and nonunital translation. Unknown physical return-state polarization absorbs the translation. Positivity and readout restrict the similarity parameter κ to an explicit interval. Every word remains indistinguishable, while leakage populations and loss rates scale by κ. A two-state classical flag also realizes this entire block model. This is a theorem about a model class, not a fitted device population.

Synthetic control on all 6,661 recorded words: κ=0.6, common conditional interval lower endpoint 0.200000, maximum all-word numerical difference 1.1e-13. Independent qutrit Kraus TP error 2.99e-16; minimum Choi eigenvalue -9.33e-17; minimum transformed return-state eigenvalue 0.488390.

| Actual recorded geometry | Maximum change when return state is incorrectly locked to I/2 |
|---|---:|
| GST row 4351, length 8195 | 0.001585122 |
| RB row 1142, length 1384 | 0.000037239 |

The locked-return perturbation is weak at the deposited exposures even with known nuisance parameters. Under independent binomial sampling, D(P_original || P_locked) is 0.01446393 jointly. Pinsker gives power ≤ α+√(D/2): a level-0.05 test has power at most 0.135041 for this synthetic pair (GST alone 0.134595, RB alone 0.058700). These are numerical evaluations of an analytic power ceiling, not an experimental fit. Nuisance refitting cannot improve a uniformly valid test past the simple-pair ceiling. Along the exact free-return orbit, power equals size for every test because the complete count laws coincide.

## Exact lifts of the compatible qubit model

The stronger construction starts from the compatible all-length general qubit map Ψ_g. For τ=I/2 it chooses Φ_g=(Ψ_g−λ_g R_τ)/(1−λ_g), return state τ_g=τ+(Ψ_g(τ)−τ)/η_g, and leakage effect z=Tr(Eτ). Complete positivity holds exactly when λ_g≤2λ_min(J_g); return positivity requires η_g≥||d_g||. The collapse Q(ρ_C,L)=ρ_C+Lτ intertwines the lifted channel with Ψ_g. Therefore the qubit, leakage and two-state-flag realizations give **identical probabilities for every word**, not only this deposit. These are sharp conditions for this construction with the qubit representation fixed, not gauge-invariant device bounds.

Unlike the original restricted-fit examples, these lifts share the prediction vector inside the full-acquisition joint count region. They require the particular unknown leakage response z=Tr(Eτ) and unknown return polarization; calibrated apparatus information could exclude them.

| Ordinary realization | Conditional leakage after recorded Gx^8192 | Maximum difference over all 6,661 probabilities |
|---|---:|---:|
| Lift fraction 0.0 | 0.000000% | 3.57e-12 |
| Lift fraction 0.25 | 9.516326% | 4.41e-12 |
| Lift fraction 1.0 | 32.571007% | 4.85e-12 |

The recorded GST row is 4273 (zero-based), with 21/50 plus counts; all models predict 0.316160321. The population differences are stipulated model resources, **not observed populations**. The all-word theorem, not more sampling of the existing binary readout, supplies the nonidentification conclusion. A calibrated leakage-sensitive monitor after one primitive gate distinguishes zero from its positive loss rate.

The GST example is `GxGxGx(Gi)^8192`. These are simulated differences on actual words, not observed residual attribution or powered detection claims. Freeing return polarization restores exact equivalence. A calibrated subspace monitor after one Gi would distinguish 0.001 from 0.0006 in this example. More of the same unknown binary readout cannot.

## Claim status after escalation

| Claim | Status |
|---|---|
| Arbitrary stationary qubit CPTP fits implemented | Physical local numerical examples; optimizer limits retained |
| Five full-class joint outer prediction bounds | Certified analytically and by exact count-tail arithmetic; conditional sampling |
| Strict contraction of complete GST-only feasible prediction region by RB | Not established; certificate removal is weaker than full-model removal |
| Nonunital/gate-dependent all-word leakage equivalence | Analytic similarity plus exact lifts of the compatible qubit fit; conditional CP budgets and Kraus controls |
| Full-acquisition joint count-region compatibility | Explicit stationary qubit point; in-sample, numerical physical certificate |
| Necessary device leakage or finite-memory cost | No positive lower bound follows from this declared region; actual populations not measured |
| All stationary qubit models excluded by acquisition | Not certified |
| Major novelty beyond original GST/RB and leakage-identification work | Not established; precise comparison in prior-art.md |

The original paper already compared GST with RB, used general Markovian gate-set estimates, and studied alternating errors. This escalation closes a model-generality gap in our implementation, supplies a valid outer bound in place of the finite-bank interpretation, and exposes return-state calibration as another anchor needed to identify leakage. It does not claim that reanalysis of this known pairing alone earns the intended major-result bar.
