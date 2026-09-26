# Relation to existing interference and ontology tests

Primary-source review: 2026-09-26 UTC. This comparison identifies the closest
reviewed argument and experiment, not an exhaustive priority result. The
contribution is formal verification and a statistical robustness analysis of
established interference physics. No new experimental mechanism is proposed.

## Closest precedents

| Source and inspected version | Material read | Comparison |
|---|---|---|
| [Hardy, arXiv:1205.1439v3 (2012)](https://arxiv.org/html/1205.1439v3) | Introduction, Section 3 interferometer argument, restricted ontic-indifference discussion | Closest reviewed conceptual argument: a phase shift on one path changes the eventual port while a localized preparation in the other path is invariant. |
| [Grangier, Roger and Aspect (1986)](https://doi.org/10.1209/0295-5075/1/4/004) | Published paper, apparatus and gated Mach–Zehnder results, pp.176–179 | Closest reviewed experimental architecture: the same triggered source supports beam-splitter anticorrelation and complementary-port fringes. Its dark-corrected visibility is not an absolute heralded contrast, and is insufficient to establish this loss-complete bound's violation. |
| [Błasiak, New J. Phys. 17, 113043 (2015)](https://arxiv.org/abs/1502.07308) | Published PDF, Sections 1–2 and model overview | Explicit adverse example to identifying particle position with all phase-carrying ontology: a detectable particle and undetectable “ghost” degrees of freedom interact locally. A phase element can affect the other-path ghost and hence the later detector outcome. |
| [Biswas, García Díaz and Winter, arXiv:1701.05051v3](https://arxiv.org/html/1701.05051v3) | Introduction, phase-unitary/POVM formulation Eqs.1–4, two-path visibility | Ordinary fringe visibility probes coherence; a dephased-null rejection is established physics. |
| [Hacker et al., arXiv:2305.03641v2](https://arxiv.org/html/2305.03641v2) | Introduction, setup and feedback scope | Count-based phase locking is an apparatus precedent, without certifying the independent heralded trials assumed here. |

Hardy's ontic-indifference premise concerns the support of a localized pure
preparation. Our null instead constrains the final response for every Q-located
ontic state in a possibly coherent preparation. Those sets are not automatically
the same. Neither Hardy's full ψ-ontology theorem nor his restricted-indifference
premise is being identified with the finite response inequality here. In
particular, a localized Bohmian particle can coexist with a phase-sensitive wave
in the other arm. This violates the null's Q-response invariance while retaining
a definite particle position. Błasiak's local construction reinforces that the
null's response condition is stronger than generic spatial locality.

The occupation equality μ(P)=w adds a separate cross-context identification.
The new analysis makes this assumption explicit, permits arbitrary occupied-arm
responses and loss, and computes finite-statistics thresholds with occupation
calibration. The argument and interferometer are known; no literature result
reviewed here establishes priority for the precise robust inequality. The design
therefore stops at a formal/statistical re-analysis, with actual apparatus
feasibility unresolved.

## Modular variables and the past of a photon

| Source and inspected version | Material read | Relationship |
|---|---|---|
| [Tollaksen et al., arXiv:0910.4227v1](https://arxiv.org/abs/0910.4227v1), New J. Phys. 12, 013023 (2010) | Introduction and modular-variable/causality discussion | A localized-particle account with interactions involving the other slit need not satisfy Q-response invariance. The equivalence of Heisenberg and Schrödinger descriptions does not give different interference predictions. |
| [Aharonov and Rohrlich, arXiv:2011.11667](https://arxiv.org/abs/2011.11667) | Primary PDF abstract and initial cavity/unitary construction | Their counterfactual-communication analysis uses a modular-angular-momentum current. It concerns additional physical degrees of freedom, not a derivation of the occupied-particle-only response premise. The complete current calculation has not been rederived here. |
| [Vaidman, arXiv:1304.7474v1](https://arxiv.org/abs/1304.7474v1) | Sections II–V on weak traces, ensemble inference and the nested interferometer | A weak-trace criterion for a pre/postselected ensemble is different from a trial's definite-region occupation. A strong region measurement changes the context; it is not a simultaneous path record for our phase runs. |

These are comparisons of explicit premises, not empirical evidence for a
backward-time influence. The current design has no intervention that distinguishes
causal directions or equivalent quantum representations.

## Contextuality and the occupation premise

[Pusey, arXiv:1409.1535v2](https://arxiv.org/html/1409.1535v2), Theorem 1 and its
proof, use measurement noncontextuality, outcome determinism for sharp
measurements, an unbiased noisy-projector measurement and a quantified
disturbance equivalence. His p₋ has a specified normalization by the undisturbed
postselection probability, not the ordinary detected-subset denominator.

[Kunjwal, Lostaglio and Pusey, arXiv:1812.06940v2](https://arxiv.org/html/1812.06940v2),
Theorem 3, Eqs.29–31 and Section IV, give a robust finite-pointer version:

\[
p_-\leq p_F(1+p_m)/2+(1-p_F)p_d.
\]

Here p₋ is a joint pointer-negative/postselection-success probability. The
measurement equivalence mixes a sharp measurement and trivial noise; the
transformation equivalence mixes identity and disturbance. Their empirical use
needs operational equivalences across a suitable preparation set and quantified
disturbance. A single measured occupation fraction supplies neither. Our equality
between ontic occupation and the which-region probability is one
noncontextuality-type premise, not a proof of generalized noncontextuality or an
application of either theorem. PH05 records the distinct unresolved instrument
and equivalence problem.

## Statistical sources

[Hoeffding (1963), Theorem 2](https://doi.org/10.1080/01621459.1963.10500830)
supplies the retained bounded-variable baseline. The exact binomial alternative
is the [Clopper–Pearson construction (1934)](https://doi.org/10.1093/biomet/26.4.404).
The coverage argument here follows directly by inverting binomial tails;
[the SciPy exact proportion-interval documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats._result_classes.BinomTestResult.proportion_ci.html)
identifies the method. `design.py` evaluates beta quantiles, and tests independently
sum binomial masses to check marginal coverage, unequal-count null risk and
small-design power. Conditional assignment coverage and the sufficient power
envelope are derived in the protocol; neither is claimed as a new theorem.

## Search record

Earlier searches used `single photon path interference local phase shift weak
measurement phase shifter contextuality experiment`, the Grangier/Roger/Aspect
1986 title, phase-cycling/coherence-witness terms and the Kunjwal/Lostaglio/Pusey
title. They identified the phase-unitary formulation, apparatus precedents and
operational-equivalence requirements. The current focused search added Hardy's
ontic indifference, local qubit interferometer models, Aharonov/Rohrlich modular
variables and Vaidman's past-of-a-particle argument. Primary versions and actual
review depth are distinguished above; failed HTML retrievals were followed by
public author/publisher PDFs where available. No source PDFs are stored here.

Further apparatus, leakage and finite-equivalence questions are PH01–PH08 in the
[question register](open-questions.md). No pilot evidence is available from which
to identify apparatus-specific mechanisms.
