# From interventions to causal supports and probability budgets

The intervention is a local control X changing a source/modulator at A. For a
fixed preparation and complete receiver record B define
Δ(A→B)=½∑b|P(b|do X=1)−P(b|do X=0)|. The chosen complete record is binary,
so Δ=|p₁−p₀|. There is no inference about unobserved worldlines.

## Three distinct claims

The finite model in `SpacetimeInfluence.Model` has preparation μ(λ), normalized
sender kernel Kₓ(a|λ), and a fixed receiver kernel R(b|λ). Its joint law is
Jₓ(a,b)=∑λ μ(λ)Kₓ(a|λ)R(b|λ). Summing *every* sender outcome gives
Pₓ(b)=∑λ μ(λ)R(b|λ), independent of x. This derives the null from channels;
zero influence is not a constructor field. The shared-bit example is nonempty
and perfectly correlated. This classical sufficient class does not purport to
represent all joint entangled measurements.

The separate Lean theorem `local_operation_fixes_effect` checks the quantum
algebraic bridge: if ∑j K†ⱼKⱼ=I and each Kⱼ commutes with the receiver effect E,
then ∑j K†ⱼ E Kⱼ=E. Taking any state expectation preserves the receiver
probability, including correlated quantum states. The identity does not require
a tensor-product factorization. It is conditional on completeness and
commutation; it does not prove local QFT algebras, state positivity, or which
laboratory operations belong to those algebras.

The physical bridge is that microcausality of local observables and a localized,
nonselective operation imply the required commutation when the complete supports
are spacelike. Tjoa's delta-coupled scalar detector analysis and Martín-Martínez's
smeared-detector analysis are reviewed in the [evidence table](evidence-table.md).
The de Ramón–Papageorgiou–Martín-Martínez causal-factorization and noncompact-tail
analyses further restrict admissible extended-detector operations and provide
model-dependent effective-cone bounds. None is a calibrated optical-apparatus theorem. Finite switching, tails,
matter response, preparation and detector approximations must be justified for
the actual system. Feynman/Wightman correlations outside a light cone are not
the commutator controlling response. A lattice propagation cone or a rotating-wave
approximation supplies no exact continuum light cone here.

Finally, evidence must establish the apparatus supports, record integrity and
calibration assumptions. No such measurements have been supplied. Thus the
result is a conditional protocol/resource study and formal verification of
known causal reasoning, not an empirical spacetime bound.

## Regions and uncertainty

Let A's entire potentially setting-bearing support have times [a₀,a₁], and B's
interaction/detection/irreversible local latch have times [b₀,b₁]. Enclose their
spatial supports in balls of radii rA,rB with center-distance lower error uD.
The triangle inequality yields d_min=D−rA−rB−uD. Include the RNG, driver cables,
modulator and any earlier electrical/radiative setting leakage in A, rather
than just the nominal optical actuation time. If one ball is inappropriate,
bound each support pair and take the worst margin.

Every pair of events is spacelike if

\[
D-r_A-r_B-u_D>c\max(b_1-a_0,a_1-b_0).
\]

The maximum bounds |tB−tA| throughout the intervals. `spacelike_of_budget`
checks this arithmetic implication. Strictness matters: equality is null
separation. The test uses vacuum c=0.299792458 m/ns. Fiber route length, refractive
index and group velocity cannot enlarge the allowed spacelike window. Distance
and timing errors must be bounds with stated coverage; an RMS jitter is not a
hard support and a Gaussian waveform has no compact support.

For the hypothetical 30 m example, rA=rB=0.5 m, uD=0.1 m, A=[−5,20] ns,
B=[35,65] ns. Then d_min=28.9 m and max time difference=70 ns. The margin is
7.91452794 m or 26.40002351 ns. An earlier setting leak at −100 ns destroys
this certificate. A 3 m tabletop arrangement with the same electronics fails
it as well. These are exact parameter checks, not achieved jitter values.

For an earlier-record test use the different condition recordLatest<choiceEarliest.
Fix the record in independent local storage before generating the later choice.
A hash made *after* the choice proves nothing about earlier fixation; a hash
alone also cannot prove the sensor's time. Bound detector-to-latch latency and
keep the RNG inaccessible to the record producer. Prevent B→X dependence as
well as X→B leakage: a later RNG that reads the earlier record creates a
forward common-data explanation of correlation. Predictable seeds defeat the
measurement-independence premise even if the announced choice time is later.

## Deriving nuisance in the witness's units

Couple each actual outcome Bₓ to an ideal outcome B*ₓ on one probability space.
For any selected event E, the indicators agree whenever the records agree, so
|P(Bₓ∈E)−P(B*ₓ∈E)|≤P(Bₓ≠B*ₓ). For binary outcomes the Lean `coupling_gap`
proves the same statement from the four joint masses. If mismatches can occur
only on timing/support failures, ordinary leakage, or recording failures, their
union bounds the mismatch probability eₓ by εtime,ₓ+εleak,ₓ+εrecord,ₓ, capped at 1.
`three_event_bound` proves the union bound without independence. With one ideal
null receiver law q, `calibrated_gap` gives |p₁−p₀|≤e₀+e₁, capped at 1.

For the memory analysis, apply this coupling argument conditionally on every
pre-choice history Hᵢ. The needed allowance bounds the realized average of
e₀,ᵢ+e₁,ᵢ on one simultaneous calibration-success event. A uniform bound over
histories suffices. Independent pilot counts under an IID calibration model
estimate marginal failure rates, not automatically these conditional bounds.
Transport requires a device model establishing a uniform envelope, or a
separately justified high-probability average bound. This is an explicit
unmeasured premise; it cannot be inferred from the synthetic runs.

Nanoseconds do not themselves add to a probability gap. Independent timing
calibration must bound the probability that a trial lies outside the declared
support or crosses the chosen window boundary. For a timing error bounded by
j, disagreement with the ideal window requires an event within j of a boundary;
bound that event mass and the error-tail probability. Pulse amplitude outside
a support is not automatically a bad-trial probability. Converting optical
tails into a mismatch bound needs a detector/channel model and energy/count
calibration; absent that bridge εtime is unbounded for this purpose.

Similarly, an RMS voltage or an RF suppression ratio is not εleak. Use controls
to bound a selected-event probability change or a channel-distance/coupling
quantity, with an explicit transport assumption from controls to primary runs.
Shielding alone does not establish the bound. Correlated timing/leakage/record
failures can be combined by union bound, but overlapping estimates must be
interpreted as bounds on specified events rather than fitted residuals.

A sharper bound max(e₀,e₁) is valid if actual laws are
pₓ=(1−eₓ)q+eₓrₓ with the *same conditional good law* q. `contamination_gap`
proves it. Knowing only a bad-trial fraction does not identify q: failures can
select different preparations in the two settings. For example q=1/2 and
opposite coupling flips can give p₀=0.3,p₁=0.7 at e₀=e₁=0.2. The general budget
is 0.4, not 0.2. The shipped analysis uses the general coupling budget.

## Premise countermodels and surviving assumptions

| Removed premise | Normalized countermodel | What survives / consequence |
|---|---|---|
| Independent preparation / assignment | A fair latent bit λ gives B=λ and observational X=λ; do(X) leaves B fair. Alternatively μₓ=δₓ with fixed identity readout gives an actual preparation-dependent table. | Fixed receiver and complete records; observational conditioning or a changed source is not the intended intervention. |
| No receiver setting channel | Fixed one-point preparation and P(B=1|x)=(1+(2x−1)g)/2, 0≤g≤1. | Independent preparation, normalization, no selection; Δ=g. A wire realizes g=1 inside the light cone. |
| Complete outcomes / no future selection | Independent fair R,Y and S=1 exactly when R=Y. | Unconditional R stays fair; among S=1, R=Y. Selection succeeds with probability 1/2 for each Y. |
| Correct supports | An ordinary wire sends X, and B=X only after its causal arrival, while a too-short path or wrong timestamp labels that record early. | Fresh settings, complete counts and standard causality; incorrect geometry creates the alleged anomaly. |
| Fixed secure record | Sensor bit R is fair; stored value is overwritten with later X. | Sensor has zero influence; the saved table has gap 1. Storage/latch is part of B's support. |
| Fresh present setting given history | Xᵢ repeats a prior known bit and the receiver copies that bit. | No contemporaneous superluminal channel; balanced marginal settings can coexist with perfect predictability. |
| Constant effect, for a single-TV interpretation of the memory interval | Alternate dᵢ=1 and dᵢ=−1. | The memory test and coverage remain valid for d̄=0, but mean |dᵢ|=1. |

Sender conditional independence is sufficient for the classical construction,
but is not necessary for operational no-signaling; the quantum algebraic
identity is a distinct route. Neither a definite path nor a hidden-variable
dimension bound is needed. The receiver-access equivalence is only for the
explicit shared-bit and independent-sender examples, not for every hidden-path
or retrocausal interpretation.

## What the binary upper bound cannot say

Coarsening keeps every trial but loses information. For three raw outcomes,
(1/2,1/2,0) and (1/2,0,1/2) have selected-first-bin gap zero while full TV=1/2.
Thus the reported upper bound concerns the specified binary record only. A
positive binary gap lower-bounds full TV, by grouping the signed differences
and applying the triangle inequality; its converse is false. There is no
multinomial-TV upper bound or universal ontology exclusion in this study.
