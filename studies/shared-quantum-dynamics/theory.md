# What shared terminal sequences cannot identify

The scientific question is whether one ordinary set of operations predicts distinct
sequence families, and which resources the observations force. A physical fit gives
an upper compatibility example; a residual from a local fit is no exclusion proof.
The endpoint here is a **restricted, exact obstruction**, accompanied by measured
source-removal witnesses. It is not a proof that all leakage or memory is invisible.

## Shared operations and scope

For every word w = g₁…gₙ over Gi, Gx, Gy, reset prepares
ρ = (I + c Z)/2 in a fixed computational subspace. The plus effect in that subspace
is diag(l,u), with the complement I−E. In the deposited convention plus is nominally
outcome 1, so the ideal empty word has probability zero. Every occurrence of g
uses the same unitary U_g and depolarizing channel
D_g(ρ) = exp(−r_g)ρ + (1−exp(−r_g))Tr(ρ)I/2, r_g ≥ 0.
U_g is a free small rotation perturbation of I, X(π/2), Y(π/2).
Thus p(+|w) = Tr[E G_gₙ … G_g₁(ρ)]. No individual circuit gets its own gate.

The numerical qubit family is deliberately restricted: unital depolarization and
unitary errors, with physical preparation contrast and asymmetric readout. It is
not all qubit CPTP maps. Parameter bounds in `protocol.json` and `sqd_models.py`
are sensitivity choices, not independently measured calibration. Any rejection of
an individual prediction vector or this restricted family is not a rejection of
stationary qubits in general, much less of ordinary quantum mechanics.

The qutrit extension uses ρ_C ⊕ L, with Tr(ρ_C)+L=1, and applies

ρ_C′ = (1−λ) D_g(U_g ρ_C U_g†) + η L I/2,

L′ = λ Tr(ρ_C) + (1−η)L.

Here λ is loss **per primitive gate**, η is return per primitive gate; L is a
population, not a rate. The binary effect is diag(l,u,z), 0≤z≤1. Loss and return
are gate independent within this declared model, and leakage is incoherent.
Coherent coupling to the third level is excluded. Reset has L=0. These assumptions
are essential to the closed form below. The stored fits do not establish them
for the physical ion. The computational subspace is stipulated; its fitted
population is conditional, not a gauge-independent experimental measurement.

Kraus operators independently realize this map: embed sqrt(1−λ) times each qubit
Kraus operator in the computational block; add sqrt(λ)|L><j| and
sqrt(η/2)|j><L| for j=0,1, and sqrt(1−η)|L><L|. Their adjoint products sum to I₃.
They destroy block coherences, so this is also a valid map on arbitrary qutrit
inputs. States are positive and normalized; effects and complements are positive.
The implementation independently checks Choi positivity and trace preservation.

The alternating model instead adds a classical bit q=±1, initially uniform and
reset each circuit. X/Y rotation error includes qθ; every primitive, including Gi,
flips q. There are exactly two memory states, no circuit-label encoding. Each
conditional channel is CPTP and the classical update deterministic. The quasistatic
model keeps q constant within a circuit and redraws it on reset: it models a
specified distribution of slow control error, not a stationary memoryless qubit.
No acquisition-time drift law can be inferred from these condensed records.
Leakage and alternating memory are different resources, not a total ordering.

The unrestricted per-circuit reference sets p_w to each empirical frequency, giving
zero row-wise deviance. This is only a saturation reference. It is not a predictive
physical explanation. No probability-law deformation is included.

## Exact all-word observable coordinates

Let a=(l+u)/2, C=c(l−u)/2, s=λ+η, assume 0<s<1, and write

A = a − B,   B = (a−z)λ/s,   γ_g = r_g − log(1−λ).

Let t(w)=e_zᵀ R_gₙ…R_g₁ e_z, with R the SO(3) action of U.
Directly solving the population recurrence and propagating the traceless block gives

**p_w = A + B(1−s)ⁿ + C exp(−Σ_g n_g(w)γ_g) t(w).**

This identity holds word by word, without Clifford averaging, twirling, a small-error
expansion, or an asymptotic limit. It covers all GST words and all RB words in the
selected files, regardless of their length or order. Taking n=0 gives a+C as required.
The population is Lₙ = λ[1−(1−s)ⁿ]/s. The rotations retain their exact sequence
order; only isotropic attenuation depends on gate counts.

Consequently even perfect knowledge of every terminal word probability cannot,
in this model, uniquely recover λ, η, r_g and z. Hold a,C,s,B,γ and all U_g fixed.
For any allowed λ′>0 set

η′ = s−λ′,   z′ = a−Bs/λ′,   r_g′ = γ_g + log(1−λ′).

These operations give **identical probabilities for every word** but generally
different Lₙ. Preparation and computational readout can remain fixed. This is
stronger than changing a GST representation by a similarity gauge: the stipulated
projector |L><L| has different expectations. It is operational equivalence only
for the available terminal binary measurement interface.

### Sharp conditional physical boundary

For fixed observable coordinates above with 0<a<1, the necessary and sufficient
conditions within this fiber (r_g′≥0, η′≥0, 0≤z′≤1) are

λ_min ≤ λ′ ≤ λ_max,

λ_min = Bs/a if B≥0, and −Bs/(1−a) if B<0,

λ_max = min{s, 1−exp(−min_g γ_g)}.

Necessity follows separately from z′≥0 or z′≤1, η′≥0, and each r_g′≥0.
Sufficiency follows by substituting any point in the interval into the Kraus
construction. Thus this is a global, analytic boundary **on this fixed-coordinate
fiber**, not a global solution over all qutrit gate sets. B=0 permits λ′=0 with
η′=s and z arbitrary; the all-word formula then has no B term. s=0 is the ordinary
no-transfer case and is treated separately, not divided by zero. The implementation
uses interior positive-rate examples for certificates.

For the numerical fit box λ,η≤.001 and r_g≤.001, additionally intersect with
λ′≤.001, λ′≥s−.001, and λ′≥max_g[1−exp(.001−γ_g)]. These box boundaries are
not hardware bounds. Floating-point endpoint checks include a 1e-9 parameter
rounding tolerance; predicted equality is checked to 2e-9. The proof, not rounding
or optimizer convergence, establishes the interval.

Intervals evaluated at a fitted coordinate vector are **conditional compatibility
examples, not confidence intervals for actual leakage**. The coordinates themselves
have uncertainty and may not even all be identifiable in this finite design.
No claim of globally necessary nonzero leakage follows from a positive displayed
λ_min. In particular compatible qubit examples are retained.

## A bounded classical-memory realization

Encode the qutrit block state as two classical flag sectors of a qubit:
active = ρ_C, dormant = L I/2. On each gate, active follows the qubit channel,
a fraction λ of its trace moves to dormant I/2, and a fraction η of the dormant
trace returns to active I/2. Readout is E on active and zI on dormant. Reset is
(ρ,0). Induction gives the same recurrences and every terminal probability.
Only **two** classical states are used. They cannot encode arbitrary circuit labels.
This restricted flagged realization can preserve memory across noncommuting words.
It proves nonidentification of a specific incoherent leakage model versus a
specific ordinary qubit-plus-memory model, not all qutrits versus all memory.

Removing GST, removing RB, or adding arbitrarily many other words using only the
same operations and binary readout leaves this equivalence intact. The joint
collection can remove other observable freedoms; it cannot remove this kernel.
Under the same sampling law, both constructions induce the same distribution of
all counts. No statistical test has power exceeding its size along that direction,
for any shot budget. This is not merely a small Jacobian singular value.

## Small added observations that break the constructed pairs

A calibrated measurement of the fixed computational-subspace projector after one
gate gives L₁=λ directly. More generally, if s is known and n≥1,
λ=s Lₙ/[1−(1−s)ⁿ]. It separates any two distinct positive λ points on the fiber.
An uncalibrated third detector output does not suffice: its response can acquire
its own compensating nuisance. A monitor with known false-positive/false-negative
rates gives L=(observed−false_positive)/(true_positive−false_positive) when its
contrast is nonzero; monitor calibration uncertainty must propagate to the bound.
No such monitor is in the selected acquisition. A leakage monitor from another
device cannot supply it. A perfectly discriminating physical leakage projector
also distinguishes this qutrit from the dormant **computational-qubit** flag
realization, conditional on the latter not changing the monitor apparatus.

There is a separate scheduling ambiguity. The two-state alternating model equals
an external classical clock that applies +θ,−θ on consecutive primitive slots,
with random initial phase. Every available record stores only a word, not slot
start times. Under a one-slot-per-gate implementation the predictions coincide
exactly. This is a **conditional clock construction**, not a claim that the hardware
actually used that timing. A calibrated idle that advances the external clock by
one slot but leaves the gate-count flag unchanged separates the pair: between
two X gates, the alternating errors normally cancel; after the idle the clock
errors add. The exact separation for the stored synthetic example is checked
independently. A real idle with unknown dynamics would not establish that distinction.

## Observability, gauges, and status

The source prototype computes weighted observable Jacobian singular values at a
generic interior point, keeping all null directions. At that point GST resolves
15 of 18 numerical coordinates at the declared relative threshold; RB only 12.
Joint diagnostics are also provided. These thresholded local ranks depend on
parameter units and finite differences; they are not global theorems. Exact
redundancies include preparation/readout contrast, a common rotation about the
preparation/measurement axis, and the leakage fiber. Comparisons use probabilities
and certified operational equalities, never differences of fitted gate matrices
as evidence for physically different devices.

The sharp fiber formula specializes familiar leakage/return and readout
confounding. The flagged-channel representation and clock dilation are familiar
state-space constructions. The contribution is an explicit all-word boundary and
independently checked constructions attached to this matched measured collection,
with removal and power limitations. No new general dimension theorem, first
observation of correlated noise, or publication-level priority is claimed.
