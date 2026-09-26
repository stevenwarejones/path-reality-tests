# Phase intervention with complete outcomes and occupation calibration

The baseline null is the quantum P/Q-dephased class. Rejecting it in a working
heralded Mach–Zehnder interferometer is expected: it establishes arm coherence.
The baseline is a calibrated, loss-complete template. The substantive conditional
ontology test below instead allows arbitrary phase action on the occupied arm.
All calculations are prospective; no apparatus or achieved power is reported.

## Quantum predictions

For a region projector P, Q=I−P, fixed preparation ρ and fixed downstream effect
E, the intervention Vθ=Q+exp(iθ)P gives

\[
p_\theta=A+2\operatorname{Re}(e^{i\theta}c),\quad
A=\operatorname{Tr}(EQ\rho Q)+\operatorname{Tr}(EP\rho P),\quad
c=\operatorname{Tr}(EP\rho Q).
\]

The contrasts are p₀−pπ=4 Re(c) and pπ/₂−p₃π/₂=−4 Im(c). The dephased state
PρP+QρQ gives A at every setting. Fixed preparation mixtures and fixed detection
loss preserve this invariance. A coherent state can also give zero contrast at
an insensitive detector.

For the balanced lossy family, the complete table is

\[
p_\theta(\pm)=\eta(1\pm v f_\theta)/2,\quad
p_\theta(\varnothing)=1-\eta,\quad f=(1,0,-1,0).
\]

Its selected-bin contrast is ηv and its quantum which-region probability is 1/2.
The [formal companion](https://github.com/stevenwarejones/ontology-separation/blob/91590c324c444c836efe2b7d7049784a81c6a475/docs/PHASE_INTERVENTION.md)
contains the density-matrix and complete-POVM realization and the exact
local-model boundary. The linked revision is under review in PR #85.

## Local-phase definite-region null

A finite ontic state λ has one region r(λ)∈{P,Q}; its preparation distribution μ
is independent of the setting. Responses Rθ(o|λ) are normalized on the complete
outcome space. For λ in Q, Rθ is independent of θ. For λ in P it is arbitrary,
including arbitrary phase-dependent loss. The occupation premise identifies
μ(P)=w with the which-region measurement probability for the same preparation.
Then, for every outcome including failure,

\[
p_s(o)-p_t(o)=\sum_{\lambda\in P}\mu(\lambda)
[R_s(o|\lambda)-R_t(o|\lambda)],\qquad |p_s(o)-p_t(o)|\leq w.
\]

This follows because each response difference is in [−1,1]. Unlike the dephased
null, the class admits nonzero fringes. The balanced quantum family exceeds its
bound exactly when ηv>1/2. At and below that boundary a two-region response model
reproduces the entire table: if η≤1/2 its P response is the lossy table (2η,v)
and its Q response always fails; if η≥1/2 the responses are (1,2ηv) and
(2η−1,0), each with preparation weight 1/2. These parameters are physical precisely
on the stated side of the boundary.

A violation rejects the conjunction of definite region, cross-context occupation
matching, invariant Q response and setting independence. Occupation calibration
estimates an operational probability; it does not prove the ontic identification.
Context-dependent occupation models, Bohmian-style models and models with
phase-sensitive fields in an empty arm can evade these premises. Q-response
invariance is stronger than ordinary spatial locality: an empty-arm field can
carry phase information locally to the later recombination point.

## Acquisition and physical correspondence

The candidate apparatus is a heralded two-mode Mach–Zehnder interferometer with
both output ports and failure recorded. The preparation, phase convention,
selected bin, nuisance budget, allocation and stopping rule are fixed using
independent calibration. Training and confirmation data have separate roles.
The [physical correspondence](physical-bridge.md) describes additional spatial
projection and propagation approximations if a spatial mask replaces two modes.

Phase settings 0, π/2, π and 3π/2 are interleaved with a fifth, which-region
calibration context on the same preparation. In that context a region readout
replaces recombination. Its P fraction estimates w. A complete calibrated P/Q
readout, or a separately justified error budget for its missing and mistaken
outcomes, is needed. Calling an undetected occupation trial Q would bias the
bound. The region measurement can disturb the system: no simultaneous path
record is claimed for the interference trials.

Every eligible source herald contributes to its assigned setting's denominator.
The outcome space includes selected detection, other detection, no detection,
and a fixed assignment for multiple detections. Detection-conditioned fringes
are not the absolute probabilities used here. Backgrounds and contamination
enter the forward model or nuisance allowance; background-subtracted fractional
pseudo-counts are not binomial observations.

Source-off, transmission, detector, settling and deliberately randomized-phase
controls quantify preparation/readout changes. A blocker is not automatically a
dephasing operation: it changes boundaries. Logs retain assignment, acquisition
order, timestamps, setting calibration and every outcome.

If Q responses vary by at most ℓ between the compared settings, the bound becomes
w+(1−w)ℓ≤w+ℓ. An entrywise behavior discrepancy δ contributes 2δ. Thus B=ℓ+2δ
is a conservative contrast allowance, with any region-readout bias included
separately or in B. Unbounded setting leakage leaves the exclusion unresolved.
Moving the element to the other arm gives the separate bound 1−w. The smaller
bound min(w,1−w) applies to a *common* observable contrast only under an additional
identification of the two placements' contrasts and complementary occupations. In the ideal quantum
model V_Q(θ)=exp(iθ)V_P(−θ), so a fixed preparation/readout gives equal absolute
opposite-phase contrasts after relabeling; physically moving a device still
needs that correspondence and its nuisance bound.

## Confidence intervals and decision rules

`design.py` retains the original Hoeffding `classify` rule and supplies
`interval_classify` with either Hoeffding or Clopper–Pearson (CP) intervals.
For Nᵢ complete independent Bernoulli trials at setting i and count Kᵢ, let
[Lᵢ,Uᵢ] be a two-sided interval with individual error α/m. Here m=4 for the
coherence baseline and m=5 when occupation calibration is included. The union
bound gives simultaneous coverage at least 1−α, with unequal Nᵢ allowed.

For CP, each tail has probability allocation α/(2m): the lower endpoint is the
Beta(Kᵢ,Nᵢ−Kᵢ+1) quantile at that probability, and the upper endpoint is the
Beta(Kᵢ+1,Nᵢ−Kᵢ) survival quantile. At Kᵢ=0 the lower endpoint is 0; at Kᵢ=Nᵢ
the upper endpoint is 1. Inversion of binomial tails gives coverage at least the
nominal level; numerical evaluation uses SciPy beta quantiles.

For each opposite-phase pair (i,j), define

\[
C_{ij}=\max(0,L_i-U_j,L_j-U_i).
\]

The baseline rejects if max(C₀₂,C₁₃)>B. The occupation-calibrated test rejects if

\[
\max(C_{02},C_{13})>U_w+B.
\]

Equality is inconclusive. On simultaneous coverage the inequalities cannot hold
under the corresponding null, so false rejection is at most α. A separate
statistical calibration failure αcal adds to this bound. Selection of additional
bins or regions changes the multiplicity problem.

The untruncated Hoeffding radius is rᵢ=√[log(2m/α)/(2Nᵢ)]. In symmetric-radius
notation the occupation rule is max contrast−2rα>ŵ+rw+B; for unequal counts
the selected pair uses rᵢ+rⱼ. The CP implementation uses asymmetric endpoints
instead of replacing them with a common radius.

## Random assignment and blocks

Conditioning on realized per-setting counts preserves these intervals when
setting assignment is independent of the trial outcomes and, conditional on
assignments, outcomes are independent Bernoulli with a fixed probability pᵢ
for each setting. Given a full assignment sequence, the joint count law factors
as ∏ᵢ Bin(Nᵢ,pᵢ). Every sequence with the same Nᵢ has this same count law, so
averaging over those sequences gives the same conditional law given the counts.
Conditional coverage therefore implies unconditional coverage. A missing setting
has no usable interval and the implementation returns an error, not rejection.

Fixed balanced blocks with an independently randomized order have fixed counts
and the same factorization under these outcome assumptions. Randomizing a block
does not remove within-block correlation or drift. The binomial guarantees here
exclude correlated/drifting blocks and outcome-adaptive stopping. Power tables
use fixed counts per setting; per-trial randomization has a conditional power
certificate at its realized counts, not the fixed-count total-trial guarantee.

## Conditional power

The Hoeffding baseline uses hγ=√[2 log(8/γ)/n] and sufficient condition

\[
d>B+h_\alpha+h_\beta,\qquad
n>\frac{(\sqrt{2\log(8/\alpha)}+\sqrt{2\log(8/\beta)})^2}{(d-B)^2}.
\]

Strict integer rounding is checked explicitly. `certified_power` supplies a
sharper *sufficient certificate*, not exact power: at specified anticipated pᵢ,
integer binomial quantiles with tail allocation β/(2m) give a simultaneous count
envelope with probability at least 1−β. Confidence endpoints are monotone in
counts. Their worst-case contrast lower bound and occupation upper bound over
that envelope certify rejection throughout it. `interval_plan` searches for
and rechecks a sufficient integer n; discreteness precludes a minimality claim.
A zero reported power lower bound means no certificate, not zero actual power.

With α=0.01, β=0.10, B=0.02 and anticipated phase probabilities
(0.5+d/2,0.5,0.5−d/2,0.5), the total trials are:

| Contrast d | Original Hoeffding power bound | Hoeffding intervals + binomial power envelope | CP intervals + binomial power envelope |
|---|---:|---:|---:|
| 0.80 | 288 | 172 | 144 |
| 0.10 | 27,364 | 21,648 | 17,272 |
| 0.05 | 194,588 | 154,660 | 123,200 |
| 0.03 | 1,751,288 | 1,390,956 | 1,108,884 |

At d=0.03, CP reduces the original total by 36.68%. With a lower-probability bin,
p=(0.115,0.10,0.085,0.10), CP needs 397,300 trials versus 998,548 for Hoeffding
intervals with the same binomial power envelope, a 60.21% reduction. CP adapts to
the binomial variance through its exact tails; an estimated variance bound is
not substituted. Contrast alone does not specify either of these power curves.

For the local-phase null, w=0.5 and B=0.02, CP uses five equally sampled contexts:

| Efficiency η | Visibility v | Total trials sufficient for at least 90.00% power |
|---:|---:|---:|
| 1.00 | 0.60 | 40,010 |
| 0.90 | 0.90 | 2,350 |
| 0.80 | 0.80 | 16,045 |
| 0.80 | 0.70 | 155,275 |

At 10,000 trials in each of the five contexts (50,000 total), the checked grid
at B=0.02 certifies at least 90.00% power for: η=0.60 with v=1.00; η=0.70 with
v≥0.90; η=0.80 with v≥0.80; η=0.90 with v≥0.70; and η=1.00 with v≥0.60,
on the grid spacing 0.10. No assertion is made between grid points.
[design.json](results/design.json) contains the full B∈{0,0.01,0.02,0.05} grid.
These are plausible target parameters, not an established realistic operating
region: a source's achieved efficiency, visibility, nuisance and independence
still need measurement. An uncertain pilot signal adds its own failure
probability to any claimed power assurance.

## Scientific scope and prior experiments

The closest reviewed experimental architecture is the heralded beam-splitter
anticorrelation and Mach–Zehnder interference experiment of
[Grangier, Roger and Aspect (1986)](https://doi.org/10.1209/0295-5075/1/4/004).
Its high fringe visibility alone does not establish ηv>1/2 in absolute
herald-denominated probabilities. [Hardy (2012), Section 3](https://arxiv.org/html/1205.1439v3)
gives a closely related phase-on-the-other-path argument. This design adds an
explicit restricted response class, complete outcomes and reproducible finite
statistics. It is a formal/statistical re-analysis of established physics,
not a new experimental mechanism. Development stops at this conditional template;
no new apparatus campaign or general trajectory exclusion follows from it.
The [literature comparison](literature-comparison.md) distinguishes this premise
set from ontic indifference, modular-variable arguments and contextuality tests.
