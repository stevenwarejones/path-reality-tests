# General CPTP escalation: proof, scope and unresolved questions

The question remains whether shared ordinary operations explain the acquisition,
and which ordinary resource costs the combined observations identify. This
retrospective escalation keeps the original length split. The earlier finite-bank
spread remains an **inner** construction, not an outer prediction region.

## General stationary qubit baseline

Each of Gi, Gx and Gy has its own arbitrary qubit CPTP map. In the Pauli basis,
write a positive process matrix as χ = LL†, with complex lower triangular L.
The raw Kraus operators are A_j = ∑_i L_ij σ_i. Put S = ∑_j A_j†A_j and
K_j = A_j S^(-1/2) U_g, where U_g is only a nominal-coordinate reference.
Then ∑ K_j†K_j = I. Every CPTP qubit map is represented: apply its inverse
nominal unitary on the right, factor its positive Pauli process matrix, and S
is already I. Singular Choi matrices are allowed; singular S is rejected.
Since a TP Pauli process matrix has trace one, bounding each unscaled factor
entry in absolute value by one does not exclude a TP map.

A common unitary gauge aligns the preparation Bloch vector with +Z and the
effect Bloch vector with the XZ plane. We fit preparation length c ∈ [0,1],
effect eigenvalues l,u ∈ [0,1], and its polar angle θ ∈ [0,π]. This covers all
physical qubit preparation/binary-effect pairs, including uninformative limits.
No interpretation attaches to differences in gauge-dependent fitted coordinates.
The redundant 52-coordinate factorization is not a claim of 52 identifiable
parameters. Independent Kraus propagation, Choi positivity and trace preservation
check saved fits; local optimization is not a global compatibility decision.

## A full-class outer bound without fitting

Let Φ be the CPTP channel of any complete word w, with the same reset state ρ
and effect E on each circuit. No reset occurs **inside** ww. Define
p_j = Tr[E Φ^j(ρ)], including p_0 for the empty circuit. Then every stationary
qubit model satisfies

\[
p_2\leq \min\{1,(2\sqrt{p_1}+\sqrt{p_0})^2\}. \tag{1}
\]

**Proof.** Diagonalize E with eigenvalues l ≤ u. For u > l let |0⟩ be its
minimum-eigenvalue eigenvector and q_j = (p_j-l)/(u-l). The Bures angle to |0⟩
is α_j = arcsin √q_j. Triangle inequality and CPTP contractivity give
α_2 ≤ A(Φ²ρ,Φρ)+α_1 ≤ A(Φρ,ρ)+α_1 ≤ 2α_1+α_0.
The elementary capped-angle inequality implies
√q_2 ≤ min(1,2√q_1+√q_0). Hence
p_2 ≤ l+(2√(p_1-l)+√(p_0-l))².
The right side is 4p_1+p_0-4l+4√((p_1-l)(p_0-l)), a decreasing function of
l ∈ [0,min(p_0,p_1)]. Its maximum is attained at l=0, proving (1).
If u=l, p_0=p_1=p_2 and (1) is immediate. The inequality concerns observable
probabilities and therefore respects every gate-set gauge.

The Bures facts are standard; see [Gutoski and Wu, section 2.3](https://arxiv.org/abs/1011.2787).
This is an application of contractivity, **not a priority claim for a new
quantum-information inequality**. The specific measured certificate selects one
RB training word by minimum SHA-256 in each of five length bins before inspecting
its count bound. Together with GST's empty word, six one-sided binomial bounds
use Bonferroni tail 1/120. `sqd_repeat.py` verifies each outward rational endpoint
by an **exact integer binomial-tail comparison**; square-root bounds also use
integer arithmetic. SciPy supplies a starting guess, not the certificate.

This produces simultaneous 95% outer bounds on five unmeasured doubled-word
probabilities, for **all** stationary qubit CPTP models consistent with the six
count constraints. It does not require a finite model bank or a successful fit.
Dropping GST frees p_0; dropping RB frees p_1. In either case this particular
relaxation returns [0,1]. This is an ablation of the certificate, **not proof
that the complete single-source prediction sets reach one, nor proof of strict
contraction of their complete feasible regions**. Full-source constraints could
already be stronger. Distinguish a valid joint outer bound from an identified
increment over every implication of GST alone.

Independence is needed within each count row for its binomial model, not between
rows for the union bound. With unrestricted within-row dependence the reported
95% endpoints are unjustified; our dependence sensitivity therefore returns the
vacuous upper bound one. A qutrit three-cycle with outputs 0,0,1 violates (1), as
does a qubit with a two-state classical clock. These are ordinary dynamics that
violate the stated dimension/stationarity assumptions, not quantum mechanics.
A sequence-dependent reset or readout also invalidates the proof.

## All-word equivalence beyond uniform and unital leakage

Fix a computational two-dimensional subspace and one leakage level. On block
states (ρ_C,L), with Tr ρ_C+L=1, let every gate have independently chosen loss
λ_g and return η_g, and **arbitrary** computational qubit CPTP channel Φ_g:

\[
ρ_C'=(1-λ_g)Φ_g(ρ_C)+η_g L τ,\qquad
L'=λ_g\operatorname{Tr}ρ_C+(1-η_g)L,\quad τ=I/2. \tag{2}
\]

The binary effect is (E_C,z), 0 ≤ E_C ≤ I and 0 ≤ z ≤ 1. Reset lies entirely in
the computational subspace. Kraus operators consist of surviving computational
Kraus operators, loss jumps from each computational basis state, return jumps
preparing τ, and survival on the leakage level. This defines a full CPTP qutrit
extension that kills inter-block coherences. Coherent leakage is outside (2).

For 0 < κ ≤ 1 define the invertible linear similarity
T_κ(ρ_C,L) = (ρ_C+(1-κ)Lτ, κL). Set

\[
\begin{aligned}
λ'_g&=κλ_g,&η'_g&=η_g+(1-κ)λ_g,\\
Φ'_g&=\frac{(1-λ_g)Φ_g+(1-κ)λ_g R_τ}{1-κλ_g},\\
τ'_g&=τ+\frac{(1-κ)(1-λ_g)}{κ η'_g}[τ-Φ_g(τ)],\\
z'&=a+(z-a)/κ,&a&=\operatorname{Tr}(E_Cτ),
\end{aligned} \tag{3}
\]

Here R_τ(X)=Tr(X)τ. The transformed model allows an unknown gate-dependent
return state τ'_g. Φ'_g is a convex mixture of CPTP channels. Direct block
multiplication gives G'_g = T_κ G_g T_κ^(-1), the reset is unchanged, and
E'T_κ = E. Therefore **every word over every ordering and length has exactly
the same terminal probabilities**. This is not a finite-data interpolation.
Meanwhile the stipulated leakage population is L'_w=κL_w for every word,
and every loss rate becomes κλ_g. A two-state classical flag coupled to a qubit
also realizes (2): its inactive branch discards the qubit state, survival stays
inactive, and return prepares τ. After transformation return prepares τ'_g.
Only two memory states are needed, independent of word length or circuit label.

The same construction extends to an arbitrary initial physical return state
τ_g in each gate: add η_g(τ_g-τ)/(κ η'_g) to τ'_g in (3). If t_g is the
Bloch vector of that initial return state, its exact positivity condition is
||η_g t_g-(1-κ)(1-λ_g)d_g|| ≤ κ η'_g. Thus a nontrivial interval near one
exists around generic strictly positive return states with interior rates and
readout, not only at the maximally mixed-return slice. The implementation and
independent Kraus tests include this extension. The closed-form endpoint below
specializes to t_g=0; it is not claimed for arbitrary polarized returns.

The similarity alone does not guarantee physicality. Write d_g for the Bloch
vector of Φ_g(I/2). Assuming η'_g>0 and κλ_g<1, the exact return-state condition is

\[
(1-κ)(1-λ_g)\|d_g\|\leq κ[η_g+(1-κ)λ_g]. \tag{4}
\]

Also require η'_g≤1 and 0≤z'≤1. Put c_g=(1-λ_g)||d_g|| and
b_g=η_g+λ_g+c_g. The lower root in (4) is
2c_g/[b_g+√(b_g²-4λ_gc_g)] for c_g>0, and zero otherwise.
The rate condition adds κ≥1-(1-η_g)/λ_g when λ_g>0. Readout adds
κ≥(a-z)/a if z<a, or κ≥(z-a)/(1-a) if z>a (zero if z=a).
The maximum of these lower bounds over gates gives the physical interval
for **this constructed orbit**. Degenerate divisions are treated separately;
the implementation rejects η'_g=0 or κλ_g=1 instead of silently extrapolating.
When rates/readout are interior and return states are positive, an open interval
near κ=1 survives even for nonzero d_g. This is not a global minimum-leakage bound.

If return is required to remain I/2, the similarity fails for generic nonunital
Φ_g. The omitted term is precisely the polarization in (3). The synthetic audit
uses gate-dependent rates and amplitude damping on the **actual recorded words**.
GST's recorded `GxGxGx(Gi)^8192` responds to locking the return state; allowing
(3) restores equivalence on all 6,661 words. Those numbers are a mechanism control,
not an estimate of the experimental device, a powered separation at its exposures,
or a claim that these synthetic models fit the counts.

A calibrated computational-subspace monitor after one occurrence of any gate
with λ_g>0 measures λ_g versus κλ_g and breaks this orbit. It must have known
contrast and a fixed physical subspace. Another unknown binary detector can
transform with the model and does not supply that anchor. More GST words, more
RB words, intermediate **uncalibrated** operations transforming under the same
similarity, or more shots of the same terminal experiment cannot break it.
Known return-state polarization can break the nonunital extension, but unital
instances retain the obstruction. None of those calibrations is present in the
selected count files.

## What has and has not been certified

Equation (1) is a full-class conditional outer prediction certificate. Equations
(2)–(4) prove a physical operational equivalence in a substantially broader
block-leakage model than the first report. Neither result certifies rejection of
all stationary qubit models for this acquisition. The general CPTP fits are local
feasible examples, and their prediction failures do not settle that question.
Nor does a valid outer bound alone prove strict complete-region source-removal
gain. We have not solved the nonlinear all-source feasibility problem or proved
that such certification is mathematically unattainable.

Thus the exact impossibility result is **identification of the stipulated leakage
scale from arbitrary binary terminal words over model class (2) with unknown
return-state polarization and readout (3)**, on the stated physical interval.
What has merely not been found is a successful held-out stationary qubit predictor,
a global exclusion certificate, or a priority claim beyond the known GST gauge
and leakage/readout literature. This distinction remains part of the result.

## Recorded observable-basis obstruction

`observable_audit.py` uses fixed fiducials {}, Gx, Gy, GxGx for preparation and
measurement, with every concatenation present in GST. Although the empirical
4×4 matrix is invertible, its entrywise count box contains rational matrices
with opposite determinant signs. The segment remains in the box and respects
all duplicate-word equalities; continuity therefore certifies a singular matrix
in the relaxation. An H⁻¹H_g realization cannot be uniformly bounded by inverting
this box. This is a concrete limitation of that method at these exposures,
not a no-go theorem for the full physical, all-source feasible set.

For the synthetic locked-return perturbation on actual exposures, the audit also
sums independent-binomial KL divergence. Pinsker's inequality bounds any level-α
test's power by α+√(D(P_original||P_locked)/2). A uniformly valid composite test
must in particular control size at this locked-return point, so nuisance fitting
cannot evade the simple-pair ceiling. This is an analytic power bound evaluated
numerically for the declared synthetic pair; it does not say these models fit
the data. On the exact free-return orbit the entire sampling distributions agree,
so power equals size for any test, at any exposures, under the same sampling law.

## An in-sample stationary qubit compatibility example

The separate all-length fit supplies a physical stationary qubit model inside a
conservative simultaneous region for the full acquisition. We allocate α=.025 to
all 6,661 two-sided cell intervals (exact outward rational binomial endpoints)
and α=.025 to ten fixed family/length mean constraints (Hoeffding). This additional
audit was declared after observing the all-length cell result and is explicitly
retrospective. It must not replace the frozen predictive test. The source route
checks every probability and independent Kraus propagation; the smallest region
slack is orders of magnitude larger than the stated numerical tolerance.

This constructive example needs neither a leakage level nor a nontrivial memory
state. Consequently this declared region cannot imply a strictly positive cost
in either resource. The claim is restricted to that conservative region: no
optimal-likelihood certificate, universal goodness-of-fit claim, device population
measurement, or prediction of untouched data is inferred. In particular the
training-only general models continue to fail their long-GST predictive checks.

## Exact leakage lift of a compatible nonunital qubit model

There is a stronger construction than the nonzero-leakage similarity orbit.
Let Ψ_g be **any** qubit CPTP gate with Choi matrix J_g in the unnormalized
column-vector convention. Let τ=I/2, let d_g be the Bloch vector of Ψ_g(τ), and
choose, separately for every gate,

\[
0\leq λ_g<1,\quad λ_g\leq 2\lambda_{\min}(J_g),\qquad
\|d_g\|\leq η_g\leq1,\quad η_g>0.
\]

Define the computational survival channel and return state by

\[
Φ_g=\frac{Ψ_g-λ_g R_τ}{1-λ_g},\qquad
τ_g=τ+\frac{Ψ_g(τ)-τ}{η_g}. \tag{5}
\]

The Choi matrix of Φ_g is (J_g-λ_g I_4/2)/(1-λ_g), so the stated λ condition
is **necessary and sufficient for this particular replacement construction**
to be CP; trace preservation follows directly. The return-state Bloch vector
is d_g/η_g, so its stated condition is also exact. Full-rank J_g gives a strictly
positive loss interval. Boundary unitary channels need not permit positive loss.

Use these channels and return states in the block-leakage map (2), now with
τ_g rather than τ, and choose leakage readout z=Tr(E_Cτ). The trace-preserving
collapse Q(ρ_C,L)=ρ_C+Lτ obeys the exact **intertwining identity**

\[
QG_g=Ψ_gQ.
\]

Indeed, the coefficient of ρ_C is (1-λ_g)Φ_g+λ_gR_τ=Ψ_g, and the coefficient
of L is η_gτ_g+(1-η_g)τ=Ψ_g(τ). Reset with L=0 and the chosen effect then give
identical probabilities for **every** word, including arbitrary mixed-gate
orderings, despite strictly positive loss and return. This is a nonminimal
physical realization, not an invertible qubit/qutrit gauge transformation.
The two-state flag instrument supplies a third, finite-memory realization.

`lift_certificate.py` attaches this construction to the **all-length qubit model
inside the full-acquisition count region**, rather than to the earlier failed
held-out fit. It uses fixed fractions of the CP budgets, not residual tuning,
and verifies every recorded word. The loss, return and population numbers are
conditional constructions in the saved representation. The Choi budget is not
gauge invariant and is not a measured leakage rate, a global leakage maximum,
or a globally necessary cost. It is a sharp interval for (5) with Ψ_g fixed.

Consequently there are zero-leakage, positive-leakage and two-state-memory
ordinary realizations with the same all-word probabilities, including the
compatible measured prediction vector. Adding arbitrary terminal GST/RB words
cannot distinguish this class even with unlimited shots. A calibrated subspace
monitor distinguishes it immediately at the probability level: after one gate
it measures zero versus λ_g. This specifies one additional **observable type**,
not a claim that one noisy shot suffices.
