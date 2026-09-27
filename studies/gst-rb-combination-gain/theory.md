# What is certified

Let M be **all** stationary qubit gate sets with a normalized positive state,
three CPTP maps G_i,G_x,G_y and one binary POVM {E,I−E}. There is no proximity to
target gates, unitality, Pauli-diagonality, gate-independent error, or leakage
assumption. For a deposited word w=(g₁,…,g_L),

\[
p_w=\operatorname{Tr}[E G_{g_L}\cdots G_{g_1}(\rho)].
\]

One preparation and final measurement are shared across both families and every
word. There are no intermediate measurements. Repeated occurrences of a primitive
mean exactly the same map. Global basis rotations are gauge: we can take the
initial Bloch vector along +z and the effect in the xz plane. This does not restrict
observable predictions. No difference in gate coordinates is itself evidence of
separation. All certificates compare probabilities.

## Confidence region and source deletion

There are m=6661 count rows. Conditional on their fixed words and denominators,
assume mutually independent Bernoulli shots with a fixed probability p_i in each
row. A row gets two-sided Clopper–Pearson error .025/m. Each upper endpoint a/10¹⁰
is rounded outward and checked using the exact integer inequality

\[
(80m)\sum_{j=0}^{k}{n\choose j}a^j(10^{10}-a)^{n-j}
\le (10^{10})^n.
\]

The lower endpoint is one minus the corresponding upper endpoint for n−k.
Boundary counts use endpoints zero/one as appropriate.

Twelve additional two-sided constraints cover each family's overall mean and its
five fixed length bins. Their target is the exposure-weighted mean
μ_A=Σ_{i∈A}n_i p_i/N_A, not an iid approximation with a common circuit probability.
Hoeffding for independent, possibly nonidentical Bernoulli variables gives

\[
\Pr(|\widehat\mu_A-\mu_A|>r_A)\le 2e^{-2N_A r_A^2}.
\]

Choose r_A=a_A/10¹⁰ such that exp(2N_A r_A²)≥960. A **rational positive Taylor
partial sum** certifies this inequality, so the bound for each group is at most
1/480=.025/12. Overlapping groups require no additional independence argument:
the union bound yields simultaneous coverage at least .95 for every cell and
group. This is a conservative confidence region, not a likelihood acceptance
region or a guarantee of a satisfactory overall residual distribution.

Define C_G,C_R from their respective cells and groups. C_J=C_G∩C_R. Source deletion
only removes constraints; it does **not** redistribute confidence budgets or
change the physical class. The construction was fixed for this retrospective
follow-up before searching the new witnesses; earlier fits and observations were
already inspected. This is not preregistered validation. Witness searches may use
all observations: they prove existence inside the specified regions, not
out-of-sample estimation. We make no post-selection novelty or p-value claim.

## Observable and complete-class outer bound

Set q=(1+μ_R−μ_G)/2. It is an ordinary binary experiment: choose R or G with equal
probability, choose a recorded word with probability n_i/N_F, run it, and score
plus for R or minus for G. The choice affects the classical scoring only; the
state, gates and terminal physical measurement remain shared. For this mixture
of recorded contexts, all models in C_J satisfy

\[
q\le U_J=(1+U_R-L_G)/2.
\]

This is an exact linear dual certificate: multiply the upper-R mean inequality
and the lower-G inequality by 1/2 and add. It remains valid over the larger set
of all probability vectors; hence it bounds **every** physical gate set in M,
including ones never visited by an optimizer. It need not be sharp over M.

The saved normalized Kraus models establish

\[
\varnothing\ne C_J\cap M,\qquad
\sup_{C_G\cap M}q\ge q_G^- > U_J\ge\sup_{C_J\cap M}q,
\]

and independently the same strict inequality with G replaced by R. Here q⁻ is a
certified lower prediction. These are nonvacuous source-removal statements about
the **full** class, not source-removed relaxations or finite-bank widths. We do
not claim to optimize either supremum. The direct aggregate bound is elementary;
the data-specific content is that explicit shared CPTP witnesses obey **every**
retained count/group constraint while crossing it.

## Exact physical witness semantics

Each saved binary64 parameter is an exact dyadic rational. For each gate form a
4×4 lower triangular L in the Pauli basis with exact scale factors 1,1/100,1/1000
as specified in `sqd_general.SLOTS`. Let A_j=Σ_i L_ij σ_i and S=Σ_j A_j†A_j.
The witness is

\[
K_j=A_j S^{-1/2}U_g,\quad G_g(X)=\sum_jK_jXK_j^\dagger.
\]

Use exact nominal U_i=I, U_x=(I−iX)/√2, U_y=(I−iY)/√2. Directed intervals verify
S is positive definite. Consequently ΣK_j†K_j=I **exactly**, and complete
positivity follows from the Kraus form. A floating eigenvalue tolerance is not
the physicality proof. General CPTP maps are representable: use their positive
process matrix factorization; an already TP process has S=I. The nominal U is
just an invertible change of coordinates, not a near-ideal constraint.

The state is (I+cZ)/2, 0≤c≤1, and
E=(l+u)I/2+(l−u)(sinθ X+cosθ Z)/2, with l,u∈[0,1]. Its eigenvalues are exactly
l,u. These inequalities are checked. Rank-deficient factors are allowed in the
class; only S must be nonsingular for this normalization representation.

The GST witness was discovered by adding a depolarizing replacement component
of nominal weight 3×10⁻⁵ to the prior all-length joint fit. The RB witness uses
nominally identical identity/depolarizing maps with replacement weight .00014,
c=.995, and a two-variable readout LP. Factorization adds a tiny positive process
regularizer. The exact saved factors, rather than these rounded construction
recipes, define the certified maps. No claim is made that the uncalibrated RB
witness implements accurate nominal X/Y rotations. Such a calibration would be
an additional constraint absent from this declared class.

## Numerical enclosure

`interval_physics.py` uses mpmath 1.3.0 directed intervals at 60 decimal digits.
For positive 2×2 S it evaluates the explicit identity

\[
S^{-1/2}=\sqrt{\operatorname{Tr}S+2\sqrt{\det S}}\,
(S+\sqrt{\det S}\,I)^{-1}.
\]

It encloses the exact affine Bloch coefficients and compares each with the
binary64 coefficient used by the C++ propagator, obtaining an outward rational
entry error ε_M. Effect coefficients receive ε_E. This also covers the floating
nominal-unitary and factor-scale conversions used by that propagator.

Assume IEEE binary64 round-to-nearest with gradual underflow and no fast-math.
The compiled kernel uses at most seven operations for each affine coordinate;
γ₇=7u/(1−7u), u=2⁻⁵³. Set
D=ε_M+γ₇(1+ε_M)+32·2⁻¹⁰²². The final term conservatively covers underflow in these
few operations. CPTP trace-distance contraction implies the Bloch linear part
has Euclidean operator norm at most one. If the previous Bloch error is e, the
new coordinate error contributed by coefficient error and rounding is bounded
by (4+3e)D. Since √3<7/4,

\[
e_{j+1}\le (1+A)e_j+B,\quad A=21D/4,\quad B=7D.
\]

With exact initial dyadic Bloch vector and LA<1,

\[
e_L\le LB(1+A)^L\le LB/(1-LA).
\]

The last inequality follows by comparing binomial coefficients with L^k.
The exact effect's Bloch coefficient norm is at most 1/2. With
D_E=ε_E+γ₇(1+ε_E)+32·2⁻¹⁰²² the probability error is at most

\[
e_L/2+(4+3e_L)D_E.
\]

All operations in this accumulation are rational. Each model's bound is below
the common outward envelope 10⁻⁹ used for **all** cells, means and q. Exposure
weights are positive and sum to one; q's two signed weight sums have total
absolute weight one. Thus the same envelope applies to q. Clipping a floating
probability to [0,1] cannot increase its error. Prefix reuse changes no word's
roundoff bound. A platform change is checked against the common envelope in the
full-source route, including the difference from saved probabilities.

The verifier converts saved probabilities to exact rational numbers, subtracts
the envelope from every constraint slack and prediction separation, and checks
all inequalities. The full-source route additionally recomputes all probabilities
with the prefix kernel and with a separate complex density-superoperator
implementation, using powers of deposited compressed germs. That independent
floating check is a cross-check; the interval/rational calculation supplies the
numerical guarantee. The offline route verifies the saved propagation artifacts
and all exact inequalities; reconstructing the words requires the external cache.

## What dependence or model changes invalidate

Under unrestricted dependence between shots the Hoeffding and binomial coverage
claims fail. Repeating a single Bernoulli draw N times is a concrete counterexample:
the sample mean remains 0 or 1, however large N becomes. With no recorded shot
chronology or justified dependence bound, this study supplies **no** .95
unrestricted-dependence separation. The deterministic inequalities about the
chosen regions still hold; their statistical interpretation does not.

Stationarity, shared preparation/readout and the physical qubit model are
conditional. A mislabeled primitive, context-dependent reset or gate, drift,
leakage, or memory can invalidate their use for this acquisition. The result
neither measures such resources nor excludes ordinary quantum dynamics outside
M. It does not repair the previous study's held-out long-sequence failures.
