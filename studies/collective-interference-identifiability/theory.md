# Physical model and preserved original certificate

The physical class and inequality below remain current. The numeric certificate,
summary witnesses and robustness budget in this document describe the original
reviewed head. The expanded 2D region, full parity models and current budget are
in [revision.md](revision.md). The original exact artifacts remain unchanged.

## Declared physical alternatives

Prepare exactly one atom in each of three distinct external input modes.
A fixed contraction T maps these modes into detected external modes; loss
modes complete T to an isometry and then a unitary. Propagation is linear and
noninteracting, acts identically on hidden states, and is shared between
singleton and three-atom preparations. Detection attenuation is incorporated
in T and must also be independent of hidden labels. This permits arbitrary
external phases, inhomogeneity and path-dependent attenuation; it does not
permit a new channel on every shot or preparation-dependent transport.

The null C2 is the convex hull of hidden density operators with a partition
of the input labels into blocks of size at most two, where different blocks
occupy mutually orthogonal hidden subspaces. Within those subspaces the state
may be mixed or entangled; even correlations between blocks are allowed if
they retain the orthogonal block labels. Partition mixtures are allowed but
must share T. This includes fully distinguishable particles, pairwise pure
or mixed interference with a labelled third particle, and correlated mixtures
of different pair choices. It excludes generic three-particle states lacking
such a block decomposition. It is **not** every possible definition of
“no three-particle interference.” It is not a claim about entanglement depth.

For any hidden density matrix rho, define J(pi)=Tr(rho P_pi). Then J(e)=1,
J(pi^-1)=conj(J(pi)), and |J(pi)|<=1. The matrix
J(sigma^-1 tau) is PSD because its quadratic form is
Tr(rho A†A) for A=sum_tau c_tau P_tau. These moments are not independent free
weights. For a component with block partition H, J(pi)=0 outside the subgroup
preserving those blocks. For C2 in three particles, |H|<=2. This argument
covers the allowed correlated states without assuming a product Gram matrix.

For a collisionless external output o, let a_sigma be the product of the
three single-particle amplitudes for assignment sigma. Its probability is

    p(o) = sum_(sigma,tau) a_sigma* a_tau J(sigma^-1 tau).

(The inverse convention can be changed consistently without changing the
bound.) The Gram property ensures positivity; a unitary dilation and the
Born rule ensure normalization including losses and collisions. Grouping
assignments by cosets of H gives

    p(o) <= sum_cosets (sum_sigma_in_coset |a_sigma|)^2
         <= |H| sum_sigma |a_sigma|^2 = |H| p_D(o).

Thus P_C2(B)<=2 P_D(B) for any collection B of collisionless outputs.
For arbitrary hidden states the same proof with H=S3 gives P(B)<=6 P_D(B).
For an n-particle block partition the analogous factor is the product of
block factorials. These are elementary interference bounds; novelty of that
inequality is not claimed.

## Actual observation channel and crop tails

The final image measures occupation parity in each 2D site. Define B as
exactly three occupied sites in one y row. With exactly three input atoms,
no false positives and passive loss, B implies three surviving atoms at
three distinct 2D sites. A same-site pair has even parity and cannot produce
B. We count B per initial-selected prepared shot, including final losses and
collisions in the denominator. We never treat a zero pixel as certainly an
empty physical site, and never divide by observed survival.

The distinguishable probability of B is bounded above by the probability
that all three independently propagated particles occupy the same row,
without demanding distinct x positions:

    P_D(B) <= D = sum_(y=8..19) q1(y) q2(y) q3(y).

Here row indices are array indices in the NC files, not absolute lattice
coordinates. For the selected singleton record, rows 9..20 are observed.
The arrays for the other inputs are copies rolled by +1 and -1 in y.
We explicitly assume that this copying correctly transfers the physical
single-particle detected probabilities. Therefore q1(y)=q(y),
q2(y)=q(y-1), q3(y)=q(y+1).

The reference calibration does not observe rows 7 or 8. Treat q(7),q(8) as
free nonnegative variables. Define r to include true loss and all other
unobserved rows. The measured empty/outside probability is
q(7)+q(8)+r. The model has 15 variables:

    [q(7), q(8), q(9), ..., q(20), r],  sum = 1.

No zero count outside the crop is treated as a measured zero probability.
All contributions to the many-particle observed rows are consecutive
triples of these variables. Allowing the tail mass to lie in rows 7 and 8
is a conservative enlargement even when the true lost mass cannot return.

The bridge assumes the same detected response after translation, including
x-dependent attenuation. Arbitrary hidden-dependent detection, extra atoms,
weak interactions and drift are not automatically covered by T.

## Finite-sample coverage and exact numerical verification

One 931-shot multinomial record supplies 12 row categories and one
empty/outside category. Marginally, each count is binomial. Construct
two-sided Clopper–Pearson intervals with tail error 1/4160 for every
category. Their union failure probability is at most 26/4160=1/160.
For 93 bunching events in 2,999 prepared shots construct a one-sided lower
interval with error 1/160. No independence across categories or across the
two families is needed for the union bound. IID observations within each
record are needed for the binomial models.

The selected n=3 statement therefore has failure probability at most 1/80.
We reserve the same budget for each of the four inspected particle numbers
2,3,4,5, so choosing n=3 after inspecting that family costs at most 0.05.
The certificate does not license an arbitrary further scan over outcomes,
times, archives or alternative witness definitions. Other exploratory
calculations are not additional confirmed discoveries.

Floating-point inverse-beta calculations only propose interval endpoints.
Endpoints are rounded outward to a rational grid and then verified by exact
integer evaluation of binomial tails. The confidence guarantee consequently
does not depend on floating-point quantile accuracy.

To bound D introduce z_i=q_i q_(i+1) and t_i=z_i q_(i+2). Four McCormick
inequalities enclose each product on a box. Together with normalization and
the interval for q(7)+q(8)+r, they form a linear outer relaxation. Every
physical model maps into it; infeasible points in the relaxation are not
needed for the argument. Binary splits partition the entire initial box.

The solver supplies nonnegative rational inequality multipliers lambda and
a free normalization multiplier mu. With c the objective coefficients,
residual v=c-A^T lambda-E^T mu, the checker evaluates

    U = lambda.b + mu + sum_j max(v_j lower_j, v_j upper_j).

This is a valid upper bound for **any** nonnegative multipliers, even if
rounded multipliers are not dual feasible. The residual term repairs them
exactly. The four leaf boxes cover the whole region; their maximum is below
0.01166. No approximate eigenvalue or LP success flag is used as a proof.

## Inner examples and the resource interpretation

The committed physical witness uses the exact empirical singleton 2D cell
probabilities p(y,x), translated for the three inputs, as squared magnitudes
of T. Row-dependent rational unit phases cancel the column overlaps.
Entries have the form sqrt(p) times a rational complex unit number.
Exact square-root intervals and Gershgorin bounds prove
I-T†T is positive definite, with slack above 0.04189. Appending
sqrt(I-T†T) as loss amplitudes gives a physical isometry.

With identical hidden states, row-output permanents are sums of positive
square roots times a common phase. Exact interval summation over distinct
x triples gives the reported bunching interval. It lies between the checked
lower endpoint and the valid upper endpoint 0.04. It fits the confidence
constraints actually used. It is not certified against every three-atom
image frequency or against a particular lattice Hamiltonian.

With mutually orthogonal hidden states, the same T reproduces all empirical
singleton cell frequencies but gives the distinguishable probability
0.00462813088179467. `sources.py --check` compares those cell frequencies
back to the original record. This is a complete-singleton-source witness,
not merely a compatible row marginal.

If rho=(1-w)rho_C2+w rho_other, using the same fixed T for both components,
then P(B)<=(2+4w)D. Hence every such decomposition in the joint region has
w>=0.007005789. A decomposition always exists at w=1. The parameter is the
minimum fraction outside the declared cluster class; it is not a Schmidt
number, triad-phase estimate or universal computational resource measure.

## Robustness, adversaries and power

If actual input-j row transport differs from its translated calibration by
TV at most delta_j, coupling independent draws gives
D_actual <= D_calibrated + sum_j delta_j. If all other effects can increase
P(B) by at most epsilon, a sufficient exclusion condition is

    2 sum_j delta_j + epsilon < 0.00032675.

This is a sufficient robustness budget, not a sharp necessary experimental
precision or proof that every larger error admits a fit. With epsilon=0,
the sum must be less than 0.000163375; if equally divided, each delta must
be less than about 0.00005446. The archive audit does not establish this.
Interactions or hidden-dependent loss may be bounded by an event-probability
allowance only after a physical derivation, not by assigning a convenient
number to epsilon.

There is a physical synthetic drift counterexample: on each shot choose
one of two rows with probability 1/2 and apply a mode permutation mapping
three distinct inputs into three distinct x sites in that row. Give the
particles orthogonal hidden labels. Singleton row probabilities are 1/2,
so multiplying averaged singletons gives D=1/4, yet P(B)=1. Every shot uses
ordinary distinguishable particles and unitary propagation. The error is
multiplying averages of a common fluctuating channel. Our fixed-channel
assumption excludes this example; no claim about arbitrary drift survives.

For the frozen calibration ceiling 0.02332, an exact binomial rejection
threshold at tail error 1/160 is 92 bunching events in 2,999 shots. Its null
tail is approximately 0.006019. Conditional power at the observed frequency
93/2999 is approximately 0.5564, and at the physical inner-example probability
it is only approximately 0.03739. These are binomial probabilities conditional
on that calibration ceiling, not end-to-end power including calibration
resampling, nor evidence of non-IID robustness. `diagnostics.py` reproduces
them. Analytic coverage above, not a Monte Carlo success rate, justifies the
stated composite-null error control.
