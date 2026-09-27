# Response to the PR20 review

The reviewed certificate at `c878cb26e6c89192d3e983dd240579d90feceb32` is
preserved. The new result strengthens conditional inference from the selected
records; it does not establish apparatus robustness or consequential novelty.

## 2D uncertainty and the exact event

Write the reference detected cell probability as `p(y,x)=q(y) r_y(x)`.
For each measured row, `r_y` is a normalized distribution over the 12 observed
x positions. When `q(y)=0`, its conditional distribution is arbitrary. The
translated inputs use reference rows `y`, `y−1`, `y+1`. For three such
conditional distributions a,b,c, the probability of distinct x outcomes is

    R(a,b,c) = 1 − sum_x (a_x b_x + a_x c_x + b_x c_x)
                + 2 sum_x a_x b_x c_x.
    D_B = sum_(y=8..19) q(y−1) q(y) q(y+1) R(r_(y−1),r_y,r_(y+1)).

This is the exact distinguishable probability of the observed event. It sums
ordered assignments once through independent draws; there is no missing or
extra factorial. The existing physical proof gives `P_C2(B) <= 2 D_B` and
`P_quantum(B) <= 6 D_B` for the declared fixed-channel null and alternatives.

Allocate calibration failure `1/160` as follows:

| Constraints | Marginal tail error | Union failure |
|---|---:|---:|
| 12 row probabilities + empty/outside | 1/8320 on each of two sides | 1/320 |
| 12 x categories in each of 12 rows, conditional on the observed row total | 1/92160 on each of two sides | 1/320 |

Conditional on a row total N_y, the cell counts have a multinomial law.
Clopper–Pearson coverage holds for every N_y, hence also after averaging over
that random total. N_y=0 imposes `[0,1]` for all its conditional cells. The
union bound requires no independence among categories. The copied records
are counted only once. All proposed endpoints are rounded outward and
checked by integer binomial tails.

Rows 7 and 8 remain unmeasured probability masses bounded only by the
empty/outside category and total normalization. Their x distributions are
arbitrary simplexes. All other unobserved mass is retained in residual loss;
unobserved x sites cannot contribute to the selected crop because the bridge
translates only y. Nothing outside the crop is inferred to have zero physical
probability from its absence in the images.

For each observed y, three pair products per x and one triple product per x
receive McCormick outer constraints. Three simplex equalities are retained.
The objective R is linear in those products. The resulting exact repaired
LP dual bounds R by a weight w_y, also capped at the physical bound one.
Shared conditional rows are allowed to optimize separately at this stage,
which enlarges the feasible set. All q constraints and normalization are then
retained in a weighted consecutive-triple relaxation for D_B. Three rational
branch leaves prove `D_B <= .0112`. Four leaves prove the comparison row-only
bound `.0119` on exactly the same revised q region. Every physical model maps
into these relaxations, so decoupling cannot remove an adversarial model.

The alternative scalar-bunching design uses the old `1/160` bunching tail and
has margin `.00124675`. Its sum-TV allowance `.000623375` is sufficient, not
necessary. The older row-only design had margin `.00032675`; its intervals
are narrower because no calibration budget went to conditional cells. These
comparisons do not prove any of the upper bounds optimal.

## One complete selected-source region

The final design keeps that same full singleton region. For the higher-order
source it allocates `1/320` to the bunching lower endpoint and `1/320` to all
parity-pattern intervals. There are

    K = sum_(m=0..3) binomial(144,m) = 497785

possible patterns under exactly three prepared particles, passive loss and
parity readout in the 12×12 crop. Every pattern, including every zero-count
pattern, receives two-sided tail error `1/(640 K)`. This is a fixed alphabet,
not an allocation over only the observed bins. Initial images are the fixed
selected input; their preparation mechanism is an assumption rather than a
new inferred source. No crop or accepted-shot denominator is changed.

The two source budgets total `1/80` at n=3. Reserving the same amount for
each n in {2,3,4,5} gives familywise .05 for this declared analysis family.
Derived event choices or optimization objectives within this simultaneous
region need no new category-selection penalty. The final bunching endpoint
is `.023025033926510948` and the C2 ceiling is `.0224`, leaving
`.000625033926510948`. The number is not a claim of familywise validity over
all exploratory archives, crops, dates and analyses ever inspected.

### Exact physical models and forward probabilities

Each new model stores nonnegative integer base amplitudes A_yx and Gaussian
integer row phases u_yj of exactly common modulus R=1105. Its detected entries
are

    T_(y,x),j = A_(y−s_j),x u_yj / (10^7 R),  s = [0,+1,−1].

There is no cyclic wrap; the master 28-row grid includes all translated
support. Translation of probabilities is exact. All three leading principal
minors of the Hermitian matrix `I−T†T` are strictly positive by integer
arithmetic, so appending loss modes completes a physical isometry. Arbitrary
row phases are allowed by the original class; these models are not asserted
to arise from the published fitted lattice Hamiltonian.

For m surviving atoms at an occupation list o, let
`V_sigma=product_(i<m) T_(o_i),sigma_i` and let
`G=I−T_crop†T_crop`, including all unobserved master modes as loss. Then

    P(o) = sum_(sigma,tau) V_sigma conj(V_tau) J_sigma,tau
             product_(i=m..2) G_tau_i,sigma_i
           / ((3−m)! product_s occupation_s!).

For identical bosons J is all ones. For the cluster example it is one only
when sigma and tau differ by exchanging input labels 0 and 1, representing
an identical pair with an orthogonal third particle. This is an explicit
member of C2, not merely arbitrary PSD moments. Tracing losses and applying
parity gives, for example,

    P(empty) = P(no survivors) + sum_i P(i,i),
    P({i}) = P(i) + P(i,i,i) + sum_(k!=i) P(i,k,k).

All probabilities have a common denominator `6 (10^7 R)^6`; summation and
comparisons to interval endpoints are exact. The checker covers every
empty-, one-, and two-click outcome and every observed three-click pattern.
For unobserved triples it uses the analytic cluster/boson multiple of the
exact distinguishable probability, refining with the exact occupation
probability if needed. Their certified probability upper bounds are
`.001630106` (quantum) and `.002625071` (cluster), below the zero-count upper
endpoint `.006507375`. No negative floating-point probability is accepted.
An independent small complex isometry test enumerates explicit loss modes
and block permanents and matches this loss-trace formula for every parity
outcome and total normalization.

| Source set | Concrete feasible model or exclusion |
|---|---|
| Joint | `full-boson.json` passes the full singleton region, every higher pattern interval and B; P(B)=.0259974396797. C2 is empty by the 2D bound. |
| Singleton only | The original `physical-witness.json`, with orthogonal labels, exactly reproduces all empirical cells; `cell_certificate.py` verifies inclusion in the new singleton region. P(B)=.00462813088179. |
| Higher only | `full-cluster.json` retains fixed linear propagation and exact translation but chooses calibration-free amplitudes. It passes all 497,785 pattern constraints and B; P(B)=.0289752778116. |

The higher-only model is not joint-compatible: e.g. its reference row 14
probability is `.27745506669436`, while the actual singleton region permits
only `[.01219991,.05337102]`. Thus the existing singleton measurement
separates this particular focusing explanation under the sharing bridge.
This is not a fitted shared-drift model of both sources. The old quantum
summary model assigned zero probability to 68 observed patterns; the new
quantum model adds small cell mass, stays inside every singleton constraint,
and resolves that support failure. The original exact empirical singleton
model remains valid for source deletion.

Coverage of every full-pattern constraint is a limited and reproducible
statement. The sparse Bonferroni region is weak; it need not capture all
powerful aggregate tests. We do not infer that the higher data are incapable
of rejecting C2 under every possible procedure, or that all released times,
particle numbers, survival curves or Hamiltonian controls are fitted.

## Matched controls

The notebook hash and four extra NC member hashes are pinned in `sources.json`.
`control_audit.py` checks every one of the 17 3NN settings, parses the actual
run maps, and verifies the selected 3NN/1D arrays agree after physical-coordinate
alignment. The selected date is **July 2, 2022**, correcting the earlier broad
date range. The 2.45/2.46 nominal-time discrepancy is duplicate packaging,
not an independently acquired calibration or a measured time uncertainty.

| Records | Actual acquisition and exposure | What is shared or constrained | What stays free |
|---|---|---|---|
| Selected singleton + three-atom images | July 2: singleton `220702/29`, 931 shots; many `220702/21,24,31`, 2,999 shots | Same nominal selected experiment; singleton estimates one reference input distribution. Common T across preparation and exact shifts to the other inputs are assumptions. | No independent selected-time test of the other inputs or shot-to-shot stationarity. |
| Nearest scan point, 2.428571428571 ms | July 7: singles `220707/50,31,52`, 98/94/94 shots; many from `220707/29,36,38,40,42`, 388 shots | Independent preparations at physical sites (y,x)=(24,24),(25,24),(26,24), at this time/date. | Each site's channel at this setting is free; no validated time/date bridge to July 2 at 2.45 ms. |
| 4.65 ms | July 2: singles `220702/39,40,43`, 230/230/233 shots; many `220702/35,37`, 1,202 shots | Same date and independent selected input sites, but different evolution time. | Separate T(4.65); relating it to T(2.45) requires a justified dynamical model. |
| Two-particle/HOM controls | Notebook coarse runs May 28, June 25/29, July 1; fine scans August 25 (including independent per-input runs) | Pair interference and timing information at their own settings. Notebook uses dip optima 1.45 and 1.50 ms to illustrate daily timing variation. | They do not give a uniform three-atom interaction, preparation or detection-error bound at the selected setting. No extra indispensable source is claimed. |

At the nearest scan time the empirical aligned row-TV differences between
input pairs are `.19626574`, `.14046895`, `.10638298`; at 4.65 ms they are
`.13478261`, `.15252846`, `.08565031`. Sampling noise and real transfer
variation are not separated by these point estimates. They are **not** lower
bounds on physical drift. As one scale diagnostic, the nearest scan's
center-displacement-row probabilities have ordinary marginal 95% intervals
`[.01123132,.10121782]`, `[.00258716,.07475177]`,
`[.00258716,.07475177]`. These are exploratory, not simultaneous transfer
bounds and not inputs to the exclusion.

The actual controls demonstrate independent site measurements exist, but
not their transfer to the selected channel with the needed precision. NC
shot order is not a timestamped interleaving schedule. The notebook's daily
HOM-time range does not provide a stochastic law, a within-run drift bound,
or a uniform event-error guarantee. No such bridge is manufactured here.

For full detected-cell distributions (including unobserved/lost outcomes),
a coupling gives `D_B(actual) <= D_B(calibrated) + sum_j delta_j`. The final
sufficient condition is

    2 sum_j delta_j + epsilon < .000625033926510948.

With no other allowance, sum delta must be below `.000312516963255474`
(about `.0001041723` per input if equal). Unlike the old row-event relaxation,
this statement needs cell-TV control; row-TV alone does not control the
collisionless event. The old row-TV tolerance can be compared on this common
full-cell metric because marginalization cannot increase TV. Larger errors
are not proved fatal. Epsilon is a symbolic event-probability bound, not an
arbitrarily assigned interaction/detection allowance. No measured epsilon,
selected-channel TV guarantee, or calibrated common-drift region has been
established. The published synthetic drift adversary remains illustrative;
fitting an apparatus-admissible fluctuating cluster model jointly to both
records and identifying a separating matched control is still open.

## Interpretation, selection history and closest comparison

For the final region, any common-T decomposition
`rho=(1−w)rho_C2+w rho_other` obeys `w >= .0139516501453`. This is specific
to the labelled-block null, not entanglement depth or a universal measure of
interference order. No new fundamental probability law is proposed.

On these same cell constraints, the fully distinguishable ceiling is `.0112`,
the identical-pair-plus-labelled-third ceiling is `.0224`, and the general
three-particle Cauchy upper bound is `.0672`. The broader allowed mixed and
correlated within-block states retain the same factor-two ceiling. Thus the
numerical event ceiling is also the one for the closest perfectly labelled
cluster class; this reanalysis does not demonstrate a stronger numerical
witness than that prior class on these observations. Young et al. already
used singleton calibration and 2D bunching. Brod et al.'s operational genuine-
indistinguishability witnesses are especially close prior art. The concrete
increment here is conservative finite-count crop/loss certification and
explicit region-compatible source-removal models, not invention of the
resource or consequential physics established by numerical exactness.
See [source-audit.md](source-audit.md) for primary literature and distinctions.
A head-to-head implementation of all prior witnesses and apparatus-robust
experimental significance remain incomplete.

Discovery was exploratory: the original pass inspected atomic n=2–5 and two
photonic archives, selected the atomic n=3 original full-row event, and proved
the row relaxation. The review explicitly prompted the conditional-cell
region, the exact distinct-x objective, and neighbouring controls. The new
cell region was fixed before solving its dual bounds. The complete-pattern
region then split the higher-source budget before the full-model search.
The search varied singleton pseudocounts, calibration-free focusing widths
and phase seeds. Witness selection within a fixed simultaneous region is
an existence search and adds no coverage penalty. All final physical models
are checked independently of that search. Earlier cyclic-crop prototypes
were discarded; final translation has padded physical support. No alternate
crop, new bunching definition or extra archive is a second confirmed finding.
The comparisons among confidence designs are retrospective, not an extra
joint confidence claim. Fresh confirmation would need an independently
specified analysis and matched calibration acquisitions.

On the development container, the complete 16-test study suite takes about
3 seconds; the 2D proof generation plus exact check takes about 0.16 seconds.
The complete-model exact check is about 2–3 seconds, without a solver. These
are observed runtimes, not guarantees. Original conditional power diagnostics
are preserved as historical frozen-ceiling benchmarks. No end-to-end power
number is added because a justified selected-channel drift/control model is
still missing; a plug-in resampling number would not resolve that gate.
