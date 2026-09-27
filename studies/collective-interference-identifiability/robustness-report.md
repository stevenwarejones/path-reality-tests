# What the archive forces when transport may fluctuate

The fixed-channel source-combination result remains intact. The next step is
now quantitative: a shared fluctuating cluster model fits the selected
singleton region, every selected parity-pattern constraint, and new measured
row-hit constraints. It uses the **same two-point acquisition law** in both
preparation families. Its fluctuation is large, not a demonstrated tiny-error
explanation of the experiment. All 16 independently calibrated 3NN settings
also admit fixed cluster models in their own declared simultaneous regions.

This establishes a limitation of inference in the broad passive-channel
class. It does **not** establish a plausible fluctuating nearest-neighbor
lattice Hamiltonian, an error actually present in the experiment, or the
minimum required fluctuation. Apparatus-supported robustness and consequential
novelty remain open.

## Decision table

| Route | Physical nuisance class | Measured controls | Certified exclusion over the class | Explicit physical feasibility | Source complementarity | Scientific change |
|---|---|---|---|---|---|---|
| Preserved baseline | One fixed hidden-independent noninteracting T, exact input translation | Selected 931 singleton and 2,999 three-atom shots | C2 ceiling .0224 below .0230250339265 | Identical-boson joint model; original source-deletion models on the original region | Original region-based criterion remains closed | Baseline retained, not re-solved |
| Preparation-dependent transport | Separate fixed singleton and many-particle channels; sum of full-output column-TV changes <= d | Same selected records | Excluded through d=.00554, using the actual event's multilinear structure | No endpoint claimed for the new row-hit region; lower-budget candidates failed complete-pattern constraints | Uses both records; no new deletion claim | Sufficient allowance is 17.7 times the previous generic .000312517, without tightening the ideal .0112 bound |
| Shared shot fluctuations | T_Z with a common IID latent law for both preparations; mean translated calibration; ess-sup sum of column-TV deviations from the mean <= d | Same selected records plus all 12 row-hit events | Excluded through d=.0397 | Equiprobable physical pair-plus-singleton channels fit at d=.89852380063206 | Joint compatibility returns in the enlarged class; no resource exclusion at this endpoint | A fitted empirical-region countermodel replaces the toy drift adversary; critical d is bracketed, not solved |
| Independent-site settings | One fixed T per setting, three independently measured columns; no translation assumption; no cross-time bridge | All 15 July 7 scan points and July 2 at 4.65 ms | Exact-event outer bounds are inconclusive | Physical fixed cluster models pass every singleton and higher-pattern constraint at all 16 settings | No binary resource exclusion or indispensable extra control family is claimed | Calculates the lower-statistics alternative instead of assuming it fails |
| Output geometry | Same selected region, expanded with the event “row y has at least one parity click” | All selected images; 12 simultaneous row-hit intervals | Rejects the old focusing example on eight rows; not a full-class exclusion by itself | Both the fixed quantum model and fitted shared-drift model pass all 12 | Additional information from the existing many-body source, not an additional source | The geometry eliminates one actual fitted nuisance explanation, but not the drifting one |
| Dynamics / other imperfections | A shared uncertain lattice Hamiltonian, contamination, detection errors and interactions would be distinct extensions | Notebook/HOM diagnostics and published apparatus analysis | No uniform selected-channel bound established | No apparatus-specific Hamiltonian countermodel claimed | Cross-date and cross-time transfer remain unvalidated | Precise remaining physical gate; no unexplained epsilon is assigned |

## Nuisance coordinates and the bracket

For a normalized detected-output distribution include every detector cell and
an additional loss/outside symbol. A finer loss decomposition can be used;
our bounds remain conservative after adding unobserved modes. Let
`p_j(z)` be column j's probability distribution under channel T_z,
`pbar_j=E p_j(Z)`, and `d_j(z)=TV(p_j(z),pbar_j)`.
The main fluctuation budget is

    d = ess sup_z sum_(j=0..2) d_j(z).

It is neither a time-jitter percentage nor an RMS fluctuation. The common law
is IID between shots and shared **after the archive's initial preparation
selection**. On each shot all present atoms experience the same latent
channel. Calibration averages probabilities, never amplitudes. The proof
only requires mean translation; the two exhibited latent channels satisfy
translation individually as well. Hidden labels form the same identical
pair {0,1} and orthogonal singleton {2} in both components. Loss modes complete
each exact contraction to an isometry. The construction allows no extra
atoms, false clicks, hidden-dependent response, or interactions.

Keep the following coordinates separate:

| Coordinate | Treatment here |
|---|---|
| Mean input-site translation mismatch | Zero in the selected-data drift construction; eliminated as an assumption at the independent-site settings. |
| Change with singleton versus three-atom preparation | Zero in the drift construction: identical law and channel family. A separate fixed-preparation TV envelope is provided. |
| Shared shot fluctuations | The sole departure in the selected fitted countermodel; explicitly bounded and constructed. |
| Preparation contamination | Zero; exactly three selected prepared particles is retained. |
| False positives / hidden-dependent detection | Zero; the actual parity map and label-independent loss are retained. |
| Interactions | Zero in the certified classes. Published simulation estimates are not inserted as a uniform event error. |

For the augmented selected region the exact checker establishes

    0.0397 < d_critical <= 0.89852380063206.

The lower endpoint is an outer proof over every channel law in the stated
class; the upper endpoint is a physical inner construction satisfying every
retained constraint. Neither is an optimizer's infeasibility report. The gap
is substantial. Each input in the symmetric construction fluctuates by about
.299508 in TV, so the feasible endpoint is **not small**. No conclusion that
.04 fluctuation suffices, or that .8985 is necessary, is supported.

### Event-specific exclusion envelopes

Write `D(p0,p1,p2)` for the exact probability that three independent draws
occupy three distinct cells in one selected row. It is a nonnegative
multilinear functional. The preserved 2D certificate gives `D(pbar)<=U=.0112`.
The calibration row intervals imply

    b = max_(y,j!=k) pbar_j(row y) pbar_k(row y) <= .055853496221322,
    g = max_(y,j) pbar_j(row y) <= .2881834.

For a fixed preparation change, write each signed difference as positive
part minus negative part. Dropping the negative parts can only increase D.
Their masses are the corresponding TV distances. One changed draw has a
coefficient at most b; two have a coefficient at most g. Consequently,

    D_many <= U + b d + g d^2/3 + d^3/27.

Here `sum d_j<=d`, `sum_(j<k) d_j d_k<=d^2/3`, and
`d0 d1 d2<=d^3/27`. Comparing twice this expression with the measured lower
endpoint excludes the entire fixed-change class through d=.00554. This uses
the exact collisionless event; the row probabilities bound only its
perturbation coefficients. It is not the generic `D_many<=U+d` coupling.

For common fluctuations the first-order terms vanish because
`E[p_j(Z)−pbar_j]=0`. For a signed pair product only the ++ and −− parts can
increase the event, each bounded by `g d_j d_k`. For the signed triple product
there are four positive-sign combinations. Thus, uniformly over all latent
laws and correlations within the common shot,

    E D(p0(Z),p1(Z),p2(Z)) <= U + 2 g d^2/3 + 4 d^3/27,
    P_C2(B) <= 2 E D.

The factors two and four matter: treating signed perturbations as only their
positive parts after canceling the linear terms would understate the bound.
The checker compares these rational expressions to the exact bunching lower
endpoint; .0397 lies strictly in the exclusion range. Higher-pattern and
row-hit constraints can only shrink the allowed set, so omitting them in
this outer proof is conservative.

## A fitted shared law, not separate family fits

`drift-minus.json` and `drift-plus.json` specify integer amplitudes and rational
unit phases. Every shot independently chooses either channel with probability
one half. Both singleton and many-particle acquisitions use that same law;
there is no preparation-dependent choice of weights or channels. Their
probability-averaged singleton distributions satisfy every original row,
conditional-cell and empty/outside interval. Each component's `I−T†T` has
strictly positive exact principal minors.

The pair-plus-singleton model's selected bunching probability is
`.023075036505304354`. Its mixture probabilities satisfy every observed and
unobserved pattern interval, including empty, one- and two-click outcomes.
All 12 row-hit constraints also pass. The mean probabilities, mixing weights,
component physicality and every comparison are checked using fractions and
integers, with no floating-point acceptance test.

The model's substantial left/right redistribution is permitted by the broad
channel class. Nothing establishes that the apparatus could switch between
these channels under its measured tunneling Hamiltonian or reported timing
variation. It is an **empirical-region ambiguity under an explicit physical
extension**, not a diagnosis of the apparatus. Compatibility with sparse
simultaneous intervals remains weaker than likelihood adequacy or agreement
with every possible aggregate statistic.

## Geometry and confidence coverage

Define H_y as at least one detected parity click in row y. Zero parity is not
zero atom occupation. To predict H_y, restrict T to that row and trace all
other modes as loss; then `P(H_y)=1−P(empty parity in that row)`. The empty
outcome includes zero surviving atoms and same-site pairs, exactly as in the
full loss/parity calculation. It is not approximated by one minus the
probability of zero physical particles in the row.

The selected-data design retains calibration failure 1/160. Higher-source
failure is also 1/160, divided into 1/320 for the bunching lower endpoint,
1/640 for all parity patterns, and 1/640 for all 12 row-hit intervals.
Every parity pattern receives two tails of `1/(1280 K)`, K=497785. Each
row-hit interval receives two tails of 1/15360. The row endpoints use exact
integer checks of the binomial Chernoff bound, rounded outward to a rational
grid. They are conservative simultaneous intervals, not plug-in error bars.
The existing n=2–5 family reserve remains .05. The fixed-channel numerical
exclusion is unchanged; no ideal-bound shaving was needed.

The old higher-only focusing model fails rows 10, 11, 13, 14, 15, 16, 17 and
18. The existing quantum model passes all rows. The new shared-drift model
passes them too. The original source-combination result remains valid on its
original region. We do not claim that fixed-channel higher-only insufficiency
has been re-established on this **augmented** region, or that one rejected
focusing example excludes every focusing mechanism.

The independent-setting analysis separately reserves .05 over all 16
settings. At each setting, calibration gets 1/640, divided equally between
all three inputs' full cell/loss intervals and row/loss intervals; the
higher-order record gets 1/640, split equally between bunching and every
parity pattern. Three-input dependence is handled by union bounds, not by
assuming separate acquisitions are independent of each other. IID sampling
within records is retained. This setting search is covered uniformly.
The selected-data and independent-setting analyses are **separate declared
confidence families**, not a single 95% guarantee over their union. No
archive-wide 95% exclusion is claimed.

## Independent calibrations at their own settings

All 16 independent 3NN settings are evaluated, including the strongest
empirical candidate at 1.571429 ms and the explicitly requested 4.65 ms
setting. Each has its own physical T with three independently constrained
columns and the same T for singleton and many-particle preparations. The
checker uses the complete crop, each recorded loss/empty outcome, and every
possible higher parity pattern. No translation relation is imposed.

The empirical independent-draw probability is calculated with the exact
distinct-cell formula. Valid full-class upper bounds take the minimum of its
entrywise-upper evaluation, the row-product upper bound, and one. These
bounds are deliberately not advertised as optimal. Physical joint C2 models,
not merely loose bounds, establish compatibility at every setting. At
1.571429 ms the point estimate gives `2D≈.022238`, below the bunching lower
endpoint `≈.0296054`, yet a calibration-compatible physical cluster model
gives `P(B)≈.0326121`. This is a concrete illustration of why the empirical
calibration frequencies cannot be treated as known probabilities.

The exact results for all settings, including bounds and witness predictions,
are in `results/robustness-check.json`. The matched NC counts are reconstructed
from the already pinned source members. Their smaller samples do not support
a resource exclusion in these regions. This does not prove that every
possible inference method lacks power.

For a joint physical extension covering these acquisitions, use one common
fair-bit latent law Z. At the selected setting use the two committed channels;
at each other setting use its committed fixed channel for both values of Z.
This uses the same law across preparation families, and fits all included
independent-site controls. Channel response is otherwise free by setting and
date. It is not a dynamical bridge or a fitted common Hamiltonian. Two-particle
HOM acquisitions at other dates, n=4/5 records and other archive experiments
are not silently added as fitted constraints.

## A concrete separating control

At the selected July 2 nominal 2.45 ms setting, prepare **two orthogonally
labelled atoms** in inputs 0 and 1, with the same channel-acquisition law.
Count the event that both survive at distinct sites within reference array
rows 9–13 and x indices 8–19. Parity detection resolves this event with two
clicks; same-site pairs are excluded. For a channel component its probability
is exactly

    p0(S) p1(S) − sum_(s in S) p0(s) p1(s).

The two exhibited explanations predict:

| Explanation | Control-event probability | Expected events in 1,000 prepared trials |
|---|---:|---:|
| Existing fixed quantum channel, with orthogonal labels in this control | .1432359739972 | 143.2 |
| Fitted shared-drift cluster law, with orthogonal labels in this control | .2000702698266 | 200.1 |

The difference is about 5.68 percentage points. The atoms must be labelled
within the same shot, with independently verified distinguishability and
label-independent transport. Time labelling by combining separate singleton
shots would erase precisely the shared-channel covariance being tested and
is not an equivalent control. A concrete implementation of suitable labels
in this apparatus still needs validation. This proposes an experimentally
specified measurement, not an acquisition already performed or authorized.
It distinguishes these two explicit models; it is not yet a uniform
resource witness against every drifting cluster model. A future analysis
would need confidence bounds on these matched pair moments and source
deletion for the added control. Shot order or non-matched HOM measurements
cannot substitute for this acquisition.

## Apparatus and prior-art gate

[Young et al., supplementary sections VI–VIII](https://arxiv.org/html/2307.06936v2)
model spatial calibration and interactions; their interaction simulations
report errors at or below 10^-4 for the examined few-atom signals. That is not
a uniform bound over the present adversarial channel class. Their quoted
HOM error budget concerns a postselected measurement and cannot be copied as
an error allowance for our per-prepared-shot event.

As a deliberately broad analytic comparison, for three particles an on-site
interaction has spectral range at most `3 |U_c|`. Subtracting a scalar energy
and applying Duhamel's bound gives an event-probability deviation no larger
than `min(1, 3 |U_c| t/(2 hbar))`, assuming that bounded Hamiltonian model.
At the paper's nominal `|U_c|/hbar=2 pi ×1.7 Hz` and nominal t=2.45 ms this is
about .0393, far above the statistical margin. The nominal input is not a
certified interaction-confidence bound. This calculation illustrates the
cost of a uniform bound; it does not assert the actual interaction error is
.0393 or refute the published simulations.

[Brod et al.](https://doi.org/10.1103/PhysRevLett.122.063602) already provide
operational genuine-indistinguishability witnesses. Our fixed-data factor-two
ceiling remains the same closest labelled-cluster ceiling, not a new witness
invented here. [Pont et al.](https://doi.org/10.1103/PhysRevX.12.031033) measure
collective indistinguishability through a cyclic-interferometer fringe; that
protocol's measured fringe is absent from these atomic records, so claiming
an equivalent-data numerical superiority would be unjustified. The targeted
literature search does not establish priority for this robustness result.

The concrete increment is an event-specific drift envelope, a physical
common-law model fitting the declared measured regions, and a specified
matched control separating it from the quantum example. Establishing a
credible bounded Hamiltonian family, optimizing the wide robustness bracket,
and obtaining an apparatus-supported exclusion remain the next gate.

## Reproduction and discovery record

The old exact artifacts remain unchanged. Run:

```sh
python studies/collective-interference-identifiability/robustness.py --check
python studies/collective-interference-identifiability/build_robustness.py --cache /tmp/collective-sources --check-source
```

`build_robustness.py --cache ... --models` regenerates proposed intervals and
witnesses; solver termination never proves physicality or exclusion. The
standard-library checker supplies that decision. Matched model search
parameters are deposited separately. These are fitted existence witnesses,
not held-out predictions for the measured records.

The search first tried fixed preparation changes with row and cell
optimization. Several lower-TV proposals satisfied bunching but failed
one-click or rare-pattern constraints, so none is reported as a feasible
endpoint. The row-hit family was introduced after inspecting the old focusing
model; all 12 rows then received simultaneous coverage, with the selected
higher-source budget reallocated before accepting the replacement. The
shared two-channel search used tighter numerical CP row intervals than the
final conservative exact Chernoff region. Both components were subsequently
rounded and checked independently. Early conditional profiles and failed
phase searches are not proofs of impossibility. The independent search
included all 16 settings and used one uniform allocation; the highest apparent
point-estimate signal is not a postselected confirmatory discovery.
