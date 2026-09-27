# Two-sided selection and a certified limitation of the transfer method

The requested two-source payoff for a fixed-input sequence remains unestablished. A subsequent [operational discrimination-resource result](contrast.md) establishes both source-deletion crossings for a different quantity, largely detector contrast. This follow-up changes the selection objective, checks a stronger shared-dynamics relaxation, and proves exactly why adding more instances of the current transfer inequality cannot answer the GST-only question. It does **not** turn a relaxation point into a physical source-deletion witness.

| Question | Result and status |
| --- | --- |
| Does the revised screen find a crossing on both complete-source sides? | No, in the declared existing physical point pool; discovery evidence only |
| Does complete GST already imply the original target lower bound over physical qubit models? | Still unresolved |
| Can the entire transfer-inequality family plus **every** GST cell/group prove a positive target lower bound? | **No. Certified explicit scalar assignments attain 0 and 1, and interpolation gives [0,1]** |
| Do those assignments establish physical GST-only feasibility? | **No. They are not quantum gate sets** |
| Does a richer Bloch-Gram outer relaxation improve the certificate? | Numerical probes do not; no new dual certificate claimed |
| Do combined data fail to improve any unmeasured prediction? | Not proved; the existing one-sided result remains valid |

## Two-sided selection objective

For a candidate lower bound `L`, the score is

\[
\min\{L-\min_{M\in B_{GST}}p_M(w),\ L-\min_{M\in B_{RB}}p_M(w)\}.
\]

For an upper bound `U`, replace each term with the corresponding largest point prediction minus `U`. The point sets are the original complete-source witnesses and the two previously verified all-GST search points. They are **discovery aids**, not approximations certified to cover the physical prediction regions. Every point's retained count constraints are independently checked by the existing verifier.

The unchanged 9,975-candidate family was rescored in both available bound directions. Zero directed bounds had positive margins on both sides. The best score was approximately **−0.65967**, for a 34-gate RB word followed by `GyGy`: its lower bound is approximately 0.11136, its lowest displayed GST-only prediction approximately 0.77103, and its RB-only prediction approximately 0.00410. The previous four-gate target scores approximately −0.66005. Ranking only the RB gap hid this large bottleneck.

This does not prove that another physical GST point or another target cannot work. It does show that changing the ranking alone within this screen is insufficient. [Machine-readable screen](results/two-sided-screen.json) and `two_sided_screen.py` preserve the objective, point pool, count and selected hashes.

## The complete transfer relaxation

Let `W` be the finite set of all 4,657 measured GST words, and let `T=GiGiGyGy`. For every finite word `w` introduce a scalar `p(w)` in [0,1]. Retain all original GST cell intervals and all six original GST group constraints. In addition impose, for **every** pair of prefixes `u,v` and **every** suffix `f`,

\[
|p(uf)-p(vf)|\leq h(p(u),p(v)),\qquad
h(a,b)=\sqrt{a(1-b)}+\sqrt{b(1-a)}.
\]

The qubit composition lemma in [the derivation](derivation.md) applies to the two states `G_u(rho), G_v(rho)` and the common suffix. Therefore every stationary physical qubit model maps into this relaxation. The converse is false and is not used.

The following construction satisfies this **infinite** inequality family, not merely the instances used in PR21 or a finite numerical sample.

Set `b=21/100`. For each measured GST word, assign

\[
p(w)=b+(1-2b)k_w/n_w.
\]

Exact rational arithmetic verifies all **4,663 two-sided GST constraints**, with minimum inequality slack **greater than 0.0028389004**, using the unchanged confidence endpoints. No new uncertainty allowance or confidence reallocation is introduced: these are exact scalar values, not numerically propagated physical predictions.

Define two extensions to the entire word alphabet:

| Assignment | `p(T)` | Every other unrecorded word |
| --- | --- | --- |
| Lower endpoint | 0 | 2/5 |
| Upper endpoint | 1 | 3/5 |

The original source geometry supplies one crucial fact: **the only recorded GST word having `T` as a prefix is `TY=GiGiGyGyGy`**, at global row 121. It has 26/50 counts, so its assigned probability is `d=1279/2500=0.5116`. The proper prefix `IIY` is row 119, with 22/50 counts and assigned probability `q=1163/2500=0.4652`. The target itself is unrecorded. The dedicated full-source route reconstructs all words and verifies this exhaustive proper-extension statement; the offline check uses its pinned row/hash metadata and exact arithmetic.

### Proof for all words

Every non-target probability is in `[b,1-b]`. On that square,

\[
h(a,c)\geq2\sqrt{b(1-b)}\geq1-b,
\]

where the last inequality is equivalent to `b>=1/5`. Thus when neither input prefix is `T`, the transfer inequality holds: the greatest possible output difference is `1-b`, even if one output word equals `T`.

If both input prefixes are `T`, the output difference is zero. It remains to consider exactly one input prefix equal to `T`.

For an empty suffix, the inequality reduces to `p<=sqrt(p)` at the lower endpoint, or `1-p<=sqrt(1-p)` at the upper endpoint.

For a nonempty suffix other than `Y`, the output from `T` is respectively 2/5 or 3/5. Its difference from any possible other output is at most 2/5, including when that other output is `T`. Meanwhile the right side is respectively `sqrt(p(v))` or `sqrt(1-p(v))`, at least `sqrt(b)`. The exact inequality `(2/5)^2<=b` proves this case.

For suffix `Y`, the output from `T` is the exceptional measured value `d`. If the other output is not `T`, its difference from `d` is at most `max(|d-b|,|d-(1-b)|)`, whose square is less than `b`. If the other output **is** `T`, its input must be `IIY` by cancellation in the free word monoid. The remaining inequalities are exactly

\[
d^2\leq q,\qquad (1-d)^2\leq1-q.
\]

Both hold by rational arithmetic. This exhausts all prefix/suffix cases and proves both endpoint assignments feasible.

Finally, `h` is jointly concave: each summand is a geometric mean of two nonnegative affine functions. Absolute output differences are convex. Consequently the set satisfying all transfer inequalities and the linear count/group constraints is convex. Interpolating the endpoint assignments attains **every** target value in [0,1].

This proves a precise methodological obstruction: **no valid inference from only these transfer inequalities, probability normalization, and the complete GST confidence constraints can certify any positive GST-only lower bound for this target.** It does not prove that GST fails to imply such a bound over shared CPTP gate sets. In particular, neither endpoint assignment is offered as the physical GST-only witness required for combination gain.

## A richer shared-dynamics outer relaxation

We also implemented a finite-prefix Bloch-Gram SDP. Write a qubit state as `(I+r_w.sigma)/2`, a gate as `r -> T_g r+t_g`, and the common effect as `a I+e.sigma`. Put all prefix vectors, `e`, and the three translations into one Gram matrix `Q`.

Every physical model satisfies:

* `Q` positive semidefinite; `||r_w||<=1` and `||t_g||<=1`;
* `p(w)=a+e.r_w`, with `0<=a<=1` and `||e||^2<=min(1/4,a/2,(1-a)/2)`;
* for each primitive, the Gram matrix of its input vectors minus the Gram matrix of the corresponding centered outputs `r_wg-t_g` is positive semidefinite;
* `||2t_g-r_wg||<=1`, since the antipode `-r_w` is also a physical qubit state.

The matrix contraction constraint follows from trace-distance contractivity: a qubit channel's real linear Bloch part satisfies `T_g^T T_g<=I`, and the centered output is `T_g r_w`. These conditions preserve a shared affine contraction across many prefixes, rather than assigning an independent map to each circuit. Dropping rank three, full effect positivity, native Choi positivity and unused data constraints makes this an **outer** relaxation. It does not narrow the declared physical class.

| Source constraints used in probe | Rows | Prefix states | Numerical minimum for target |
| --- | --- | --- | --- |
| GST, lengths at most 4 | 73 | 78 | Approximately 0 |
| GST, lengths at most 5 | 139 | 148 | Approximately 0 |
| Joint, lengths at most 4 | 82 | 78 | Approximately 0.08455 |

These are **uncertified solver outputs**, not global physical minima, not verified dual bounds, and not certified feasible physical points. They suggest this outer relaxation does not improve the existing analytic 0.11136 joint certificate. Small negative values near machine tolerance are not negative physical probabilities. No exclusion follows from solver status.

An explicit strictness check matters here: Bloch transposition `diag(1,-1,1)` is a ball isometry, so it satisfies the contraction conditions, yet its Choi matrix is the swap and has eigenvalue −1. Thus contraction alone does not enforce complete positivity. The exact all-word transfer obstruction and the numerical Gram probes concern different relaxations; the former is not presented as a proof about the latter.

Initial small probes used Clarabel; the 148-state Clarabel process exceeded available memory. SCS completed all three declared configurations. Final artifacts were regenerated in an isolated environment with CVXPY 1.7.3 and SCS 3.3.1; solver versions are recorded. Neither memory exhaustion nor an optimizer result is a physics exclusion.

## What is retired and what remains open

Retire **transfer inequalities alone** as a route to the GST-only assessment of this target. Do not retire the target on the claim that GST already implies the current joint bound: that has not been proved. Also do not launch another tiny resource-boundary crossing; the existing operational obstruction is preserved unchanged.

The unresolved mathematical step is a certificate that retains enough of the actual shared **CPTP** dynamics. One concrete formulation uses fixed 2x2 state/effect matrices, three Choi-positive trace-preserving native channels, and the repeated-prefix equations `r_wg=T_g r_w+t_g`. A verified polynomial/moment relaxation would need those bilinear composition equations and the qubit dimension restriction, with independently checked dual residual bounds. Long repeated-germ records are promising constraints on persistent channel modes, but no new inference from them has been certified here. The current Gram formulation deliberately omits information essential to that proposed next step; merely running its solver harder does not supply it.

A future candidate is successful only after a full-class joint bound and **both physical complete-source crossing witnesses** are verified under one declared confidence construction. No finite-bank spread, scalar relaxation assignment, or numerical SDP minimum meets that condition.

## Reproduction

```sh
python studies/unmeasured-dynamics-gain/transfer_barrier.py --check
python studies/unmeasured-dynamics-gain/transfer_barrier.py --cache /tmp/unmeasured-source --check
python studies/unmeasured-dynamics-gain/two_sided_screen.py --cache /tmp/unmeasured-source --check
python -m venv /tmp/dynamics-gram-venv
/tmp/dynamics-gram-venv/bin/python -m pip install -r studies/unmeasured-dynamics-gain/discovery-requirements.txt
OPENBLAS_NUM_THREADS=1 /tmp/dynamics-gram-venv/bin/python studies/unmeasured-dynamics-gain/gram_probe.py --cache /tmp/unmeasured-source --check
```

Populate the external cache with the original README's hash-checked downloader first. Exact barrier verification needs no SDP solver. The Gram reproduction checks its numerical snapshot to a declared tolerance, not a dual certificate. [Protocol](followup-protocol.json), [barrier certificate](results/transfer-barrier.json), [two-sided screen](results/two-sided-screen.json), and [Gram probes](results/gram-probes.json) retain those distinctions. The transfer lemma, Gram construction and convexity arguments use established geometry and convex optimization; no novelty claim is added.
