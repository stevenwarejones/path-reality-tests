# Sharp boundary for the complete reference table

The model and representation assumptions are exactly those of the
[path-contextuality study](../path-contextuality/derivation.md): a shared
preparation μ, shared stochastic final response r, normalized instrument kernel,
negative-response cap q, and total channel (1−d) identity + dD. For a probability
interpretation restrict q,d to [0,1]. This is the reduced observable interface,
not reproduction of every tomographic calibration context.

Write a=P(−,success), b=P(+,success), f=P(success|bypass), and g=P(−).
For every finite ontic cardinality:

1. g≤q by averaging the pointwise negative cap.
2. a+b≤f+d(1−f), since the identity part retains r and the residual channel's
   final success is at most one.
3. a≤qf+d(1−f), the existing finite-pointer bound.
4. b≥(1−d−q)f. At state λ the diagonal total mass is at least 1−d. The negative
   diagonal mass is at most the whole negative row, hence at most q. Therefore
   the positive diagonal mass is at least 1−d−q. Multiply by r(λ)≥0, add the
   nonnegative off-diagonal terms, and average with μ. The inequality remains
   valid when its right side is negative.

The reference data are f=49/625, g=337/625, a=1369/15625, b=144/15625.
Consequently feasibility requires g≤q and

    d ≥ max(1/50, (a−qf)/(1−f), 1−q−b/f).

## Attainment: no unproved cardinality truncation

Let u=min(q,1081/1225), x=(a−fu)/(1−f), and s=1−u−144/1225.
Prepare deterministic-success S with mass f and deterministic-failure F with
mass 1−f. Use the following complete kernel:

| Starting state | −→S | −→F | +→S | +→F |
|---|---:|---:|---:|---:|
| S | u | 0 | 144/1225 | s |
| F | x | 49/100 | 0 | 51/100−x |

For g≤q≤1, every entry is nonnegative and each row sums to one. The negative
response at F is x+49/100=(g−fu)/(1−f)≤u because g≤u. The complete joint table
is the reference table, including both failure cells; bypass success is f.

The S→F rate is s and the F→S rate is x. Thus any d≥max(s,x) supplies a
normalized residual kernel with flip probabilities s/d and x/d, respectively.
The reference region guarantees d≥1/50>0. If q≤1081/1225, these are the two
nonconstant lower bounds. If q is larger, the fixed choice u=1081/1225 gives
s=0 and x=1/50. This establishes sufficiency for every point of the region.
Universal necessity plus this construction proves sharpness across **all finite
cardinalities for this particular table**; it does not claim a general two-state
reduction or cover arbitrary measurable ontologies.

The two sloped branches cross at q=22223/25823. The horizontal branch starts at
q=1081/1225. No model exists below q=337/625 even at d=1.

At q=16/25 the exact minimum is 297/1225≈0.242449. Inverting the negative-success
inequality alone gives only 13/320≈0.040625, which is not sufficient for this
complete table. At d=1/50 the exact minimum cap is 1081/1225≈0.882449. These
numbers quantify this ontic representation class, not measured apparatus disturbance.

`boundary.py` verifies every entry with Fractions. An independent eight-variable
transition LP in the tests checks several boundary points; its floating-point
output is corroboration, not the proof certificate. The analytic construction
is checked at exact breakpoints and a dense rational grid. The Lean companion
is [PR #87](https://github.com/stevenwarejones/ontology-separation/pull/87).

## Conditional probability and shared-procedure errors

Let a model satisfy the cap and disturbance premises. If
|b_obs−b_model|≤ε_b and |f_obs−f_model|≤ε_f, then

    b_obs ≥ (1−d−q) f_obs − ε_b − |1−d−q| ε_f.

A separately assumed total-variation distance τ between probe and bypass
preparations contributes at most τ to ε_f for a response in [0,1]. A separately
assumed supremum readout discrepancy ρ between those procedures contributes at
most ρ. Errors affecting the **probe** positive-success probability must also be
included in ε_b. These are bounds on specified model/observable discrepancies,
not consequences of a small operational-equivalence residual.

`positive_floor` first expands the bypass interval by the supplied τ and ρ,
intersects with [0,1], and minimizes the linear expression over its two endpoints.
It subtracts ε_b and uses positivity. This is a sound conditional bound, not a
claim that every error budget is jointly attainable or statistically certified.
No Wen SD column is converted into any of these budgets.

## Literature scope

The parent theorem is Kunjwal, Lostaglio and Pusey (2019),
[Theorem 3 and Appendix B](https://arxiv.org/html/1812.06940v2).
The present increment uses more cells of the same reduced interface and gives
an exact reference-data compatibility characterization. No novelty claim is made;
this is not a new assumption-free contextuality test or a trajectory exclusion.
