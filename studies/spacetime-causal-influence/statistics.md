# Inference for a complete local binary record

Let X be the intervention and B the predeclared binary receiver record. The
stationary target is d=P(B=1|do X=1)−P(B=1|do X=0); binary TV is |d|.
All initiated trials contribute, including no-click and invalid-detection trials.
These are established bounded-variable and exact-binomial methods applied to
this protocol. Probability coverage below is a mathematical derivation; the
Lean companion checks algebra, interval transport and the finite union bound,
not the exponential or binomial tail theorem.

## Scope of the earlier-record variant

The derivation below is for the selected optical protocol. For a receiver bit
fixed before the current choice, fresh assignment conditional on that bit gives
a valid independence-null test: the conditional score mean is zero. A detected
association rejects the conjunction of fresh assignment, secure records and
the no-backward-influence null. It does not by itself identify a backward
`do`-effect. Under a proposed backward-causal alternative the earlier bit could
be a descendant of the later choice; conditioning on it cannot simultaneously
serve as an unexamined randomization premise. Such an extension must specify
its causal or potential-outcome model, exogenous assignment and filtration.
The optical memory interval and causal power claims are not automatically
backward-effect intervals. Under an explicit IID sampling model, binomial
intervals can still describe the observed conditional association.

## Fixed horizon with device memory

Fix N before acquisition. Let Hᵢ contain all relevant pre-choice history,
including device memory. Require a fresh fair setting conditional on Hᵢ and
the causal randomization/consistency premise identifying the conditional
responses pₓ,ᵢ with interventions. Do not infer this premise from balanced totals.
Define Zᵢ=(2Xᵢ−1)(2Bᵢ−1), equal to +1 for a match and −1 otherwise. Directly,

\[
E[Z_i\mid H_i]=\tfrac12(2p_{1,i}-1)-\tfrac12(2p_{0,i}-1)=d_i.
\]

The relevant filtration contains each completed trial before the next setting.
Receiver probabilities may depend arbitrarily on the earlier history. Past
settings may reach the receiver between trials; fresh present settings must
remain independent of that history. This is why randomly permuted balanced
blocks, whose last settings are predictable, do not satisfy this particular test.

Here is the concentration argument with its constant made explicit. Under any
conditional law of Z∈{−1,1}, write f(t)=log E[exp(t(Z−EZ))]. Exponential tilting
gives f″(t)=Var_t(Z)≤1, and f(0)=f′(0)=0, so f(t)≤t²/2 for either sign of t.
Applying this conditional inequality successively by the tower property yields
E exp(t∑(Zᵢ−dᵢ))≤exp(Nt²/2). Markov's inequality and t=r give

\[
P\left(\left|\bar Z-\bar d\right|\ge r\right)
\le2e^{-Nr^2/2},\quad
r_\alpha=\sqrt{2\log(2/\alpha)/N},\qquad
\bar d=N^{-1}\sum_i d_i.
\]

This is the conditional bounded-range exponential argument associated with
[Hoeffding](https://doi.org/10.1080/01621459.1963.10500830) and
[Azuma](https://doi.org/10.2748/tmj/1178243286). It requires no IID receiver
outcomes. Clip [Z̄−rα,Z̄+rα] to [−1,1]. The interval covers the *realized
average conditional signed effect*. It does not bound N⁻¹∑|dᵢ|: alternating
dᵢ=+1,−1 cancel. With a constant dᵢ=d it is an interval for the original binary
effect. Under the causal null each dᵢ=0, so cancellation is not a null loophole;
it is a limitation of what a null result can exclude.

For a calibrated predictable bound |πᵢ−1/2|≤ζ on πᵢ=P(Xᵢ=1|Hᵢ),

\[
E[Z_i|H_i]=d_i+(2\pi_i-1)(p_{1,i}+p_{0,i}-1).
\]

The extra bias has magnitude at most 2ζ. Add 2ζ to the interval radius, as
`score_interval` does. An observed frequency imbalance is not a bound on this
conditional predictability. A reused random seed known to the receiver can
make ζ=1/2 even with exactly balanced settings. No time-uniform confidence
sequence or optional stopping guarantee is claimed. Changing N after looking
at these results invalidates the stated coverage.

## Independent-binomial comparison

The optional `iid_difference` assumes, conditional on the assignment sequence,
independent Bernoulli outcomes with fixed p₀,p₁. The assignments are independent
of the relevant preparation and potential responses. Given any assignment
sequence, K₀,K₁ have product Bin(n₀,p₀)×Bin(n₁,p₁) law. Every sequence with
the same counts has that same law; averaging over those sequences establishes
the binomial law conditional on the random counts. This also establishes
unconditional coverage. The memory test above does not require this stronger law.

Construct each Clopper–Pearson interval with error α/2, hence α/4 in each tail:
Lₓ=Beta⁻¹(α/4;Kₓ,nₓ−Kₓ+1) and
Uₓ=Beta⁻¹(1−α/4;Kₓ+1,nₓ−Kₓ). Use Lₓ=0 for Kₓ=0, Uₓ=1 for Kₓ=nₓ;
an empty setting gets [0,1]. Inverting monotone binomial tails gives individual
coverage ≥1−α/2. A union bound gives simultaneous coverage ≥1−α and therefore

\[
d\in[L_1-U_0,\ U_1-L_0].
\]

The [SciPy exact-interval documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats._result_classes.BinomTestResult.proportion_ci.html)
identifies this construction. The implementation uses pinned SciPy beta quantiles;
independent tests evaluate binomial tail sums and small-count coverage, including
unequal counts, empty settings and endpoints. The existing phase-study code was
inspected: its conditional-IID assumptions fit this *optional* comparison but
not the memory analysis. No phase-specific threshold or all-successes bound is
reused as a two-distribution theorem.

## Calibration, rejection and an upper bound after a null result

For either signed interval [L,U], simultaneous coverage implies

\[
\max(0,L,-U)\le|d|\le\min(1,\max(-L,U)).
\]

If nuisance permits |d|≤B under the causal null, reject only if the lower bound
is strictly greater than B; equality is inconclusive. The [physical bridge](physical-bridge.md)
derives B from couplings. For the memory analysis these are conditional on Hᵢ:
require |dᵢ|≤e₀,ᵢ+e₁,ᵢ and a simultaneous bound
N⁻¹∑ᵢ(e₀,ᵢ+e₁,ᵢ)≤B on the calibration-success event. Then |d̄|≤B.
A uniform per-history bound is sufficient; pooled IID calibration of an
unconditional failure rate alone is not. A weaker high-probability bound on
this realized average may also suffice, if separately proved and its failure
risk included. Without such transport, the memory guarantee does not apply
to the calibrated physical null. To bound a nuisance-free effect d* when |d−d*|≤B,
the upper limit is min(1,U_abs+B), **not** U_abs−B. This remains useful after
a null result. Setting-predictability correction is already in the score
interval; do not count it a second time in B.

For the purely synthetic example N=100,000 and 50,000 matches, α=0.01,
ζ=0, the absolute observed-effect upper bound is 0.0102939956932. With B=0.002
the clean-effect upper bound is 0.0122939956932 (1.23%). With memory this bounds
|d̄*|; with a constant effect it bounds binary TV. This is an example calculation,
not a measured limit. Calibration failure probability αcal adds by union bound,
even when calibration records share a source. For M predeclared windows/bins,
allocate α/M to each and report family coverage; do not pick the best window.
The shipped design uses exactly one primary window and binary outcome.

## Prospective power

For the shipped power grid assume ζ=0. With a nonzero predictability allowance,
add 2ζ to the threshold (or replace B by B+2ζ for this calculation only).
Assume the *observed* conditional signed effect is at least d>B on every
history. A one-sided version of the same argument gives missed detection ≤β
whenever d>B+rα+sqrt(2 log(1/β)/N). Thus

\[
N>\frac{(\sqrt{2\log(2/\alpha)}+\sqrt{2\log(1/\beta)})^2}{(d-B)^2}.
\]

`sufficient_trials` uses strict integer rounding and rechecks this inequality.
These input d values are observed probability gaps, including losses; they are
not hidden influence parameters or predicted new-physics effects. If instead
d* describes the clean alternative and nuisance can cancel up to B, substitute
d=d*−B, so worst-case separation requires d*>2B before sampling error.
When d≤B no shot count can certify uniform power against every admitted null:
the specified alternative's table is already inside the nuisance envelope.

For IID stationary score trials K_matches is Bin(N,(1+d)/2). The code computes
the strict rejection-region tail sum and separately simulates it with PCG64,
seed 20260926. The numerical binomial power applies to that simulation model;
the preceding lower bound covers the stated history-uniform alternative.
Monte Carlo frequencies of 0 or 1 are not claims of zero risk or certain power.
Acquisition-time estimates omit calibration, resets, storage and control runs.
