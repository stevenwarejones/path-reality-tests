# Retrospective click-only analysis (v2)

This replaces the headline centered-score bound, while preserving v1 in
[statistics.md](statistics.md) and `results/analysis.json`. The data, windows and
v1 results were already inspected when v2 was specified. This is a new
retrospective analysis, not preregistered or independently confirmatory evidence.
The finite grid and simultaneous family below are explicit. They cover choosing
within that family; they do not cover an arbitrary search over statistics,
pairing algorithms or outcome rules. Each version has its own 99% family
coverage under its premises; a joint claim across both versions would require
an additional error allocation (the union bound gives at least 98%).

## Target and click-only identity

Use the same potential-response probabilities and exogenous joint assignment
premise as v1. For fixed receiver setting y, under ideal uniform joint assignment,

    W'_i = 4 1{Y_i=y} (2X_i−1) B_i
    E[W'_i | H_i] = p_1y,i − p_0y,i = d_y,i.

No-click trials have score zero but remain in the full n-row denominator.
Thus setting imbalance among no-click trials cannot dominate the estimate.
The full-row click score is not exactly the descriptive difference of two
empirical arm rates: those use separate random denominators. Both are reported.
The target remains |mean(d_y,i)|. Changing signs can cancel; this does not bound
mean(|d_y,i|) or a maximum per-trial influence.

## A count-sensitive confidence sequence

For each receiver side, setting context (x,y), pulse group and phase radius,
define C_i = 1{X_i=x,Y_i=y} B_i in {0,1}, with predictable mean
mu_i = E[C_i | H_i]. For an interval starting at the deterministic row s, let
K_n = sum C_i and M_n = sum mu_i over its next n rows. For every real lambda,

    E[exp(lambda C_i) | H_i] = 1 + mu_i (exp(lambda)−1)
                              ≤ exp(mu_i (exp(lambda)−1)).

Consequently exp(lambda K_n − (exp(lambda)−1) M_n) is a nonnegative
supermartingale starting at 1. Ville's inequality bounds the probability that
it ever reaches 1/delta by delta, without an IID or fixed stopping-time premise.
The variance sensitivity comes from binary event counts: the resulting width
scales roughly with sqrt(K), rather than sqrt(n) from a range-only score bound.
This elementary Poisson-style envelope is conservative and not an exact
Clopper–Pearson interval or the empirical-Bernstein construction of Howard et al.

Use the fixed positive grid
L = {0.005, 0.01, 0.02, 0.04, 0.08, 0.16, 0.32, 0.64, 1, 2} and its negatives.
There are J = 2 receiver sides × 4 setting contexts × 3 phase radii × 3 pulse
groups = 72 event labels. For each signed grid element and start s≥0 allocate

    delta_s = alpha / [2 |L| J (s+1)(s+2)],   b_s = log(1/delta_s).

The start weights telescope to 1. The two signs, grid and J labels consume
alpha total. Ville already covers every endpoint, so no separate length or
reported-block multiplicity penalty is needed. Simultaneously, except on a set
of probability at most alpha,

    max(0, max_{lambda in L} (lambda K_n−b_s)/(exp(lambda)−1)) ≤ M_n,
    M_n ≤ min(n, min_{lambda in L} (−lambda K_n−b_s)/(exp(−lambda)−1)).

The denominators' signs explain the reversed upper inequality. Intersecting
multiple grid bounds is valid because every grid element was allocated error.
Empty intersections indicate failure of the conjunction of premises/confidence
events, not a zero upper limit. The implementation exposes an empty flag.

History may include device memory and drift. The sequence and each outcome rule
must be adapted: retrospectively fitting pairing, phase centers or undeclared
windows to these same outcomes is not justified by endpoint uniformity. This
remains an explicit ledger blocker, not a property inferred from good test results.

## Unknown outcomes: two envelopes

For each event label the trusted count is k. Timestamp-excursion and ambiguous
setting rows remain in n. Stored unambiguous settings restrict compatible labels;
codes 0/3 allow either latent binary setting. A latent binary intervention is
assumed even on ambiguous rows, as in v1.

The unrestricted envelope allows a click on every compatible uncertain row.
The event-supported envelope allows it only when that receiver has **any raw
channel-0 event in its sync interval**, regardless of phase or pulse number.
This needs reliable record ordering, complete detector-event recording within
retained trials, and settings consistent with unambiguous stored codes. It
allows arbitrary phase/window eligibility on supported rows. It does not cover
an event moved between sync intervals, a detector event lost by acquisition,
or additional unenumerated physical trials. Use the unrestricted envelope when
that event-support premise is not accepted. Neither is an apparatus certificate.

If m compatible rows can contain a click, pathwise k ≤ K_n ≤ k+m.
The confidence endpoints above are increasing in K_n. Substitute k in the
lower endpoint and k+m in the upper endpoint. This remains valid for
outcome-dependent missingness and for simultaneously completing all labels;
it never deletes a no-click or uncertain row. Aggregating event support across
all detector phases is intentionally conservative.

## Assignment and ordinary leakage sensitivity

Under exogeneity, mu_i = q_xy,i p_xy,i. Exact uniform assignment gives
mean(p_xy,i) = 4 M_n/n, and subtracting the two arm probability intervals bounds
the signed effect. A useful nonideal alternative is a **per-history** joint-TV
cap epsilon < 1/4, which implies q_xy,i in [1/4−epsilon,1/4+epsilon]. Then

    L_xy / [n(1/4+epsilon)] ≤ mean(p_xy,i)
                            ≤ U_xy / [n(1/4−epsilon)],

intersected with [0,1]. Subtract arm 0 from arm 1, widen by the assumed ordinary
leakage gap B, and intersect with [-1,1]. This cap must hold for every relevant
history; v1's average epsilon cap alone is insufficient for this inversion.
The observed setting imbalance is not such a calibration. Nor may an observed
click frequency be substituted for a uniform conditional p_max when asserting
a setting-bias term of order epsilon p_max. The inversion avoids that substitution.

The v2 leakage grid is B = 0, 1, 10 and 50 per million, on the scale of its
primary limits. The v1 reference retains its original 0, 1,000 and 10,000 per
million grid and numerical results.

All numerical epsilon and B values, including zero, are assumptions. Any
calibration failure probabilities must be added to alpha. Exogeneity, physical
alignment, adapted outcome rules, acquisition completeness and transported
leakage evidence remain uncertified. Therefore these are conditional signed-effect
bounds, not experimental exclusions of faster-than-light or retrocausal models
that preserve operational no-signaling.
