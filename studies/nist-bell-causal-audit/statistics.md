# Conservative centered-score reference (v1)

The primary retrospective v2 analysis is derived in [click-statistics.md](click-statistics.md).
This reference is retained because its broader unknown-score allowance remains valid
without the additional detector-event support restriction. At this run’s click rates
it is too imprecise to exclude a setting completely suppressing receiver clicks.

This is a retrospective audit of one documented NIST run. The primary record
is any receiver click in zero-based pulse bits 4, 5, 6. No click is outcome 0;
multiple clicks in the group are outcome 1. Receiver settings remain separate.
The two other pulse groups and three phase radii are declared in protocol.json.
The ten blocks are contiguous row-index partitions, not equal wall-clock times.

## Target and assignment premise

For Alice → Bob at fixed receiver setting y, write
p_xy,i = P(B_i=1 | do(X_i=x,Y_i=y), H_i). H_i includes all relevant pre-choice
history and device memory. Require consistency and exogenous joint assignment
relative to the potential receiver responses. This is a causal assumption;
balanced observed settings do not establish it. Let d_y,i=p_1y,i−p_0y,i.
The reverse direction exchanges parties. The target in an interval is the
realized mean of these conditional signed effects.

For an exactly uniform joint assignment conditional on history, define

    W_y,i = 2 1{Y_i=y} (2X_i−1)(2B_i−1).

Then E[W_y,i | H_i]=d_y,i, by summing over the four settings. W is in [-2,2].
The factor 2 accounts for the receiver-setting probability. We use the full
paired-row denominator; selecting the receiver-setting rows and pretending
that their random sample size was fixed would require another argument.

More generally, if the joint setting law has history-conditional TV distance
at most epsilon_i from uniform, the expectation differs from d_y,i by at most
4 epsilon_i: the range of the conditional response score is at most 4.
A bound on the realized average epsilon_i over the reported interval suffices.
Observed setting frequencies are not that bound. This premise is stronger and
different from simply importing the paper's Bell-test excess-predictability
number. The latter is not used as this audit's calibration.

## Simultaneous coverage over memory and retrospective endpoints

For D_i=W_i−E[W_i|H_i], the conditional bounded-range exponential inequality
is E[exp(t D_i)|H_i] ≤ exp(2t²). Iteration and Markov's inequality give

    P(|mean(W−E W)| ≥ r) ≤ 2 exp(−n r²/8)

for any deterministic contiguous interval of length n. No IID receiver model
is assumed. This follows the same bounded-range argument derived in the
[protocol study](../spacetime-causal-influence/statistics.md), with doubled range.

To avoid assuming that archival stopping and the chosen block endpoints were
fixed in advance, allocate error to EVERY start s≥0 and length n≥1:

    alpha_s,n = alpha / [M (s+1)(s+2) n(n+1)]
    r_s,n = sqrt(8 log(2M(s+1)(s+2)n(n+1)/alpha) / n).

Both reciprocal sequences telescope to sum 1. A union bound across starts,
lengths and the M declared score labels gives simultaneous error at most alpha,
even if endpoints are chosen after inspecting the sequence. M=396 conservatively
allocates a label to each of 3 phase radii × 3 pulse groups × 4 directed receiver
contrasts × 11 displayed intervals. This is deliberately conservative; repeated
labels could be allocated more efficiently. It is not a newly claimed theorem.

This argument requires the underlying ordered paired sequence and its outcome
rules to satisfy the premise at every row. It does NOT justify postselection,
unrecorded interior trial deletion, a data-dependent reassignment of trial pairs,
or an arbitrary search over undeclared window rules. Raw reconciliation supports
record fidelity, but does not independently establish physical alignment or RNG
exogeneity. The phase calibration is inherited from the archive. The confidence premise
requires its rule to be fixed independently of the scored future outcomes (or
covered by a larger simultaneous family). The archive does not establish this
adaptedness; it is another explicit blocker, even with uniform endpoints. Its transfer
to a secure physical record remains a separate premise.

## Unknown rows and ordinary leakage

Two setting codes are ambiguous. Completing their scores assumes that a binary
assignment/response record exists behind the ambiguity; if the apparatus instead
performed an undefined intervention, an additional contamination model is needed. We also treat all rows spanning each paired
huge timestamp excursion (through the following interval) as unknown. Their
union has 19,968 rows including the ambiguous settings. None is deleted.
For m unknown scores, every completion changes the score mean by at most 2m/n
relative to a sum using zero for the unknown entries. This is a pathwise
statement, so missingness may depend on settings or outcomes. It covers only
unknown scores whose trial slots and denominators are known, not uncounted
missing physical trials. `missing_outcome_difference` separately gives exact
identification limits for known per-arm missing-outcome counts; tests enumerate
all completions.

For average joint-setting TV allowance epsilon and a coupling-derived bound B
on the shift from observed to nuisance-free signed effects, the interval is

    [known_score_sum/n − r_s,n − 2m/n − 4epsilon − B,
     known_score_sum/n + r_s,n + 2m/n + 4epsilon + B] ∩ [-1,1].

Calibration failure risks must be added to alpha. The code widens by B; it does
not subtract B from an upper bound. An empty intersection is a failed conjunction
of premises, not an upper limit. Source-backed history-conditional calibration
and ordinary-leakage transport are unavailable here. The numerical epsilon and
B grid is explicitly assumed. Even epsilon=B=0 is an assumption, not an estimate.

The descriptive difference of arm rates and the centered score need not agree
numerically: the latter also reflects finite-sample setting imbalance and uses
worst-case treatment of ambiguous records. Do not interpret that difference as
a detected physical effect. The conservative bounds trade precision for coverage
under memory, uncertain scores and retrospective endpoints.

## Limits of a null result

These intervals cover |mean(d_y,i)|, not mean(|d_y,i|). Alternating effects can
cancel. Ten blocks provide a diagnostic at one temporal scale, not a proof of
zero effect on every trial. A constant-effect assumption would identify the
interval with a binary total-variation influence bound for that receiver setting.

No rejection of hidden superluminal or retrocausal models preserving operational
no-signaling follows. No backward do-effect is identified. This audit does not
reanalyze the published Bell-test rejection or claim that a reconstructed local
marginal difference is a violation of relativity.
