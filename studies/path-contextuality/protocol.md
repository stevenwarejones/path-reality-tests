# Conditional protocol and statistical analysis

## Procedures and denominator

An eligible trial is a source herald recorded **before** independent context
selection. Keep its record even if the probe is lost, the final detector fails,
postselection fails, or more than one output fires. Fix a deterministic handling
of multiple clicks before acquisition (for example, classify as an explicit
invalid/no-detection outcome, and include its probability in the device model).
No background-subtracted counts are passed as binomial observations.

Required procedures are the two-mode source |ψ⟩, probe K, bypass I, phase reversal
Z, success/failure readout F, path readout E, trivial fair coin, and the four
calibration preparations |Q⟩, |P⟩, |+⟩, |+i⟩. X,Y,Z analyzers span all Hermitian
qubit effects. The preparations span the real vector space of Hermitian 2×2
matrices. On a trusted linear two-mode operational description, equality on these
spanning preparations and analyzers extends by linearity to all states/effects.
It does not prove that uncharacterized extra modes or generalized operational
states are absent. Preparation and measurement characterization are prerequisites,
not results of fitting the same witness data.

The ancilla is measured in a binary basis; no continuous pointer is truncated.
Record P(m,f|s,y,g) with m∈{−,+,lost}, f∈{success,failure,no detection}, preparation
s, final setting y and fixed g=(a,b). Loss categories of the complete table remain
in every denominator. The Lean bound applies to its binary coarse graining:
negative versus all other pointer records, and success versus all other final
records. The nine-cell table is still retained for calibration and diagnostics. Calibration contexts without a probe can use an explicit
“not applied” pointer label. All rates use their context's eligible-trial count.

## Calibration-to-premise map

| Procedure / measured quantity | Purpose | Additional premise still required |
|---|---|---|
| Source herald, randomized context and timestamp | Common preparation and stable sampling | No context-dependent source drift or selection; IID model for the stated statistics |
| Probe minus rate on characterized Q input | q (including post-instrument record erasure β) | This state reaches the response cap under the declared noisy-path equivalence |
| Four calibration preparations × probe / noisy sharp-path simulator | 8 contexts, compare three-outcome pointer distributions | Measurement equivalence for all allowed preparations; noncontextual representation of it |
| Four preparations × X/Y/Z analyzers × M/I/Z | 36 contexts, compare complete output distributions | Linear two-mode completeness and transformation equivalence M=(1−d)I+dZ, represented noncontextually |
| + input, M, characterized X readout | d via probability of the − result | Complete calibrated X detector, known readout error; unknown loss cannot be identified with d |
| Source, bypass, final F | f | Same preparation and final response as in probe trials |
| Source, probe, final F | a from negative pointer AND success | Complete trial table; no division by postselection counts |

The noisy sharp-path simulator uses a coin with p_m=2q/β−1 when β>0, followed by
E or a fair binary output, and the same record erasure channel. Operational
comparison of the full pointer distribution also checks the lost outcome.
Characterize the classical randomizer independently or include its uncertainty.
The I/Z comparison uses the *same fitted interval for d* across all 36 contexts,
not 36 independent choices. A joint confidence region propagates shared parameters.

The 44 audit contexts are acquisition settings, each with three outcomes. Two
coordinates per context determine the third, giving 88 marginal intervals.
Reported audit radius 0.01 is per coordinate; a difference of two unrelated
coordinates can have radius 0.02, and a inferred third coordinate radius 0.02.
No claim uses 0.01 as a uniform equivalence residual. Covariance between these
coordinates does not invalidate a Bonferroni bound, but independence cannot be
assumed for parameter fitting.

**Decision boundary:** A small residual or passing an interval consistency check
cannot establish an exact operational equality. This design proceeds only
conditionally on the two representation premises (or separately justified
near-model allowances). A fully empirical generalized-noncontextuality claim
would require a validated implementation of the exact equivalences, for example
a suitable secondary-procedure construction for the *instrument and channel*
family, including its uncertainty and accessible supplementary transformations.
A preparation-only mixing LP does not solve that channel problem. None is
claimed here. Without that input, report the witness and conditional conclusion.

## Confirmatory rule

Choose one (g,source,F), risk budget and sample allocation using theoretical or
separate pilot data. The supplied grid is wholly synthetic. For each of a,f,q,d,
retain successes k_i and eligible counts n_i. With main error α, construct four
two-sided Clopper–Pearson intervals at error α/4. Empty contexts get [0,1].
Let a_L,f_L,f_U,q_U,d_U denote endpoints. Reject only if

a_L > min(q_U, max(q_U f_L+d_U(1−f_L), q_U f_U+d_U(1−f_U))) + B.

Here B is a predeclared, independently justified bound on the actual witness
relative to a model with the representation premises. B=0 in the ideal reference.
One may set B=ε_a+|q−d|ε_f only with the explicit near-model premise and appropriate
uncertainty bounds; measured tomography residuals alone do not supply it.
The maximum over both f endpoints handles q<d. All four confidence intervals
cover simultaneously with probability ≥1−α by a union bound; this suffices for
the false-rejection guarantee without assuming their independence.

For the power certificate the implementation also supplies Hoeffding intervals
with radius r_α(n)=sqrt(log(8/α)/(2n)). The full capped witness is 1-Lipschitz
in each of its four coordinates on [0,1]^4. If the true gap is Δ, four IID
estimation errors each at most r_β occur with probability ≥1−β. Thus
Δ>4(r_α+r_β) guarantees rejection by the Hoeffding rule. The CP rule's Monte Carlo
power is reported separately; it is not the theorem certified by that calculation.

Audit intervals reserve a separate α_cal=0.005; main α=0.005. Combined statistical
coverage is ≥99.00%. This does not include an unquantified chance that the device
model or equivalences are false. The 90.00% power certificate is conditional on
valid device premises; it is not a guarantee that the apparatus will pass an
unimplemented equivalence-validation procedure.

## Trial allocation, selection and stopping

Fixed counts: draw a randomized permutation of the predeclared context labels,
with n trials in each of the four main contexts and n_cal in each of the 44 audit
contexts. Complete the schedule independently of observed outcomes. The budget
is 4n+44n_cal. The main q and d strata repeat settings in the audit family
using separate trials; 48 labels denote acquisition strata, not 48 distinct
physical configurations. This includes calibration audit collection but excludes unspecified
pilot or initial apparatus characterization; their unknown cost prevents an
end-to-end experiment cost claim. For n_cal=ceil(log(176/α_cal)/(2t²)), all 88
coordinate radii are at most t with simultaneous coverage ≥1−α_cal.

Random counts: an alternative fixes total N and assigns each of K=48 acquisition strata
independently and uniformly. Conditioned on n_i, IID outcome sampling gives the
same coverage. The probability that any n_i<N/(2K) is at most
K exp(−N/(8K)); `random_context_budget` adds that count risk to the power budget.
It retains a logarithmic bound to avoid representing numerical underflow as exact
zero. This deliberately conservative alternative has a much larger total budget.

Do not stop early, retry blocks until rejection, search among postselections, or
remove failed/lost trials. If several candidate witnesses are tested, preallocate
familywise α or use held-out confirmatory data. Memory/drift require a different
conditional-probability analysis; IID CP intervals are not justified by marginal
rates alone. Source/herald efficiency is needed to convert eligible trials to
emitted trials. Calibration detector specifications and pilot stability data are
not available, so no clock time or empirical feasibility is asserted.

## Joint decision using both success cells

The extended predeclared rule is `dual_decision`, with order (a,b,f,q,d).
Here b is **nonnegative-pointer success**, including pointer erasures in the
nonnegative group. All eligible trials stay in the denominator. At unit efficiency
this is ordinary positive-pointer success. a and b are disjoint cells of the same
probe records; bypass, q-calibration and d-calibration supply the other three
acquisition groups. Thus five estimated probabilities still use four main strata
and the existing 44 audit strata.

Use two-sided Clopper–Pearson intervals with error alpha/5 each. Dependence of a
and b does not invalidate the Bonferroni coverage. Reject if either

- a_L > min(q_U, max over f endpoints of q_U f + d_U(1−f)) + allowance;
- min over f endpoints of (1−d_U−q_U)f > b_U + allowance.

The endpoint minimum remains necessary if 1−d_U−q_U is negative. The allowance
must cover any separately justified model discrepancy. The rule is a valid
familywise test of two necessary inequalities, not an exact finite-sample test
of all four compatibility constraints.

The simultaneous Hoeffding alternative uses radius sqrt(log(10/alpha)/(2n)).
Both witnesses change by at most one per unit change of each participating input.
For equal counts, max(gaps)>4(r_alpha+r_beta) certifies power at least 1−beta.
`dual_certified_budget` includes all four main and 44 calibration groups. Snapshot
comparisons retain the single-witness budgets; adding a witness incurs a
multiplicity penalty. At the ideal point its conservative budget is slightly
worse; the declared lossy scenarios show a smaller sufficient budget, including
a scenario where only the positive witness survives. The shared
probe simulation uses multinomial records, not independent a/b binomial draws.

At the ideal reference point, the positive gap is 109/6250=0.01744, versus
297/15625=0.019008 for the negative witness. The rare positive cell is
144/15625. Its advantage is robustness to disturbance: thresholds 297/1225
versus 13/320 at q=16/25. The old grid and random-allocation outputs remain
explicitly single-witness baselines; the new dual budget is a fixed-count design.


### Loss-dependent power comparison and calibration ordering

These are prospective scenarios under the declared erasure/flip model, not
measured apparatus performance or universal statements about realistic loss.
All counts include the same calibration audit.

| Final efficiency / probe efficiency / flip | Negative gap | Positive gap | Negative-only trials | Joint-rule trials |
| --- | ---: | ---: | ---: | ---: |
| 1 / 1 / 0 | 0.019008 | 0.017440 | 4,351,904 | 4,432,252 |
| .9 / .95 / 0 | 0.013422 | 0.014011 | 6,411,796 | 6,221,616 |
| .8 / .9 / .01 | 0.007081 | 0.009852 | 17,065,004 | 10,228,608 |
| .5 / .8 / .02 | −0.005477 | 0.003771 | No violation | 56,409,916 |

Thus the second witness is slightly more expensive at the ideal point, improves
this sufficient trial budget in the specified lossy cases, and is the only
surviving witness in the last case. Probe loss lowers the calibrated negative
cap q′, raising the positive floor; b must simultaneously include lost-pointer
successes. The table evaluates both changes together. A smaller sufficient
Hoeffding budget is not a proof of sample-optimality.

The q-calibration preparation must attain the **largest eigenvalue of the actual
lossy negative effect**; a generic expectation value would underestimate the
required cap. With path-dependent erasure the effect may remain diagonal, but
its ordering must still be checked: E−=diag(η_Q·16/25,η_P·9/25), so Q maximizes
only if 16η_Q≥9η_P. Diagonality alone does not imply this. For example η_Q=1/10,
η_P=1 reverses the ordering. The characterization/tomography audit must support
both the eigenbasis and ordering (with uncertainty included). If it cannot,
characterize and prepare a maximizing eigenstate or use a conservative upper
bound on the largest eigenvalue; do not reuse the Q count as the cap. The existing
single-q Bernoulli calibration and power budget are conditional on that
characterization; extra characterization trials are not silently included.
