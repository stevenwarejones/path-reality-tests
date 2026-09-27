# Actual-setting timing feasibility protocol

This is a retrospective exploratory protocol, fixed before extracting the
comparisons below. The earlier audit and Bob marginal timing distribution have
already been inspected; this is not a preregistration or an independent discovery
sample. It extends PR #10's artificial-label recovery experiments.

## Fixed population and outcomes

Use the pinned run4 HDF overlap, all 107,109,596 archived rows. Reconcile each
receiver's raw settings, legacy click words and timestamp excursions against the
HDF. Remote and local settings use the archive's row pairing. In each receiver,
select baseline pulse bits 4, 5, 6 without a phase-radius cut. Select the first
such event in original raw-record order. The categorical outcome is no event,
early (phase minus configured peak < 0), core (0 <= residual < 2 bins), or late
(residual >= 2 bins). All phase tails are included. No event is a real outcome,
not a discarded row. These coarse bins test a recorded distribution; they do
not estimate a continuous-time delay or a full point-process distribution.

Compare actual remote setting 1 minus 0 at each local setting, for Bob and
Alice, on the full overlap, first half and second half. Halves are fixed by
source row index (floor(n/2)); no window, gate, phase or outcome optimization.
The four binary features are any event, early, core and late. No-event contrasts
are the negative of any-event contrasts.

A row is uncertain if either setting is ambiguous or it falls in either side's
paired timestamp-excursion ranges, including the decoder's guard. Keep all these
rows. Lower event counts use trusted records only. Upper completions include
all compatible uncertain rows (unrestricted), or only compatible uncertain rows
containing any receiver detector record, regardless of pulse/phase
(event-supported). The latter requires complete, correctly ordered detector
records. Known settings restrict compatibility; ambiguous codes allow either
setting. Neither envelope covers physical trials absent from the archive or
arbitrary pairing failures outside the marked ranges.

## Simultaneous inference and assumptions

For each binary event-and-setting-context indicator C_i, let mu_i be its
conditional mean given the past. For fixed signed lambda,
exp(lambda sum C_i - (exp(lambda)-1) sum mu_i) is a nonnegative
supermartingale. This follows from 1 + mu*(exp(lambda)-1) <=
exp(mu*(exp(lambda)-1)); it does not assume independence across trials.
Use the fixed magnitudes (.005,.01,.02,.04,.08,.16,.32,.64,1,2), both signs,
32 labels (2 receivers * 4 setting contexts * 4 features), and two deterministic
starts (0 and floor(n/2)). Union allocation alpha=.01 gives log threshold
log(2*10*32*2/.01). Ville's inequality covers both first-half/full endpoints
from start 0 and the second half from the other start. Interval-censored counts
widen the same event, so both completion envelopes and nuisance grids share
this simultaneous coverage. No IID permutation p-value is used.

To infer setting-conditioned response laws, ASSUME exogenous current setting
assignment with per-history joint assignment probabilities in
[.25-epsilon, .25+epsilon], even conditional on relevant device state before
assignment. This is a causal/randomness assumption, not a consequence of
observed setting frequencies. With factorization mu_i=q_xy,i*p_xy,i, invert the
count bounds using n*(.25 +/- epsilon). Differences cover time-averaged response
probabilities at a fixed local setting, averaging over the realized pretrial
histories. Arbitrary temporal memory is allowed under those assumptions.
Effects may cancel across time; these are not bounds on average absolute
per-trial effects, long-term interventions, or worst-case effects.

Report feature contrasts and a conservative TV interval for the four-category
averaged record laws. TV = half the sum of absolute category contrasts.
Use epsilon = 0, .0001, .001, .01, .05, .1 as assumed sensitivity values.
Separately allow B = 0, 1e-6, 1e-5, 5e-5 for the SUM of the two arm TV distances
between observed and target laws; triangle inequality widens TV by B once.
No calibration certificate for epsilon or B is supplied by this study. If a
future certificate can fail with probability eta, add eta to the .01 error
budget. Event-supported completeness and archived row pairing are additional
assumptions. Confidence is conditional on all these assumptions.

## Decision rules and replication

Publish every fixed comparison, including nulls and both missingness envelopes.
Proceed as an archive-methods/sensitivity project if reconstruction is exact
and conditional bounds remain interpretable. A nonzero lower bound is only an
assumption-conditional recorded association requiring nuisance investigation;
it is not a violation of relativity. A zero lower bound is not proof of no
signalling. Alice is a second-receiver cross-check within the SAME run, sharing
photons, apparatus and timing; it is not an independent-run replication.

Before any physics claim: obtain or justify per-history assignment and causal
nuisance calibration, verify physical timing/pairing/completeness, and freeze
an independently sourced second-run analysis. Without those, stop at
conditional archive-level feasibility. Do not select a second run on the basis
of which result it produces.
