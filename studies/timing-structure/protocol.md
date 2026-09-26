# Protocol: timing structure and limits in an existing Bell archive

This protocol fixes new comparisons before their extraction. It is retrospective:
PR #10's marginal timing results and the existing causal audit have been seen.
There is no untouched discovery sample or preregistration. All analyses use the
same pinned 03_43 run4 originals. No new experiment or separately selected run
is required for this study. This is not independent replication of PR #10.

## Shared records

Reconcile both raw streams to every local setting, legacy click word and paired
timestamp excursion in the HDF overlap. Cache only outside the repository.
Retain every interior source row i in [64, n-64), including no-event rows.
Exactly 128 boundary rows are outside the declared population because every lag
and history must exist. This is fixed before inspecting outcomes. Neither end
is wrapped around. Full, first-half and second-half summaries use the original
source index floor(n/2), not a data-dependent cut.

The outcome is PR #10's first event in baseline pulse bits 4,5,6, with no phase
radius cut: no event, early residual<0, core 0<=residual<2, late residual>=2.
Residual is decoded phase minus configured local peak. Preserve all phase tails.
The four binary features are any event, early, core, late; no-event contrasts
are minus any-event contrasts. Raw multiplicity is counted even though only the
first selected event defines the feature. The source decoder remains unchanged.

Uncertain rows include either side's timestamp excursion plus decoder guards,
or ambiguous local/remote settings. Keep them as interval-censored data.
Unrestricted completion allows an event in any compatible uncertain row.
Event-supported completion allows it only if the receiver has any detector
record, across ALL pulses and phases. The latter assumes detector-record
completeness/order; neither envelope covers physical trials absent from the
archive. Known setting codes restrict possible completions; ambiguous codes
allow either setting. No remote outcome enters the features or selection.

## A. Cancellation hidden by averaging

Two fixed, binary, local-history conditions:

1. Clock: duration from receiver sync i-2 to i-1 exceeds 129102 timestamp tags.
   This interval has ended before row i; it is not a fit to current outcomes.
2. Recovery proxy: at least one receiver channel-0 detector record in rows
   i-64,...,i-1, across all pulses. This is a history proxy, not a calibrated
   detector dead-time model.

Report pooled and both states of each condition, each local setting and receiver,
on full/first-half/second-half populations. An uncertain clock-history row makes
clock state unknown; any uncertain row in the recovery window makes recovery
state unknown. Such rows are compatible with either state, not deleted.

Bound weighted stratum contributions on the ALL-interior-row scale, rather than
dividing by a possibly sparse/uncertain stratum occupancy. Infer current-setting
responses under the explicit per-history exogenous joint-assignment assumption
q_xy in [.25-epsilon,.25+epsilon]. Use epsilon=0,.001,.01. No certificate is
inferred from empirical frequencies. The event/stratum/context indicators are
Bernoulli adapted variables. A finite-lambda supermartingale inversion uses
alpha=.005, 160 labels (2 receivers*5 groups*4 contexts*4 features), two starts
and both signs of 10 fixed lambdas. This covers all reported periods jointly.
See statistics.md for the complete estimand and proof.

For each condition, bound the sum of stratum TVs and its excess above pooled TV.
A positive lower bound on that excess would establish cancellation within this
coarse conditional model. Otherwise report a bound, not an inferred hidden effect.
Also report the signed state-1 minus state-0 feature interaction. No fitted
classifier, chosen sign, window search or conditioning on a current outcome.

## B. Temporal fingerprints

Fixed lags l=-64,-16,-4,-1,0,1,4,16,64, defined by receiver O_i against remote
X_(i+l). Positive means the remote setting is later in archived row order.
At every lag publish exact trusted counts, exposures and missing-row completion
ranges for the assignment-normalized signed score 4*(k1-k0)/N. Nonzero-lag
scores are descriptive, not causal effect estimates or confidence intervals.
Do not apply the current-setting q_xy inversion to past settings. No circular
shifts, permutation p-values, IID bootstrap or uncorrected lag maximization.

For the four POSITIVE lags only, implement a separate assignment-time betting
test. The earlier local feature and local setting must already be measurable
before the later remote bit is assigned. Under the null that each remote bit
has conditional probability in [.5-eta,.5+eta] given the entire past, the fixed
stake factors have conditional expectation <=1. This is a setting-freshness
/predictability diagnostic, not a test of retrocausality. The archive cannot
certify the required cross-station physical ordering; state it as an assumption.

Use both signs of stakes .005,.01,.02,.04,.08,.16,.32,.64 in an equal mixture.
Test only the full interior population. Allocate alpha=.005 over 64 tests
(2 receivers*2 local settings*4 features*4 positive lags). Unknown contributions
receive a pathwise worst-case factor; report both completion envelopes. Show
eta=0,.0001,.001,.01 as nested sensitivity hypotheses. Do not optimize a lag,
stake, period or eta for an unadjusted discovery. Sum family error with part A
for a .01 conditional error budget. A zero-lag confidence bound and a
future-setting test address different nulls and must not be merged into one
causal statement.

## C. Identifiability and adversarial witnesses

Provide executable exact-rational witnesses for:
- opposite analog offsets inside each timestamp cell that leave the entire
  recorded tag distribution identical;
- distinct fine timing laws with identical coarse timing-category laws;
- opposite effects in two local states that cancel in the pooled distribution;
- ordinary delayed response with persistent settings producing current and
  future associations, plus two causal models with identical observed tables;
- a current-outcome/setting-based state selection that manufactures a conditional
  association, demonstrating why part A requires predictable conditions.

Attach the first witnesses to the archive's measured event frequencies without
claiming that their hypothetical latent laws describe the apparatus. Distinguish
full-tag invisibility from coarse-summary invisibility. Quantization witnesses
are conditional on the stated ideal quantizer; physical timing calibration is
not inferred. Show explicit assumptions that would remove each ambiguity and
whether this archive supplies them. Proving an information-loss limitation is a
completed result, not an unfinished empirical search.

## Validation and stopping rule

Use small exact adaptive-memory enumerations, independently computed toy counts,
chunk-boundary/lag tests, source reconciliation with two chunk sizes, and frozen
synthetic recovery experiments including null, reversing, memory and quantization
controls. Report all seeds, grids and outcomes; Monte Carlo checks do not prove
tail calibration. Full-source CI must re-extract both receivers and reproduce
committed aggregate results. Stop after this declared family; do not add tests
in response to interesting fluctuations. Existing methods and their scope will
be credited in sources.md; novelty is not presumed.
