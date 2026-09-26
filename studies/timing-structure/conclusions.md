# Completed findings and limits of identifiability

## What the existing records answered

The fixed analysis was completed on both receivers and every declared lag,
local setting, history state, chronological period and uncertainty envelope.
All 107,109,596 overlap rows reconcile against the raw streams; the analysis
uses the same 107,109,468 interior rows for every lag. No-clicks remain outcomes.
Both raw extraction and aggregation reproduce with independently changed chunk
sizes. The analytic error guarantees and adversarial tests are supplied; neither
ordinary least-squares standard errors nor an IID permutation is substituted.

Within the stated assignment model, no ideal-assignment cancellation interval
has a positive lower bound, and no signed state interaction excludes zero.
Full-record event-supported hidden-TV upper bounds range from about 82 to 119
ppm of all interior rows. These are coarse, weighted, averaged-record bounds.
They do not exclude all effects within the less common state or all effects
within a category. Enlarging assignment or missingness uncertainty only weakens
the interpretation; both sensitivities are retained in the machine-readable output.

No future-setting freshness test rejects in the declared 64-test family.
Worst-case treatment of uncertain rows makes these diagnostics conservative.
The minimum known-event imbalance needed for each test to reject is included
in results.json, making a non-rejection's sensitivity explicit. Some descriptive
lag completion ranges exclude zero; these contain no sampling uncertainty and
MUST NOT be reported as significant lag effects. All lags were retained, including
those fluctuations. There is no fitted preferred lag or claimed propagation speed.

## What the recovery study establishes

Balanced-state synthetic reversing responses can be detected by the interaction
analysis while the pooled analysis misses them. Resolving a strictly positive
hidden-TV lower bound is harder than detecting a signed interaction. The actual
clock/recovery occupancies are uneven. In the occupancy-matched toy designs,
none of the four cases resolves a positive hidden-TV lower bound in 400 trials,
even with a strong rare-state redistribution chosen to cancel exactly in the
population average. Interaction detection ranges from 0/400 to 400/400. This is
a useful demonstrated sensitivity limit, not grounds for tuning a replacement
state definition on the same observed associations.

A process with fresh bits and ordinary outcome memory does not produce a
rejection in 400 synthetic repetitions. Persistent bits with ordinary delayed
response, and a deliberately misaligned record, each reject the ideal-freshness
null in 400/400. Their rejection therefore does not identify backwards causation.
The probabilities and sample sizes are declared toy models; these frequencies
are not calibrated hardware power or a proof of exact tail probabilities.

## Which ambiguity can this archive resolve?

| Ambiguity | Completed result | What would remove it | Supplied by this study? |
|---|---|---|---|
| Opposite effects cancel in averages | Conditional stratum and interaction bounds; exact pooled-null witness; recovery checks | Sufficient counts in an externally justified predictable state, with valid assignment assumptions | Bounds for the two stated states; no universal coverage of every possible state |
| Earlier/current/later setting association | Fixed descriptive lag map and a separately derived future-freshness test | Valid filtration/physical ordering and per-history setting calibration; a specified causal model | Neither certified by the archived bit frequencies; conclusions stay conditional |
| Different analog laws give identical integer tags | Exact within-cell construction; 0.98-bin illustrative gap | A calibrated within-bin response law or a finer/differently dithered measurement | No; digital data alone cannot determine sub-bin latent distributions |
| Different fine tag laws give identical coarse categories | Exact distinct-tag witness inside one category | Analyze finer categories or full tags with a new justified test family | The witness identifies this study's resolution limit; it does not claim the full archive lacks those tags |
| Direct effect versus latent common cause | Two structural causal models with identical observed tables and different interventions | An exogeneity/intervention premise or independent physical evidence | Not inferred from the observational table |
| Condition selected using current setting/outcome | Exact collider counterexample | Predictable state definition, fixed before looking at current records | Yes for the stated mathematical record-order model; physical timing assumptions remain explicit |

The quantization witness is a theorem about the stated ideal rounding map. For
any common integer-tag law Q, lift each tag K to K-a in one arm and K+a in the
other, with 0<a<1/2. Both push forward to Q while the event-conditioned supports
are disjoint. More generally the same lifting preserves *each* observed arm law
separately even if the laws differ: changing the lifts changes the latent mean
gap by 2a without changing any observed integer tags. The construction is
pointwise in K and can be applied to an entire recorded tag sequence; adding
no-event mass preserves the equality. It is not a claim that the apparatus
actually uses those latent point masses. Full physical plausibility requires
assumptions about jitter, synchronization and response which the equality alone
does not supply. Measured event frequencies merely scale hypothetical witnesses.

## Status

All planned code paths, mathematical arguments, counterexamples, synthetic
controls, raw-data checks and generated outputs are included. No unimplemented
analysis is needed to interpret this PR. The limitations above are explicit
non-identifiability or conditionality results, not promises that the archive
will answer them after one more regression. No new data collection was required.
This is a completed archive-methods study, not evidence for a law-of-physics
violation and not a certification of the experiment's physical assumptions.
