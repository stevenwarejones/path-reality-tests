# Timing structure and what the archive can hide

One existing run; both receivers. Retrospective, assumption-conditional analysis.
The declared family separates current-setting effects, lag diagnostics and exact
information-loss witnesses. No physics violation or RNG certification follows automatically.

**Nothing was detected; these null results are not strong evidence against structured effects.**
None of the four occupancy-matched toy designs resolves positive hidden TV in
400 repetitions per case (see recovery-report.md). The displayed future-setting
examples need 873–1,247 events of imbalance to reject, versus 11–133 observed:
lag +1, any event, eta=0, event-supported completion, holding known-event total
and uncertain rows fixed. These limits do not establish that hidden effects exist.

## Source reconciliation

| Receiver | Reconciled rows | Interior rows | Full-overlap selected events |
|---|---:|---:|---:|
| alice | 107,109,596 | 107,109,468 | 33,608 |
| bob | 107,109,596 | 107,109,468 | 31,534 |

The same 128 endpoint rows are outside every lag/history population. Interior
uncertain rows and no-event rows are retained. Every local raw setting, legacy
click word and timestamp excursion reconciles. No remote outcome is used.

## Cancellation: fixed local-history conditions

Units are ppm of ALL interior rows. Bounds are weighted stratum contributions,
not probabilities conditional on a detected photon. The table uses epsilon=0
and the event-supported completion assumption. All other cases are in results.json.

| Receiver | Local setting | Condition | Sum of stratum TVs | Hidden TV beyond pooled |
|---|---:|---|---:|---:|
| alice | 0 | clock | [0.0, 84.7] | [0.0, 84.7] |
| alice | 0 | recovery | [0.0, 86.6] | [0.0, 86.6] |
| alice | 1 | clock | [0.0, 116.8] | [0.0, 116.8] |
| alice | 1 | recovery | [0.0, 119.1] | [0.0, 119.1] |
| bob | 0 | clock | [0.0, 82.2] | [0.0, 82.2] |
| bob | 0 | recovery | [0.0, 83.5] | [0.0, 83.5] |
| bob | 1 | clock | [0.0, 113.2] | [0.0, 113.2] |
| bob | 1 | recovery | [0.0, 113.5] | [0.0, 113.5] |

Across all 48 ideal-assignment condition/receiver/setting/period/envelope comparisons:
**0 positive hidden-TV lower bounds; 0 weighted state-contribution contrasts exclude zero.**
Overlapping comparisons are covered jointly, not treated as independent replications.
The weighted state-contribution contrast is d_1-d_0, including each state prevalence.
If conditional effects are equal to delta, this contrast is (w_1-w_0)*delta:
unequal state sizes alone can make it nonzero. It establishes neither different
conditional effects nor cancellation by itself. Cancellation refers to hidden TV.
The JSON key interactions retains this weighted meaning; see statistics.md.
Failure to resolve cancellation is not proof that instantaneous effects vanish.

### Observed condition occupancy

| Receiver | Clock state 1, known rows | Recovery state 1, known rows | Recovery state unknown |
|---|---:|---:|---:|
| alice | 103,581,980 | 100,755,823 | 20,661 |
| bob | 100,451,378 | 99,077,782 | 20,661 |

These fixed states were not balanced or optimized after seeing their occupancy.
Recovery means a detector record somewhere in the previous 64 rows; it is not a
measured dead-time classification. Clock uses the completed i-2 to i-1 interval.

## Temporal fingerprints

Lag means O_i against X_(i+lag). Positive lags use later archived settings.
All lag scores are descriptive 4*(k1-k0)/N; their missingness ranges are not
confidence intervals. Unequal assignment frequencies can change these scores.

| Receiver | Local setting | Lag | Any-event score (ppm) | Event-supported completion range (ppm) |
|---|---:|---:|---:|---:|
| alice | 0 | -64 | 0.71 | [-7.3, 8.8] |
| alice | 1 | -64 | 4.74 | [-4.0, 13.6] |
| alice | 0 | -16 | 0.93 | [-6.6, 9.2] |
| alice | 1 | -16 | 0.34 | [-7.7, 9.4] |
| alice | 0 | -4 | 2.80 | [-5.5, 10.3] |
| alice | 1 | -4 | -9.37 | [-17.9, -0.8] |
| alice | 0 | -1 | 1.46 | [-6.4, 9.3] |
| alice | 1 | -1 | -0.56 | [-9.5, 7.6] |
| alice | 0 | +0 | 1.76 | [-6.5, 9.2] |
| alice | 1 | +0 | -6.54 | [-15.8, 1.3] |
| alice | 0 | +1 | -2.13 | [-9.9, 5.9] |
| alice | 1 | +1 | 4.97 | [-3.4, 13.7] |
| alice | 0 | +4 | 4.52 | [-3.7, 12.1] |
| alice | 1 | +4 | -6.09 | [-14.6, 2.5] |
| alice | 0 | +16 | -3.02 | [-10.9, 5.0] |
| alice | 1 | +16 | -8.70 | [-16.6, 0.6] |
| alice | 0 | +64 | -3.25 | [-11.9, 4.2] |
| alice | 1 | +64 | -4.82 | [-13.6, 4.1] |
| bob | 0 | -64 | -2.13 | [-10.3, 5.0] |
| bob | 1 | -64 | 12.06 | [4.3, 19.5] |
| bob | 0 | -16 | -0.41 | [-7.5, 7.6] |
| bob | 1 | -16 | 0.34 | [-7.4, 7.4] |
| bob | 0 | -4 | 1.08 | [-6.5, 8.5] |
| bob | 1 | -4 | -5.79 | [-13.8, 1.0] |
| bob | 0 | -1 | 3.32 | [-4.9, 10.2] |
| bob | 1 | -1 | 7.13 | [0.4, 15.2] |
| bob | 0 | +0 | -6.76 | [-14.2, 0.9] |
| bob | 1 | +0 | -7.58 | [-14.4, 0.3] |
| bob | 0 | +1 | 2.50 | [-5.3, 9.7] |
| bob | 1 | +1 | -0.41 | [-7.7, 7.0] |
| bob | 0 | +4 | -2.95 | [-10.6, 4.4] |
| bob | 1 | +4 | 5.49 | [-1.2, 13.6] |
| bob | 0 | +16 | 2.88 | [-4.6, 10.5] |
| bob | 1 | +16 | 4.82 | [-2.5, 12.4] |
| bob | 0 | +64 | 0.56 | [-7.2, 8.3] |
| bob | 1 | +64 | -1.23 | [-8.8, 6.6] |

### Separate future-setting freshness test

This asks whether earlier local records predict a later remote bit under the
stated per-history freshness bound. It does not test retrocausality. The stake
mixture and all 64 tests are paid for; unknown rows receive worst-case factors.

| Assumed eta | Envelope | Rejections / 64 | Largest log E lower bound |
|---:|---|---:|---:|
| 0 | unrestricted | 0 / 64 | -51.740 |
| 0 | event_supported | 0 / 64 | -3.860 |
| 0.0001 | unrestricted | 0 / 64 | -51.758 |
| 0.0001 | event_supported | 0 / 64 | -3.875 |
| 0.001 | unrestricted | 0 / 64 | -51.890 |
| 0.001 | event_supported | 0 / 64 | -3.950 |
| 0.01 | unrestricted | 0 / 64 | -53.216 |
| 0.01 | event_supported | 0 / 64 | -4.219 |

### How much imbalance would the future test need?

Fixed example: lag +1, any-event feature, eta=0, event-supported completion.
The threshold holds known-event total and uncertain rows fixed. It is a
sensitivity calculation, not an additional test or a fitted effect estimate.

| Receiver | Local setting | Known events | Uncertain candidate rows | Observed absolute imbalance | Minimum rejecting imbalance |
|---|---:|---:|---:|---:|---:|
| alice | 0 | 9535 | 421 | 57 | 925 |
| alice | 1 | 24065 | 457 | 133 | 1247 |
| bob | 0 | 8015 | 403 | 67 | 873 |
| bob | 1 | 23509 | 394 | 11 | 1173 |

## Exact limits on interpretation

- **Full-tag invisibility:** opposite offsets of 0.49 bins inside each rounding cell
  give a 0.98-bin event-conditioned analog shift with identical recorded tags.
  Under the adopted 78.125 ps/bin unit this illustrative shift is 76.5625 ps.
  The event-conditioned latent TV is one; full-trial latent TV is the event probability.
  This is a quantizer-model witness, not a measured apparatus delay or sensitivity limit.
- **Coarse-summary invisibility:** moving between distinct recorded tags in the same
  early/core/late category is invisible to this study even though the full tags differ.
- **Cancellation:** the exact toy model has pooled TV zero but state-averaged TV
  200 ppm. A pooled null therefore does not exclude structured effects.
- **Lag ambiguity:** an ordinary one-row delayed response with persistent settings
  produces current AND future associations. Equal observational tables can also
  arise from direct influence or a latent common cause with different interventions.
- **Selection artifact:** conditioning on O=X manufactures a perfect association
  from independent variables. Our chosen states exclude current O and X.

Exact fractions and four measured-frequency witness anchors are in results.json.
No latent law was fitted to establish these equalities.

## Decision

This completes the declared three-way interrogation of existing records. Its
deliverables are reproducible conditional bounds, separately justified temporal
diagnostics, and constructive identifiability limits. The archive alone supplies
no certified assignment/calibration premise that converts them to a new physics claim.
See conclusions.md for the precise scope and which additional assumptions break
each ambiguity; those are limits of the result, not unimplemented analysis items.

See protocol.md, statistics.md, sources.md and validation.md. All fixed results,
including halves and uncertainty sensitivities, are committed in results.json.
