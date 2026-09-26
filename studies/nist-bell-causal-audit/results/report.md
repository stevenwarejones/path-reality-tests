# Recorded marginals and conditional sensitivity

The source contains 107,109,596 aligned stored rows. Both raw setting streams
and the archived click-update behavior reconcile exactly over that interval.
The analysis below uses bitwise-OR reconstructed clicks, includes no-click rows,
and keeps ambiguous settings in the population. No spacelike causal claim is enabled.

## Primary three-pulse record

The marginal differences describe valid-setting rows. The intervals instead use
the full-row click-only score and count-sensitive confidence sequence.
19,968 rows remain uncertain; 877 Alice and 797 Bob rows contain any raw
detector event. The tighter envelope requires the explicit detector-record
completeness/order premise in [click-statistics.md](../click-statistics.md).
Here joint-setting TV allowance and ordinary-leakage allowance are **assumed zero**.
Bounds are on absolute average signed effects, not average absolute influences.

All numbers below are probabilities per million trials. Baseline comparisons
are descriptive scale comparisons, not confidence bounds on a relative effect.

| Direction | Receiver setting | Arm click rates | Descriptive gap | Click-only score | Event-supported upper | Unrestricted upper | v1 upper |
|---|---:|---:|---:|---:|---:|---:|---:|
| A → B | 0 | 140.0 / 131.8 | -8.19 | -8.29 | 38.64 | 226.92 | 2652.44 |
| A → B | 1 | 428.3 / 420.9 | -7.45 | -7.43 | 54.82 | 237.07 | 2314.13 |
| B → A | 0 | 129.8 / 130.9 | 1.04 | 1.05 | 30.95 | 217.34 | 2546.20 |
| B → A | 1 | 399.0 / 395.8 | -3.17 | -3.14 | 51.91 | 233.18 | 2425.42 |

The event-supported upper limits are 13.0–29.3% of the smaller observed arm click rate.
The largest is 0.005482%, compared with the old 0.2652% reference.
The old limit is 6–20 times the click rate and cannot exclude complete
suppression of a detector at these rates. The unrestricted v2 envelope
remains wider than the lower-rate baselines. The useful tighter comparison
therefore depends materially on the detector-event support premise.
This revised analysis was specified after inspecting the original results;
it is retrospective, not preregistered or independently confirmatory.

The simultaneous interval construction allocates total error 0.01 across the
72 event labels, a fixed lambda grid and every contiguous interval; it tolerates arbitrary device
memory under the stated conditional assignment premise. Every reported interval
contains zero. This does not validate the premise or certify no signaling.

## Assumed calibration sensitivity

Largest primary upper limit across both directions and receiver settings:

| Assumed per-history joint-setting TV cap | Assumed leakage gap (per million) | Event-supported upper limit (per million) |
|---:|---:|---:|
| 0 | 0 | 54.82 |
| 0 | 1 | 55.82 |
| 0 | 10 | 64.82 |
| 0 | 50 | 104.82 |
| 0.0001 | 0 | 55.16 |
| 0.0001 | 1 | 56.16 |
| 0.0001 | 10 | 65.16 |
| 0.0001 | 50 | 105.16 |
| 0.001 | 0 | 58.25 |
| 0.001 | 1 | 59.25 |
| 0.001 | 10 | 68.25 |
| 0.001 | 50 | 108.25 |

These allowances are sensitivity parameters, not measured calibration results.
Lack of a valid history-conditional calibration prevents promoting these
numbers to an apparatus-certified causal bound.

## Reconstruction and window sensitivity

| Side | Nominal words differing | Narrow-radius words changed | Wide-radius words changed |
|---|---:|---:|---:|
| alice | 105 | 5036 | 1093 |
| bob | 94 | 1097 | 337 |

## Selection diagnostic

Conditioning on a sender click yields the following receiver differences.
These are selected correlations, not causal effects; the primary analysis
does not make this selection. Their size illustrates why complete trials matter.

| Direction | Receiver setting | Sender-click-selected difference |
|---|---:|---:|
| A → B | 0 | -0.4561 |
| A → B | 1 | -0.6887 |
| B → A | 0 | -0.4154 |
| B → A | 1 | -0.6620 |

The bitwise-OR words differ from the archived words because of buffered
fancy-index semantics. Effect on published analyses has not been assessed.
The reconstructed record has these additional primary-window
receiver clicks across all setting codes (counts, not inferred effects):

- Alice: 31.
- Bob: 27.

The original diagnostic table is reproduced after assigning the two A=3
no-click rows to A=2. The causal audit does not adopt that ambiguous-setting
assignment. This establishes count reconciliation, not a reproduction of the
published Bell-test p-value or a challenge to that result.

![Drift and conditional sensitivity](diagnostics.svg)

The left panel shows descriptive contrasts for ten contiguous index blocks.
The right panel shows full-run conditional limits under ideal joint assignment
and zero ordinary leakage for all pulse groups and phase radii. Narrow/wide
change each phase radius by one timetagger bin (78.125 ps), keeping its center
fixed. Pulse grouping is distinct from phase-radius variation. Neither is a
new spacelike-separation calibration. No best window is selected.
