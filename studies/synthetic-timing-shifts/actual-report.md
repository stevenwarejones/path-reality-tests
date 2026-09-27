# Actual-setting timing feasibility

Retrospective exploration of the pinned run4 overlap. These are 99% simultaneous
**assumption-conditional** confidence bounds, allowing arbitrary temporal memory.
Assignment, record completeness, pairing and physical nuisance calibration remain unverified.
No result below proves or disproves a law of physics.

## Reconciliation and population

| Receiver | Archived rows | Uncertain rows retained | Selected detector events | Multiclick rows |
|---|---:|---:|---:|---:|
| bob | 107,109,596 | 19,968 | 31,534 | 0 |
| alice | 107,109,596 | 19,968 | 33,608 | 0 |

Every local raw setting, legacy click word and timestamp excursion reconciles exactly.
No-event rows remain in the denominator. The fixed outcome is no event / early / core / late
for the first detector record in baseline pulse bits 4, 5, 6, without a phase-radius cut.

## Full-record bounds under ideal assignment (epsilon = 0)

TV is total variation between the two **time-averaged four-category record laws**.
Units are parts per million of all archived trials, not of detected photons.
These do not bound fine-scale timing shifts, instantaneous effects or effects that cancel over time.

| Receiver | Local setting | Unrestricted TV interval (ppm) | Event-supported TV interval (ppm) |
|---|---:|---:|---:|
| bob | 0 | [0.0, 432.2] | [0.0, 53.8] |
| bob | 1 | [0.0, 443.7] | [0.0, 76.9] |
| alice | 0 | [0.0, 427.2] | [0.0, 54.7] |
| alice | 1 | [0.0, 445.9] | [0.0, 80.1] |

The event-supported envelope assumes complete receiver detector records; the unrestricted
envelope allows an event in every compatible uncertain row. Neither covers missing physical trials.

## Fixed chronological cross-check

| Receiver | Local setting | Period | Event-supported TV interval, epsilon=0 (ppm) |
|---|---:|---|---:|
| bob | 0 | first_half | [0.0, 72.6] |
| bob | 1 | first_half | [0.0, 104.1] |
| bob | 0 | second_half | [0.0, 66.9] |
| bob | 1 | second_half | [0.0, 94.8] |
| alice | 0 | first_half | [0.0, 71.4] |
| alice | 1 | first_half | [0.0, 111.7] |
| alice | 0 | second_half | [0.0, 68.7] |
| alice | 1 | second_half | [0.0, 104.2] |

Alice and Bob share this run and apparatus. This is a second-receiver cross-check,
**not an independent-run replication**. The halves also share apparatus and are not independent replications.

## Assignment sensitivity

| Assumed per-history epsilon | Largest full-record event-supported TV upper bound (ppm) |
|---:|---:|
| 0 | 80.1 |
| 0.0001 | 80.4 |
| 0.001 | 83.7 |
| 0.01 | 117.0 |
| 0.05 | 274.9 |
| 0.1 | 533.0 |

At epsilon=0: 0 of 24 fixed receiver/setting/period/envelope TV intervals
have a positive lower bound; 0 have incompatible confidence/model constraints.
These are overlapping comparisons with joint error control, not independent tests.
Observed setting frequencies cannot certify the required per-history bound.

For a sum-of-arm nuisance TV budget B, subtract B from the lower TV bound (floor zero)
and add B to the upper bound (cap one). The JSON includes B=0, 1, 10 and 50 ppm.
Those are sensitivity assumptions, not measured apparatus tolerances.

## Feasibility decision and next gate

- **Proceed with conditional archive methods:** both receivers reconstruct and the same fixed
  analysis runs with no-clicks, temporal-memory bounds and missing-record sensitivity.
- **Do not advance to a physics claim yet:** no certified per-history assignment bound,
  causal nuisance budget or independent-run replication has been supplied.
- Next obtain calibration evidence and an independently selected run, freeze this protocol,
  then rerun it without tuning gates or categories. Decide the scientifically useful effect
  threshold before examining that run; compare it with the calibrated upper bounds.

All signed feature intervals and sensitivity comparisons are in `actual-results.json`.
See `actual-protocol.md` for the derivation, estimand and limitations.
