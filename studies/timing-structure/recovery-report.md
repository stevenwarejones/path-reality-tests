# Frozen recovery experiments

These are independent synthetic trials, not new physical experiments and not
estimates of power against a calibrated picosecond shift. Rates and 95% exact
binomial intervals are in recovery-results.json. No seeds or amplitudes were
changed after seeing the results. Statistical validity comes from the proofs.

Seed 20260928; 400 repetitions per case.

Detection of a weighted state-contribution contrast compares d_1-d_0, including
state prevalence. It alone is not evidence of different conditional effects or
cancellation; cancellation is assessed by the hidden-TV result. The JSON key
interaction_detected retains this weighted meaning.

## Sparse categorical responses

107,109,468 trials; two balanced fixed local states; fair independent settings.
Local-setting event probabilities .00015/.00045; conditional early/core/late
probabilities .25/.5/.25. Amplitude moves mass between early and late in
opposite directions for the two remote settings. Reversing also flips by state.
The same conservative inference constants as the actual study are used.

| Model | Amplitude | Pooled detected | Weighted state-contribution contrast detected | Hidden TV resolved |
|---|---:|---:|---:|---:|
| null | 0.0 | 0 | 0 | 0 |
| drift_null | 0.0 | 0 | 0 | 0 |
| fixed | 0.2 | 400 | 0 | 0 |
| reversing | 0.1 | 0 | 400 | 38 |
| reversing | 0.2 | 0 | 400 | 400 |
| reversing | 0.25 | 0 | 400 | 400 |

## Occupancy-matched cancellation stress test

Additional toy designs use each observed marginal condition occupancy, with
the same declared sparse event law. The rare-state amplitude is .25 and the
common-state amplitude is scaled by N_rare/N_common, forcing exact cancellation
in expectation. This panel uses descriptive occupancies, not fitted remote effects.
It is not a calibrated physical power curve. Strong rare-state effects can be
hard to detect when their weighted contribution is small.

| Receiver | Condition | Rare-state rows | Weighted state-contribution contrast detected | Hidden TV resolved |
|---|---|---:|---:|---:|
| alice | clock | 3,508,163 | 0 | 0 |
| alice | recovery | 6,334,206 | 350 | 0 |
| bob | clock | 6,639,349 | 378 | 0 |
| bob | recovery | 8,012,571 | 400 | 0 |

None of these four occupancy-matched cases resolves positive hidden TV in
400 repetitions. Consequently, the archive non-detections are not strong evidence
against structured effects. These toy results are not calibrated physical power.

## Temporal controls

200,000 trials, a drifting .01 baseline event rate, and ordinary outcome
dependence on an earlier setting. One fixed future-lag/feature/local-setting
test is evaluated using the full 64-test family threshold.

| Model | Freshness rejection at eta=0 | Rejection allowing eta=.4 |
|---|---:|---:|
| fresh_with_outcome_memory | 0 | 0 |
| persistent_delayed | 400 | 0 |
| one_row_misalignment | 400 | 0 |

Persistent settings have stay probability .9, so eta=.4 is necessary.
Their future association is produced by ordinary memory, with no influence
from the future. The misalignment case deliberately pairs the current causal
setting as if it were the next row: it is an indexing/order failure, not
retrocausality. The large eta=.4 row is a sensitivity illustration, not an
assertion that it certifies the misaligned process. Exact quantization and
coarse-category invisibility witnesses are verified separately in results.json.
