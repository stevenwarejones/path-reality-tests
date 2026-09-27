# Empirical-background injection recovery

**This is an injection study on one fixed local record, not a test of
the actual remote settings and not experimental evidence of signaling.**

See [protocol and limitations](empirical-protocol.md). Fresh independent
artificial setting bits are assigned to trials. All repetitions reuse the
same measured background; their intervals concern conditional recovery.

Receiver: bob. Verified overlap: 107,109,596 rows.
Retained: 107,107,360; exclusions before artificial assignment: 2,236.
Selected baseline-pulse events: 31,532 in 31,532 trials; 0 multiclick trials.
No-event trials retained: 107,075,828. No phase-radius cut.

All selected events remain in the external cache. The tested feature is
the first selected event per trial in original record order. All timing
bins and overflow tails survive injection. Other events in that trial
do not become independent observations.

## Actual local timing background

Residual phases are in 78.125 ps tag-bin units relative to the configured
local peak. Quantiles are descriptive; they do not establish clock accuracy.

| Segment | Local setting | Trials | Event trials | 1% | 25% | Median | 75% | 99% |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| training | 0 | 13,388,121 | 1,999 | -104.961 | -0.012 | 0.988 | 1.714 | 22.722 |
| training | 1 | 13,389,095 | 5,929 | -78.939 | 0.369 | 1.040 | 1.768 | 4.343 |
| test | 0 | 40,168,923 | 6,018 | -108.624 | -0.012 | 0.988 | 1.714 | 19.411 |
| test | 1 | 40,161,221 | 17,586 | -78.001 | 0.318 | 1.040 | 1.714 | 4.660 |

## Conditional detection frequency

400 artificial relabelings per injection; identical labels are reused
across injection sizes within each repetition. Training uses the first
quarter; all four tests use the remaining three quarters. Each rejects
at 0.01/4. Trained features use local-setting × block-parity strata.

| Injection | Arm shift (ps) | Count | Pooled timing | Trained timing | Trained histogram | Any test (95% MC interval) |
|---|---:|---:|---:|---:|---:|---|
| null | 0.00 | 0.00% | 0.50% | 0.25% | 0.00% | 0.75% (0.15–2.18%) |
| local_clock_control | 0 (local-only shift) | 0.00% | 0.00% | 0.00% | 0.00% | 0.00% (0.00–0.92%) |
| fixed_0.125 | 9.77 | 0.00% | 100.00% | 99.00% | 100.00% | 100.00% (99.08–100.00%) |
| fixed_0.25 | 19.53 | 0.00% | 100.00% | 100.00% | 100.00% | 100.00% (99.08–100.00%) |
| fixed_0.5 | 39.06 | 0.00% | 100.00% | 100.00% | 100.00% | 100.00% (99.08–100.00%) |
| fixed_1 | 78.12 | 0.00% | 100.00% | 100.00% | 100.00% | 100.00% (99.08–100.00%) |
| reversing_0.125 | 9.77 | 0.00% | 0.75% | 98.75% | 99.75% | 100.00% (99.08–100.00%) |
| reversing_0.25 | 19.53 | 0.00% | 0.75% | 100.00% | 100.00% | 100.00% (99.08–100.00%) |
| reversing_0.5 | 39.06 | 0.00% | 0.75% | 100.00% | 100.00% | 100.00% (99.08–100.00%) |
| reversing_1 | 78.12 | 0.00% | 0.00% | 100.00% | 100.00% | 100.00% (99.08–100.00%) |
| regime_transfer_failure | 78.12 | 0.00% | 0.00% | 0.00% | 0.00% | 0.00% (0.00–0.92%) |
| digital_fixed_2 | 156.25 | 0.00% | 100.00% | 100.00% | 100.00% | 100.00% (99.08–100.00%) |
| digital_fixed_4 | 312.50 | 0.00% | 100.00% | 100.00% | 100.00% | 100.00% (99.08–100.00%) |
| digital_reversing_2 | 156.25 | 0.00% | 0.00% | 100.00% | 100.00% | 100.00% (99.08–100.00%) |
| digital_reversing_4 | 312.50 | 0.00% | 0.00% | 100.00% | 100.00% | 100.00% (99.08–100.00%) |
| redigitized_fixed_0.125 | 9.77 | 0.00% | 74.25% | 51.00% | 47.50% | 84.00% (80.03–87.45%) |
| redigitized_fixed_0.25 | 19.53 | 0.00% | 100.00% | 100.00% | 100.00% | 100.00% (99.08–100.00%) |
| redigitized_fixed_0.5 | 39.06 | 0.00% | 100.00% | 100.00% | 100.00% | 100.00% (99.08–100.00%) |
| redigitized_fixed_1 | 78.12 | 0.00% | 100.00% | 100.00% | 100.00% | 100.00% (99.08–100.00%) |
| redigitized_reversing_0.125 | 9.77 | 0.00% | 0.75% | 51.25% | 45.50% | 68.25% (63.44–72.79%) |
| redigitized_reversing_0.25 | 19.53 | 0.00% | 0.50% | 99.50% | 100.00% | 100.00% (99.08–100.00%) |
| redigitized_reversing_0.5 | 39.06 | 0.00% | 0.25% | 100.00% | 100.00% | 100.00% (99.08–100.00%) |
| redigitized_reversing_1 | 78.12 | 0.00% | 0.00% | 100.00% | 100.00% | 100.00% (99.08–100.00%) |
| centered_latent_reversing_0.25 | 19.53 | 0.00% | 0.50% | 0.25% | 0.00% | 0.75% (0.15–2.18%) |

The local-clock control adds a regime-dependent 78.125 ps offset to
both artificial assignment arms equally. It is a null control.
The transfer-failure injection changes the relation between regime and
shift direction after training. Unlike the ideal Gaussian simulation,
real regime weights need not balance, so pooled cancellation need not
be exact and this case need not be completely invisible.

The fixed/reversing rows use symmetric floating-point perturbations of
decoded phases. Digital rows instead add symmetric whole-tag offsets,
which are exact integer offsets to copied detector time tags. Redigitized
rows assume uniform latent positions inside each original tag bin and
round after the shift; this is an explicit, uncalibrated rounding model.
The centered-latent counterexample places true arrivals at tag-bin centers;
its small analog shift changes no stored tag. Sub-bin power therefore
depends on an unobserved latent assumption. None establishes
instrumental picosecond resolution or arbitrary-memory real-label
significance. Null rates validate artificial randomization
on this fixed background only. Full per-test Monte Carlo intervals and
source hashes are in [empirical-results.json](empirical-results.json).
