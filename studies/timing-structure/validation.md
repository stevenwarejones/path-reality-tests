# Validation record

## Mathematical and implementation checks

- Independent row-by-row oracle for every group and lag, including odd midpoint,
  both ambiguous setting codes, timestamp guards, unknown states, all event
  categories and no-events. Streaming aggregation agrees at chunk sizes 1, 17
  and 1,000 on the synthetic fixture.
- Exact enumeration of adaptive, history-dependent setting processes verifies
  the fixed-stake mixture expectation bound; all three-valued completions of
  missing rows verify the pathwise conservative factor.
- Count inversion checked against every paid exponential boundary; assignment
  and missingness sensitivity checked for monotone widening.
- Exact integer rejection threshold for future-event imbalance, including the
  immediately preceding non-rejecting count; mixture and invalid-input checks.
- Exact rational quantization, coarse-graining, cancellation, latent-confounding,
  ordinary-memory and collider witnesses; negative and positive quantization cells.
- Cache corruption rejected before aggregation. Unknown remote bits are counted
  once per physical row in the future test, even when compatible with both arms.

## Actual sources

Both pinned archives were independently decoded at 1,000,000 and 200,003 raw
records per chunk. Aggregation used 1,000,000 and 173,011 source rows per chunk.
The outputs agree byte-for-byte, including the external per-row cache hashes.
Every archived local setting, legacy click word and timestamp excursion agrees
with the HDF in each extraction. The 128 boundary rows are the same across all
lag/history comparisons. No interior row is deleted because of its outcome.

## Reproduction gates

The regular CI workflow runs all unit tests and regenerates the aggregate
analysis and fixed-seed recovery results. The full-source workflow downloads
hash-pinned originals and rebuilds both timing-structure aggregates with --check,
then verifies the analysis. Existing baseline workflows remain intact.
The plots are generated from committed results as a CI artifact, rather than
being a second manually maintained source of numerical values.

The simulation seed is 20260928, with 400 repetitions per scenario and exact
pointwise binomial intervals. The balanced designs were specified in the
simulation code; the additional occupancy stress panel uses only observed
marginal state sizes. It does not fit remote-setting effects or alter the
empirical test family. The null runs alone cannot certify a .005 tail; the
supermartingale arguments supply the error control.
