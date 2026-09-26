# Data and tooling validation

The included deposit is version 3 / 450489 of Dryad DOI
10.5061/dryad.x0k6djj14. All 19 original files match the unchanged byte counts and
SHA-256 digests in `../manifest.json`. The source copies are CC0; paper and
supplement PDFs and complete cell inventories are excluded from git.

## Checks implemented

- Verification fails on missing files, wrong sizes, same-size corruption or
  unlisted dataset entries. Verification-only mode performs no download or repair.
- Mocked access-control and invalid-HTML responses stop subsequent download
  requests and never become saved originals. CI does not depend on Dryad access.
- A generated test workbook checks hidden-sheet inclusion, cell coordinates,
  formula text and a distinct known cached result. These are synthetic fixtures,
  not invented experimental observations.
- All requested workbooks must pass integrity checks before any extract is
  written. A missing later input leaves no partial output directory.
- The numerical tests check complex path/matrix composition, coherent versus
  incoherent probabilities, zero-valued SMAPE cases, and rejection of changed,
  truncated or non-finite numerical snapshots.
- The full synthetic calculation checks the committed JSON with relative
  tolerance 1e-10 and absolute tolerance 1e-12 for floats; integer counts, strings,
  keys and array lengths must match exactly. These tolerances accommodate floating
  arithmetic, not empirical uncertainty. Material intentional changes require
  reviewed updates to both code and snapshot.
- The tracked-data allowlist permits exactly the manifest-listed originals plus
  `raw/README.md`. Local Markdown links are checked; PDFs, ZIPs, unapproved
  workbooks and complete extracted inventories are rejected.

The 14 behavioral tests pass locally under Python 3.12 with the study's pinned
dependencies. Verification has been exercised from a different working directory.
GitHub Actions repeats those tests, verifies all 19 files, reads all 18 workbooks,
and regenerates the synthetic JSON and plot into a temporary artifact directory.
The plot is a review artifact; CI does not compare PNG bytes. Extracted inventories
are temporary and are not uploaded as artifacts.

The workflow uses read-only repository permissions and pinned action commits.
Neither the scientific scripts nor the workflow updates the manifest or commits
generated files automatically. Inputs are verified again after processing and
`git diff --exit-code` checks that tracked contents remain unchanged.

## Scope

This validates integrity, extraction behavior and declared synthetic calculations.
It does not validate the real scientific column mappings, reproduce Figures 3/4,
assess empirical agreement, or establish the completeness of per-repeat records.
Those are separate replication obligations. A disagreement with the paper is not
itself a failing CI test. The phase-intervention study and its tests remain in
their own later PR.
