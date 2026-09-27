# Verification and reproducibility

The frozen protocol/constants were committed locally before fine extraction,
and published at `2ddbc65c65aa29d3c8d2b59a142ffd325ae63889` before the new
association analysis. Prior #13 coarse outputs were known. No bin boundaries,
thresholds, periods, family budgets or simulation amplitudes were selected from
the new setting-dependent outputs. Later changes implement/document the frozen
choices; they do not retroactively label this analysis preregistered.

## Measured sources

- Both receivers: all 107,109,596 overlapping local settings, stored legacy words
  and timestamp-excursion lists reconcile. Identical 128 endpoint exclusions.
- Fine histograms reproduce #13's pooled half-by-context coarse event counts,
  trusted exposures and both compatible-uncertainty budgets exactly.
- First selected event remains first in raw order. No phase-radius cut. Selected
  full-overlap events: Alice 33,608, Bob 31,534; no selected multi-event rows in
  this run. Tests nevertheless cover multiplicity and ties explicitly.
- Default 1,000,000-record and independent 200,003-record full raw extractions
  reproduce byte-identical aggregates AND external lossless-cache SHA-256s.
  No per-row source extracts or raw archives are committed.
- The retained float64 phase comes directly from the audited decoder. For Alice
  the selected phase range is approximately 0.0375..161.3775, Bob
  0.0925..161.145; these observed ranges do not choose or shrink the legal grid.
  Floating modulo near an exact period/integer boundary is left unchanged,
  so the comparison reproduces the decoder's law rather than silently replacing
  it with rational re-rounding. The enclosing [0,162) domain is fixed in advance.

## Independent tests

All 248 repository tests pass locally after merging main, including 16 new
fine-timing tests.
Analysis, recovery and four SVG figures regenerate exactly; tracked data/link
and whitespace checks pass.

The new suite checks legal support/family cardinality and all refinement maps;
first-event raw ordering, pulse limits, tied events, phase boundaries, tails and
multiplicity; a separately coded scalar row oracle for ambiguous labels, guarded
rows, no-events and block/half endpoints; and chunk sizes 1,17,1000.

Statistical checks include exact enumeration of adapted Bernoulli exponential
processes, all completions of small uncertain counts, inversion boundaries,
propensity/missingness monotonicity, batch/scalar agreement, impossible-tree
flags, and additive-tree projections compared against independent linear
programs. Exact within-category changes preserve the coarse law; sufficiently
large fine changes give positive fine/J lower bounds. #13's exact full-tag
quantizer invisibility witness remains unchanged and is invoked directly.

The committed recovery includes 216 paired empirical-background cases with
400 repetitions each and six row-model sensitivity cases. Null controls, fresh
settings with outcome memory, persistent assignments, intentional misalignment,
and outcome-dependent missingness distinguish statistical assumptions from
implementation correctness. Monte Carlo is a power/implementation diagnostic,
not a substitute for the coverage proof.

## CPU-independent serialization

The first full-source CI run verified both raw extractions and their exact
coarse reproduction, then failed byte-for-byte measured JSON regeneration.
NumPy SIMD dispatch changed the final bits of scalar `expm1` denominators
across CPU configurations. Count inversion now uses scalar `math.expm1`,
consistent with the earlier studies. The regression compares serialized CDF
and partition bounds in fresh processes with default dispatch and with
AVX512F, AVX2 and FMA3 disabled. Full measured outputs and recovery also
reproduce exactly under the disabled-feature configuration.

Regeneration changes 1,292 measured floating-point entries by at most
1.70e-21 and one recovery floating-point entry by 2.72e-20. Every integer
(including every detection count), report and figure remains unchanged.
This fixes reproducibility; it does not change the statistical procedure
or any scientific conclusion.

## Automated gates

Offline CI runs every repository test, repository data/link checks, analysis and
recovery `--check`, and deterministic regeneration of all four SVG figures.
Full-source CI downloads the pinned originals and runs the fine extractor on
both receivers with `--check`, including exact coarse reproduction, then checks
the measured analysis again. Figures are generated solely from committed outputs.
The measured source pipeline and prior study files are preserved.
