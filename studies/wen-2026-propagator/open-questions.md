# Question register

Updated 2026-09-26 UTC. IDs preserve continuity with the earlier review.
Closed means the stated narrow question has evidence; conditional or blocked
means it does not. [Method comparison](method-comparison.md) explains why a
numerical difference need not overturn a scientific comparison.

| ID | Question | Evidence / disposition | Next action and closure criterion |
|---|---|---|---|
| Q01 | Can every deposited original be acquired and verified? | **Closed:** all 19 files match the unchanged version-3 manifest. | Reopen only on a new version or a digest mismatch; never silently revise the manifest. |
| Q02 | Is the supplement available? | **Closed:** primary supplement read; its source and hash are in provenance. | Recheck relevant sections when a revised source appears. |
| Q03 | Which quantities are acquired versus reconstructed? | **Closed at method level:** pointer images and E are acquisitions; K, path amplitudes, Q and C are reconstructions. | End-to-end calibration verification depends on Q04/Q08. |
| Q04 | Does this deposit contain the complete propagator/repeat records? | **Inventory closed; full replication blocked:** 18 books / 19 visible sheets, no complete labelled K[slice,out,in,repeat] tensor. Images are one 321×101 sample set. | Need all complex K estimates by repeat, raw frames/counts, references, exposures, gains, backgrounds, acquisition order and cross-channel covariance. Do not replace these with means. |
| Q05 | Are millions of path labels independent trials? | **Closed: no.** Products share propagator factors; the five 17×17 matrices have 1,445 entries before repeats. | Any uncertainty calculation must preserve shared factors and calibration dependence. |
| Q06 | Can the reported 4.45% be reproduced? | **Partly resolved:** it compares Q to theory, not E to Q. Deposited-curve calculation gives 5.571440%; E–Q is 4.730009%. Tested alignment/normalization/metric choices do not recover the quote. | Need original 17-point theory reference, exact normalization, averaging order and evaluation code. A rounding or plotting convention explanation is acceptable if demonstrated; do not select fits or subsets to hit the target. |
| Q07 | Which postulate-I fidelity is intended: 94.9% or 94.4%? | **Open:** introduction and Results give different values. | Need the evaluated path vectors and code/version; distinguish a textual discrepancy from a demonstrated computational error. |
| Q08 | How are a.u. scales, exposures and losses related? | **Open:** Q/C and the dense theory curves have matching numerical means on different grids. Shape-normalized and deposited-unit results preserve the E–Q/E–C ordering. | Need normalization steps and calibration linking E, K and all success/failure outcomes. Cannot infer absolute detection probabilities from a.u. labels. |
| Q09 | Is the path-vector fidelity normalization implemented as printed? | **Printed-expression question:** extra powers of self-inner-products matter for unnormalized vectors; unit-normalized vectors remove that distinction. | Need actual vector normalization and fidelity code. Do not replace path fidelity with probability-distribution overlap. |
| Q10 | Did text extraction omit mathematical structure? | **Closed for the checked square-root example:** MathML contains the radical; plain-text omission is not a paper error. | Inspect MathML or rendered equations before raising further notation questions. |
| Q11 | What are the physical center coordinates? | **Partly resolved:** data use −8…8. Reported 5.73 μm gives half-bin edges ±48.705 μm, plausibly explaining the stated ±48.72 μm range. Conditional 6.09 μm sensitivity retained, not preferred. | Need camera-to-object calibration, mode/bin boundaries and quadrature weights to close exact physical mapping. |
| Q12 | Which action bins and phase reference were used? | **Open beyond replots:** 100 display coordinates include 0 and 2; deposited theoretical phase is not simply that display coordinate. No raw memberships or phase-reference code. | Need action convention, endpoint-specific reference, wrapping, bin edges/counts and arithmetic/circular averaging; then recompute groups. Preserve stored coordinates meanwhile. |
| Q13 | Do matching path sums force literal simultaneous trajectories? | **Closed as identifiability statement:** finite path and transfer descriptions predict the same observations with the same detector. | For a new claim, specify a different observable prediction and the tested model-class restrictions. |
| Q14 | Does this test faster-than-light/backward-time intervention? | **Closed for this study: no such intervention is implemented by the ordered paraxial slices.** | Any causal test needs a separate operational protocol and relativistic apparatus correspondence. |
| Q15 | Can the noise simulation and 17.4% path-level error be reconstructed? | **Limited check complete:** independent shared-entry noise simulation exists. Group means/SDs do not determine the empirical path error or author's exact seeded run. | Need full K repeats and simulation seed/correlation/normalization code; no path-label bootstrap. |
| Q16 | Are corrections, code or clarifications available? | **Search-limited:** staged two-engine searches and fresh Crossref metadata identified no usable correction/code release. Some search results were irrelevant. Absence is not proved. | Refresh primary sources when new data, code, author clarification or a correction appears, and before stronger public claims. |
| Q17 | What does the formal companion certify? | **Closed at declared scope:** exact ideal tables, restricted incoherent exclusion, enlarged contextual response and phase-access distinction. | Empirical calibration remains outside that certificate; pinned links in README. |
| Q18 | Does the finite intermediate sum model the endpoint apparatus? | **Open physical correspondence:** projections/quadrature and unfiltered propagation need not coincide. | Need actual apertures, preparation and mode profiles, retained losses, omitted-amplitude and discretization bounds. See physical-bridge.md. |
| Q20 | Are repository checks complete? | Local integrity, mapping, residual and synthetic checks are implemented; current run status belongs to the PR's Actions checks. | Require the full workflow to pass for the submitted head; a green run is not empirical proof. |
| Q21–Q24 | What follow-up can discriminate declared models, and with what calibration, contextuality assumptions and novelty? | **Separate study:** phase-intervention protocol and supporting work are preserved for their own review. | Complete apparatus feasibility, nuisance budget, operational controls, focused literature search and finite-sample design before execution. No experimental success or novelty is claimed. |
| Q25 | How should the S2 histogram tail be interpreted? | **Recorded, not repaired:** S2A totals 1,445; S2C first 99 rows total 1,445, with one extra tail count at R104. Edge/center conventions could explain 119/99 row counts, not by themselves the nonmonotone tail. | Need histogram-generation code or source clarification to identify intended rows. Tail does not enter the dedicated Fig3/4 replot. |
| Q26 | Are length labels individual lengths or bin display positions? | **Open:** varying-endpoint labels 0…71; fixed labels 0,2,…,62. Last-bin inclusion is a possible explanation. | Need grouping rules/counts before path-level inference. Do not reinterpret SDs or drop endpoints. |
| Q27 | What explains the S3B phase-column correspondence? | **Recorded:** its simulation-labelled phase equals Fig4C's theory-labelled phase. Labels retained as deposited. | Need the figure-generation arrays/code to determine intent; no automatic relabelling. |

## Diligence and stopping criteria

Pursue every new question to either a reproducible answer or a documented
dependency with a concrete closure criterion. Do not quietly drop a discrepancy,
upgrade a plausible explanation into a fact, or call an externally blocked item
complete. Record adverse evidence and update dependent prose, calculations and
tests together. Searches must be repeated at source acquisition, after numerical
ambiguities emerge, and before any stronger conclusion or experimental design.

The remaining external questions require specific data/calibration/code rather
than more numerical fitting to the same summaries. No author or Dryad contact
is authorized; Steven decides whether to request those materials. A lack of
those inputs is a reproducibility boundary, not evidence of misconduct.
