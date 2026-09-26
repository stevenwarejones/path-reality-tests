# Data dictionary — source-described schema, not yet inspected workbooks

These roles come from the author README rendered on Dryad and the main/supplement captions. Exact sheet names, row offsets, formula caches, missing values and normalization must be verified from originals. Do not treat this dictionary as an extraction result.

| File family | Variables / units | Origin and uncertainty | Permitted analysis once acquired |
|---|---|---|---|
| Fig3.xlsx | Theory and experiment x_f/δx; theoretical Q,C; observed E; reconstructed Q,C in a.u.; SE(Q), SE(C) | E is endpoint camera acquisition; Q,C derive from shared K. Caption says SE of ten measurements; no E uncertainty field listed. | Separate E–Q, E–C and Q–Q_theory residuals; document scale choices. No covariance-aware significance without covariance/trials. |
| Fig4A.xlsx | L/δx, mean path weight, SD | L=Σ|Δx| is transverse discrete variation, not full optical distance; dispersion across grouped paths. | Replot mean/SD. Group means cannot recover per-path MAPE. |
| Fig4B.xlsx | Action bins, mean path weight, SD | 100 action intervals; Record the source action-unit convention when mapping columns. Shared K factors induce dependence. | Replot and inspect counts/wrapping. No treating bins as independent shots. |
| Fig4C.xlsx | Action bins, experimental mean phase/π, SD, theoretical phase/π | Phase/action summary; reference, wrapping and circular versus arithmetic averaging require clarification. | Circular residuals with documented convention; no per-path fidelity from bin means. |
| FigS1B-L/Minus/Plus/R/PSI.xlsx | First row/column: transverse coordinates in δx; matrix: grayscale | One illustrated propagator input, including reference image. Intensities not automatically raw photon counts. | Check Eq.8 reconstruction only if gains/exposures/backgrounds and reference scale are available. |
| FigS1C.xlsx | x/δx; real/imaginary K; theory and experiment; component SE | Sample propagator for x_input=4δx, slices2→3. | Verify sample reconstruction and residuals; does not establish full K tensor. |
| FigS2A/C.xlsx | Phase/amplitude histograms by displacement | Aggregate frequencies; A uses 120 bins, C uses 100; Record the exact phase-residual definition when mapping columns. | Recompute aggregate statistics only with exact bin edges/counts; do not infer trial-level covariance. |
| FigS2B/D.xlsx | Displacement; means, SD; theory | Summary of K phase and magnitude variability. | Check noise parameter consistency, not independent data for another significance calculation. |
| FigS3A/B.xlsx | Simulated length/action/phase means and SD | Numerical noise model; not observations. | Compare with explicitly seeded synthetic runs; no statistical power inferred from a single plotted simulation. |
| FigS4A/B.xlsx | Fixed versus variable endpoint versions of path summaries | B reportedly has Sheet1 variable and Sheet2 fixed endpoint. Exact column labels await workbook inspection. | Resolve headers against supplement before parsing; preserve both selections. |
| README.md | Author metadata and column descriptions | Documentation; CC0 deposit | Check the pinned hash and compare descriptions with actual sheets. |

This preliminary dictionary records source descriptions only. Actual sheet names, column mappings and available repeat-level records will be documented during replication.
