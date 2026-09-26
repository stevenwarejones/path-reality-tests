# What the numerical differences do—and do not—change

This is a conditional calculation from the deposited figure summaries. The
reproducible numbers are in [summary.json](results/replication/summary.json),
with input filenames and SHA-256 digests. The paper is
[Wen et al. (2026)](https://doi.org/10.1126/sciadv.aeh1011); the dataset is
[Dryad version 3](https://doi.org/10.5061/dryad.x0k6djj14).
No statement here establishes an author error or an invalid experiment.

## Two comparisons with different purposes

E is the separately acquired endpoint distribution; Q is computed coherently
from the inferred propagators; C discards the cross terms of the same path
decomposition. Comparing E with Q and C probes agreement of those displayed
predictions with the endpoint acquisition. Q versus theoretical Q instead checks
agreement of the reconstruction with its theoretical reference. The paper's
4.45% refers to the latter comparison. Neither is the path-space fidelity.

Our primary calculation keeps all 17 endpoints and the deposited a.u. scales.
Theory is interpolated at the experimental coordinates, because only the two
outer coordinates are exact matches on its 850-point grid. We use the paper's
symmetric relative error, `100 mean(2 |a−b| / (|a|+|b|))`; a pair of zeros
contributes zero. No fitted gain, translation, endpoint exclusion, or phase
offset enters the primary comparison.

| Comparison | Symmetric relative error | RMSE (a.u.) |
|---|---:|---:|
| E versus Q | 4.730009% | 0.064499 |
| E versus C | 12.308281% | 0.140730 |
| Q versus deposited theory, linear interpolation | 5.571440% | 0.067284 |
| C versus deposited theory | 4.185665% | 0.051874 |

E is closer to Q than to C at 13 of the 17 coordinates by absolute difference.
This is descriptive, not a binomial test: coordinates share reconstruction and
calibration errors. E has no deposited uncertainty column. Q's SE alone is
insufficient to assign a significance to E−Q, and the many path products are
not independent observations.

## Why might 5.57% differ from 4.45%?

The difference is about **1.12 percentage points**. Its relative size is about
25% of the quoted error, but that does not mean a 25% change in the predicted
probabilities. We checked plausible processing choices before interpreting it.

| Explanation | Check and outcome | What remains possible |
|---|---|---|
| Alignment of a dense plotting curve | Linear 5.571440%; nearest point 5.574579%; local four-point cubic 5.571440%. | These ordinary interpolation choices cannot explain the difference. An underlying analytic/reference array not identical to the plotted curve could. |
| Rounded plotting values | Under a piecewise-linear slope bound and 5-decimal x / 9-decimal y rounding, the computed change in the metric is below 0.00005 percentage points. | This rules out simple rounding of these displayed numbers as an explanation of 1.12 points, not an unknown upstream transformation. |
| Normalization / plotting grid | Unit-sum normalization of Q and the interpolated theory gives 5.338614%. | It changes the comparison, but does not recover 4.45%. The exact normalization pipeline is absent. |
| Different definition of “MAPE” | Ordinary theory-denominator MAPE gives 5.644976%; Q-denominator 5.519757%; pooled symmetric error 5.395314%. | None reproduces the quote; the stated symmetric definition remains the primary calculation. We do not choose a metric to force agreement. |
| Computing an error before versus after averaging | The operations do not commute; a small exact counterexample is tested in `tests/test_replication.py`. Path multiplication and squaring introduce further nonlinearities. | The ten individual records and calculation order are unavailable. This is a possible explanation, not a demonstrated one or a correction to the paper. |
| A different data/code revision or theory model | Searched the DOI/title and repository records for released code, clarification and corrections; no usable additional implementation identified. | The actual 17-point theoretical reference, normalization, and author evaluation code would settle this more directly than further arbitrary fits. |

There is a concrete normalization clue. Q and C over their 17 deposited points,
and each theoretical curve over its own 850 points, have mean approximately
1.030184403. The 17 interpolated theory-Q values instead have mean
1.012591817. Equal means on different grids suggest a plotting-normalization
step, but the table alone does not identify how it was done. The deposit's a.u.
values must not be presented as calibrated absolute detection probabilities.

An additional, explicitly conditional free finite-grid calculation uses 17
centers, five intervals, wavelength 795 nm and interval 15 mm. It compares
unit-sum shapes for 5.73 μm and 6.09 μm spacings; it does not fit an apparatus
model. Its source and results are retained in `replicate.py` and `summary.json`.
Neither grid replaces the deposited theory in the primary comparison. Finite
windowing, mode shapes and quadrature weights need their own physical evidence.

## Does the difference matter scientifically?

**For exact numerical reproducibility:** yes. We cannot claim to have reproduced
the printed 4.45%. The appropriate status is an unresolved calculation/pipeline
question, with tested explanations and a specific missing-input list.

**For the displayed E–Q versus E–C comparison:** the tested choices preserve the
ordering. With each series normalized to unit sum, the errors are 4.491211%
and 12.463748% respectively, compared with 4.730009% and 12.308281% in deposited
units. This supports a limited descriptive statement of closer agreement with
Q. It neither provides a significance level nor certifies the physical meaning
of every reconstruction assumption.

**For path reality:** neither 4.45% nor 5.57% distinguishes mathematically
equivalent path-sum and transfer descriptions. Such a separation requires a
specified alternative with different observable predictions and experimentally
checked assumptions. A narrower incoherent null is not all definite-trajectory
theories. The finite formal companion makes this distinction explicit.

## Other apparent differences need proportionate treatment

- **Spatial range:** the stored centers are labelled −8…8. With 5.73 μm spacing,
  outer centers lie at ±45.84 μm and half-bin edges at ±48.705 μm, close to the
  quoted ±48.72 μm. This is a plausible center/edge explanation, not a reason
  to replace 5.73 by 6.09 μm. Exact calibration remains unavailable.
- **Histogram counts and bin labels:** 120 edges give 119 bins, and 100 edges
  give 99 bins. That could explain some documentation/count differences. FigS2C
  also has a nonmonotone tail with one additional count; the code records it
  without deletion. Whether that tail was used by the authors is unknown. It
  does not enter our Fig3/4 calculations, which use their dedicated tables.
- **Figure 4 phase:** the deposited mean-phase minus theory RMSE is 0.021815π;
  shortest-arc and ordinary residuals agree here. The mean reported SD is
  0.068134π. Those SDs are not errors on independent means, so their ratio is
  not a significance test. Display labels, action bin membership, phase reference
  and arithmetic versus circular averaging remain distinct questions.
- **Reported fidelity:** the introduction gives 94.9% for postulate I while the
  Results give 94.4%. Also, the printed path-vector normalization raises a
  printed-expression question when read for unnormalized vectors. Unit-normalized
  vectors remove that particular distinction. Neither observation determines
  the number the authors' code actually computed; the full complex amplitudes
  are needed. We do not substitute a classical overlap of plotted probabilities
  for their path-space fidelity.

The 17.4% path-level error likewise cannot be recovered from means and SDs of
groups. Unweighted averages of the group means would change the estimand and
cannot be offered as a replication. Preserving that boundary matters more than
turning every formatting discrepancy into a scientific objection.
