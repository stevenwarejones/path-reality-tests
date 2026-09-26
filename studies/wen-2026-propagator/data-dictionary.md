# Verified workbook dictionary

Inspected all 18 workbooks and 19 worksheets in Dryad version 3 / 450489,
against the unchanged [manifest](manifest.json). All sheets are visible; no
formula cells were found. Row 1 contains headings or an image coordinate axis.
Ranges below are inclusive Excel coordinates. The machine-readable
[column map](column-map.json) checks sheet sets, dimensions, headers and numeric
blocks before analysis. Formatting dimensions are not data dimensions.

| File / sheet | Actual columns and rows | Meaning and limits |
|---|---|---|
| Fig3 / Sheet1 | A2:C851: 850 theory points, x, Q, C. D2:I18: 17 points, x, E, Q, SE(Q), C, SE(C). | x in δx units; ordinates in a.u. Match curves by x, not row number. No E error column. SEs are described as ten-measurement errors by the main caption. |
| Fig4A / Sheet1 | A2:C73: 72 rows, length, mean weight, SD. | Deposited length labels 0…71 in δx units. Preserve labels; exact grouping/edge convention is unavailable. |
| Fig4B / Sheet1 | A2:C101: 100 rows, action coordinate, mean weight, SD. | Action labels span 0…2, inclusive, in the paper's πℏ units. They are plotting coordinates, not a deposited array of bin boundaries. |
| Fig4C / Sheet1 | A2:D101: action coordinate, measured mean phase, SD, theoretical mean phase. | Phase in π units. Same action coordinates as Fig4B. No phase offset or unwrapping is fitted. |
| FigS1B-L, Minus, Plus, R, PSI / Sheet1 | A2:A322: 321 row coordinates; B1:CX1: 101 column coordinates; B2:CX322: grayscale. | A1 is empty. Worksheet formatting extends to 322 columns, but columns after CX contain no values. Actual images are **321 × 101**, not 321 × 321. First-column coordinates −40…40; first-row coordinates −12.25…12.75, both step 0.25 in δx units. README assigns these to x and y respectively. These are one sample image set, not a complete acquisition tensor. |
| FigS1C / Sheet1 | A2:C341: 340 theory points, x, Im(K), Re(K). D2:H18: 17 observed points, x, Im(K), SE(Im), Re(K), SE(Re). | Sample slice 2→3, input position 4δx, per supplement. Do not swap real/imaginary columns. |
| FigS2A / Sheet1 | A2:A120: **119** phase coordinates. B1:R1: displacements 0…16. B2:R120: histogram counts. | 1,445 total counts. Source description says 120 bins; 120 edges giving 119 centers is a possible explanation, not confirmed processing history. Histogram populations are not measurement repeats. |
| FigS2B / Sheet1 | A2:C18: 17 displacement, phase mean, SD values. D2:E201: 200 theory displacement/phase values. | Phase in π units. The caption's residual-language and the nonzero displacement-dependent theory curve require care; do not subtract theory twice. |
| FigS2C / Sheet1 | A2:R120: 119 numeric histogram rows. Rows 2:100 contain 99 increasing amplitude coordinates, 0.50505…1.49495. Rows 101:120 return to 0.80336…0.99496. | Preserve both segments. First segment totals 1,445 counts; the tail contains one nonzero count at R104, giving 1,446 across all rows. No rows are deleted or repaired. |
| FigS2D / Sheet1 | A2:C18: 17 displacement, amplitude mean, SD values. D2:E201: 200 theory displacement/amplitude values. | Amplitude divided by the reported mean amplitude. These aggregates cannot restore each propagator entry. |
| FigS3A / Sheet1 | A2:C73: 72 length, simulated mean weight, SD rows. | Simulation, not experimental observations. |
| FigS3B / Sheet1 | A2:F101: action, simulated mean weight, SD, phase labelled simulation, phase SD, phase labelled theory. | Retain the literal labels. D equals Fig4C's theoretical phase column; this numerical equality does not establish which source labels or arrays should change. |
| FigS4A / Sheet1 | A2:C73: 72 varying-endpoint length/mean/SD rows. **D2:F33: 32 fixed-endpoint length/mean/SD rows.** | Fixed-endpoint x has its own column D, labels 0,2,…,62. The author README omits that second x column in its variable list. |
| FigS4B / Sheet1 | A2:F101: action, varying-endpoint weight, SD(weight), phase, SD(phase), theory phase. | Exactly reproduces the corresponding Fig4B/C columns. |
| FigS4B / Sheet2 | A2:F101: same fields for fixed endpoint. | Distinct selection; never merge its rows with Sheet1 as independent repeats. |

`.xlsx` is understood for each workbook name. Every mapping is based on the
verified file's actual headers plus the author README/main/supplement captions.
The table records ambiguities instead of rewriting the originals.

## Units, bins and uncertainty

The action headers literally use forms such as `S(pi/\\hbar)`. We label the
replots as action coordinates in πℏ units, consistent with the paper's range
0…2πℏ, while preserving all numerical x values. The 100 stored labels include
both endpoints; they are neither the left edges 0,.02,…,1.98 nor the centers
.01,.03,…,1.99 of the stated 100 equal intervals. They may simply be display
positions. Rebinning requires path actions, membership/counts and author code.

Length labels can also denote bins rather than individual integer lengths.
Endpoint inclusion in the final histogram bin is a plausible reason the stored
maximum label differs from a finite-grid maximum. This is an unresolved mapping
question, not evidence that paths were discarded.

SD columns are not converted to SEs. Ten repeats are reported for Fig3, but the
individual records and correlations are not deposited in these sheets. Fig4's
SDs describe spread, not independent errors on each mean. No effective sample
size is inferred from the number of path labels or bins.

## Availability conclusion

The inspected sheets contain plotting summaries, histograms, one example image
set and one example propagator slice. They do **not** supply a labelled complete
complex tensor `K[slice, output, input, repeat]`, endpoint-image repeats, or the
calibration/phase/binning code. There are no additional hidden worksheets.
We can replot the supplied means and compute descriptive differences; we cannot
recreate every path amplitude, the reported path fidelities, or covariance-aware
uncertainty from these tables. This conclusion concerns this pinned deposit,
not all data the authors may hold elsewhere.
