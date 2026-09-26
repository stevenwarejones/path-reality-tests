# Replication validation

This records a conditional reanalysis of Dryad version 3 / 450489. It reproduces plotted summaries, not the full experiment. All source files and SHA-256 digests appear below. Code: `replicate.py`; exact numerical outputs: [summary.json](replication/summary.json).

## Inputs and mapping

- All 19 originals verified before calculation; no originals modified.
- All 18 workbooks / 19 visible sheets inspected, including nonempty regions inside formatted bounds. No formula cells were found.
- Every sheet has an explicit checked entry in [column-map.json](../column-map.json); [data-dictionary.md](../data-dictionary.md) records meanings and ambiguities.
- Figures preserve deposited coordinates and units. Fig3 Q/C errors are supplied SEs; Fig4 bands are supplied SDs. E has no deposited uncertainty.

## Descriptive results

| Comparison | n | Symmetric error (%) | MAE (a.u.) | RMSE (a.u.) |
|---|---:|---:|---:|---:|
| E vs Q | 17 | 4.730009420 | 0.048313611 | 0.064498849 |
| Q vs deposited theory | 17 | 5.571439833 | 0.055107096 | 0.067283643 |
| E vs C | 17 | 12.308280654 | 0.121669368 | 0.140730057 |
| C vs deposited theory | 17 | 4.185665074 | 0.043230722 | 0.051874336 |

Q–theory uses coordinate-based linear interpolation of the deposited 850-point curve at all 17 observed positions. It is a separate comparison from E–Q. No scale is fitted or endpoint omitted. The paper quotes 4.45%; the table-based value above does not recover it. See [method comparison](../method-comparison.md) for tested explanations and why this does not by itself overturn the E–Q/E–C ordering.

## Sensitivity checks

| Alignment | Q–theory error, deposited units (%) | Unit-sum shapes (%) |
|---|---:|---:|
| linear | 5.571439833 | 5.338613730 |
| nearest | 5.574578913 | 5.338224838 |
| local_cubic | 5.571439596 | 5.338617878 |

Figure 4 contains 72 length-group rows and 100 action-group rows. Its phase RMSE is 0.021814792π; shortest-arc RMSE is 0.021814792π. The mean supplied phase SD is 0.068133595π. This is not a significance ratio or a reconstructed phase covariance.

The Fig4 summaries match their varying-endpoint FigS4 counterparts exactly (maximum absolute difference zero). FigS3B’s simulation-labelled phase also equals Fig4C’s theory-labelled phase; no labels were changed. Histogram row/count details are retained in the numerical summary, including the FigS2C tail. Those histograms are not inputs to the Figure 3/4 residual calculations.

## What remains unavailable

The inspected tables do not contain the complete per-repeat propagator tensor, endpoint repeats, or acquisition/calibration/phase/binning code. The reported path-level error and path-space fidelities cannot be reproduced from group means, SDs and one sample propagator. No probability-distribution overlap is substituted for path fidelity. No fabricated repeats, significance levels, confidence intervals or author random seeds are used. [Open questions](../open-questions.md) gives closure criteria.

## Software and visual validation

Twenty-four unit tests pass locally. Checks include independent raw-OOXML/40-digit Decimal Figure 3 calculations and hand-calculated residuals, irregular coordinate alignment, sparse blocks, forbidden extrapolation, phase branch cuts, normalization and nonlinear averaging counterexamples. The main workflow checks all inputs again after processing, regenerates the numerical summary/residual tables into temporary storage and checks reviewed snapshots. Original and generated Figures 3/4 were inspected at readable resolution; no clipped data or labels were found. CI status for a submitted revision is recorded by its Actions run, not inferred from a local run.

The committed derived outputs are two replots, two residual CSVs and a metadata/numerical summary. Each PNG includes source filename/digest in visible footnotes and metadata; each CSV starts with a source/digest comment; the summary identifies all inputs. No full cell extracts or paper images are included.

## Input provenance

| Original file | SHA-256 |
|---|---|
| Fig3.xlsx | `dccdde7960394bca2df8a28cd667f4ab54a39d15ce8ef9dce997bbb962ebe092` |
| Fig4A.xlsx | `5897f9f31b51fe387c32e329581066d4276a6d27df7c9012435cb21e85469d82` |
| Fig4B.xlsx | `da16adb40e929cc2d4a0aea59f0feeafe37a6d19c6e220f051315902b332203a` |
| Fig4C.xlsx | `139739e91ec5f4d9a2aafaabb7f65e8658f9c9b73ebb8bf5b3bfc08ced9d6f6b` |
| FigS1B-L.xlsx | `bc1f311ac63766cec60f852ca2e3893e6e719dee98a0a58dac6e3f7ab1a82da6` |
| FigS1B-Minus.xlsx | `a8a2e4eafe8412929d573a10ccf6294789e603f578b447c1ece053d4159cb234` |
| FigS1B-Plus.xlsx | `3edf88999c89e89572e9ade8ccc65959def984a516f9bb4d4e086722bcccb285` |
| FigS1B-PSI.xlsx | `c11590c741d5eb8b8d7cc62fd7c39084fbde31f5bdddd8bc2bc70e97c50fbef6` |
| FigS1B-R.xlsx | `9e16dc59ae14c4ded6e3a14c6b87ee5186f8e4390016190dffb2302a285fb4f7` |
| FigS1C.xlsx | `5d4ca2ba8f27aa0fefbe45b8f467bc0f63469b82e76fdef1588d6e9e44a412c1` |
| FigS2A.xlsx | `b7348b277afba3280d024e4ca9945b6b69a46333adeee909ae4807d4e01781d8` |
| FigS2B.xlsx | `482e635f6644c9d58c13989c395cf8e74ae7ffcbd88704196d02986bf2fe300e` |
| FigS2C.xlsx | `7412bc1f0977a280b4a7028671cb0f307f961f02c8ab991ba1c250952c1e1806` |
| FigS2D.xlsx | `a60e0ec1f82db52734801c206d136b959bb37ad7067eed6631c626d10cd71ab4` |
| FigS3A.xlsx | `d922b67d7fc2e700f44e5ac551e15ae6e73d05316d3e578dfcc9d7274f464393` |
| FigS3B.xlsx | `e042db4a8a4fbaa632d65a6680da913f24e936b57e44526b8bed39d9f5bb944d` |
| FigS4A.xlsx | `b42c2adf09d162d7967857f8a4a45fa8c95b92877757fef083ab18c0b97f3d5c` |
| FigS4B.xlsx | `49edd6c8faeba3f6ab65cd5523d21c363a6d08b11e94e556b707b45e3ea51b9f` |
| README.md | `e523bc79c4d882ea20c27c189fe5cea73c604fc8166a7737603cd34bd890875b` |
