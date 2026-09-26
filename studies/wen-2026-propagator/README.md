# Wen 2026: propagator and path reconstruction

Independent reproducible checks of Wen et al., *Direct experimental test of
Feynman's path integral postulates with single photons*, Science Advances 12,
eaeh1011 (2026), [paper DOI](https://doi.org/10.1126/sciadv.aeh1011).
This study reproduces the deposited Figure 3/4 summaries and separates observable
comparisons from claims about individual paths. It does not presume that the
paper is wrong.

## Getting the data

The 19 verified CC0 originals are included unchanged in `raw/dataset/`.
[Download and verification instructions](raw/README.md) identify
[Dryad DOI 10.5061/dryad.x0k6djj14](https://doi.org/10.5061/dryad.x0k6djj14),
version 3 / 450489, the authoritative source. All 19 files must match
[manifest.json](manifest.json); never change a digest to accommodate a local file.
The dataset's own README stays separate from repository instructions.

## Reproduce

From this folder, using Python 3.12:

```sh
python -m pip install -r requirements.txt
python fetch_data.py --verify-only
python extract_workbooks.py
python replicate.py
python synthetic_checks.py --check
```

The extractor writes complete cell inventories to ignored `results/workbooks/`.
Those extracts are not committed or uploaded as CI artifacts. `replicate.py`
reads the verified workbooks directly using [column-map.json](column-map.json)
and writes only derived figures, residual tables and a numerical summary.
Paths resolve relative to the scripts, so invocation from another directory
works. `--directory` chooses another copy of the same originals;
`--output-dir` chooses a separate destination for derived files.

To check reviewed results without replacing them:

```sh
python replicate.py --check
python replicate.py --check --output-dir /tmp/wen-review
```

The second command also regenerates plots for inspection. Snapshot tolerances
describe floating-point reproducibility, not experimental uncertainty. Agreement
with the paper's quoted numbers is not a CI pass condition.

## What has been reproduced

![Figure 3 replot](results/replication/figure3-reproduced.png)

![Figure 4 replot](results/replication/figure4-reproduced.png)

The figures preserve the deposited coordinates, arbitrary units and supplied
SE/SD columns. Figure 3's endpoint acquisition E is descriptively closer to
coherent reconstruction Q than to incoherent sum C: symmetric relative errors
are 4.73% and 12.31%, respectively. That ordering persists after unit-sum
normalization. This is not a significance claim.

The separate Q-versus-theory calculation gives 5.57% using the deposited theory
curve, compared with the paper's reported 4.45%. Interpolation, simple rounding,
unit-sum normalization and alternative common error definitions were checked.
They do not reproduce that exact number. The underlying reference array,
processing order and normalization code are unavailable, so this remains a
conditional calculation and pipeline question, not an established author error.

Figure 4's phase residual RMSE is 0.021815π relative to its deposited theoretical
mean phase. Group means and SDs do not supply the full path amplitudes needed to
recompute the reported 17.4% path-level error or path-space fidelities.

- [Validation and source hashes](results/validation.md)
- [Numerical results and all-sheet inventory summary](results/replication/summary.json)
- [Figure 3 residuals](results/replication/figure3-residuals.csv)
- [Figure 4 phase residuals](results/replication/figure4-phase-residuals.csv)
- [Actual workbook mappings and ambiguities](data-dictionary.md)
- [Why results can differ and whether it matters](method-comparison.md)

## From measurements to path claims

The inference chain has distinct steps:

1. Polarization-resolved camera signals and a reference image are acquired.
2. An optical measurement model maps calibrated signal differences to complex
   propagator estimates K, with reference, normalization and phase assumptions.
3. Products of shared K entries define path amplitudes; coherent sums give Q
   and sums of squared magnitudes give the restricted incoherent comparator C.
4. A separately acquired endpoint distribution E is compared with those predictions.
5. An ontological conclusion requires additional assumptions about which models
   the experiment distinguishes.

The separate E acquisition makes the comparison useful; it is not made vacuous
merely because Q is a reconstruction. Conversely, reconstructed paths are not
independent per-photon trajectory records. The 17⁵ products reuse entries of
five 17×17 propagators. Even those entries can share calibration and source noise.
Resampling path labels cannot replace acquisition-level uncertainty propagation.

For finite matrices the path sum equals the corresponding matrix-product entry:

\[
\sum_{x_1,\ldots,x_{M-1}}\prod_{k=1}^{M}K_k(x_k,x_{k-1})
=(K_M\cdots K_1)_{x_M,x_0}.
\]

The same detector rule therefore predicts the same observations under either
description. Excluding C excludes that defined class, not every theory with
definite trajectories or contextual responses. This paraxial experiment has
ordered propagation slices and does not implement an intervention testing
superluminal or backward-time signaling.

## Uncertainty and remaining inputs

All 19 worksheets were inspected, including the five sample camera images.
There is no complete labelled per-repeat propagator tensor, endpoint repeat
dataset, acquisition ordering or calibration/phase/binning implementation in
these sheets. Processed means cannot reconstruct those records. No p-values or
confidence intervals are inferred from the number of paths or plotted bins.
The supplied Fig4 SDs must not be treated as errors on independent means.

[Physical correspondence](physical-bridge.md) details preparation, finite-window,
projection, loss and detector obligations. [Open questions](open-questions.md)
records what is closed, what is conditional, and exactly which inputs would
settle the rest. [Literature comparisons](literature-comparison.md) and the
[staged search log](literature-search-log.md) document primary sources, review
depth, adverse checks and refresh triggers. Original paper/supplement PDFs remain
outside git; [provenance](provenance.json) records their retrieval identifiers.

The existing synthetic calculations remain separate from empirical replication:
[synthetic results](results/synthetic.json) check finite sums, coordinate
sensitivity and shared-entry noise under stated assumptions. They are neither
a fit to the measured tables nor a reconstruction of the authors' random seed.

## Formal companion and follow-up

The ideal finite case study lives in `ontology-separation`:
[pinned Lean source](https://github.com/stevenwarejones/ontology-separation/blob/a4d294ddb3847134ac2bdbb509c435c2eaff7a98/OntologySeparation/Experiments/PathInterference.lean)
and [physicist-readable case study](https://github.com/stevenwarejones/ontology-separation/blob/a4d294ddb3847134ac2bdbb509c435c2eaff7a98/docs/PATH_INTERFERENCE_CASE_STUDY.md).
It certifies finite tables and access-relative model distinctions, not this
apparatus or its empirical numbers.

A phase-intervention experiment is a separate proposed study, with its own
explicit dephased null, calibration obligations and literature review. Stronger
contextuality tests would require additional operational-equivalence and
disturbance controls; neither the current tables nor the Lean example supplies
those controls.

## Automated checks

[GitHub Actions](../../.github/workflows/ci.yml) verifies the original hashes,
tests integrity failures and mapping/numerical counterexamples, inventories all
workbooks in temporary storage, checks numerical snapshots and residual tables,
and regenerates review figures. It verifies inputs again after processing.
[CI details](results/ci-validation.md) distinguish those software checks from
scientific validation. No workflow updates the data manifest or commits results.
