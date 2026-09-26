# The deposited Wen example: recoverable signals and missing calibration

Sources: Wen et al., [main article](https://doi.org/10.1126/sciadv.aeh1011),
Eq. (20)/(23); supplementary Section 6 (PDF page 3) and Fig. S1 (page 9).
The primary XML and supplement were retrieved on 2026-09-26 and their SHA-256
hashes agree with the existing study's provenance:

- XML: a89dd991b2af69fc04b72171c69487adc758485f03fb864abe9c7fc156169693
- Supplement: 47f4a063550529c47979c216fb0b4162d051b38c6f9c1f78f147772d908cd567

The source PDFs are not copied into the repository. The implementation reads the
six already deposited workbooks and verifies their existing manifest hashes.

## Executed pipeline

1. Read all five 321×101 images with the deposited x/y axes. Retain grayscale as
   grayscale: exposure, gain, trigger counts and dark calibration are unavailable.
2. Select the exact x coordinates in the 17-point FigS1C slice (−8,…,8). There is
   no coordinate fitting, interpolation, row deletion, clipping or sign choice.
3. Declare the spatial reduction. We report all-y averaging, the y=0 row, and
   the strip |y|≤1. The first is the displayed baseline, not a recovered author ROI.
4. Under equal channel gains and zero background, form Sx=I_plus−I_minus and
   Sy=I_R−I_L. The R−L sign is given by supplementary Section 6.
5. Compute C=(Sx+i Sy)/sqrt(I_PSI). This follows the intensity-based pointer
   bilinear formula. It is **not** the normalized Stokes expression obtained by
   additionally dividing each difference by its pair sum.
6. Retain the unknown common normalization and reference phase explicitly.
   Eq. (20) contains the complex reference ψ*, whereas Eq. (23) replaces its
   spatial phase by the paper's near-collimation approximation. We cannot
   independently validate that phase approximation from intensity-only images.
7. Compare C with the deposited complex slice after projecting onto one common
   complex scale c=<C,K>/<C,C>. This fit is explicitly a descriptive shape
   diagnostic, not a measurement of the missing normalization or phase.

Equal-weight intensity averaging before the square root assumes a compatible
separable y profile if interpreted as a one-dimensional propagator. The deposit
alone does not establish that reduction. All alternative ROIs are retained;
we do not select the one giving the smallest residual.

The relative L2 residual is about 10.31% for all-y averaging, 36.82% for the
central row, and 10.75% for the narrow strip. These compare one example image
set with a reported repeated-measurement mean; the original repeat pairing is
missing. They are not the paper's Figure 3 percentage error, a fidelity estimate,
a hypothesis test, or evidence of experimental inconsistency.

`results.json` contains every input reduced signal, reconstructed complex value,
pair-total discrepancy, scale fit, and the deposited real/imaginary SEs. The SEs
are preserved and not used to manufacture an independent-point likelihood.

## Conditional sensitivity, not confidence

The declared sensitivity cases vary each channel's common gain over 1±1% or
1±5% and its background over 0…1 grayscale unit. These are **chosen sensitivity
ranges**, not estimated detector errors. `conditional_box` propagates them to
rigorous coordinate enclosures before common scale and spatial reference phase.
The boxes are outer bounds: shared calibration parameters mean their endpoints
need not be jointly attainable. Random shared-parameter tests independently
check enclosure, while the analytic endpoint argument supplies the guarantee.

If the reference denominator can vanish, the pipeline reports the failure rather
than clipping it. Without a bounded global reference normalization, absolute
propagator scale is unbounded. Without a spatial phase calibration, relative
phase remains unidentified. Neither omission is repaired by fitting FigS1C.

## Explicit indistinguishable completions

For any spatial phase θ(x), replace ψ(x) and K(x) by exp(iθ(x)) times themselves.
Then ψ* K, |ψ|², |K|² and every local polarization intensity are unchanged.
This also corresponds to left-multiplying an output operator by a diagonal
unitary, so it is not merely a numerical re-pairing. Requiring an independently
known free-propagation operator can rule out such changes—but that is an extra
model constraint, not information in these local images alone.

The output contains two concrete completions θ(x)=±πx/16 with identical existing
signals and different propagator phases. They do not purport to reproduce the
unavailable whole apparatus history. Repeating the same local pointer/intensity
measurements cannot distinguish them. A phase-sensitive measurement relative to
an independently calibrated reference can, under the restricted hypotheses in
[measurement-design.md](measurement-design.md).
