# Candidate combinations and selection

The named datasets were treated as starting points. Four combinations were
evaluated; two have executable source-backed feasibility calculations.

| Combination | Shared modification and responses | Added information | Access / nuisance risk | Decision |
|---|---|---|---|---|
| Sodium + layered sensor | Gaussian mass-density coupling; Talbot-Lau coherence and force noise | Direct superposition survival versus induced diffusion | Both originals downloaded; optics/phase calibration and colored moving-path map incomplete | Preserve complete baselines and numerical white-response diagnostic; do not invent joint confidence |
| Layered sensor + LPF torque | Spatial Gaussian mass coupling; kHz load response and mHz cube torque | Separates parts of a declared temporal spectrum over ~six frequency decades | Sensor original plus published torque envelope; full science product not reconstructed; noise-frame motion matters | Selected **conditional envelope feasibility** calculation |
| Layered sensor + XENONnT ER | Same mass coupling; mechanical diffusion plus spontaneous radiation with charge cancellation | Tests whether the strong white radiation constraint remains useful with spectral freedom | 433 events downloaded, efficiency/background audited; incomplete migration/2026 energy range | Executable runner-up, no-background count bound and approximate response |
| LPF translation + rotation + outgassing/housekeeping | Force/torque geometry and gas versus CSL noise | Potential cancellation of common spectral amplitude at matched frequency | Ratio idea already in 2025 paper; unmeasured pressure and electrode noise covariance limit subtraction; NASA sample housekeeping only | Reject new detection/exclusion claim; no arbitrary cross-epoch ratio |

## Quantitative comparison

The generated [report](results/report.md) includes separate A, separate B and joint
bounds for the selected candidate. Removing either source is exactly the relevant
single-source column. At white noise the layer load wins near rc=100 nm. At
correlation time 1 s the space envelope becomes much stronger. This is
complementary frequency response, not added statistics.

The radiation runner-up calculates a Poisson upper bound on signal counts without
subtracting any background, then integrates the full xenon charge-cancellation
factor over the released efficiency. This uses an explicitly approximate
identity energy-response map. It is compared at alpha=1 and 1.5 and under a 20%
efficiency loss. It cannot reproduce the 2026 official 1–140 keV analysis.

All frequency comparisons use angular frequency internally: omega=2 pi f for
mechanics, omega=E/hbar for photons. A keV photon samples a frequency many orders
of magnitude above either mechanical instrument.

## What survives nuisance freedom?

- Multiplying the selected spectral envelopes by 0.5 or 2 does not remove the
  asymptotic long-memory escape. These are stress factors, not confidence coverage.
- Changing layer thickness over 366–374 nm is calculated separately. Full assembly
  geometry and cross-correlations are still required for a physical certificate.
- Radiation sensitivity to the atomic distance approximation is computed; missing
  detector migration is a more fundamental obstacle to a certified recast.
- The stationary-body factorization itself is a **material assumption**, not a
  nuisance that can be ignored. The two experiments cannot both be assumed at
  rest in an arbitrary common colored-noise frame. The moving-overlap diagnostic
  shows how strongly this can change temporal sensitivity.
- The two spatial benchmarks are illustrative and chosen before the final
  validation run. They specify geometry, density, separation, time and suppression.
  They are not a universal definition of definite outcomes.

## Why the selected calculation is bounded in scope

A new exclusion is not justified by the available audited response/calibration
records. The completed selected analysis is therefore a deterministic envelope
feasibility study, including explicit survival sequences; it is not labeled a
full joint empirical fit. The original pairing is not declared infeasible in
principle, and the uncalibrated alternatives are not used to claim existing data
can never answer the question.

The most consequential next analysis input is a validated moving-body,
noise-frame-dependent colored response, together with calibrated finite-duration
instrument kernels. For the original experiment this means the colored
Talbot-Lau map. For a new measurement, a slow or held massive superposition with
several holding times and independent optical/noise calibration would directly
test the zero-frequency dephasing that finite-band PSD summaries leave open.
