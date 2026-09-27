# Field-level sufficiency and interpretation

| Required quantity | Actual field / source | Uncertainty or missingness | Consequence |
|---|---|---|---|
| Sodium counts | 95 `*_s1.dat` second columns | 41 integer bins per scan | Observable fits enabled |
| Grating position | First column, nm | Nonuniform measured positions; workbook 15 nm is nominal | Preserve actual coordinates |
| Integration time | Workbook column C | 1 s | Rate/count distinction retained |
| Incident rate | Workbook B, kcounts/s | Rounded empirical source rate, no trial-by-trial denominator | Author binomial posterior reproduced, not certified |
| Grating powers | Workbook F:H, mW | Calibration covariance absent | Fixed optics only |
| Velocity distribution | Author code: 158 ± 9 m/s Gaussian width | Width is not uncertainty on mean; no per-scan velocity records | Author quadrature reproduced |
| Mass distribution | QMS scan 16; mass-to-charge times 2 | Doubly charged assumption; acceptance .85–1.15 times 170 kDa | Conditional interpolated mass weights |
| Cluster geometry | Author code: bulk Na density 971 kg/m³ | Homogeneous sphere assumption | Finite-size author damping reproduced |
| Phase | Author macroscopicity code: -2.48 rad | Estimated from same scans | No independent calibration claim |
| Optical functions | Author code: polarizability, absorption, waists, perfect reflection/loss | Joint nuisance likelihood absent | Model comparisons descriptive |
| PSD frequency | PSD first column, Hz | Finite Blackman FFT | No invented time series |
| Flux PSD | Second column, Phi0²/Hz | Third column = PSD/sqrt(nav) | Infer 40–80 averages; marginal gamma motivation |
| Leakage mask | Deposit says remove six central points | Exact indices not deposited | Store chosen six bins and ±1-bin sensitivity |
| Apparent Qprime | Required by PSD denominator | Not supplied in final workbook | Profile illustrative values; no exact raw reproduction |
| Intrinsic Q | Final workbook B3:B16 | Paper says resolution better than 1%; covariance absent | T/Q fit; do not substitute Q for Qprime |
| B vs T/Q | Workbook A:E | B fits and errors; correlated with source PSD | Use summary or spectra, never both as independent data |
| Stiffness | Workbook L20=0.43 N/m, M20=0.013 | Paper rounded differently | Literal workbook reference |
| Force upper limit | Workbook N29 | Published 95% construction; no joint coverage | Published envelope, not a density |
| Layer geometry | Paper supplement | Thickness 370 ± 4 nm; orientation/full assembly positions not fully numeric | Isolated-load coefficient; no full-assembly certificate |
| XENON events | `data_unbinned_1to30kev.txt` | 433 reconstructed energies | Conservative count upper bound enabled |
| XENON efficiency | `efficiency.txt` | Reconstructed-energy curve, no migration matrix | Identity-energy response is feasibility only |
| XENON background | `bkg_model.txt` | Best-fit rate, not full nuisance likelihood | No subtraction in count bound |
| Xenon charge cancellation | 2026 supplement S1/Table S1 | alpha in [1,1.5], beta=1.04; rc >= 1.15e-8 m | Explicit model sensitivity |
| LPF torque | 2025 paper Sec. III | Minimum inferred from published figure; no local confidence/covariance | Deterministic envelope only |
| LPF raw sample | NASA FITS HOUSEKEEPING | No science HDU in chosen object | Access audit only |
| Common noise rest frame | None of these releases | Apparatus trajectories/attitudes relative to hypothetical frame unspecified | Long-memory universal-field inference disabled |

Measured values, published constraints, predictions, simulations and exact
conditional certificates are tagged separately. Source-derived predictions are
not claimed to be observed collapses. The published sensor summary fit is
orthogonal; the independent baseline fit weights only B errors and reports the
difference explicitly.

The Blackman window correlates adjacent spectral estimates. A product of gamma
marginals is not asserted as an exact likelihood. Gaussian Monte Carlo in this
study validates only the separately named synthetic procedure.

## Additional motion artifacts

`results/motion.json` associates its calculation with all source and code hashes,
a hash of the baseline data and a hash of the auxiliary summary. Its `inputs`
contain derived normalized mass quadrature, optical settings, frequency/averaging
metadata and SI-converted JPL vectors. The original transport responses remain
external. LPF vectors are explicitly tagged as predicted, not reconstructed.
`data` separates the benchmark, sensor assembly, window calculations, sodium
nonlinear intervals, error budgets, one-channel/joint bounds and sampled common
frames. `results/moving-report.md` and `results/moving-diagnostics.svg` regenerate
from those rounded numerical values. None supplies a calibrated joint likelihood.
