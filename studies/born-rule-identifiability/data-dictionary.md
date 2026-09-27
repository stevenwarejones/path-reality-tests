# HDF5 audit and processing map

| File | PD shape (V) | Shutter shape | Time shape | Temperature |
|---|---|---|---|---|
| 23degrees measurement.h5 | (3360,40,1) | (3360,1,1) | (3360,40,6) | Absent |
| 30degrees measurement.h5 | (3600,50,1) | (3600,1,1) | (3600,50,6) | (3600,1,1), housing voltage; one NaN at zero-based row 1166 |
| interference contrast_21to35.h5 | (400,5,1) | (400,1,1) | (400,5,6) | (100,1,1), attribute V but values in thousands; unresolved possible resistance scale |

Numerical arrays are float64; masks are validated integer-valued. Main cycles
are eight consecutive entries, a permutation of 0–7. Contrast cycles are four
entries with codes 0,2,4,6. Bit 0=A, bit 1=B, bit 2=C. Actual 23°C metadata do
not contain the set-temperature array mentioned by the deposit's general prose.
Do not substitute a nominal setpoint for measured housing temperature.

Timestamps are [year,month,day,hour,minute,fractional second], checked strictly
increasing across all samples and date rollovers. No timezone is invented.
The missing temperature value is recorded but not used in intensity inference.
Contrast data are inventoried, not assumed to calibrate the proposed A phase.

Average PD samples within each setting, group cycle means, and apply the union
of published exclusions to whole cycles with one-based indexing, retaining
328/420 and 433/450 cycles. Zero-based and no-exclusion sensitivities are retained.
Zero-based cuts leave gross outliers and disagree with published summaries;
no additional cuts are introduced to improve an anomaly.

Peres correlations use per-cycle background subtraction. I3 uses absolute means
and coefficients (−1,1,1,−1,1,−1,−1,1). κ is mean(I3)/mean(sum absolute pair
interference), as paper B.3, not mean of ratios. Covariance stays at cycle level;
HAC retains original cycle gaps. No per-setting renormalization, photon-trial
interpretation, inferred V-to-W gain or clipping of invalid denominators is used.

## Measured table extension

`table-audit.json` stores integer count pairs `[0,1]`, indexed by job, preparation,
measurement and outcome. Scan has 9 preparations, Viviani and the Aria control 5;
all have 4 measurements. Explicit job indices preserve the missing scan job 10.
`joint-tables.json` stores physical parameters, predictions, residuals, spectral
bounds, frozen-job results and count simulations. Angles are radians; probabilities,
contrasts and leakage fractions are dimensionless. `theta` is universal only in
the single-system family in [joint-tables.md](joint-tables.md). `leakage_cap` is
a sensitivity setting, never an observed calibration.
