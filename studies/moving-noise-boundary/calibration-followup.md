# Calibration follow-up after review of `bea7bf9`

**The review's scientific assessment stands.** The combination gain is still
1, the original moving SSB witness is unresolved, and this PR does not meet the
mission's minimum measured complementarity requirement. This follow-up tried
to acquire the missing bridge, rather than tighten the redundant certificate.
It adds a reproducible physical-quantity reconstruction from CAL source data
and a second, independently calibrated mechanical-source investigation. Neither
is promoted to a joint limit. No spectrum is claimed to be empirically rejected
by their combination.

The machine-readable results are [calibration-gate.json](results/calibration-gate.json).
The small extracted observations and nominal reconstructions are in
[calibration-inputs.json](results/calibration-inputs.json); their original
files, URLs, sizes and SHA-256 hashes are in [gate-manifest.json](gate-manifest.json).
The implemented extraction and response checks are in [calibration_gate.py](calibration_gate.py).

## CAL: a numerical physical baseline, with a normalization qualification

The public [2020 Figure 4 source workbook](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41586-020-2346-1/MediaObjects/41586_2020_2346_MOESM1_ESM.xlsx)
contains six expansion-width observations in micrometres at 168–768 ms, as well
as 23 centre-of-mass positions. These are measured, already reduced data,
not synthetic expansion curves or the subsequently deposited 2019 lens images.
The [public figure caption](https://www.nature.com/articles/s41586-020-2346-1/figures/4)
gives energies 231(9) and 720(79) pK. Its error bars are labelled standard
deviations with at least three observations; the covariance and exact counts
per point are not supplied in the workbook.

Using all six points, with the supplied error columns as diagonal weights,
fit `w(t) = a + b t`. The calculation independently reproduces both energies
and their quoted uncertainties to within 1 pK **if** the source-width slope is
divided by sqrt(2) before applying `E/k_B = m v_rms²/k_B`:

| Axis | Source-width slope (µm/s) | Slope / sqrt(2) (µm/s) | Reconstructed energy (pK) | Propagated formal error (pK) | Caption |
|---|---:|---:|---:|---:|---:|
| x | 210.2572 | 148.6743 | 231.0485 | 8.6407 | 231(9) |
| z | 371.1319 | 262.4299 | 719.8760 | 78.8858 | 720(79) |

There is a real qualification: **the workbook itself calls these columns rms
widths**. Taking that label literally yields 462.097 and 1439.752 pK under the
caption's velocity-to-energy convention. The factor sqrt(2) is inferred from
the numerical match, not established here by a recovered imaging fit function.
It is compatible with conversion of a Gaussian 1/e width to its standard
deviation, but that explanation is not verified for these reduced data.
The [published correction](https://www.nature.com/articles/s41586-020-2517-0)
addresses millimetres versus micrometres on the figure axis, not this convention.
The script retains both interpretations instead of silently correcting the
source labels. Its formal errors are a reproduction of diagonal weighted
regression, **not a new confidence region for correlated observations**.

This resolves part of the baseline question but does **not** reproduce the
2019 52 pK lens result from images. A 2020 reduced-data length column cannot
calibrate a different image sequence by itself.

### What was checked in the image archive

In addition to the pinned nine-shot day-197 sweep, the public PSI-32 version-4
inventory was searched for matching free expansion and late-2019 lens sequences.
An additional science ZIP was actually downloaded and decoded:

* file ID `66eb934dc6ad4bdc12fa5edd`, version 1;
* `2019_271_190927_V24_3652462534_Bigelow_STA_2p4_hold_24_free_exp_TOF_2_L1.zip`;
* 71,713,667 bytes; SHA-256
  `f46eaadd8a7f2298cf936903962529eb81faad7a8d64e9745f00244906de6a7c`;
* three 87 ms and three 97 ms image/reference pairs, plus the mapping CSV.

The archive name describes a promising sequence; it does not identify the
published selected shots. The old day-197 region of interest does not transfer
to these frames. Whole-frame inspection shows strong optical structure; no
validated condensate/spin fit or physical-width measurement was obtained from
these six pairs. They are discovery evidence, not inputs to the regression
above, and are not added as another fitted dataset to the joint study.

Inspected PNG metadata specifies 2048×2048 storage, a 16-bit container with
12 significant bits, and NI image-format metadata. It does not supply the
object-plane magnification. The 2018 instrument paper's sensor pixel pitch is
not a measured object-plane scale for a 2019 run. Later post-upgrade CAL recoil
calibrations likewise cannot be transferred without an instrument-state join.
The PSI-33 “calibration” workbook and expansion CSV are commanded trap ramps,
not an image-scale or measured expansion calibration.

| Missing bridge | What existing records recover | What remains unbounded or unjoined |
|---|---|---|
| Length/width definition | 2020 reduced widths, units, reported scatter; numerical energy reproduction | 2019 magnification and uncertainty; actual width convention behind the 2020 workbook |
| Imaging response | Published 2022 PSF-floor and spin/thermal/TF fitting description | Shot-level PSF, background model, accepted spin populations, fit covariance and selection |
| Lens dynamics | Published timing/model and sensitivity to effective duration | Effective waveform and initial correlations for the specific downloaded shots |
| Acquisition/motion | Image labels joined to CSV; no borrowing of sequence-creation dates | Release-to-image clock mapping, time scale, ISS orbit and orientation/error join |
| Statistical inference | Six measured reduced widths and their reported marginal errors | Shared calibration covariance, sample counts/correlations and likelihood |

The earlier centroid-subtracted finite-path kernel in [cal-response.md](cal-response.md)
still specifies what to propagate. It cannot turn an uncalibrated image width
or an unjoined lens Green function into a constraint. Neither image download
nor reproduction of a different reduced-data energy is a pass of the full gate.

## Pivot: a levitated 26.7 Hz mechanical source

[Fuchs et al., Science Advances 10, eadk2949 (2024)](https://doi.org/10.1126/sciadv.adk2949)
explicitly discuss collapse-model applications. The choice here concerns its
mechanical force response; it is not an anomaly search in the gravitational
drive. [Zenodo record 10300430](https://zenodo.org/records/10300430) provides eight
raw lock-in exports. The [versioned author manuscript](https://arxiv.org/pdf/2303.03545v2)
includes the calibration and analysis supplements.

Two records were downloaded and numerically processed:

| Lateral / vertical source-wheel position | Export date label | Samples in each of 9 channels | Duration | Fitted amplitude decay time |
|---|---|---:|---:|---:|
| −3 cm / 48.1 cm | 2022-05-14 09:56:06 | 151,085 | 43,000.54 s | 210,397.7 s |
| +3.5 cm / 48.1 cm | 2022-05-13 05:03:23 | 75,429 | 21,467.83 s | 111,170.0 s |

The second geometry is explicitly illustrated in Supplement Fig. S6. This
selection was for a documented baseline, not selection of the quietest record.
The export time scale is not established as UTC, and the numerical clock is
relative to the export. No ephemeris is attached to those labels.

The parser verifies the actual nine sequential `(time,value)` blocks and their
identical time grids. They are three amplitudes in volts, three phases in
degrees and three oscillator frequencies in hertz, sampled about every
0.2846135 s. Demodulator 2 is identified provisionally with the decaying 26.7 Hz
mode, and demodulator 3 with the detuned wheel signal. The frequency and amplitude
behaviour support this assignment, but a wiring/configuration record is absent.

Supplement C supplies the nominal displacement gain `0.16 V/µm`, a 7% relative
uncertainty, and a SQUID gain `0.43 V/Φ0`. Supplement E describes ringdown
subtraction, a phase rotation and selection of integer wheel cycles. Under
that flux conversion the +3.5 cm fit evolves from 15.827 to 13.461 mΦ0 over
five hours, matching the scale and direction of Fig. S6. Its fitted decay time
is about 2% above the published 109,000 s lower bound. **This is not a
reproduction within the quoted 14.7 s fit error**: the paper explicitly reports
amplitude-dependent and coupled-mode drift, and its exact selected sample
interval is not supplied in the export.

### Nominal force reconstruction and failed baseline

The code subtracts the fitted decaying envelope, removes a linear demodulated
phase, and transforms the complex residual. With the RMS convention
`x(t)=sqrt(2) Re[z(t) exp(iω_c t)]`, the two-sided complex baseband PSD maps to
the positive-frequency one-sided physical PSD without another factor two.
Independent real-signal and Parseval tests include off-grid lines.

The nominal inverse susceptibility is

\[
\chi^{-1}(\omega)=m(\omega_0^2-\omega^2+i\gamma\omega),\qquad
m=0.43\ {\rm mg},\quad\gamma=1.84\times10^{-5}\ {\rm s}^{-1}.
\]

The dimensional `Q/f` bandwidth typo in Supplement E is not propagated;
the implemented susceptibility uses the angular damping convention above.
The fixed diagnostic band is 2–4 mHz on both sides of the estimated mode
frequency. Hann and rectangular reconstructions are both retained; there is
no optimized frequency/lag/calendar search or empirically calibrated CI.

The Hann nominal ASDs are approximately **4.91 and 2.77 fN/√Hz**, respectively.
They are 9.83 and 5.54 times the published best floor of 0.5 fN/√Hz. That best
floor need not characterize these two whole records. The mismatch could involve
record/interval selection, acquisition settings, phase treatment, nonlinear
mode response or a gain join; it is not evidence against the publication.
No multiplier was fitted to force agreement. A 7% gain nuisance cannot explain
these factors by itself. This is a failed noise-baseline reproduction, with the
discrepancy preserved as a regression output.

### The acquisition freedom that blocks an exclusion

The export does not contain the demodulator low-pass order/time constant,
external gain, input range or complete analysis history. The stated voltage
calibration does not establish a lower response throughout a baseband interval.
For illustration, a conventional n-pole low-pass response gives

\[
|H(\delta f)|^2=[1+(2\pi\delta f\tau)^2]^{-n}.
\]

This is an explicit sensitivity family, **not recovered instrument settings**.
At a 4 mHz offset an eight-pole filter has power gain 0.613 for τ=10 s,
1.22×10⁻⁷ for τ=100 s and 3.90×10⁻²³ for τ=1000 s. A justified bound
τ≤4.58 s would instead guarantee at least 90% power gain there for n≤8.
No such bound is established by this archive examination. If τ is left
unbounded, its response infimum at nonzero offset is zero; a finite scan of
assumed settings would not supply the missing lower certificate.

The injection test sends a known harmonic force through the declared oscillator,
gain, low-pass filter and finite Fourier acquisition. Ignoring the filter
recovers the attenuated force, exactly as it should. It does not test a fully
calibrated archive pipeline or infer an empirical filter setting. This control
prevents interpreting attenuation as an exclusion of an added physical field.

## Common-field relevance and why no new combined limit is reported

For a scalar mass-coupled field the additional source must constrain the same
measure, not independently normalized force noise. For any mechanical mode
with displacement shape `u(r)`, its plane-wave force factor is

\[
A(\mathbf k)={\hbar\over m_0}
 \int\rho(\mathbf r)\,i\mathbf k\!\cdot\!\mathbf u(\mathbf r)
 e^{i\mathbf k\cdot\mathbf r}\,d^3r.
\]

The finite-observation response must include the common frame phase along the
apparatus trajectory, this factor, the mechanical Green function, the
acquisition filter and the fitted-nuisance projection. Correlated terms from
the three magnets, attached glass bead and glue must remain inside the squared
sum. The publication's nominal total mass and outline do not establish a
uniform lower mode-force factor, particularly with its documented nonlinear
coupling. Treating the assembly as a point mass would be inappropriate at
`r_c=100 nm`.

The 26.7 Hz mode offers a distinct response from the 3.5 kHz layered sensor
in an expanded spectral class. **It does not by itself close this PR's
0.1–1 Hz Earth-following blind band.** A narrow lock-in archive is not a
calibrated broadband low-frequency displacement record. Its stationary
terrestrial trajectory in that common Earth-following law cannot be replaced
by a separately chosen frame to shift the band into resonance. There is
therefore no reason to claim a global combination gain merely from obtaining
this lower-frequency instrument. For inertial fields, its own epoch/motion
uncertainty and the sodium quantum-response problem also remain.

## Other serious alternatives and exact stopping obstruction

| Candidate | Freedom it could remove | Verified finding | Why not promoted |
|---|---|---|---|
| CAL raw and reduced expansion | Atomic relative acceleration, free-flight/lens response, ISS motion | Science arrays; reduced-data energy reconstruction above | Specific lens/imaging/motion join still missing |
| Fuchs levitated mass + sensor | Different mechanical mode and geometry, tens-of-hertz spectral weight | Two actual multichannel records and published voltage calibration | Noise baseline, acquisition lower response and assembly mode not validated; slow band persists |
| ESA LPF calibrated acceleration | Low-frequency space response with different velocity and baseline geometry | Public frontend's documented synchronous metadata and file routes recovered | A documented metadata example returns HTTP 600 with a null Hibernate-session error; document retrieval returns HTTP 500; no calibrated science file obtained by this route |
| Vinante 2020 Meissner sphere | Simpler translating body with low-frequency modes | DataCite resolves DOI [10.5258/SOTON/D1402](https://doi.org/10.5258/SOTON/D1402) to the Southampton archive | Public landing/file requests returned HTTP 401; no raw dataset obtained, no bypass attempted; damping alone does not bound random-unitary heating |
| Westphal 2021 torsion balance | Millihertz differential force geometry | Publication and data-availability statement inspected | Raw data available on author request; not obtained and no author contacted |
| Ren 2026 levitated magnetometer | Low-frequency torsional response | Public paper and Science Data Bank landing found | Main cylindrical magnet is axisymmetric about the torsion axis; magnetic torque sensitivity is not its scalar mass-density torque sensitivity. Ancillary assembly geometry would be essential |

The ESA failure is an access/service observation on 2026-09-27, not evidence
that the calibrated products do not exist. The documented failing example was
`/lpfsa-sl/metadata-action?RESOURCE_CLASS=EXPERIMENT&SELECTED_FIELDS=EXPERIMENT.TITLE,EXPERIMENT.COLOR_CODE&PAGE_SIZE=100&RETURN_TYPE=JSON`.
The response reported `TransactionManager.unbindAndReleaseSession(): Input
hibernate session is null`. Metadata requests for analysis-object IDs and the
document download also failed; housekeeping was not substituted for science.

**Next requirement:** recover existing shot-level CAL calibration/selection and
clock/orbit records, or an independently documented low-frequency science
bundle with a bounded acquisition response and known scalar-force geometry.
The Fuchs route specifically needs demodulator settings, run gain/wiring,
selected analysis intervals and validated mode response, and would still need
a source covering the slow band. These are existing-record and analysis
requirements; this work does not establish that new measurements are necessary.
If those records cannot be recovered, a new measurement must publish those
quantities together, rather than just a quoted best noise floor.

No authors were contacted, access controls bypassed, or experiments initiated.
The chosen endpoint is **not achieved**. The remaining obstruction is a
calibrated lower/likelihood response and the sodium moving quantum bridge,
not computation time or insufficient spectrum-grid density.

## Reproduce and interpret

```sh
python studies/moving-noise-boundary/calibration_gate.py --check
python studies/moving-noise-boundary/calibration_gate.py --data-dir /tmp/gate-sources --fetch --check
python -m unittest discover -s tests -p test_noise_calibration_gate.py -v
```

Offline checking recomputes the CAL regressions and filter sensitivities from
the committed small summaries. Full-source checking verifies all hashes,
reopens the workbook and both time series, checks the clocks and recalculates
the nominal spectra and ringdowns. Neither route claims the absent calibration
products have been supplied. This is retrospective exploratory source work;
the numerical regression tolerances and diagnostic band are reproducibility
checks, not preregistered tests or selection-adjusted statistical claims.
