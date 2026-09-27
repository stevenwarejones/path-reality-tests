# CAL source audit and the response that is still needed

The initial literature route understated availability. NASA's current public PSI
catalogue **does contain downloadable CAL absorption/reference images**. This
follow-up retrieved and decoded a real nine-shot sweep rather than stopping at
the paper's 2022 author-request/future-release statement. Lack of an established
calibration join, not lack of public images, is the remaining obstruction.

## Actual acquisition and reproducible audit

The public UI at <https://psi.nasa.gov/physci/repo/search> uses these read-only
services (no account or restricted file was accessed):

- GET `https://psi.nasa.gov/geode-py/ws/repo/investigations/PSI-32` gives record
  version 4, public and unrestricted, DOI `10.60555/6cnc-re60`.
- POST `https://psi.nasa.gov/geode-py/ws/files/v2` with
  `{"folder":"genelab-data/","fileType":"study","studyId":"PSI-32","version":4}`
  lists files. The record's release status, not the search index's zero-valued
  `PublicReleaseDate`, establishes release in this audit.
- The versioned download resolvers and exact SHA-256/lengths of the three retained
  source objects are in [cal-manifest.json](cal-manifest.json). Expiring transport
  URLs are not stored. Files stay outside git.

The selected 107,889,947-byte L1 ZIP contains 18 2048×2048 16-bit PNGs and a
measurement mapping CSV. It is a nominal 1.8 ms lens sequence with nine images at
17, 17, 37, 57, 77, 97, 117, 137 and 157 ms labels. These are sequence labels,
not a verified release-to-exposure timing calibration. The CSV joins the images
to **2019 day 197 (July 16), 00:52:09–01:06:14**, with no clock scale/offset
established. The `190622` sequence name must not be used as an acquisition date.
Four additional mapping rows for the 157 ms path have `not_applicable` image
times; the executable records these and requires exactly one image-bearing row.

The PSI-33 data-layout document identifies image0 as containing atoms and image1
as the reference. The PSI-32 record independently describes the same destructive
absorption/reference acquisition. This convention is verified, not guessed from
an attractive cloud. The PSI-33 `Calibration_test1.xlsx` is **commanded ramp
currents/fields**, not camera magnification, shot noise, or measured cloud radii.
Its text calls the total 0.3 s while its tabulated endpoint is 250 ms; it cannot
silently calibrate a different investigation's 1.8 ms lens. The separately
inspected PSI-33 "Expansion Profile" CSV is also a command table, not measured
expansion data. A 2018 instrument paper specifies physical camera pixels of
5.5 μm; that does not establish the object-plane magnification of the 2019 shots.

`cal_audit.py` hash-checks all three objects, decodes every pair, joins the CSV,
and forms log(reference/atoms) in one fixed pixel ROI. The ROI has no zero or
saturated intensity pixels. A Gaussian plus planar background gives descriptive
pixel diagnostics, **not** the paper's spin-resolved Thomas–Fermi/thermal/PSF fit.
The late-shot peak amplitudes (about 0.02–0.05) are comparable to or below the
roughly 0.054 pixel residual RMS. Even the two 17 ms diagnostic fits disagree
substantially. These fits cannot reproduce the reported 52 pK or justify a
confidence region. They are retained to make the source inspection executable,
not to extract an opportunistic noise association. ROI selection followed visual
discovery; this is not a preregistered selection or a modulation test.

The [2022 paper](https://www.nature.com/articles/s41467-022-35274-6) and its
[Supplementary Figure 6](https://static-content.springer-cdn.com/esm/art%3A10.1038%2Fs41467-022-35274-6/MediaObjects/41467_2022_35274_MOESM1_ESM.pdf)
were inspected. They explicitly treat magnetic-lens switching, interacting cloud
expansion, split spin clouds and rejected low-density points. The supplement
compares an effective lens shortened by 0.4 ms; this changes inferred energies
substantially. The main methods use a roughly 42 μm resolution term and a
Thomas–Fermi radius-to-energy map with factor 7. Those published model parameters
are leads for a reproduction, not independent measured uncertainty envelopes for
this selected archive segment.

## Same-field finite-observation kernel and an exact blind spot

Use the covariance measure, mass-normalized coupling and **single frame law** of
[theory.md](theory.md). Write its full real spectral measure as μ(dk,dω), with
phase k·y−ωt and reflection symmetry. A point particle's acceleration in a fixed
measurement direction e is −(ℏ/m₀)e·∇φ, independent of its mass. For prescribed
free paths y_j(t)=vt+r_j and weights w_j summing to one, define

\[
 I_T(\Omega)=\int_0^T(T-t)e^{-i\Omega t}dt
 =\frac{1-e^{-i\Omega T}-i\Omega T}{\Omega^2},\qquad I_T(0)=T^2/2,
 \quad F(k)=\sum_jw_je^{ik\cdot r_j}.
\]

To first order in force-induced displacements, the added **centroid-subtracted
cloud width variance** is the integral of

\[
 K_{\rm width}(k,\omega)=\left(\frac{\hbar}{m_0}\right)^2
 (k\cdot e)^2|I_T(\omega-k\cdot v)|^2\,[1-|F(k)|^2].
\]

The centroid kernel has |F|² in place of 1−|F|²; the single-particle kernel has
1. This follows by averaging squared displacements and subtracting the squared
weighted mean. **Uniform force produces centroid motion but no width increase.**
An individual atom's velocity variance cannot automatically be called cloud
heating. No spatial independence approximation was used. The line kernel is
positive, handles off-grid atoms and has the frozen observed-force limit T⁴/4.
The free-flight transfer samples low observed frequencies missed by the sensor's
retained kHz band. Two different flight durations give nonproportional frequency
kernels. This establishes *potential* response diversity, not measured added
identification; flight time, geometry and motion all still need a calibrated map.

For expanding, accelerated or rotating paths and a lens, replace I exp(ik·r_j)
by A_j=∫G_j(T,t)(k·e(t)) exp[i k·y_j(t)−iωt]dt. The width kernel is
(ℏ/m₀)²[Σw_j|A_j|²−|Σw_jA_j|²], retaining all cross and delayed correlations.
G_j is the displacement Green function including the physical lens/interaction
linearization; simply putting G=T−t through the lens is not justified. This
formula is only a controlled approximation when force-induced trajectory changes,
cloud interactions and imaging can be bounded. It is not a new exact quantum
CAL propagator. The same Earth-following law requires ISS trajectories in those
coordinates; assigning ISS its own comoving stationary field would change the
hypothesis. We did not substitute a nominal ISS speed or a borrowed epoch.

`cal_response.py` implements the free-path diagnostic with a small-argument
series, avoiding cancellation at arbitrarily narrow/zero-frequency atoms. Tests
compare independent quadrature, frozen limits, nonproportional durations, Doppler
sign and Gaussian-quadrature injected force realizations with centroid refitting.
A common acceleration/pointing adversary cancels exactly. These are **kernel
injections**, not recovery through the unestablished CAL calibration chain.

## The bridge that failed and the next records needed

The subsequent [calibration follow-up](calibration-followup.md) obtains a
numerical physical baseline from the separate 2020 reduced-data workbook,
with a documented width-convention qualification. It does not establish the
calibration of these 2019 lens images or replace the missing entries below.

| Needed item | Obtained status | Why it matters |
|---|---|---|
| Science image pairs and relative sequence labels | Obtained, pinned and decoded | Establishes that an actual measured alternative exists |
| Acquisition labels | Image-to-CSV join obtained | Prevents borrowing the sequence creation date for orbit matching |
| Object-plane scale, dark/flat/saturation response and PSF uncertainty | Not joined for these shots | Width in pixels is not width in metres; an unconstrained scale can reverse a physical bound |
| Spin-resolved selection/fit covariance, initial correlations and effective lens | Paper supplies model description, no reproduced calibrated shot fit here | Ordinary lens/interaction changes can mimic a change in expansion |
| Clock convention, release/exposure mapping and ISS position/attitude errors | Not established | Sets the common-field Doppler kernel; metadata labels alone are insufficient |
| Published physical baseline | **52 pK not reproduced** | CAL does not pass Gate A and is not used in a joint limit |

These are existing-data/calibration/analysis requirements, not evidence that new
measurements are necessary or that the records do not exist elsewhere. The
implemented kernel specifies which calibration and trajectory map a follow-up
must validate. Until then no CAL upper bound, likelihood, source-removal rejection
or global spectrum optimum is reported. This is the stopping obstruction after
credible source retrieval, rather than an excuse to scan an uncalibrated model.
