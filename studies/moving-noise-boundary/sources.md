# Candidate audit and source-to-response feasibility

Audit: 2026-09-27 UTC. Current main was 95d209e; PR #12 is merged as 96b7c5c.
Its reviewed head is 77e9ba9bc2ab847b8b037f83839ffda29fb68b56. The existing
moving calculation was rerun unchanged and reproduced. All eight objects in
[the existing pinned manifest](../collapse-compatibility/manifest.json) were
freshly acquired outside the repository and passed byte/hash checks. No contact,
account creation, experimental request, or purchase was made.

## Materially different combinations

| Candidate | Freedom the added source could remove | Executed evidence | Bridge and decision |
|---|---|---|---|
| Sodium + layered sensor, full-quantum perturbation bound + sampled response | Interferometry tests coherence; mechanics probes force derivatives and temporal bands | 95 scans, selected mass distribution, 14 actual spectra and workbook loaded; all fringe fits and thermal slope reproduced; finite observation calculation executed | New norm theorem avoids optical propagation approximation only in a low-variance class. Source-removal gives gain 1 in the selected band. **Fails complementary identification**, retained as quantitative feasibility boundary |
| Sodium + LPF calibrated translational/rotational science + trajectory | Low-frequency mechanical response and different motion/geometry could close terrestrial slow-field freedom | ESA PreProcessed description and live archive frontend retrieved; matching NASA February summary retrieved; prior eight-file baseline reproduced | Examined NASA object has housekeeping only. No calibrated angular science/control reconstruction or reconstructed trajectory uncertainty joined. Uniform orientation bounds already exist, so attitude alone is not the blocker. **Not promoted to joint inference** |
| Sensor + CAL cloud expansion / lensing on ISS | Atomic momentum diffusion/expansion probes a free-flight kernel and space motion, unlike a 3.5 kHz mechanical band | PSI-32 v4 public L1 ZIP downloaded; 18 science/reference PNGs decoded, nine shot labels joined to CSV; PSI-33 layout and command workbook inspected | Real images are available. Object-plane scale/PSF, effective lens/spin fit and clock/ISS trajectory error join remain unestablished; 52 pK is not reproduced. **Measured source obtained; calibrated response gate blocked** |
| Sensor + sodium + radiation | High-frequency acceleration/radiation response could constrain tails of a broad spectral class | Pinned XENONnT ZIP verified as part of the prior full-source acquisition; prior radiation audit inspected | Existing 433-event subset lacks the required energy-migration bridge for the later radiation test. Radiation does not directly close the selected sub-Hz blind band. **Reject as substitute for low-frequency identification** |

These are alternative combinations, not independent successful analyses. Failure
to obtain a product through the inspected routes is not proof that it is absent
from all archives. No retrospective scan for epoch/lag associations was conducted.

## Verified source objects used by the executable calculation

Exact bytes, SHA-256 hashes, and URLs are repeated in
[the generated source summary](results/sources.json). The original immutable
record versions are Zenodo 17502163 v1 and 3956096 v1.

| Source | Measured observable, units, epoch | Response/calibration | Unresolved input |
|---|---|---|---|
| Na-Cluster-Interference.zip, 88,626 bytes | 3,895 integer counts at 95×41 nonuniform measured positions; coordinates nm converted to m; powers mW converted to W; April 22, 2025 directory; rate kHz and integration s | OLS mean/cos/sin in actual coordinates; author mass selection and optical baseline; nominal 158±9 m/s author velocity model; 0.98375 m grating spacing | True mass/exposure tail probabilities, independent preparation, correlated-event likelihood, optical calibration and residuals; model velocity distribution is not event velocity |
| PSD.zip, 290,283 bytes | 14 measured flux spectra, normalized Phi0²/Hz, f in Hz, supplied error column; temperatures 30–1000 mK; acquisition timestamps not deposited | 100 kHz, 2^22 samples, Blackman, averaging inferred independently from (PSD/error)²; 30 mK has a nonuniform retained frequency selection; no invented time series | True mode/attachment response, exact Qprime/mask, epoch and control transfer; differing temperature is not differing geometry |
| Final_data.xlsx, 12,910 bytes | T, Q, T/Q, fitted B and error, cached calibration formulas; published fit values | Reproduce weighted thermal line and use published orthogonal slope to calculate C_phiF; stiffness 0.43±0.013 N/m; resonance 3532.752872... Hz | Calibration transfer to low-frequency driven response is assumed; 10× multiplier is a stress, not measured coverage |

The output also records first/last rows, minima, sample counts, averaging counts,
calibrated slope and the source-model support. Full-source reproduction asserts
fringe and regression agreement at 2e-10 relative/absolute tolerances as specified
in the code. These are numerical reproduction tolerances, not physical accuracy.
The previous residuals and published-vs-refit discrepancy are not hidden.

## Preserved same-frame regression

The unchanged `moving_analysis.py --check` reproduces the original SSB witness:
τ=10⁶ s, λ=35.6535 s⁻¹ and λ/τ=3.56535×10⁻⁵ s⁻², with benchmark exponent 10.
Its consistently stationary benchmark requires λ/τ≈4.62016×10⁻¹⁴ s⁻², a ratio
about 7.72×10⁸. No new input here overturns that sufficient response-envelope
result. The new full-quantum norm bound at its variance (λ/2τ) saturates at 1;
therefore it **cannot repair its sodium approximation bridge**. The new low-band
Earth-following calculation has its own consistently held benchmark and is not
an inertial-frame refit or a demonstration that motion makes the old witness safe.
The two frame laws are kept distinct. No new tracking accuracy is inferred.

## LPF routes actually checked

The [post-review follow-up](calibration-followup.md) additionally recovers the
ESA frontend's documented synchronous machine interface. Its metadata example
fails with a server-side database-session error, and the documentation-product
download fails with HTTP 500. This is an observed service failure, not proof of
missing products. That follow-up also records the CAL workbook reconstruction,
two actual levitated-mass time-series analyses, and other low-frequency pivots.

- [ESA PreProcessed](https://esdcdoi.esac.esa.int/doi/html/data/astronomy/lisa-pathfinder/PreProcessed.html),
  DOI 10.5270/esa-9z4vm0u, v1.0, covers 2016-02-29 through 2017-07-18. It explicitly
  describes corrected timestamps. The metadata is not a calibrated force product.
- [Archive frontend](https://lpf.esac.esa.int/lpfsa/) returned HTTP 200, 3,088 bytes
  of GWT application shell. This establishes an accessible frontend, not a
  downloaded matched science bundle. The prior endpoint audit is preserved;
  no unsupported query URL was treated as proof of data availability.
- [NASA matching summary](https://heasarc.gsfc.nasa.gov/FTP/lpf/data/summ/drs_20170205_235441__20170225_235341.sum)
  returned HTTP 200, 144 bytes, and reports no science HDU with 28,800 housekeeping
  rows. Its displayed MJD numbers do not themselves match the filename year; they
  are not silently used for a trajectory join. The filename overlap is not a
  calibrated science-time overlap.
- [Dynamics calibration](https://arxiv.org/abs/1806.08581) concerns physical force,
  actuation and stiffness calibration. Raw displacement or commanded acceleration
  is not the desired independently calibrated residual force/torque.
- [Mission-wide noise analysis](https://arxiv.org/abs/2405.05207) and the previous
  run audit remain the specific run/filter reference. We did not retrieve a new
  matched angular-science/control bundle or promote the earlier nominal JPL
  prediction to reconstructed tracking. PR #12's 0/2/10 km/s radii remain choices.

The uniform-orientation result in PR #12 can remove a demand for measured attitude
for a sufficient covariance bound. It does not remove the need to know the
observable's calibration or justify a motion-error set for a physical claim.
The new study therefore does not impose unavailable attitude as an unnecessary
prerequisite and does not reuse the torque envelope as raw data.

## Cold-atom alternative beyond the original pairing

[The detailed CAL audit](cal-response.md) corrects the initial availability
assessment: current NASA PSI-32 v4 releases real absorption/reference arrays.
A 107.9 MB nine-shot sweep, its embedded timestamp table, and two PSI-33
interpretation documents are pinned in [cal-manifest.json](cal-manifest.json).
[The executable audit](cal_audit.py) loads all 18 images and reproduces diagnostic
pixel fits; [its output](results/cal-audit.json) records the exact missing bridges.
The supplement and primary instrument documentation were inspected. A nominal
"Calibration" workbook contains trap commands, not imaging calibration.

The additional [free-cloud kernel](cal_response.py) retains centroid subtraction,
spatial correlation, Doppler shift and finite observation, and is tested by
independent quadrature and Gaussian injections. It shows why cloud expansion
could supply a different response, and why common force cannot simply be called
heating. It does **not** supply the missing calibrated CAL quantum/lens/imaging
response or a new empirical noise limit. No ISS ephemeris is substituted without
the release/image clock and calibration join. [Bilardello et al.](https://arxiv.org/abs/1605.01891)
supplies prior cold-atom collapse theory, not that archived apparatus map.

## Next step: analysis versus records versus a new observation

1. **New analysis on existing archives:** prove a useful moving-field quantum
   bound (or implement validated finite-grating propagation) including independent
   preparation and mass/exposure tails. The new conservative bound is too weak
   at the old SSB strength. A better theorem could change that without new data.
2. **Existing-but-unobtained records:** matched LPF calibrated science/control
   output and reconstructed trajectory uncertainties; or CAL imaging/PSF, effective-lens, clock and ISS trajectory calibration joined
   to the **obtained** shot records. The required calibrated bundle remains unestablished.
   These must pass calibration and source-removal gates before a broad scan.
3. **If those cannot close the slow band:** a held massive superposition in the
   same declared frame at several holding times, or calibrated force observations
   covering the field band, would interrogate the explicit surviving kernels.
   For the 100 nm sphere and the constructed variance, the analytic lower bound
   reaches D=1 at no more than about 3.2 ms (with the same spatial separation and
   frame law). This is a prospective response requirement, not an experiment
   request, cost estimate, or claimed practical feasibility.

No raw archive, workbook, paper, or external source code is redistributed.
