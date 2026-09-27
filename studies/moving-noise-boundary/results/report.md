# Moving-noise feasibility boundary

**Gate A is blocked for the requested complementary empirical result.** This is
an implemented response-bound diagnostic, not completion of that research mission.
The selected slow terrestrial field band gives no measured-source identification
gain when the sensor is added. No large scan or joint confidence claim follows.

The [post-review calibration follow-up](../calibration-followup.md) reproduces
the separate 2020 CAL reduced-data energies with an inferred width conversion
and processes two levitated-mass time series. Their calibration gates remain
blocked; they do not change the source-removal result below.

| Executed quantity | Result | Status |
|---|---:|---|
| Shared field variance | 1.5807119e-13 s^-2 | Constructed |
| Held-sphere dephasing lower bound | 10 | Analytic continuum certificate |
| Uniform-band dephasing | 68.418346 | Numerical evaluation |
| Sodium probability perturbation | ≤ 3.7759645e-06 | Conditional exact quantum bound |
| Largest visibility perturbation / measured SE | ≤ 0.010928485 | Descriptive sensitivity, not coverage |
| Sensor response / minimum measured PSD | ≤ 2.1370589e-06 | Conditional assembly/calibration bound, 10× gain stress |
| Field-induced displacement / correlation length | ≤ 1.7449452e-11 | Linear-mode self-consistency diagnostic |
| Unbounded input tail allowed before 1-SE certificate fails | 0.00034055187 | Needed bound; not measured |
| Prior SSB witness under new quantum bound | 1.0 | Vacuous; prior gap remains |

The field law is stationary in one Earth-following rigid coordinate system, with
Gaussian spatial covariance at 100 nm and any positive temporal variance measure
on 0.1–1 Hz. All atoms and between-grid spectra in this declared band are covered;
other frequencies and inertial frame laws are not. The benchmark is recomputed
in that same law. This does not reinterpret PR #12's SSB field as terrestrial.

## What the measured sources do and do not add

All 95 sodium scans (3,895 count bins) and all 14 force-sensor spectra are loaded
from hash-verified originals in the full-source route. The deposited fringe fits
and thermal-regression calibration are reproduced independently. The sensor bound
uses the actual discrete Blackman sample count, both sidebands and both window
conventions. The primary sensor bound additionally covers any fixed orthogonal
sample-space detrending, at a quantified sensitivity cost; no frequency grid
supplies the uniform certificate.

| Source subset | Variance ceiling of sufficient certificate (s^-2) |
|---|---:|
| Sodium | 1.4409452e-11 |
| Sensor | 7.3966699e-08 |
| Both | 1.4409452e-11 |

The gain is exactly 1.0. Leaving this
band or exceeding either ceiling is **not exclusion**. An upper-response
certificate cannot establish a rejected spectrum. Thus these results do not meet
the mission's minimum complementarity requirement. They locate a useful boundary:
a full quantum norm estimate can avoid the sodium Markov approximation, but in
the class where it is informative the existing sensor band adds no identification.

## Why this is not an empirical survivor

The mass/exposure support, independent preparation and apparatus response are
conditional inputs. The archived velocity model does not measure its tails;
preparation history and correlated repeated events are not a calibrated joint
likelihood. A 95% single-realization trace-distance estimate is only
0.0086902, much weaker than the
ensemble estimate. Existing optical residuals remain. No simulated recovery is
misrepresented as recovery through an unavailable calibrated experimental chain.

See [theory](../theory.md), [source/candidate audit](../sources.md),
[prior art](../prior-art.md), and [claim status](../README.md).
Numerical records, nuisance stresses, off-grid spectra, and required next inputs
are in [boundary.json](boundary.json). Original PR #12 files remain unchanged.
