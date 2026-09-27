# Moving-field response and joint dephasing compatibility

## Claim, dynamics and quantifiers

This extension is a **dephasing-compatibility study**. The historical directory
name is retained. Take a real, homogeneous Gaussian potential with covariance

\[
\mathbb E[V_\rho(t)V_{\rho'}(s)]/\hbar^2
={\lambda\over2\tau m_0^2}e^{-|t-s|/\tau}
\int\rho(x)\rho'(y)e^{-|x-y|^2/(4r_c^2)}dx\,dy.
\]

Coordinates here are in **one declared field frame**. Evolution in each
realization is unitary, with Hamiltonian H0+V. The physical measure is the stated
Gaussian probability law, with no state-dependent reweighting. The white limit
reproduces the ensemble double-commutator generator in [derivation.md](derivation.md).
At long memory it does not select outcomes in individual realizations. No claim
about solving the measurement problem follows from the benchmark.

The main witness uses the solar-system barycentric (SSB/ICRF) field frame, rc=100
nm and tau=1e6 s. The benchmark is a 100 nm radius silica sphere, density 2200
kg/m³, with two centers separated by 100 nm perpendicular to its motion, held
for 10 ms in an Earth laboratory at the sodium experiment's date. It translates
in the SSB frame. Its required ensemble coherence exponent is 10. A separate
Earth-corotating field is also defined: translate by Earth's center and rotate
about ICRF z with angular speed 7.292115e-5 rad/s. This is a declared idealized
terrestrial frame, not a calibrated ITRF transform; it is a non-inertial model,
not an inertial frame assigned separately to each instrument.

The assertion is existential: the same spectrum, amplitude and frame enter all
implemented responses. A response **upper bound** below an envelope proves a
sufficient allowed region under the response assumptions. An upper bound above
an envelope does **not** exclude a model. The 13 common-velocity samples around
SSB are a sensitivity profile, not a search over all possible frames. The sampled
0/5/10 km/s neighbourhood is not a physical velocity cutoff. Unsearched frames
remain unclassified. No general outcome-selection or joint-confidence claim is
enabled.

## Distributional limit, actual channels and zero velocity

Let a=lambda/tau. With a fixed as tau grows,

\[
{a\tau\over1+\tau^2 z^2}\to \pi a\delta(z).
\]

This is convergence of positive measures against bounded continuous functions
vanishing at infinity: put y=tau*z, dominate with 1/(1+y²), and apply dominated
convergence. For nonzero velocity v, integrate the spatial kernel along the
coordinate parallel to v. If that marginalized kernel is integrable and
continuous at omega/|v|, then

\[
S^{(2)}(\omega)\to {\pi a\over |v|}
\int_{k_\parallel=\omega/|v|} W(k)\,d^2k_\perp.
\]

The Gaussian factor times the rigid finite-body force/torque form factors meets
these conditions. A delta slice can hit a form-factor zero, and a nonzero slice
need not exceed a measured bound. At v=0 and omega!=0 the point response instead
vanishes; in a finite observation a zero-frequency atom can still leak. These
limits must not be interchanged.

For force along z, W is hbar² rc³/(pi^(3/2)m0²) times
exp(-rc²k²) k_z² |mu(k)|². For torque along z, replace k_z mu with
(k_x partial_k_y - k_y partial_k_x)mu. These follow by differentiating the
potential with respect to displacement and angle, respectively. In the
separable cases implemented for independent limiting checks—layer translation
along a lateral side, or cube translation along its torque axis—the remaining
parallel factor is |H(k)|², H(k)=2sin(kL/2)/k. Consequently

\[
\lim_{\tau\to\infty}\tau F(\omega,\tau,v)
={\pi e^{-r_c^2(\omega/v)^2}H(\omega/v)^2\over |v| I_0(L)},
\quad
F={1\over\tau I_0}\int_0^\infty e^{-t/\tau}\cos(\omega t)I_0(L,vt)dt.
\]

I0(L,x) is the Fourier transform of exp(-rc²k²)|H(k)|². It is evaluated in
real space with error functions. This resolves narrow Lorentzians without a
frequency grid. Tests compare it to independent Fourier quadrature and the
surface limit, including the distinct zero-speed case.

### Two bodies can remain correlated at a delay

For identical bodies sharing a frozen field, the response of (F1-F2)/2 includes
sin²(k.D/2), not an automatic factor 1/2. Large equal-time separation does not
remove correlations when one body reaches the other's field location later.
If D is parallel to v, the frozen slice fixes k.D=omega D/v and can make the
differential response tiny. The implemented LPF bound allows the maximum pair
factor 1, so it covers correlations and anticorrelations; the old independent
mass assumption is not reused for the new moving calculation. The saved aligned
example is a geometry diagnostic, not a claim about LPF's unknown attitude.

## A bound that includes missing attitude and finite windows

Write a rigid channel's generalized mass source as a finite signed measure g:
g=-partial_z rho for force; g=(y partial_x-x partial_y)rho for torque. Define

\[
Q(x)={\hbar^2\over m_0^2}\iint g(ds)g(ds')
 e^{-|x+s-s'|^2/(4r_c^2)}.
\]

Let Q0=Q(0), B=(hbar ||g||TV/m0)² and D be the support diameter. Cauchy–Schwarz
and the support bound give

\[
|Q(x)|\le\min\{Q_0,B\exp[-((|x|-D)_+)^2/(4r_c^2)]\}.
\]

This also applies between rotated versions of the source, because isotropy keeps
the norm Q0 unchanged. Put z=sqrt(log(B/Q0)). Integrating this upper envelope
along a straight trajectory gives

\[
L_{\rm eff}=D+2r_c z+r_c\sqrt\pi\,\operatorname{erfcx}(z),\qquad
S^{(1)}(f)\le 2a Q_0 L_{\rm eff}/|v|.
\]

The same quadratic-form bound holds for a curved trajectory if its projection
onto one fixed direction has derivative at least v_min throughout the segment.
Then |x(t)-x(s)|>=v_min|t-s| and the covariance operator bound follows from
Schur's test. We use the minimum **projected** nominal LPF velocity over the
padded run, not just the minimum speed. Velocity-error radii are subtracted
from that bound. For the unknown sensor epochs, a 28–32 km/s barycentric speed
budget with a 28 km/s chord/projection lower bound over each acquisition block
is an explicit terrestrial-orbit/rotation assumption, not an inferred epoch.
The nominal ephemeris sampling/interpolation and its error budgets are not a
validated spacecraft tracking covariance. This keeps the empirical claim
conditional even though the analytic bound is uniform over orientation.

For any fixed linear observation/filter vector b, the covariance bound controls
E|b*x|² by its white-input equivalent. It therefore survives normalized windows,
differencing and fixed detrending; no favorable leakage direction is assumed.
For an actually stationary process this is the usual positive kernel integral
of its PSD. A fitted or adaptive environmental regression is not covered without
its signal transfer being derived. None is performed here.

### Assembly instead of an isolated load

The load uses all 47 density interfaces. The Si flexure is bounded by a rigid
450 by 57 by 2.5 micrometre box of density 2330 kg/m³, with a vertical flexural
mode in [0,1] normalized at the load. Positivity of the lateral Gaussian kernel
makes the rigid-box norm an upper bound in that mode class. The magnetic sphere
uses radius 15.5 micrometres and density 7430 kg/m³. These geometry values follow
the sensor papers; exact elastic mode and attachment gains are not deposited.
Epoxy is bounded by 1 pL at density at most 4000 kg/m³ (a declared density upper
assumption); its shape is unrestricted. |mu(k)|<=M gives the uniform one-sided
point-mass force bound per a, sqrt(pi)(hbar M/m0)²/(rc v).

The total response is bounded by (sum sqrt(S_i_upper))². All cross terms are
therefore retained regardless of component placement. The result does not
silently call an isolated load conservative. A tenfold additional PSD gain is
stress-tested; the exact headroom factor is saved. This is a conditional
assembly bound, not a recovered mode calibration.

For a cubic torque source ||g||TV=M and D=sqrt(3)L. The single-body Q0 is the
independently checked Fourier torque coefficient. The same upper bound applies
to half the differential torque using the triangle inequality for the two
filtered outputs. Its normalization matches I² S_Delta_gamma/4. Throughout,
S1 per Hz=2 S2(2 pi f); no published lambda label supplies this conversion.

## Benchmark recomputed in the same frame

Let g_R(s) be the Gaussian-smoothed autocorrelation of a homogeneous sphere,
including (M/m0)². We compute it by radial integration of the exact overlap
volume pi(4R+r)(2R-r)²/12. For parallel translating centers with fixed separation d,

\[
{D_{\rm coh}\over a}=\int_0^T(T-t)e^{-t/\tau}
\{g_R(vt)-[g_R(|vt+d|)+g_R(|vt-d|)]/2\}\,dt.
\]

It reduces to G T²/2 for stationary frozen noise, and to G*tau*H(T,tau) for
stationary finite OU noise. The static spatial result is tested against the
independent Fourier sphere calculation. The finite moving integral uses the
spatial crossing time as its integration scale. For the chosen transverse arms,
its value is positive and its fast-motion limit scales as 1/v. Consequently the
required a grows approximately as v: the benchmark cannot retain its old
stationary normalization. Parallel separation gives a different answer, also
tested. Spatial tails beyond the source diameter plus 12rc are exponentially
small; numerical quadrature is convergence-tested, not interval-certified.

## Sodium bridge: controlled local approximation and its limits

In the paraxial description, separate the measured classical beam translation
from transverse quantum dynamics. In the moving frame the spatial covariance
along the center trajectory decays on tc=rc/|v_beam+v_lab-u|, even for a frozen
field. If transverse motion, optical changes and recoil are negligible over tc,
the transverse master equation is the translation-covariant Markov generator
obtained by integrating this correlation. For frozen noise its spectral weights
are the positive k.v=0 slice. The effective separation loss is

\[
a\sqrt\pi {r_c\over |v|}(m/m_0)^2
[1-\exp(-|d_\perp|^2/(4r_c^2))].
\]

The Talbot–Lau characteristic-function solution of this **declared local
Markov generator**, with instantaneous loss gratings, gives the triangular
separation response multiplied by sqrt(pi)rc/|v|. This is not obtained by
identifying two prescribed paths with the exact non-Markov quantum map.
The independent two-time prescribed-path calculation only checks the local
limit. The model is valid where tc divided by the optical transit/flight time,
hbar tc/(m rc²), transverse-drift/relative-speed, and weak-scattering parameters
are small. The first two and the point-cluster size ratio are saved. Finite OU
corrections at tau=1e6 s are negligible relative to spatial crossing times,
but no rigorous operator-norm error theorem is claimed. Unknown transverse beam
emittance and full finite-grating transfer remain empirical limitations. This
extension deliberately does not extrapolate the short-*temporal*-memory formula
of Toroš et al. to large tau.

For unknown orientation, d_perp<=d and the triangle inequality supplies a lower
relative speed from the terrestrial trajectory budget and archived beam speeds.
This yields a damping interval for each mass and velocity. The signed S1
harmonics are summed with their proper interval endpoints, followed by absolute
value and division by the S0 sum. Neither log of an average nor an unsigned
visibility mixture is used. Mass nodes/weights and powers are source-derived;
158±9 m/s is the author's velocity-model input, supported by the paper's TOF
method, not event-resolved data. The 5–7% width range and a width-removal result
are saved as sensitivities. No individual mass/velocity-resolved fringes are
fabricated.

The reported one-standard-error scale is **distinguishability from the audited
optical baseline**, not a confidence limit or a new goodness-of-fit acceptance.
Existing baseline residuals persist. In the SSB frame the field sweeps the
apparatus many times during each count integration, supporting ensemble sampling.
A static Earth-corotating field does not guarantee repeated-event ergodicity;
that frame's ensemble map is not asserted to be the realized fringe experiment.
This distinction is an additional reason to retain the dephasing scope.
For LPF the terrestrial velocity includes v_LPF-v_Earth-Omega cross r_relative;
simply subtracting Earth's translational velocity would be wrong. Its instantaneous
speed range is saved. A single whole-window monotonic chord bound has not been
validated for this rotating trajectory, so the terrestrial space-envelope entry
is explicitly unavailable. No arbitrary minimum speed or point-frequency
substitution turns that frame into a certified joint survivor.

## Acquisition audit and what was actually acquired

- [Sensor paper and supplement, 2002.09782v2](https://arxiv.org/abs/2002.09782v2):
  100 kHz sampling, 2^22 samples, Blackman weights, 40–80 independent-block
  averages over 28–56 minutes; six central bins removed. The deposited
  frequency spacing verifies the block duration and the error columns verify
  averaging counts. Symmetric/periodic Blackman convention, measured Qprime,
  exact epoch, orientation and raw time streams are missing. We calculate
  retained-bin mechanical/window convolution over Qprime=1e5,1e6,2.83e6 and both
  DC-leakage conventions. The omitted frequency tail is bounded using the exact rational Blackman
  Fourier envelope and the maximum mechanical susceptibility. The convolution
  is a continuous-window approximation
  to fine sampling; it is not a replay of unpublished time series. Its sampled
  gain stress is not claimed to bound every possible calibration error.
- [Mission-wide LPF analysis, 2405.05207](https://arxiv.org/abs/2405.05207):
  Table I identifies run 10 as February 14, 2017, duration 13.3 days. Appendix B
  gives Blackman–Harris windows and 50% overlap; Appendix C uses frequency-
  dependent stretch lengths and approximately decorrelated selected bins.
  Angular acceleration subtracts calibrated control torques and stiffness
  terms. Glitches were removed. These facts prevent substituting uncorrected
  DRS telemetry for the published torque envelope.
- [ESA Full dataset](https://esdcdoi.esac.esa.int/doi/html/data/astronomy/lisa-pathfinder/Full.html)
  warns of timestamp inconsistencies and recommends PreProcessed. The archive
  landing/configuration and public frontend were retrieved. A directly usable
  calibrated acceleration/torque/attitude bundle was not resolved. Two attempted
  document retrieval URLs returned HTTP 500; these unsuccessful parameterized
  probes do not prove that the archive is unavailable. No login restriction was
  bypassed and no asynchronous/email data request was submitted.
- [NASA matching February summary](https://heasarc.gsfc.nasa.gov/FTP/lpf/data/summ/drs_20170205_235441__20170225_235341.sum)
  reports **no science HDU**, only 28,800 housekeeping rows. This is more specific
  than the original January probe. Other NASA science products exist but are
  not a matched calibrated reconstruction of this torque result. No unrelated
  DRS run is aligned to the published averaged spectrum.
- [JPL Horizons API](https://ssd-api.jpl.nasa.gov/doc/horizons.html) resolves LPF
  to **-141043**. The fetched 6-hour ICRF vectors span February 14–28, padding
  the run's unspecified start time; Earth vectors share those epochs. Earth
  vectors for April 22–23, 2025 cover the sodium date at hourly cadence. Units
  are km and km/s in transport, converted to SI; times are explicitly UT/JDUT,
  and geometric states have no light-time correction. The LPF header identifies
  source **lpf_20160413**, provided by ESA April 14, 2016, and explicitly says
  subsequent states are nominal predictions. Earth uses DE441. This is **not**
  a retrieved reconstructed mission ephemeris or measured attitude product.
  The 2 and 10 km/s error radii quantify conditional sensitivity, with no
  probabilistic coverage assigned. Current source versions are hash checked;
  only the query-generation timestamp is removed before external caching.

The published LPF envelope is kept distinct from trajectories and PSD estimates.
No averaged PSD is converted into a time-modulation observation. Its actual
segment lengths and calibration are not reverse engineered from a diagram.
The uniform covariance bound is deliberately chosen to make an exact LPF window
unnecessary *for this sufficient conditional response bound*.

## Statistical and spectral scope

The expected PSD is a quadratic form in the Gaussian covariance. Single proper
complex Fourier ordinates are exponential; averaged independent ordinates are
gamma. Overlapping windows and adjacent bins are correlated, and detrending can
make the complex coefficient improper. The code tests covariance propagation
and DC removal by direct Gaussian injection; it does not manufacture a product
of independent gamma likelihoods for these released summaries. The existing
source errors/likelihood audit is retained.

A frozen temporal measure S(omega)=pi*a*delta(omega) supplies a constructive
nonnegative spectrum in the fixed Gaussian spatial model. It is the finite-
variance limit of the OU family, not a grid LP certificate. All finite OU
mechanical bounds above hold uniformly in tau; the benchmark is normalized at
its actual tau. No physical interpolation or infinite-frequency tail certificate
is claimed from the older affine checker. General spectral exclusion remains
disabled because there is already a conditional constructive survivor and the
full nonlinear instrumental likelihood is not established.

## Added value and remaining shortfall

This implements a connection between distinct physics records (sensor and
sodium, with the LPF torque summary) and their auxiliary sampling, geometry and
velocity information. Motion creates a finite mechanical response, changes the
benchmark normalization by nearly nine orders of magnitude, and yields a tiny
sodium perturbation after its signed mass/velocity average. The assembly bound
and delayed cross-body correlations change the old model substantially.
Without trajectory information the same amplitude allows millimetre/second
relative drift, which can yield a much larger torque response; the saved
counterfactual quantifies that missing-information range. With the nominal
trajectory and stated error domain, a response-level survivor is certified by
upper bounds. Adding space data materially reduces the sufficient allowed
amplitude interval, but still leaves the witness. Window and nuisance removal,
source removal, and common-frame samples are executed and saved.

The remaining shortfall against a **fully empirical** joint survivor is precise:
there is no reconstructed LPF torque/attitude likelihood, and the sodium local
Markov approximation lacks an experimentally calibrated error envelope and a
new baseline goodness-of-fit analysis. Neither is disguised as missing generic
metadata. The largest allowed unmodeled response gain is quantified. No novelty
is claimed: frame dependence (Adler 1807.11450), colored mechanics (1805.10100),
colored/dissipative interferometry (1601.03672), rotational geometry (2501.08971)
and windowed LPF inference (2405.05207) are established prior art. The new
repository contribution is the particular executed, source-associated joint
response comparison and its conditional survivor, not the underlying methods.
