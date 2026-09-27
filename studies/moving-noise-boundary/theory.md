# Model, proofs, and boundary of the claim

## One field, one frame law

Let y = R(t)^T[x-X(t)] be one declared rigid Earth-following coordinate law.
R and X define the **physical noise model**, not a separate fitted frame for each
experiment. An instrument attached to that idealized rigid Earth has fixed y;
a sodium beam still travels through it. We do not use this law for LPF without
a matched trajectory. Actual vibration, Earth deformation, and mode motion are
not set to zero experimentally: the mechanics below is a conditional linear-mode
response. An inertial SSB-stationary field is a different model class.

Define a zero-mean **Gaussian** scalar field with units s^-1 by

    E[phi(y,t) phi(y',s)]
      = exp(-|y-y'|²/(4 rc²)) integral_[0.1,1] cos(2 pi f (t-s)) mu(df),
    sigma² = integral mu(df),          mu >= 0.

The spatial Fourier measure is rc³ pi^(-3/2) exp(-rc² |k|²) d³k, normalized to
one. With the convention exp[i(k.y - omega t)], the two-sided joint measure is
the product of this spatial probability measure and the symmetric temporal
measure with half the weight at each of +/- 2 pi f. It is nonnegative and even
under (k,omega) -> (-k,-omega), giving a real homogeneous stationary field in y.
No one-/two-sided PSD factor is hidden in sigma²: it is the equal-time variance.
For a temporal density, the one-sided PSD per Hz is dmu/df.

The mass coupling is H_noise(t) = (hbar/m0) integral rho(x-q) phi(y(x,t),t) dx,
where m0 = 1 atomic mass unit, mass density is kg/m³, and rc = 100 nm. The Gaussian
covariance already provides the spatial smearing; it is not applied a second
time. Each realization evolves unitarily with H0 + H_noise, and the probability
law is the stated Gaussian measure. There is no state-dependent reweighting or
outcome selection. A rigid positive density has Fourier form factor bounded in
modulus by its total mass. The quantum perturbation bound below thus also covers
finite cluster size, arbitrary orientation, recoil, and deterministic optical
operations within its exposure interval.

For general trajectories the covariance uses y(t)-y(s), not an instantaneous
speed. For straight motion y(t)=y0+v t, the phase is -(omega-k.v)(t-s), so
omega_observed=omega_field-k.v. The implemented finite covariance and path
quadratic forms retain curved trajectories and delayed inter-body correlations.
Those tests validate conventions; they are **not** a replacement for sodium
quantum propagation or an empirical LPF response.

## Exact ensemble quantum perturbation bound

For a noise-independent initial state, arbitrary deterministic H0(t), and a mass
at most m=m0 M exposed for at most T, let Phi and Phi0 be the ensemble and noiseless
channels. Then

    (1/2) ||Phi-Phi0||_diamond <= min(1, (exp(2 sigma² M² T²)-1)/2) = epsilon.

Proof: first take finitely many Fourier pairs and finite-dimensional truncations.
In the interaction picture each Fourier coupling is a unitary conjugate of a
position kick, of norm at most M (rigid extended form factors only reduce it).
The commutator superoperator has completely bounded trace norm at most 2M.
Expand the time-ordered channel Dyson series. Odd Gaussian moments vanish.
At order 2n, Wick's rule has (2n-1)!! pairings. Each contracted Fourier pair has
absolute integrated covariance at most sigma²: phases have modulus one and the
spectral measure is positive. Thus the order-2n diamond norm is at most

    (2M)^(2n) (sigma²)^n (2n-1)!! T^(2n)/(2n)!
      = (2 sigma² M² T²)^n/n!.

Summation gives the bound. This bound is uniform in mode number and Hilbert-space
truncation; the absolutely convergent expected Dyson series defines the ensemble
map on trace-class states as these approximations converge. It does not require
an operator-norm bound on an infinite-volume Gaussian field realization.
A reference ancilla does not change any operator norm, proving the diamond form.
Deterministic loss operations can be dilated to unitaries with retained loss
flags, so they are covered *before conditioning on detection*. Any environment
introduced this way must also be independent of the field initially. Arbitrary
noise-correlated preparation is not covered.

This is a deliberately conservative non-Markov bound. It makes no spatial
advection, paraxial, weak-coupling, or prescribed-two-path approximation. A
non-Gaussian field with the same covariance need not obey the Wick bound.
White temporal noise has infinite equal-time variance and is not in this class;
only its finite integrated prescribed-path phase limit is tested separately.

For a convex population whose bounded mass/exposure part has weight 1-eta, the
bound becomes eta+(1-eta)epsilon. Unknown tails therefore have an explicit cost.
An author velocity model is not an empirical tail certificate. The protocol's
195500 u / 25 ms domain includes the deposited selected mass nodes and nominal
12.45 ms flight, but these are **assumed** support/initialization conditions.
The bound at PR #12's SSB witness saturates at 1: it does not validate that witness.

The unconditional probability of any event changes by at most epsilon. If
noiseless survival q has a known lower bound, conditional total variation is at
most min(1,2 epsilon/(q-epsilon)) when q>epsilon. Never identify unconditional
channel closeness with postselected fringe closeness.

For the actually used OLS [1,cos,sin] fit, let l0,lc,ls be the absolute row sums
of the design pseudoinverse and let q,V be the noiseless mean probability and
visibility. A uniform per-bin error epsilon gives

    |Delta V| <= epsilon (sqrt(lc²+ls²)+|V| l0)/(q-epsilon l0).

This retains the nonuniform measured positions and signed cosine/sine components.
The saved optical q,V are fixed-calibration predictions; measured standard errors
are descriptive scales. No baseline goodness-of-fit or confidence acceptance is
inferred. The existing residuals remain. The 1-SE ceiling is a **certificate
threshold**, not a bound on the true field strength.

A single realization is different. Duhamel, unitarity, and Cauchy–Schwarz give
E||psi_phi-psi0||² <= sigma² M² T² for independent input. Markov's inequality
only gives trace distance <= sqrt(sigma² M² T²/alpha) with probability 1-alpha.
There is no union over archived particles here. Correlated events and shared
noise mean that multiplying individual probabilities or nominal intervals does
not produce a joint likelihood. Temporal atoms in particular need not be ergodic.
The continuous uniform-band example avoids an atom, but no calibrated finite-run
mixing likelihood is claimed for it either.

## Held benchmark: a continuum certificate, not a sampled maximum

A sphere of radius 100 nm and density 2200 kg/m³ has its two centers held 100 nm
apart for 10 ms at rest in the same y frame. For prescribed rigid centers,

    D = G integral mu(df) T² sinc²(f T),
    G = integral spatial_measure(dk) |rho_hat(k)/m0|² [1-cos(k.d)].

This is exact Gaussian ensemble phase dephasing for the held-branch model. Traps
are assumed to impose those paths; it is not free propagation of the sphere.
The numerical G uses the exact sphere autocorrelation. A separate **analytic
lower bound** restricts the positive integral to |k|<=1/R:

- rho_hat/M >= 1-(kR)²/10 >= 0.9, by cos z >= 1-z²/2;
- 1-cos(k.d) >= (11/24)(k.d)², for d<=R;
- exp(-rc²k²) >= exp(-rc²/R²).

Angular/radial integration gives

    G >= (M/m0)² rc³/pi^(3/2) exp(-rc²/R²)
         * (0.9)² * (11/24) * d² * 4pi/(15 R^5) = G_lower.

Also sinc²(fT) >= [1-(pi fmax T)²/6]² on this band. Normalizing sigma² with
these two lower bounds proves D>=10 for **every** positive measure of the
specified total variance on [0.1,1], including atoms and arbitrarily narrow
features. No unobserved tail is silently permitted. The nominal uniform-band
witness produces more than 10 because the spatial certificate is conservative.
Only elementary floating-point evaluations of explicit inequalities remain;
this is not an interval-arithmetic certificate of arbitrary numerical integrals.

## Mechanics and finite observation

For a generalized force source g=-partial_q rho, define

    Q0=(hbar/m0)² integral g(dx)g(dx') exp(-|x-x'|²/(4rc²)).

Then a stationary source has force covariance Q0 integral cos(2pi f lag) mu(df).
We reuse PR #12's complete assembly upper bound Q0_upper, including load, flexure,
magnet, bounded glue, and all cross terms by the triangle inequality. Its mode
and glue assumptions remain assumptions. It supplies an upper bound, not a lower
calibrated force response that could exclude a field.

For n samples spaced dt with Blackman window w, a real line of variance A at nu
has expected one-sided periodogram

    A dt/sum(w²) [|W(f-nu)|²+|W(f+nu)|²].

The code evaluates the finite geometric sums exactly (to floating precision),
retaining aliases, the negative-frequency image, and symmetric/periodic windows.
Since sum|Blackman Fourier coefficients|=1, and every shifted sideband is separated
from zero and Nyquist, |W|<=1/sin(pi dt delta). This provides the saved conservative
bound valid between grid points, for any positive measure in the entire band.

A stronger robustness requirement is also reported. For **any fixed orthogonal
sample-space projection** P (including ordinary least-squares detrending), use
||P(w exp(-i2pi f t))||_1² <= n ||w||_2². The expected one-sided PSD per input
variance is at most 2 n dt. This bound is much looser, but still leaves the sensor
redundant here. An adaptive fitted transfer or arbitrary nonorthogonal control
filter is not covered without its gain. The analysis uses this projection bound
as its primary conservative sensor comparison.

For stationary LTI mechanics chi(f)/chi(0)=[1-(f/f0)²+i f/(f0 Q)]^-1. Over 0–1 Hz
its squared magnitude is <=[1-(1/f0)²]^-2 for all Q. Thermal calibration gives

    C_phiF = B1 omega0/(4 kB k),

where B1 is the measured slope versus T/Q in the workbook. This translates the
force response to the **measured raw flux PSD**, avoiding a pointwise comparison
to the earlier subtracted force envelope. We apply a tenfold PSD gain stress.
This still assumes the thermal slope calibrates the low-frequency driven mode,
with stationary initial conditions and no unmodeled control response. The archive
does not establish a universal DC gain bound. The 14 spectra include differing
temperatures and averaging counts, not 14 independently calibrated geometries.

The added field's predicted rms displacement/k is tiny relative to rc, checking
its own linear-response backaction. This does not certify all ordinary thermal
motion or the massive benchmark trap. An uncalibrated trap's backaction cannot be
inferred from this sensor calculation.

## What is and is not identified

In the selected band the joint sufficient region equals the sodium sufficient
region. Removing the sensor changes nothing; removing sodium removes the tight
norm certificate. These are not acceptance sets. No spectrum is proven excluded
by either source, and no dual majorant of the benchmark using calibrated lower
experimental responses has been obtained. Calling this complementarity would
repeat exactly the error the mission sought to avoid.

Without an independence/nonnegative-background premise even a measured total
PSD need not upper-bound a component: F_field=Z, F_background=-Z+epsilon X leaves
arbitrarily little total noise. The synthetic covariance check explicitly
constructs this PSD-positive ordinary correlated-background adversary. The
paper's noise budget may justify a restricted background model, but the released
spectral marginals alone do not establish it.

The useful boundary is therefore specific: the full quantum bound is adequate
for a broad *low-frequency Earth-following ensemble* class, but does not transfer
to the moving SSB witness and does not furnish an empirical joint likelihood.
A certified moving quantum response, or a calibrated complementary low-frequency
mechanical/atomic record with independently bounded preparation and response, is
needed to progress. This is not a theorem that no existing archive could work.
