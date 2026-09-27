# Physical derivation, statistical scope and escape construction

## One fixed reference model

Use SI units, reference mass m0 = atomic mass constant, spatial Fourier transform
mu_tilde(k) = integral rho(x) exp(-i k.x) dx, and the white CSL generator

D[rho] = -lambda rc^3/(2 pi^(3/2) m0^2)
         integral exp(-rc^2 k^2) [mu(k),[mu(-k),rho]] d^3 k.

A rigid point mass at two stationary positions separated by d has
-log(coherence ratio) = lambda (m/m0)^2 T [1-exp(-d^2/(4rc^2))].
The Gaussian stochastic-potential representation reproduces ensemble
decoherence and diffusion; it alone does not prove individual definite outcomes.

We define two-sided angular-frequency PSD by
C_F(t) = integral S_F^(2)(omega) exp(-i omega t) d omega/(2 pi).
The corresponding one-sided PSD per Hz is S_F^(1)(f)=2 S_F^(2)(2 pi f).
For a rigid body,
S_F^(2)/lambda = hbar^2 rc^3/(pi^(3/2) m0^2)
 integral exp(-rc^2 k^2) k_x^2 |mu_tilde(k)|^2 d^3 k.

The 2020 paper labels the same printed Fourier coefficient as one-sided.
Our derivation keeps the two-sided/one-sided conversion explicit. We do not
silently identify a printed rate normalization with the common generator.
The workbook force limit itself is reproduced independently of that conversion.
For a neutron reference instead of the atomic mass constant, multiply the
coefficient by (m_u/m_n)^2. This small mass convention does not remove a factor 2.

## Interference baseline

The deposited optics use three perfect reflecting loss gratings. At each mass
and velocity, the mean and first harmonic are S0 and S1=B_-1 B_2 B_-1.
Our vectorized implementation evaluates the same Bessel functions independently.
The author-module oracle verifies both optics and finite-sphere damping.

The observed fundamental contrast is
V = 2 |sum_m,v w_m,v S1(m,v) R(m,v)| / sum_m,v w_m,v S0(m,v).
It is not the average of per-particle visibilities. The first harmonic is signed.

For a symmetric point-particle Talbot-Lau geometry, each flight lasts t=L/v.
The separation grows linearly to d_max=h t/(m grating_period), then closes.
Integrating the point-mass decay rate gives
-log R/lambda = 2 t (m/m0)^2 [1-sqrt(pi) erf(x)/(2x)],
x=d_max/(2rc).
The small-x series in the code avoids cancellation. Time-domain quadrature is
an independent check. The numerical white diagnostic uses rc >=100 nm, much
larger than the sodium-cluster radius; it is explicitly a point-particle
approximation. It is not used to claim a confidence limit.

The deposited macroscopicity calculation instead uses a homogeneous sphere at
density 971 kg/m³ and its 200-point Fourier integral. We reproduce that
calculation, the 400-point tau_e grid, empirical phase -2.48 rad, approximate
Jeffreys prior and author quantile convention. The mapping is
rc=(hbar/sigma_q)/sqrt(2), lambda=(m0/m_e)^2/tau_e.
Its posterior is conditional on incident rates, phase and optics, and is not
a simultaneous frequentist confidence region.

Even a positive mixture is nonlinear: for weights (1/2,1/2) and response
coefficients (1,3), V(lambda)/V(0)=(exp(-lambda)+exp(-3lambda))/2.
The negative logarithm is not homogeneous in lambda. Signed mixtures may also
exhibit cancellation changes. No visibility growth was found on the committed
white grid; the warning concerns the proposed general spectral reduction.
A single linear positive kernel for log of the averaged visibility has therefore
not been derived.

## Static sphere and layered geometry

For a homogeneous sphere of radius R and mass M,
G(rc,d)=4/sqrt(pi) (M/m0)^2 integral_0^infinity
 z^2 exp(-z^2) [3 j1(z R/rc)/(z R/rc)]^2 [1-sinc(z d/rc)] dz.
The stationary white exponent is lambda T G. The code integrates to z=10.
The discarded absolute tail is at most
8/sqrt(pi) (M/m0)^2 (5+1/40) exp(-100).
Quadrature is numerical; the actual physical response is not interval-certified.

For a rectangular transverse side a,
I0(a)=integral exp(-rc² k²) |2sin(ka/2)/k|² dk
=2pi[a erf(a/(2rc))+2rc/sqrt(pi)(exp(-a²/(4rc²))-1)].

For layers normal to z, integrate k_z² |mu_z(k_z)|² using the density jumps
Delta rho_i at positions z_i:
Iz=sqrt(pi)/rc sum_i,j Delta rho_i Delta rho_j
exp(-(z_i-z_j)²/(4rc²)).
The result times I0(a) I0(b) and the CSL prefactor is the isolated-load response.
The cantilever, sphere, elastic mode shape and cross terms are not included.
Dropping them is not automatically a rigorous conservative bound.

For a cubic LPF mass, the force and torque integrals are evaluated through I0,
I2 and derivatives of the top-hat Fourier transform. In particular,
I2=2sqrt(pi)/rc [1-exp(-a²/(4rc²))],
J=integral exp(-rc² k²) k H H' dk = -I0/2+rc² I2.
The derivative-square integral follows independently from the real-space
autocorrelation:
Ider=2sqrt(pi)/rc integral_0^a exp(-u²/(4rc²))
[a³/12-a²u/4+u³/6] du.
Thus torque/force = 2Ider/I0-2J²/(I2 I0), tending to a²/6 for rc/a -> 0.

For uncorrelated identical masses, the one-sided PSD of half their force
difference equals the single-body two-sided PSD. The same holds for the
half torque difference, consistent with I² S_Delta_gamma/4.
Spatial correlations at 0.38 m separation are negligible at the tested rc,
but the motion through a preferred colored-noise frame remains an open premise.

## Colored temporal responses and the missing bridge

Declare f_tau(u)=exp(-|u|/tau)/(2tau), with spectrum
S_tau(omega)=1/(1+omega² tau²), and fix lambda as its zero-frequency amplitude.
There is no redundant second free overall spectrum amplitude.

For a stationary spatial superposition the temporal integral is exactly
H(T,tau)=T-tau[1-exp(-T/tau)], with H(T,0)=T.
This supplies the finite-time exponent lambda G H, not simply lambda G T.

A stationary mechanical spectral model gives the white force coefficient times
S_tau(omega). The selected feasibility calculation uses the quoted positive
frequencies as spectral-envelope constraints, not as a substitute for the
instrument's complete finite-time likelihood.

For prescribed paths xa(t), xb(t), the exact Gaussian phase-average kernel is
K(omega)=(m/m0)^2/2 integral dt ds cos(omega(t-s))
[C(xa(t)-xa(s))+C(xb(t)-xb(s))-C(xa(t)-xb(s))-C(xb(t)-xa(s))],
where C(x)=exp(-|x|²/(4rc²)). Positivity follows by writing it as the
spatial Gaussian-weighted squared modulus of a temporal Fourier amplitude.
The code checks the static-path limit independently. It does not identify these
two prescribed paths with the full Talbot-Lau quantum evolution.

Motion matters: for locally parallel straight paths at speed v, spatial overlap
contains exp[-v²(t-s)²/(4rc²)]. Its exact OU overlap integral is
sqrt(pi) rc/(v tau) erfcx(rc/(v tau)).
At v=158 m/s, rc=100 nm and tau=1 microsecond it is about 0.00112,
while H(10 ms,tau)/(10 ms) is about 0.9999.
This large difference rules out using a stationary holding-time response as a
colored sodium response. This local calculation is not itself the missing map.

More generally, a translating body's spectral response contains
S_tau(omega-k.v) inside the spatial Fourier integral, rather than a common
S_tau(omega) factor. The Earth laboratory and LPF cannot both be presumed at
rest in an arbitrary universal noise frame. Consequently, the long-memory
sequence below is a theorem about stationary-response spectral envelopes,
not an established escape from the full moving-apparatus experiment.

## Constructive limitation of finite-positive-frequency combinations

Fix any stationary benchmark with G>0, duration T>0, required exponent C>0.
Choose lambda_tau=C/[G H(T,tau)]. The benchmark is met exactly.
Since H(T,tau) ~ T²/(2tau), lambda_tau/tau -> 2C/(G T²).
For each omega_j>0,
lambda_tau/(1+omega_j² tau²) -> 0.
Any finite set of nonzero-frequency upper envelopes with positive right-hand
sides is therefore eventually satisfied. The covariance variance
lambda_tau/(2tau) remains finite, tending to C/(G T²).

The statement extends to finite weighted spectral constraints whose nonnegative
kernels are supported above a positive omega_min and have finite integral:
bound S_tau there by 1/(omega_min² tau²).
It does not cover an instrument kernel with substantial zero-frequency response,
an observed static-superposition constraint, or moving-frame Doppler response.

This explains exactly what adding mHz data buys, and what it does not buy
under the stated approximation. It is not a universal no-go theorem for public
datasets or a viable complete ontology. A cutoff on tau, a bound on static noise,
or a calibrated coherence experiment can remove this escape.

## Radiation runner-up

The code implements the 2026 XENON supplement's charge-cancellation expression
and orbital radii. Frequencies are omega=E/hbar, not E/h or energy in keV without
conversion. Alpha=1 and 1.5 quantify the quoted intratomic-distance approximation.
The old public release has 433 selected events under 30 keV; under a Poisson
model, 0.5 chi²_(2(n+1),.95) is an upper bound on total mean counts. Nonnegative
background makes it conservative for signal counts without subtraction.

Turning that into lambda uses 1.16 tonne year, efficiency versus reconstructed
energy, and an identity energy-response map. Migration and calibration are not
provided, so this is a quantitative feasibility calculation, not a certified
parameter exclusion. Markovian limits are never reused unchanged for colored
noise. A full detector response would fold the energy-dependent signal before
integration.

## Certificates and statistical separation

The finite LP returns primal/dual vectors and numerical residuals. The exact
rational checker proves domination only for explicitly declared piecewise-affine
functions: check every common knot, plus the slope of the infinite tail.
A missing tail, a bad middle knot or negative weights is rejected.
Nothing promotes a sampled physical kernel to a continuum certificate.

Combining two individually quoted 95% limits does not yield 95% joint coverage.
The selected result is a deterministic compatibility set. The independently
named Gaussian simulation uses Bonferroni bounds at overall alpha=.05 and checks
coverage/recovery under its own known-variance model. It does not validate the
experimental likelihoods. No claim is based on treating fitted summaries and
their source records as independent observations.
