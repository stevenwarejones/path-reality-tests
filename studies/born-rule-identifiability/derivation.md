# Exact optical observation map and added measurement

All intensities are ideal linear detector-equivalent volts. The measured-summary
instantiation treats estimated means as exact: its intervals are not confidence
regions. Detector calibration is required for a physical exclusion.

## Hidden coherence and an ordinary phase explanation

Let H be the real positive definite coherence matrix from background b and six
lower settings: H_ii=I_i−b, H_ij=(I_ij−I_i−I_j+b)/2. Define
G(u)=H+i u(e_A e_Bᵀ−e_B e_Aᵀ). In this restricted completion family,
G is PSD exactly when |u|≤U=sqrt(det(H)/H_CC). Schur-complement out positive
H_CC: the remaining Hermitian 2×2 matrix has fixed positive diagonals and
determinant det(H)/H_CC−u². This proves both necessity and sufficiency, including
endpoints. The implementation requires H positive definite.

Any PSD G has an ordinary optical realization as an incoherent mixture of
coherent path fields using its spectral decomposition (with conjugation chosen
to match the v†Gv convention). Partial coherence is a substantive nuisance
premise, not a statement that this is the published apparatus state.

For mask s use its bit vector, except in all-open context where
v=(exp(iδ),1,1). This is a one-parameter shutter-context interaction: A acquires
a phase when B and C are both open. It is **not** the paper's symmetry-reduced
first-order crosstalk model. Other seven intensities are b+sᵀHs. No parameter is
fitted separately to every cycle or outcome.

The all-open deviation from q0=1ᵀH1 is exactly
D(u,δ)=2R(cosδ−1)+2u sinδ, R=H_AB+H_AC.
A diagnostic adds θ S to the same cell, S=tr(H). Thus (θ=0,u,δ) and
(θ=D/S,u,0) have identical eight mean readings. θ is limited by positivity in
every used measurement; it is not a complete probability-law modification.
With the same stipulated additive noise law these model distributions could
also coincide, but no such temporal noise fit is asserted for the real data.

## Sharp calibration region and local check

On |u|≤U, |δ|≤d≤π/2, D is continuous on a connected compact domain, so its
image is a closed interval. Extremize over u at ±U and, by symmetry, δ in [0,d].
For each sign evaluate 0,d and every stationary angle satisfying
−R sinδ±U cosδ=0. This is the analytic envelope in `phase_envelope`; continuity
proves every intermediate value is attained. It is sharp only for this family.

If ε=I_ABC−b−q0 and envelope=[L_d,U_d], the allowed diagnostic interval is
[(ε−U_d)/S,(ε−L_d)/S]. Shared θ across runs with independent u,δ gives interval
intersection. Opposite residual signs can empty the zero-phase intersection;
that rejects the restricted shared diagnostic, not quantum mechanics. Radii
0,.01,.03,.1 rad are exploratory scenarios, not measured calibration or priors.
H uncertainty is not included, so no empirical exclusion follows.

The explicit archive pairs use u=U/2 and the smallest-magnitude root δ in
[−.5,.5]. Complex matrix evaluation independently checks the resulting table.
At δ=0, θ and δ derivatives are S e_ABC and 2u e_ABC; whitening by any common
positive definite covariance preserves their dependence. The code uses measured
cycle covariance only as a metric. The nonlinear equality, not the rank check,
establishes ambiguity. Peres F<1 and allowed mixed coherence are known facts.

## One extra phase setting

Keep preparation, transmission and δ fixed and apply additional A phase φ.
The diagnostic is assumed phase-independent. For
W=[I_ABC(φ)+I_ABC(φ+π)]/2−I_A−I_BC+I_empty,
antipodal averaging cancels all cross terms involving A. The remainder is
G_AA+G_BB+G_CC+2 Re G_BC, hence W=0 ordinarily and W=θ S for the diagnostic.
Background cancels because the coefficients sum to zero. The original pair is
identical and this ninth setting differs: one extra setting suffices and zero
cannot distinguish this pair. This is not a universal minimum or tomography test.
Amplitude-changing crosstalk, BC changes or actuator side effects need not cancel.

A relative increment π+e, |e|≤e_max≤π gives
|W−θ S|≤2 sqrt(G_AA q_BC) sin(e_max/2), from PSD Cauchy–Schwarz.
Additional independently certified bias ≤r per reading adds 4r, the coefficient
L1 norm. r must include detector, leakage, gain/exposure and uncancelled drift;
it cannot be inferred by fitting away the signal.

For synthetic controls let a=|G| times maximum common power, ell bound fixed
leakage amplitude and q bound quadratic detector response. Leakage affects
A, BC and dark only; triangle inequality bounds W by
4 ell(a_AB+a_AC)+ell²(sum(a_BC block)+a_AA+sum(a)). Each field coefficient has
modulus ≤1, so intensity ≤sum(a), and response adds 4q(sum(a))². This is an
analytic synthetic envelope, not an archive calibration.

## Conditional acquisition cost

With n independent readings in each of five settings, common known range width M,
and coefficients (.5,.5,−1,−1,1), sum of squared coefficients is 3.5. Hoeffding:
Pr(|W_hat−E W_hat|≥t)≤2 exp(−2nt²/(3.5M²)). Under null |E W_hat|≤B, reject at
|W_hat|>B+M sqrt(3.5 log(2/α)/(2n)). For signal s and desired power 1−β, a
sufficient n is ceil[3.5M²(sqrt(log(2/α))+sqrt(log(1/β)))²/(2(s−2B)²)] when
s>2B. Opposite systematic shifts of null and alternative cost **2B**, not B.
Total acquisitions=5n. This is pointwise for a prespecified design. Add calibration
failure probabilities if bounds themselves are statistical.

The observed 40/50 repeated photodiode readings are not assumed independent
trials. Gaussian simulations have separately specified known noise and common
block drift; they do not validate arbitrary temporal dependence or HAC coverage.
