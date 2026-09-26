# Prospective protocol and conditional cost

No apparatus measurements enter this design. Values in the generated map use
L=2π, m=hbar=1. A dimensionless time t corresponds to physical time
`t*m*L_physical²/(4π²*hbar)`. Neither these times nor the synthetic efficiencies
are measured experimental performance.

## Fixed acquisition and null

Use the same 24 strata for every candidate integer N: pairs (0,j), j=1,2,3;
times 0.5,2,8,32; and phases 0,π/2 relative to the *continuum-calibrated* phase.
No post hoc choice of best setting is used for allocation or the rejection rule.
Prepare a coherent mode pair, evolve on a circle, apply a calibrated balanced
mode recombiner, and record plus, minus or failure on every eligible trial.
A successful preparation herald defines eligibility before readout; failed
preparations are separately counted for source cost. Mean-field interactions,
trap inhomogeneity and rotation require bounds before applying the free-particle
model. Preparation contamination can be charged as an extra TV radius in code;
the map assumes it is zero and makes that idealization explicit here.

For each N, the composite null consists of its complete outcome tables under
one common calibrated physical scale and phase convention, efficiency eta and
contrast w=eta*v. Independent interval relaxation produces a larger null box.
The alternative is the continuum under the same nuisance certificates.
Phase uncertainty is bounded by offset_radius + fractional_scale*|E_j t/hbar|.
The fractional scale bounds the combined t*hbar/(mL²), not independently chosen
kinetic energies at every setting. For fixed mode labels, geometry and mass
uncertainties enter through that scale; misidentified mode labels are excluded
by a separate preparation certificate. Bounds on visibility and efficiency
apply uniformly across all main strata; their transport is a premise.

## Rejection and coverage

For M=72 probability coordinates (3 outcomes times 24 settings), n independent
eligible trials per stratum, use r(n,M,a)=sqrt(log(2M/a)/(2n)). Reject a candidate
if at least one empirical probability lies outside its null interval enlarged
by r(n,M,0.025). Hoeffding and a union bound give simultaneous coverage without
requiring the three multinomial coordinates to be independent.

The same simultaneous box can reject any number of N values without an extra
candidate-count penalty: when one candidate is true its whole table lies in
that box. This is a family exclusion rule, not independent pairwise tests.
A robust gap g ensures rejection whenever r_alpha+r_beta<g. The code chooses the
smallest positive integer satisfying that strict inequality, with beta=0.1.
Zero g yields no promised budget, not a claim of non-identifiability.

Three separate calibration strata measure detection probability and plus
probability at known phases 0 and π. Each gets
ceil(log(6/0.025)/(2*r_cal²)) independent eligible trials. Then eta has radius
r_cal and w=p_plus(0)-p_plus(π) has radius 2*r_cal. This requires calibrated
baseline phases and the same eta,w at all main settings. Baseline phase error
or transport drift needs additional radii; it is not absorbed by a statistical
confidence interval. The synthetic map assumes exact calibration-phase settings.

For prospective power, calibration centers are random. The map doubles the
calibration radii around the stipulated true centers, to include both center
error and the interval itself. Its g is conservative for both model boxes.
Total type-I risk is ≤0.025+0.025=0.05. Unconditional power is at least
1-0.025-0.1=0.875, conditional on the deterministic external certificates.
This does **not** promise 90% unconditional power. The reported total eligible
count includes all 24 main strata and all 3 calibration strata.

## Full acquisition cost remains conditional

A real source budget additionally needs a certified preparation/heralding
success probability h, a fixed attempt cap with a binomial shortfall bound,
calibration acquisition time, clock/length/mass uncertainty and phase-reference
characterization. Dividing eligible trials by h gives an expectation, not a
high-confidence acquisition guarantee. `total_source_attempts` is therefore null.
The map deliberately supplies no invented source efficiency or photon mass.
Current output is a conditional eligible-trial guarantee and a missing-input
list, not an experimental exclusion or a complete apparatus budget.

Report the full disconnected candidate set passing/failing this rule. Do not
assume monotonicity in N or time; phase wrapping gives blind regions. Before
acquisition fix the candidate set, controls, calibration certificates, eligible
trial rule, multiplicity and stopping counts. A sequential/adaptive extension
requires different coverage arguments.
