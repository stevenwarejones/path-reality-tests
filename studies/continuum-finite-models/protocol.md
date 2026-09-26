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
confidence interval. `calibration_phase_radius=kappa` adds a contrast-bias
bound min(2,kappa²/2), since |1−cos(kappa)|≤kappa²/2. The 660-case
map uses its stated phase radius for calibration too, including twice this bias
in its prospective power envelope. The targeted shared-scale comparison records
its separate, exact calibration-phase premise explicitly.

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

## Executable count-to-exclusion rule

`decision.py` implements the rule above. The input JSON has `calibration` entries
`detection`, `phase0`, and `phase_pi`, each `[successes, total]`; `main_counts` is
a 24-row list of `[plus, minus, failure]`; `candidates` contains distinct integer
site counts. Optional `certificates` supplies `phase_offset`, `fractional_scale`,
`contamination`, `efficiency_transport`, `contrast_transport`, and
`calibration_phase_radius`. Defaults are
zero certificate radii and the risk allocation above, **not evidence that real
errors vanish**. All 24 settings occur in the order returned by
`decision.DEFAULT_SETTINGS`. The detection calibration's success is any click;
the phase calibrations' success is plus, with all trials in the denominator.

```sh
python studies/continuum-finite-models/decision.py counts.json --output exclusions.json
```

The function intersects simultaneous calibration intervals with 0≤w≤eta≤1,
rather than clipping noisy point estimates and assigning them spuriously small
uncertainty. It encloses that physical intersection by a rectangle. If the
intersection is empty it reports `inconclusive-calibration` and excludes no
candidate. This conservative refusal can reduce power off the calibration
coverage event; the 87.50% guarantee already charges that event's risk.

Each candidate is rejected only for a **strict** violation of an enlarged null
coordinate interval. A value exactly on the threshold is retained. The output
includes every rejected/retained N and witnesses naming setting, outcome,
frequency, interval and statistical radius. Failure is tested like either click
outcome. Unequal stratum totals use their own radii with the same 72-coordinate
penalty; they must still come from a fixed acquisition plan. Candidate count
does not change multiplicity. End-to-end tests drive this actual function with
adverse calibration/main fluctuations within the promised coverage events at
the advertised integer budgets, and also check a command-line round trip.


## Shared scale and phase across momenta

`joint.py` keeps a single offset and kinetic-scale multiplier across all settings.
It covers the calibrated parameter rectangle by interval cells. Every discarded
cell has a setting/outcome whose entire predicted interval misses the simultaneous
observed confidence interval. Rejection requires covering the whole rectangle;
a cell-budget limit or unresolved point returns inconclusive. This strengthens
the independent-box rule using the same coverage event and the same 5.00% risk.
No penalty is needed for searching cells or testing additional candidate N.
Efficiency/contrast rectangles are still conservative enlargements.

The prospective certificate covers **both** continuum and lattice nuisance
families in four dimensions. A cell is discharged only when one coordinate has
a gap larger than the proposed bound throughout that cell. This is continuous
interval coverage, not a grid minimum or an optimized gap. Endpoints are padded
outward and gaps below 1e-12 are not certified; floating results remain numerical
certificates rather than Lean proofs.

The certificate also saves a finite null-parameter partition. Taking the union
of its cut locations partitions the null rectangle more finely than every
certificate cell. For the true continuum parameter, every null subcell therefore
inherits a separating witness. If r_alpha+r_beta<g and calibration coverage
holds, the count decision rejects every subcell. Use this predeclared partition
for the advertised power; a capped adaptive search alone promises only sound
rejection, not termination within its budget. Tests exercise this distinction.

```sh
python studies/continuum-finite-models/joint_design.py --check
python studies/continuum-finite-models/decision.py counts.json --shared --output exclusions.json
```

For guaranteed power, the count JSON also includes `partitions`, mapping candidate
N (JSON string keys) to the two sorted cut lists in `decision_partition` from the
matching design certificate. The code checks exact coverage of the declared
parameter bounds. The design conditions, including calibration phase and transport
radii, must match the acquisition certificates.

At eta=0.8, visibility=0.9, calibration radius 0.001 and fractional scale radius
0.001, the targeted comparison certifies gaps 0.01 for N=128 and 0.001 for
N=192,256 at zero offset radius. With offset radius 0.005 it certifies N=128,192;
N=256 remains inconclusive at the searched thresholds/budget. All use the same
24 settings; the largest grid N=512 remains inconclusive. These are improved
sufficient boundaries, not exhaustive identifiability limits.

For one momentum, exact scale overlap is possible: writing d=E_a−E_c, choose
s_c=d/(2E_c), s_a=−d/(2E_a). Then (1+s_c)E_c=(1+s_a)E_a.
For j=3,N=128 both shifts are within 0.001. The tests verify equality for all
listed times/readout phases and show those same scales fail at j=1. Thus a
single-momentum degeneracy and its multi-momentum resolution are both checked.

## Calibration and main costs

Both maps report `total_calibration_trials` and `total_main_trials` separately.
For r_cal=0.001 the three calibration strata cost 8,220,960 eligible trials.
The shared-scale N=128 gap of 0.01 needs another 3,816,216 main trials; a gap of
0.001 needs 381,621,024 main trials. The calibration floor must not be interpreted
as dispersion sensitivity. Power is computed from the declared risk parameters,
1−alpha_calibration−beta_main, rather than stored as an unrelated constant.
