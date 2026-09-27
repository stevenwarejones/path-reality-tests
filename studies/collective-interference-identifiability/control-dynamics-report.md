# Shared dynamics and the measurement needed to close the drift loophole

The proposed two-atom regional control does **not** separate the competing
classes in the expanded measured region. An exact physical four-component
cluster mixture has precisely the same control probability as the committed
quantum example, while passing all retained singleton, parity-pattern and
row-hit constraints. This is a stronger negative result than a failed search
for separation. The worst-case absolute prediction gap is exactly zero.

A matched three-atom acquisition with mutually orthogonal labels provides a
sufficient class-wide control. It measures the third-order transport event
directly and gives a factor-two cluster ceiling even for arbitrarily large
shared fluctuations. This is a conditional experimental design, not a new
measured exclusion. The apparatus-constrained inference remains unresolved.

## Decision table

| Route | Shared physical restriction and data | Result | Remaining gate |
|---|---|---|---|
| July 7 joint dynamics | One separable nearest-neighbor Hamiltonian, Gaussian common scale, time offset and survival, fitted to all 4,266 singleton shots across 15 times and three inputs | Local fit reconstructed; six singleton-cell intervals fail at three settings | This fitted point cannot support an inference about clusters. No certified exclusion of the dynamics family, or adequate joint fit, is claimed |
| Published detailed apparatus model | Pinned spectroscopy, Hamiltonian, unitary and resampling artifacts are available | Inventory and calibration-frame provenance checked; later calibration cannot simply be assigned to July | July confidence region or validated transfer of alignment, timing, potential and noise |
| Hamiltonian-scale robustness | Finite spectral bandwidth W, Gaussian scale SD sigma, fixed label-independent detection, selected mean calibration | Exact sufficient certificate for t W sigma / 2 <= .00925 | Product bound is conditional, not measured; mean input translation remains assumed |
| Existing pair-control proposal | Full expanded selected region and arbitrary shared physical channel laws | Exact zero worst-case gap to the quantum example | This binary pair statistic is not decisive, even with unlimited precision |
| Full low-order information | Explicit loss-dilated isometries under common latent laws | All singleton and pair distributions can agree while the three-atom event differs | Structural result, not a claim that every pair experiment fails on this archive |
| Matched labelled triple | Same three inputs and channel law, three mutually orthogonal labels | P_C2(B) <= 2 P_labelled(B), uniformly over channel laws and pair partitions | Validate label-independent transport and common acquisition law; collect the control |
| Source removal | Existing expanded source region plus a prospective labelled-control upper bound | Physical witnesses after removing either the original many-body source or the new control | Prospective complementarity only; singleton data are redundant for this new inequality |

## 1. What can legitimately be shared

The primary paper describes a spatial lattice Hamiltonian and a Gaussian
scale-noise model. The detailed public implementation additionally supplies
spectroscopic fit parameters and uncertainties, tunneling conversion data,
Hamiltonian arrays, and calibrated unitaries. Those are real apparatus
inputs, not missing in their entirety. Their existence does not establish
an uncertainty region for the selected July acquisitions.

The hash-pinned `inference/H_resampling.py` identifies an August 5 reference
coordinate frame and August 22–23 alignment handling. The notebook explicitly
records changes in lattice alignment, including mirror adjustments on August
12 and 17. The July 2 selected run, July 2 4.65 ms run, and July 7 scan are
therefore kept as separate date groups. Sharing July 7 parameters across its
scan is a declared within-date stationarity hypothesis, not a verified fact.
No August fit uncertainty is reinterpreted as July shot-to-shot noise.

`results/apparatus-audit.json` inventories the six relevant archive members,
their hashes and array dimensions, and independently checks the three input
locations in every July 7 singleton setting. Source reconstruction checks
these against the already pinned archive and notebook. The potentially
pickled spline file is inventoried without unpickling it.

| Quantity | Candidate sharing | What is currently justified |
|---|---|---|
| Lattice geometry and Hamiltonian form | All dates | Model structure is documented; numerical coefficients and alignment are not thereby fixed |
| Jx, Jy, effective time offset, survival | July 7 scan | Fitted jointly as an explicit diagnostic hypothesis |
| Shot-scale Gaussian law | Preparations within a date | A physical model to test, not a measured noise-confidence interval |
| Spectroscopic potential and coordinate registration | August to July | No validated transfer bound established |
| Selected mean translation | Three selected inputs | Retained conditional assumption; the new Gaussian envelope does not eliminate it |
| Imperfect labels, interactions, false clicks, preparation selection | Across a new control and original data | Must be validated separately; not absorbed into a free epsilon |

Primary reference: [Young et al., Methods I.7 and Supplement II/VI/IX](https://arxiv.org/html/2307.06936v2).
The local archive provenance, rather than the paper alone, determines the
calibration-transfer limitation above.

### Joint July 7 fit: a concrete failed candidate

`dynamics.py` fits all 45 independent singleton distributions with

    H/h = -Jx A_x - Jy A_y,
    effective time = (nominal time + offset) (1 + sigma Z), Z ~ N(0,1).

The infinite-lattice transition probabilities use squared Bessel functions;
outside-crop probability is included in loss. The same five parameters apply
at all 15 times and all three input sites. Loss is fixed, independent per
atom and label independent. The Gaussian integration is common within each
shot; products are averaged after multiplying, not products of averages.
The fit uses singleton counts only. Many-body bunching predictions are then
computed for comparison, without fitting the many-body outcomes.

| Parameter | Local fitted proposal |
|---|---:|
| Jx/h, with Jx interpreted as the energy coefficient | 121.9462 Hz |
| Jy/h | 106.5812 Hz |
| Additive nominal-time correction | -0.374282 ms |
| Relative Gaussian scale SD | .0811700 |
| Survival | .962720 |

The code variables Jx/Jy are already in Hz; the displayed h notation makes
the physical energy convention explicit. The fitted scale SD is **not** a
measured fluctuation level. It can compensate for omitted confinement,
longer-range hopping, alignment or model misspecification.

The singleton multinomial deviance is 3369.897879, reported descriptively
without a chi-square p-value. More directly, six existing simultaneous cell
intervals fail at .714286, .928571 and 1.142857 ms. This candidate is therefore
not compatible even with the declared conservative calibration region. Its
C2 bunching deficit cannot be attributed to particle indistinguishability:
the underlying transport proposal already fails. Local fitting failure does
not exclude every parameter choice or a more faithful Hamiltonian model.
The 61- versus 81-point quadrature discrepancy is below 1e-10; this numerical
convergence check is not an exact physical-feasibility certificate.

## 2. Gaussian noise without an invalid ess-sup conversion

A Gaussian has unbounded support, so its SD must not be substituted for the
earlier essential-supremum TV budget. A moment bound avoids that error.
Let H0/hbar have spectral width W on the full relevant mode space, and
H(S)=S H0 with S Gaussian of SD sigma. Initial states and the detection/loss
map are fixed and label independent. Set kappa=t W sigma/2.

For two independent draws S,S', removal of a global phase and the unitary
perturbation bound give

    TV(p_j(S), p_j(S')) <= (t W/2) |S-S'|.

Convexity relative to pbar_j=E p_j(S), followed by Jensen, gives for
`d(S)=sum_j TV(p_j(S),pbar_j)`:

    E d^2 <= 9 (t W/2)^2 E |S-S'|^2 = 18 kappa^2,
    E d^3 <= 27 (t W/2)^3 E |S-S'|^3 < 135 kappa^3.

Here the Gaussian absolute third moment is `8 sigma^3/sqrt(pi) < 5 sigma^3`.
Substituting into the previously proved signed multilinear envelope gives

    P_C2(B) <= 2 [U + 12 gamma kappa^2 + 20 kappa^3].

With U=.0112 and gamma=.2881834, the exact rational checker obtains
`.0230234427369 < .023025033926510948` at kappa=.00925. This is sufficient
uniform exclusion under the listed assumptions. It is not an optimized
threshold or an estimate of the experimental error. A truncated crop cannot
be used to bound W unless coupling to omitted modes is controlled. Mean
translation, preparation selection, interactions and detection errors remain
separate assumptions. The earlier TV bracket is unchanged.

## 3. Exact evasion of the proposed pair statistic

The four `pair-{above,below}-{minus,plus}.json` components are exact physical
contractions with individually translated probability columns. Each of two
fair mixtures passes the full expanded selected region. Their pair-event
probabilities bracket the quantum target:

    .1432000624053198 < .1432359739972339 < .1432360689338087.

For endpoint probabilities p_minus and p_plus and target q, mix the two laws
with the exact rational weight `(q-p_minus)/(p_plus-p_minus)` on the upper
law. Every constraint is convex in the acquisition law, including the
conditional singleton-cell constraints after multiplying by row mass.
The checker also directly rechecks the mixed calibration and every parity
and row-hit constraint; it does not rely solely on that argument.

| Quantity | Exact-matching cluster example |
|---|---:|
| Original three-atom bunch probability | .0230750278052 |
| Proposed pair-control probability | Exactly equal to the quantum target |
| Shared fluctuation budget, sum of three column TVs | .936055521147 |
| Labelled three-atom bunch probability | .0118215439756 |
| Complete-pattern deviance at the fitted point | 17944.060078 |

This model has substantial fluctuations and no demonstrated lattice
Hamiltonian realization. The deviance has no calibrated reference
distribution here; it underscores that the deliverable is **compatibility
with this region**, not aggregate statistical adequacy. The new model does
not improve the old minimum-drift upper bracket .898524.

For the proposed binary pair statistic, the optimization
`inf_(surviving C2 laws) |P(control)-q|` has value exactly zero: nonnegativity
is the lower bound and the physical mixture attains it. No global search for
the interval endpoints is needed to settle this question. Since the quantum
class includes the target example, this control cannot separate the complete
classes. Other pair statistics or their joint use may still constrain them.

## 4. Why pair information need not determine triple information

Take three detected sites in one row, one assigned to each input, plus three
input-specific loss modes. A latent bit z_j routes input j deterministically
to its detected site if z_j=0 and to its loss mode if z_j=1. Every component
has orthonormal columns and an isometric physical realization.

Compare a uniform law on the four even-parity bit strings with a uniform law
on the four odd-parity strings. All one-input and two-input routing
distributions agree, including loss, so every corresponding singleton and
pair observable agrees. But the probability that all three detected sites
are occupied is 1/4 versus 0. The exact checker enumerates the full laws.
Both examples are within the distinguishable subclass of C2. This is a
structural obstruction, not a fit to the atomic archive, and does not prove
that every pair-based witness is powerless for the selected region.

## 5. One sufficient added statistic, uniformly over fluctuations

Repeat the selected experiment with all three atoms made mutually
orthogonally labelled, while preserving their visible inputs, the channel
law after initial selection, and label-independent transport/detection.
Measure the **same** event B: three distinct detected cells in one selected
row. Do not condition on survival. A time-labelled construction from separate
shots is not equivalent.

For every channel and each C2 partition, the two-term Cauchy-Schwarz bound
on assignment amplitudes gives `P_C2(B|Z) <= 2 D(B|Z)`. The fully labelled
control measures `E D(B|Z)` directly, so

    P_C2(B) <= 2 P_labelled(B).

The partition may itself correlate with Z. Averaging retains the inequality.
This is the established cluster factor applied to a matched third-order
transport measurement, not a claimed new resource inequality. It needs no
bound on fluctuation size, no singleton translation symmetry and no common
Hamiltonian fit. It does need the physical matching conditions above.

The quantum channel predicts `.0046670348373` for this control. Every C2
model passing the original bunch lower bound requires at least
`.0115125169633`. The new pair-evading example predicts `.0118215439756`.
Thus the relevant separation is now class-wide, conditional on an observed
control upper bound, rather than separation of two fitted examples.

A concrete fixed-size protocol uses **3,000 prepared control trials**. If
there are at most **19** B events, `.0115` is an exact conservative one-sided
binomial upper endpoint with tail 1/320. Twice this endpoint is below the
existing lower bound. The two event bounds alone have total error at most
1/160 under IID sampling and the declared event/source selection conditions.
This is a separate prospective statement; it does not silently extend the
old n=2–5 or 16-setting confidence families. Broader exploratory selection
still needs its own coverage treatment. No acquisition has been performed.
A practical scheme for three independent labels with unchanged transport
remains to be validated in this apparatus.

A future acquisition is not automatically matched to the archived July 2022
law. Without a defensible transfer guarantee, the operational protocol must
randomly interleave fresh original and labelled preparations and re-estimate
both bounds. The 3,000/19 calculation is conditional on retaining the stated
original lower bound; it is not permission to compare a later control to an
unverified historical channel.

This is the smallest **sufficient control exhibited here**: one binary
statistic from one additional preparation family. It is not a proof of global
experimental minimality or a statement that three atoms are necessary for
every possible class-separating protocol.

### Source removal for the prospective expanded design

Treat the prospective control region as `[0,.0115]`, not invented observed
counts. Exact physical witnesses establish:

* **Remove the new control:** the four-component C2 mixture passes every
  existing selected singleton, complete-pattern and row-hit constraint.
* **Remove the original many-body record:** use the existing quantum channel
  with pair-plus-singleton labels in the original preparation and fully
  orthogonal labels in the control. It passes the retained singleton region
  and the new control upper bound. The original many-body constraints have
  actually been removed, not replaced by bunching alone.
* **Keep both:** the class-wide inequality excludes C2. The original bosonic
  model passes the existing region and the prospective control upper bound.

Singleton data are **not indispensable** for this new inequality. The payoff
would combine the original many-body source with a new matched control, not
prove that every available source is necessary. For the existing expanded
fixed-channel region alone, the missing higher-only deletion witness remains
open; failed restricted searches are not presented as a class exclusion.

## Reproduction and next physical gate

```sh
python studies/collective-interference-identifiability/control_design.py --check
python studies/collective-interference-identifiability/dynamics.py --check --cache /tmp/collective-sources
```

The optional `build_control_design.py` and `dynamics.py --fit` regenerate
numerical proposals. Exact certificate checking never requires the optimizer
to reproduce its search path. The old robustness artifacts remain unchanged.

The next empirical target is still to determine whether a justified July
Hamiltonian/noise confidence set contains a cluster explanation. Required
inputs are date-matched potential/alignment and timing/noise information, or
an explicitly validated transfer from the available later calibration. If
that bridge cannot be established, the matched labelled-triple protocol is a
concrete alternative that measures the missing transport moment directly.
Neither an apparatus-supported exclusion nor an adequate shared-Hamiltonian
countermodel has yet been established. Priority and broader significance
remain unclaimed.
