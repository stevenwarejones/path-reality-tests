# Physical transport and a fresh matched-control protocol

This follow-up separates three questions: whether a deposited Hamiltonian explains the calibration records, whether three labels can be implemented without changing the experiment, and what a fresh two-arm experiment would cost. It does **not** establish an apparatus-supported cluster exclusion or a validated laboratory protocol.

## Existing records: keep the dates separate

`extended_dynamics.py` adds linear alignment gradients, quadratic confinement and axial next-nearest-neighbor hopping to the nearest-neighbor lattice. `diagonal_dynamics.py` additionally allows diagonal hopping:

\[
H/h=H_y\otimes I+I\otimes H_x-J_d A_y\otimes A_x,
\quad H_a=\operatorname{diag}(F_a c+K_a c^2)-J_a A_a-J_{2a}B_a.
\]

A common Gaussian time-scale fluctuation multiplies each shot's evolution time; a fitted time offset and independent survival probability complete the model. Coordinates are centered at the physical reference site. July 2 uses **1,624 unique singleton shots in four records**, including the 931-shot record only once, and July 7 uses **4,266 shots in 45 records**. No parameters transfer across dates. The source audit reconstructs the July 2 record from the pinned archive.

The terms are motivated by the confinement, alignment and non-nearest-neighbor effects discussed in the [source experiment](https://arxiv.org/html/2307.06936v2). Their search bounds are exploratory, **not measured parameter confidence regions**. The axes use quadratic local confinement rather than a fully calibrated optical potential. Interactions, spatially varying hopping and detection backgrounds are not fitted. The separable search uses three starts; the diagonal search starts from that result. July 2's diagonal optimizer reached its iteration limit. Deposited points are reproducible, but neither optimizer convergence nor global optimality is claimed.

| Date / candidate | Singleton deviance | Simulations at least as discrepant / 20,000 | 99% upper bound on point-null tail |
|---|---:|---:|---:|
| July 2, separable extended | 916.615 | 38 | 0.002748 |
| July 2, with diagonal hopping | 875.608 | 1,520 | 0.080469 |
| July 7, separable extended | 3253.009 | 0 | 0.0002303 |
| July 7, with diagonal hopping | 3223.789 | 0 | 0.0002303 |

Each simulation draws complete multinomial records at the deposited point and original sample sizes. This avoids assigning speculative fitted degrees of freedom. It is a **point-null diagnostic**, not a calibrated composite-family goodness-of-fit test after fitting. The Monte Carlo upper bound covers simulation uncertainty conditional on the numerical probabilities. Domain/quadrature refinements change probabilities by less than 1e-8; they are convergence checks, not interval-arithmetic certification.

The July 2 diagonal point passes the retained calibration intervals and is not rejected by this aggregate diagnostic at 5%. The July 7 point still fails three calibration intervals and has poor aggregate agreement. Thus a more faithful model improves July 2 substantially, but does not explain all July records.

At the July 2 diagonal point, the same-event bunching predictions are:

| Evolution time | Fully bosonic | Fully labelled | Maximum C2 over partitions at each fluctuation |
|---|---:|---:|---:|
| 2.45 ms | 0.0306991 | 0.00550284 | 0.0106312 |
| 4.65 ms | 0.0134125 | 0.00262008 | 0.00486230 |

The observed rates are 93/2999 and 14/1202. The C2 maximum permits the partition to change with the fluctuation. It is a calculation **at one Hamiltonian point**, not an exclusion after propagating transport uncertainty. Full many-body pattern adequacy remains untested for this point. Neither another failed fit nor a small fitted-point C2 probability excludes the physical family. Establishing a calibrated parameter region and optimizing the many-body predictions over it remains necessary.

## A candidate label implementation, with explicit validation gates

A candidate is to put the three atoms in axial motional states n=0,1,2, all ending in the same electronic ground state. In an ideal resolved-sideband clock sequence, leave one atom in g0; use g0 → e1 → g1 for the second; and repeat a blue-sideband/carrier cycle to reach g2 for the third. Here g and e denote the clock ground and excited states. This requires coherent addressing, suitable trapping during the pulses, measured crosstalk and verified return to the sampling lattice. It is a proposed sequence, **not a demonstrated capability of this archive**.

[Shaw et al., *Erasure-cooling, control, and hyper-entanglement of motion in optical tweezers*](https://arxiv.org/html/2311.15580v2) demonstrates coherent control of the two lowest motional levels and electronic/motional transduction in strontium tweezers. That provides a relevant physical building block, but does not validate three matched labels in this lattice. The source experiment treats axial excitation as hidden distinguishability; this does not establish transport equivalence across deliberately prepared n=0,1,2 states.

Even in a perfectly separable harmonic approximation, different axial states change contact interactions. The exact diagonal contact integrals relative to I00 are

\[
\frac{I_{nm}}{I_{00}}=
\begin{pmatrix}1&1/2&3/8\\1/2&3/4&7/16\\3/8&7/16&41/64\end{pmatrix}.
\]

`fresh_protocol.py` checks these rational integrals. They are not bounds on off-diagonal collision-induced transitions, motional leakage, or real-apparatus errors. Ending in the same electronic state does not by itself preserve interactions or detection.

Before using the control inequality, the experiment must establish label overlap, label-dependent transport, interactions, readout and selection errors at the required scale. Avoid cooling that erases the prepared labels. Apply duration-matched sham operations to the original arm, randomize arm assignment after common initial selection, retain preparation failures as outcomes where appropriate, and interleave arms and label permutations. Any postselection must be shown to preserve the required common fluctuation law. Interleaving alone does not prove IID sampling.

Spatial labels that change the visible inputs are not automatically a control for this full two-dimensional event. Singletons in different shots do not reproduce shared shot fluctuations. Other electronic labels require their own transport and response validation. No deposited measurement currently closes these gates.

## Robust inequality for imperfect matching

For product pure labels, suppose every pair-overlap amplitude is at most a, uniformly conditional on the relevant latent variables. The six-assignment Gram matrix has diagonal 1 and off-diagonal absolute row sum at most r=3a²+2a³. Its eigenvalues lie between 1−r and 1+r. Write λ=1−r>0. Therefore the control event obeys pC,ideal ≥ λ DC, where DC is the fully distinguishable probability under the control transport.

Let δ bound |DO−DC| for the fully distinguishable **three-particle event**, and let eO,eC bound the absolute event-probability errors in the two arms from separately justified preparation, interaction and response effects. Then every C2 model satisfies

\[
p_{O,\mathrm{obs}}\le A p_{C,\mathrm{obs}}+b,
\quad A=\frac{2}{\lambda},\qquad
b=2\delta+e_O+\frac{2e_C}{\lambda}.
\]

This follows from pO,ideal ≤ 2DO and holds under shared transport fluctuations when the bounds apply conditionally or to a justified common coupling. A sufficient δ can come from the expected sum of full-column TV deviations under such a coupling; TV between **averaged** singleton distributions alone does not supply it. Similarly, average pair interference visibility is not automatically a uniform label-overlap bound. Mixed or correlated label preparations need an appropriate Gram/response bound rather than an unverified substitution into a.

These are event-level allowances, not percentages of timing error. Their calibration failure probabilities must be added to the statistical error budget. The ideal factor two is recovered only at a=δ=eO=eC=0.

## Fresh acquisition and power

Without validated historical transfer, collect fresh original and labelled-control shots under a common acquisition law. Predeclare NO, NC and reject C2 only when

\[
L_{\rm CP}(K_O,N_O)>A U_{\rm CP}(K_C,N_C)+b.
\]

Each one-sided Clopper–Pearson tail is 1/320. Under the declared IID model, the two acquisition bounds have joint failure probability at most 2/320; add all nuisance-calibration failure allocations. `fresh_protocol.py` verifies an exact **hypothetical** example: 130/5000 original and 23/5000 control admit conservative rational bounds L=.0201 and U=.0094, so L>2U. These are not acquired observations.

For the displayed physical bosonic model (pO=.02599744, pC=.004667035), direct binomial integration gives:

| Fresh original + control | Total science shots | Ideal power |
|---|---:|---:|
| 1500 + 1500 | 3,000 | 11.40% |
| 3000 + 3000 | 6,000 | 53.36% |
| 5000 + 5000 | 10,000 | 91.26% |
| 7500 + 7500 | 15,000 | 99.55% |
| 10000 + 10000 | 20,000 | 99.99% |

The earlier 3,000-control-shot / ≤19-events rule has 92.39% power only conditional on the historical original bound applying. It is not the cost or power of a fresh two-arm experiment. `fresh_power.py` integrates both binomial counts; reported powers use numerical SciPy evaluation, not interval-certified arithmetic.

Six rational amplitude-scaled fully bosonic models all pass the expanded selected-data region, including complete parity patterns and twelve row observables. Their 10,000-shot powers range from **86.58% to 93.69%**. Their physical convex mixtures remain compatible with those intervals. Monotonicity using the enclosing rate rectangle gives conservative numerical powers of **69.54% at 10,000 shots** and **94.16% at 15,000 shots** for this declared finite convex hull. This is not optimization over all compatible bosonic models and does not establish aggregate statistical adequacy of these examples.

A uniform high-power statement over *all* compatible quantum models is impossible if that class includes the surviving C2 models: type-I error control limits rejection at those null members. For the narrower fully indistinguishable class, a global minimum has not been established.

| Mismatch bounds | 10,000-shot power with displayed rates fixed | Conservative adverse-rate power |
|---|---:|---:|
| Ideal | 91.26% | 91.26% |
| a=.1 | 87.94% | 85.18% |
| a=.2 | 70.06% | 49.75% |
| δ=.001 | 76.91% | 47.75% |
| δ=.003 | 30.00% | 0.49% |
| a=.1, δ=.001 | 71.42% | 35.86% |
| eO=eC=.001 | 66.51% | 26.48% |

Unlisted allowances are zero. The first column changes the rejection threshold only. The adverse calculation also uses pO,min=pO−eO and pC,max=(1+r)(pC+δ)+eC, truncated to [0,1], with the original ideal transport as reference. It is conservative under that explicit alternative envelope, not a universal power guarantee.

Total acquisition cost is

\[
\left\lceil N_O/y_O\right\rceil+\left\lceil N_C/y_C\right\rceil+N_{\rm validation}.
\]

Usable yields and validation costs are not established. As a scale illustration only, a validated binary error assay with zero failures needs **6,905 trials** to upper-bound its error probability by .001 at failure probability .001. Three such assays plus 10,000 science shots cost **30,715 trials before yield losses or further calibration**. This calculation does not certify transport TV, unseen failures, or coherent label-state errors. A realistic protocol must price those measurements separately.

## Scope and reproducibility

The exact mixture defeats the particular proposed pair statistic. The separate structural example proves that complete pair information need not determine triple information; neither result rules out every pair-based strategy for this archive. Archive mixture compatibility remains compatibility with the declared interval region, not aggregate adequacy. The earlier source-removal witnesses establish prospective complementarity; no labelled-triple dataset has been observed.

Run the following after installing the pinned study requirements and reconstructing the source cache:

```sh
python studies/collective-interference-identifiability/fresh_protocol.py --check
OPENBLAS_NUM_THREADS=1 python studies/collective-interference-identifiability/extended_dynamics.py --check --cache /tmp/collective-sources
OPENBLAS_NUM_THREADS=1 python studies/collective-interference-identifiability/diagonal_dynamics.py --check
python studies/collective-interference-identifiability/fresh_power.py --check
```

`--fit` is optional numerical exploration and is not part of certificate verification. Results and proposal parameters are deposited separately. The next unresolved steps are a physically calibrated transport region and a measured label/response validation budget. The archive alone presently supplies neither a family-wide apparatus exclusion nor a fully validated control implementation.
