# Joint statistics: exclusion of four earlier candidate points

**Follow-up:** [Classical-first dissipative models](classical-fit.md) now pass
both components of this unchanged region, with separately verified quantum
counterparts. The four exclusions below remain valid. The implementation now
checks the 97.5% simultaneous cell component as well as the aggregate component
for every candidate.

**The four tested physical pairs are outside a valid aggregate 95% confidence
region. The classical-memory class has not been excluded.** The previous
cell-interval certificate remains correct, but is insufficient evidence of
aggregate statistical adequacy. We now have a finite-sample calculation of
that limitation, not just a large descriptive deviance.

The calculation uses all 2,916 four-outcome IBM rows, including both middle
flags, separately at all nine delays. It conditions on the selected historical
mapping and assumes independent multinomial rows, predetermined exposures,
and stable probabilities within each row. **Cross-row independence is an
additional assumption:** the earlier cell-wise union bound did not require it.
The archive cannot establish independent trial acquisition, drift control, or
the historical mapping. Literal rows remain a conditional sampling model,
not a verified reconstruction of independently acquired trials. No conclusion
here repairs those acquisition uncertainties.

## A finite-sample confidence region without fitted degrees of freedom

For row r, let N_r be its total, n_ri its four counts and p_ri>0 a complete
candidate probability vector. Define the Pearson statistic

\[
 X(n,p)=\sum_{r,i}\frac{(n_{ri}-N_rp_{ri})^2}{N_rp_{ri}}.
\]

Under that fixed p, exact multinomial moments give

\[
 \mu=E_pX=3R=8748,\qquad
 V(p)=\operatorname{Var}_pX
 =\sum_r\left[6+\frac{\sum_i1/p_{ri}-22}{N_r}\right].
\]

Boundary probability points are included: if p_ri=0 but n_ri>0, the point
assigns the observed event probability zero and is excluded. Otherwise omit
zero-probability categories and replace each row's mean by K_r−1 and variance
by 2(K_r−1)+[sum_(p_i>0)1/p_i−K_r^2−2(K_r−1)]/N_r, where K_r is its support
size. A deterministic row has zero statistic and variance. The verifier tests
these boundary cases; positivity of the fitted points does not silently remove
boundary models from the full confidence-region definition.

These are finite-N identities, not asymptotic chi-square moments. For one
row, expand X as sum_i n_i^2/(Np_i)−N and use multinomial factorial moments
E[(n_i)_a(n_j)_b]=(N)_(a+b)p_i^a p_j^b. This gives mean 3 and the displayed
variance. Independence adds variances across rows. The tests independently
enumerate small multinomial experiments and recover both identities exactly.

For X>mu, Cantelli's inequality gives

\[
 \Pr_p\{X(N,p)\ge X(n,p)\}
 \le u(n,p):=\frac{V(p)}{V(p)+(X(n,p)-\mu)^2}.
\]

Set u=1 for X<=mu. The inverted aggregate region C_A(n)={p:u(n,p)>0.025}
has coverage at least 97.5% for every fixed true p in the declared sampling
model. To see this, u<=alpha is equivalent to
X−mu>=sqrt(V(1−alpha)/alpha), whose probability is at most alpha by Cantelli.
The bound is conservative, but already decisive for the tested points.

This is a **confidence-region membership audit**, not a fitted-model deviance
test. The region is defined over complete row probabilities. Fitting, selecting
or revisiting candidate points after seeing the counts does not change its
pointwise coverage for the true fixed probability vector. We do not subtract
estimated parameters, select a chi-square degree count, or turn a selected
candidate's u into a p-value for its whole model class. Excluding a composite
class would require proving that **every** probability vector in that class
lies outside C(n); no such optimization or proof is supplied here.

The test's retrospective choice does not supply a generic correction for an
arbitrary search over statistical tests. The defensible claim is the coverage
of this explicitly specified region and its calculated membership, not a
preregistered confirmatory finding or an omnibus selective-inference guarantee.
The two mapping regions are conditional alternatives; their union covers a
correct interpretation under its sampling assumptions.

To account for choosing between cell-wise and aggregate summaries, define the
new joint region as C_A intersected with a simultaneous **97.5% cell box**,
using tail allocation 0.025/(2M), M=11664. Bonferroni coverage is at least 95%,
without independence between the two summaries. The earlier saved 95% cell
certificate is unchanged; its box is not silently reused with the new budget.
All four points fail C_A, so checking the other component cannot rescue them.
This new joint region is not claimed to be a subset of the earlier 95% box.
All rows, delays and flags are retained; no selected residual or favorable
mapping supplies the decision. Exploratory simulation diagnostics do not
enter this certified decision. A claim selected from further tests or record
subsets would need its own simultaneous or selective-inference control.

## Results

| Conditional interpretation / candidate | Pearson X | Exact null mean | Cantelli tail upper bound | Aggregate region |
|---|---:|---:|---:|---|
| Regrouped, earlier identity-delay pair | 16478.592367 | 8748 | <0.000295 | Outside |
| Literal rows, earlier identity-delay pair | 15361.778458 | 8748 | <0.000402 | Outside |
| Regrouped, separate delay unitaries and fitted shared SPAM | 10577.602748 | 8748 | <0.005232 | Outside |
| Literal rows, separate delay unitaries and fitted shared SPAM | 9594.818295 | 8748 | <0.023961 | Outside |

The first two pairs retain their earlier exact physicality and cell-interval
certificates. Those claims and the new aggregate exclusions are compatible:
they concern different confidence constructions and different independence
assumptions. Since each pair has exactly equal probabilities, its classical
and quantum explanations receive exactly the same statistical decision.
The opposite first-bit regrouping is a label permutation of the corresponding
regrouped model and has the same statistics.

The [saved result](results/ibm-joint-statistics.json) contains outward bounds,
not only the rounded numbers above. For each rational candidate probability,
the verifier accumulates lower/upper Pearson values and exact-moment variance
bounds using directed Decimal arithmetic. It accepts an exclusion only if the
upper bound on u is <=0.025. No simulation quantile or optimizer status enters
that decision. Approximate multinomial deviances remain display diagnostics.

## New delay-dependent physical candidates

The expanded search assigns independent pre- and post-instrument qubit
unitaries to every delay pair. It also fits six shared preparation states and
three shared two-outcome readout effects. The 18 two-flag instruments remain
shared across delays, with first effects shared across re-preparation labels.
These controls are nuisance choices; no new measured calibration is asserted.

A bounded eight-round alternating search minimizes full multinomial negative
log likelihood: convex shared-instrument updates and local unitary/SPAM
updates. The initial seed was fitted to these same counts. Some SDP solves
report `optimal_inaccurate`; they are discovery only. The saved points are
rationalized and independently verified. The search is not a global optimum
and does not exhaust delay-dependent classical memory. Its failure to find
an accepted point did not exclude that larger class. The later
[classical-first search](classical-fit.md) supplies accepted physical points.

The new instruments have an exact quantum counterpart with w=1/10000 using
the previously proved B_b=(1−w)K_b+w q_b Id identity. Rational LDL tests verify
both instruments' CP and TP. The pre/post unitaries are represented by rational
Cayley coordinates; exact orthogonality and determinant checks verify the
rotations. Shared SPAM states and effects obey exact Bloch-ball inequalities.
The complete rational probability tables are independently compared with
complex Choi contractions and a Kraus/SWAP circuit check.

For each delay, wrapping the original process with these pre/post unitaries
is a local unitary transformation across the temporal AB|OC cut. It preserves
the partial-transpose witness: transform its test vector by the corresponding
partially conjugated local unitary. The expectation remains −w/8=−1/80000.
Thus **both new candidate pairs are physical and differ in process resource**;
they are not claimed to satisfy the new statistical region. A better classical
fit was not automatically labeled quantum: its separate counterpart was
constructed and checked.

The instruments differ between each pair's two explanations. Their w values
and CP endpoints concern those selected instruments, never the maximum memory
allowed by all data-compatible models. None of these rejected candidates
supports a practical ambiguity claim under the stronger region. We therefore
do not inflate the small-w probe separation into a new experimental remedy.
The earlier isolated-instrument proposal remains a discriminator of its
specified pair only.

## Decision and reproduction

**Continue the statistical/physical feasibility question; reject the use of
these four points as adequate full-record explanations.** Preserve their valid narrower
certificates and their exact model-family identities. No adequate paired fit,
classical-class exclusion, or empirical complementary-source gain has been
established. More general delay-dependent channels and classical records,
acquisition dependence and stronger independently justified calibration remain
possible directions; no local search failure resolves them.

```bash
python studies/calibrated-quantum-memory/joint_statistics.py --source-dir /tmp/nmn --check
python -m unittest discover -s tests -p 'test_calibrated_quantum_memory.py' -v
```

The committed [delay candidates](results/ibm-delay-candidates.json) contain
model coefficients, not raw source counts. Optional discovery uses the
separate solver requirements and writes a rational candidate JSON:

```bash
python -m pip install -r studies/calibrated-quantum-memory/requirements-fit.txt
python studies/calibrated-quantum-memory/fit_delay_candidates.py --source-dir /tmp/nmn --mapping regroup_first_1 --rounds 8 --output /tmp/delay-candidate.json
python studies/calibrated-quantum-memory/joint_statistics.py --source-dir /tmp/nmn --candidate-file /tmp/delay-candidate.json --output /tmp/delay-candidate-statistics.json
```

Solver output is not a certificate; the committed rational candidates are the
objects checked by the independent verifier. The empirical high-bar goal
remains open. The PR stays draft; no merge or author contact is authorized.
