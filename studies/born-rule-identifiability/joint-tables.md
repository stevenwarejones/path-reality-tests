# What the measured combination can identify

The main extension combines the **December 2022 Nairobi preparation-angle scan**
with the **October 2023 Nairobi Viviani acquisition**. Their geometries differ:
the first has four fixed anchor preparations plus five scanned preparations;
the second independently varies five nonplanar preparations and four measurements.
They are distinct measured acquisitions, not independent processing of the same events.
The [generated report](results/joint-table-report.md) contains actual-data results,
full-table profiles, removal checks, ordinary physical constructions and validation.

“The rotation dataset cannot distinguish a probability deformation from an angle
warp because its whole circuit family samples one effective projection. A
preparation–measurement table supplies relations between independently varied
preparations and measurements; under a stable two-dimensional, context-independent
instrument model, extra rank can distinguish those explanations.” That first pairing
is informative about the **scope of the gauge**, but the unrestricted Lima angle warp
still contributes no theta constraint. It is retained as a failed identification
pairing, not counted as joint gain.

“The Viviani dataset cannot distinguish the deformation from preparation-dependent
leakage because both can add a fifth independent table component. The earlier scan
supplies a different geometry of deformation responses across five additional
preparations; under the specified pure-ensemble/readout/control family, it constrains
the same theta differently. With independently allowed leakage, both complete
tables instead admit ordinary qubit-plus-one-level models.” This is the implemented
replacement combination. It tests a universal law, **not** transfer of calibration
between two epochs bearing the same device name.

## Source audit and experimental bridge

| Source and immutable input | Actually available | Used / limitations |
|---|---|---|
| [7470893](https://zenodo.org/records/7470893), `wyniki-nairobi.zip` | Count rows with `0 0`, `1 0`, preparation `n`, measurement `j`, 115 complete jobs, 8 repeats/cell/job, 100000 shots/circuit | 9×4 full tables. IDs 0–115 except 10; include job 115. Metadata says 20 jobs and job list has only 20 entries. No guessed missing jobs. |
| Same record, `Job.py`, `dim-wit-par.py`, `dim-par-anp.py` | Source circuits, shuffled index generation, analyzer | Read as text; never execute acquisition scripts. Analyzer's fixed job loop is not copied. The circuit makes `n=0` the first anchor, `n=1..3` anchors 2..4, `n=4..8` the scan. |
| [21775462 v7](https://zenodo.org/records/21775462), `results_nairobi.zip` | Counts/settings in nine 5×4 job tables, 15 repeats/cell/job, 100000 shots/circuit | Metadata says one job; job list has one entry. Counts and row indices define all nine jobs. |
| Same record, `Job.py`, `viv.py`, `dim-an-viv.py`, `viviani-nairobi.txt` | Circuit order, angle lists, published processed determinant/table | Reconstruct gates independently. Count-derived table/determinant agrees with supplied Nairobi summary. |
| Same record, `ibm_nairobi_calibrations.csv` | T1/T2, frequencies, anharmonicity, assignment/gate errors | No preparation-resolved leakage, error bars, per-job timestamps or calibration-to-job mapping. Does not certify our axis box or leakage limit. No transfer to 2022. |
| Same record, `results-aria-viv.zip`, IonQ scripts and summary | 20 single-qubit hardware jobs, 5 repeats/cell/job, 1000 shots/circuit | Separate measured drift control, not a calibration of either Nairobi acquisition. Large job variation already discussed in the paper. |
| Same v7 deposit: Brisbane, Torino, Pittsburgh, Garnet files | Downloaded listings/bytes and scripts; counts for simultaneous multi-qubit measurements; simulations separately named | Not pooled as independent single-qubit experiments here. v7 explicitly repairs Pittsburgh indices and adds missing Garnet job; earlier versions cannot be silently substituted. |
| [7470893](https://zenodo.org/records/7470893), Perth and `*-s.zip` | Additional hardware and simulator files inspected | Simulator records not independent measured replication. Perth has a different parameter start and row schema; not silently fed through Nairobi decoder. |
| [GST experimental repository](https://github.com/pyGSTio/supplemental-info-arXiv-1605.07674) | Public experimental sequences; richer alternative if table route fails | Not analyzed as measured evidence here: the table route supplies the explicit joint construction. No claim of GST-based result. |

All selected originals, including analysis/acquisition scripts and calibration CSV,
are byte/MD5/SHA-256 pinned in [table-manifest.json](table-manifest.json). They remain
external. [table-audit.json](results/table-audit.json) contains newly reduced job-cell
counts, order hashes, source metadata discrepancies, circuit-order lag diagnostics
and incomplete job-list IDs. It contains no shot event stream. Order hashes differ
across jobs; each setting has the expected repetitions. Random ordering is recorded
at the **circuit** level; actual shot chronology is not supplied. Archive job numbers
order the holdout split; complete physical time ordering is not independently verified.

The 2022 table uses the four-gate sequence S_alpha, S_beta, S_phi, S_tau;
S_t = Z_t^dagger SX Z_t. In the code, matrices multiply in execution order, and
the probability is the squared |0> amplitude. The decimal eta=1.23095942 is retained
as deposited, not replaced by a more exact optimized value. The scan uses
alpha=p, beta=p+pi/2 with p=0,2pi/5,...,8pi/5. The Viviani sequence has
alpha=tau=0, beta=(pi/4,-pi/4,3pi/4,-3pi/4,0), phi=-beta[:4].
Some Bloch-vector sign conventions in the paper describe complementary outcome
tables; four-column complementation preserves the appended determinant. Our code
fixes the outcome by matrix multiplication and the deposited count header.

## Operational alternative and ordinary competitors

For pure preparation direction r and binary projective axis m, define

\[
p_\theta(0|r,m)=\frac{1+h_\theta(r\cdot m)}2,\quad
h_\theta(x)=x+\frac\theta2x(1-x^2),\quad -1\leq\theta\leq1.
\]

The other outcome is complementary. Unitary rotations act on pure directions;
preparation mixtures average these **pure** probabilities. This is a restricted
single-system operational theory with a richer operational mixed-state space;
it is not an ordinary density matrix with a nonlinear trace rule. No composite
or arbitrary-POVM extension is assumed; Bell observations are not attached to it.

In each acquisition d the fitted preparation is an antipodal ensemble with weights
(1+c_di)/2 and (1-c_di)/2 on ±r_di. After a projective measurement, independent
classical assignment parameters l_dj and u_dj give

\[
P^{(d)}_{ij}=l_{dj}+(u_{dj}-l_{dj})
\frac{1+c_{di}h_\theta(r_{di}\cdot m_{dj})}{2}.
\]

Only theta is shared. Each acquisition has its own directions, contrasts and
assignment parameters. Readout parameters may vary by measurement setting;
their stochastic validity is enforced, not a common fidelity prior. Directions
are displaced in a deterministic two-coordinate tangent frame and normalized.
Each coordinate lies in ±b (b=.15 in main fits, .05 and .5 sensitivity). The maximum
axis angle is at most arctan(sqrt(2)b); the actual maximum is saved. These are
sensitivity choices, **not measured calibration bounds**. For fixed theta all fits
are nonconvex feasible upper bounds from two starts. The grid objective weights
acquisitions equally in probability space. It is not a likelihood ratio or a
global confidence region; a worse local fit cannot exclude a parameter globally.

At theta=0, this model supplies valid ordinary qubit states and effects,

\[
\rho_i=(I+c_i r_i\cdot\sigma)/2,\quad
E_j=(l_j+u_j)I/2+(u_j-l_j)m_j\cdot\sigma/2.
\]

Our decisive ordinary competitor adds **one** orthogonal leakage level:

\[
\tilde\rho_i=(1-\lambda_i)\rho_i\oplus\lambda_i,
\quad \tilde E_j=E_j\oplus z_j,
\quad \tilde P_{ij}=(1-\lambda_i)P^{\rm qubit}_{ij}+\lambda_i z_j.
\]

All c,l,u,z,lambda are bounded in [0,1]. Thus states are normalized PSD and effects
lie between zero and identity by construction. This is ordinary Born probability;
the difference is preparation dynamics/state space and measurement response.
There is one population per preparation and one response per measurement, **not
one independent correction per outcome cell**. Parameters and predictions are saved
so the certificate can be checked without rerunning an optimizer. These states can
be prepared by replacement channels; the construction does not claim a single
gate-level CPTP map for every repeated S gate. That stronger sequence constraint
requires extra data and is outside this table-level compatibility result.

The implementation also constructs one ordinary qutrit instrument per recorded job,
with a tested leakage cap .005 and independent apparatus parameters across jobs.
It reports every residual, including failures to attain the numerical tolerance.
This explicitly permits drift but does not fit a dynamical drift law; using a
separate instrument for every job is an assumption relaxation, not a predictive
explanation of the hardware. The frozen-job check retains that distinction.

The same algebra describes a rare classical flag retained from preparation until
measurement. Constant leakage is insufficient to produce the extra dimension;
preparation dependence is essential. This is not an assertion that leakage is the
actual hardware fault. Gate-error summaries do not bound its population dependence
or leakage readout across these preparations. Both preparations and effects remain
independent choices in the constructed qutrit model.

For comparison a fully context-dependent binary-reset channel can turn any fitted
qubit table Q into observed P: use reset-to-0 probability (P-Q)/(1-Q) when P≥Q,
and reset-to-1 probability (Q-P)/Q otherwise. Its maximum is reported as a secondary
upper bound. It has more freedom than the qutrit construction and supplies no
unique physical explanation.

## Global algebra and scope of the angle gauge

The Legendre expansion is exact:

\[
h_\theta(x)=(1+\theta/5)P_1(x)-(\theta/5)P_3(x).
\]

The spherical-harmonic addition theorem expresses this kernel in 1+3+7=11
constant/degree-1/degree-3 features. For nonzero theta in the declared domain all
three coefficients are nonzero. Harmonics of different degrees are linearly
independent functions. There exist 11 evaluation directions giving an invertible
feature matrix: otherwise the span of all evaluation vectors would have dimension
<11 and a nonzero linear combination would vanish everywhere, a contradiction.
Choose such directions on each side. The kernel matrix is the product of two
invertible feature matrices and a nonsingular diagonal coefficient matrix, so rank
11 is attained. Analyticity makes full rank generic for sufficiently rich designs.
This proves existence, **not** rank 11 for every finite or noisy design. Ordinary
theta=0 has rank at most 4 and generically attains it. The supplied deterministic
16×16 example numerically illustrates, rather than proves, attainability.

For the actual 5×4 Viviani design, append a column of ones to P. Independent
symbolic reduction gives

\[
\det[P_\theta,\mathbf1]=-\frac{3\theta(3\theta+16)}{1024}.
\]

It is nonzero for **both** signs of theta≠0 in [-1,1]. Therefore no stable ordinary
qubit preparation/measurement reparameterization can match the ideal full table.
This globally breaks the earlier one-projection angle gauge for this design.
The determinant is -.00046962890625 at +.01 and +.00046787109375 at -.01.
An overall preparation contrast c and four fixed assignment contrasts v_j multiply
the determinant by c^4 times the product of v_j: moderate contrast reduces
sensitivity continuously, unlike the earlier unit-disk witness's contrast threshold.
Unequal ensemble contrasts and freely varying controls require the full model.

This does **not** identify theta for unrestricted ensembles. For example, the
nonnegative spherical ensemble density (1+3c r·u)/(4pi), |c|≤1/3, has no degree-3
component. Averaging h over it gives (1+theta/5)c r·m, an ordinary contrast change.
The equivalence is global over all measurement axes for that ensemble class;
it does not cover nearly pure high-contrast preparations automatically. An inference
that simply substitutes a mixed-state Born probability into f would miss this issue.

## Certified relaxation and finite-sample boundaries

Let H subtract the mean over preparation rows, and let P have n×4 entries. Any
ordinary qubit table Q has HQ of rank ≤3, irrespective of pure/mixed states,
unknown directions, common channels or asymmetric measurement effects. The
Eckart–Young bound and ||H||_2=1 imply

\[
\inf_{Q\in\mathcal Q_2}\|P-Q\|_\infty
\geq \frac{[\sum_{k>3}\sigma_k(HP)^2]^{1/2}}{\sqrt{4n}}.
\]

This is a certified **relaxation lower bound**: arbitrary rank-three centered
matrices are not necessarily physical qubit tables. Physical fitted qubit tables
give upper bounds. Together they bracket a clear metric without claiming to solve
the nonconvex global distance problem. The same lower bound applies to a maximum
leakage fraction L in the specified ordinary model, because |tilde P-Q|≤L.
The smallest successfully tried leakage cap is an upper bound, not a minimum.

For counts N_ij of mutually independent bounded Bernoulli trials (not necessarily
identical means), the empirical error e_ij is sub-Gaussian with proxy
v_ij=1/(4N_ij). Independence across cells is also required. Gaussian integration
of the scalar sub-Gaussian mgf gives, for 0<s<1/(2 max v),

\[
E\exp(s\|e\|_F^2)\leq\prod_{ij}(1-2sv_{ij})^{-1/2}.
\]

The usual exponential minimization yields with probability ≥1-alpha

\[
\|e\|_F\leq R_\alpha=
\sqrt{\sum v+2\sqrt{\sum v^2\log(1/\alpha)}+2\max v\log(1/\alpha)}.
\]

Distance to the rank-three set is 1-Lipschitz in Frobenius norm. Subtract R_alpha
from the empirical tail norm, truncate at zero and divide by sqrt(4n). Each
acquisition gets alpha=.025; a union bound gives simultaneous 95% coverage without
assuming independence **between acquisitions**. The bound concerns the average
probability table if trials drift; it does not show each instantaneous device
violates the qubit model. It is not a Born-rule confidence interval.

If shots within a job can be arbitrarily dependent, treat each complete job table
as a bounded observation. Equal cell exposures make the pooled table the job mean.
Assuming independent jobs but allowing dependence between cells, Hoeffding plus
union bound gives R=sqrt(4n log(8n/alpha)/(2J)). This is deliberately conservative
and vacuous here. With dependence between jobs even that coverage is unavailable.
Empirical job standard errors and circuit-order lag diagnostics are descriptive;
large counts do not establish independence.

Pooling changing qubit preparations and effects can increase rank:
E[r_t m_t^T] is not E[r_t]E[m_t]^T. We preserve within-job singular values and
pooled versus mean determinants, and use the independent Aria acquisition to
stress this issue. Its drift was already documented by the original authors.
No published determinant significance is relabeled as our discovery.

## Model choice, holdout and simulation interpretation

The [protocol](table-protocol.json) fixes the theta grid, tangent boxes, leakage
caps, optimizer starts, split and simulation seed before held-out fitting. The
study is retrospective and the public papers were inspected; it is not preregistered.
The first floor(2J/3) archive-numbered jobs train the models. The remaining block
is assessed with frozen predictions and no re-selection. All twenty/36 cell
predictions contribute. Binomial deviances are descriptive, with no Wilks
calibration at singular rank models. Sensitivity analyses are not independent
confirmations or multiple-comparison-adjusted discoveries.

Simulations generate binomial **counts** at the observed cell/job totals, equivalent
to aggregating the recorded circuit repeats when their probabilities coincide.
They include fitted contrast/readout, theta=.01, the explicit ordinary qutrit,
two correlated drifting ordinary qubit regimes, and a preparation-memory response.
The latter two deliberately violate pooled stationarity/independence. Their
rank-test rejections measure vulnerability to misattributing apparatus to a new
rule, not false rejections of the stated stable-qubit null. A known-nuisance linear
theta projection checks injection recovery; it is explicitly not nuisance-profiled
identification. This limitation is also visible in the real held-out records.

## Prior art and precise delta

| Claim | Closest prior work | What this PR adds / does not add |
|---|---|---|
| Qubit table rank bound, null determinant | [Precise certification](https://doi.org/10.1140/epjqt/s40507-024-00230-4), [2301.03296](https://arxiv.org/abs/2301.03296) | Reuses the known bound; no new dimension theorem. |
| Viviani circuits, cross-platform anomalies, pooled drift | [Cross-platform paper](https://www.nature.com/articles/s41598-025-27248-7), [2404.06792](https://arxiv.org/abs/2404.06792) | Reconstructs deposited settings and preserves jobs; their anomalies and drift are already published. |
| Leakage as a dimension-witness explanation | [Strikis–Datta–Knee](https://arxiv.org/abs/1811.05220) | Supplies explicit small, preparation-dependent qubit-plus-one-level states/effects for these **two full measured tables**, not a claim that leakage was discovered here. |
| Self-consistent calibration and gauges | [Nielsen et al., Quantum 5,557](https://quantum-journal.org/papers/q-2021-10-05-557/) | Calibration-free GST still uses Born probabilities (review §2) and its model class has stationarity/Markov assumptions (§6). We do not use a Born-based fitted gate set as independent proof of Born's rule. |
| Higher harmonic operational dimension | Generalized probability/operational sources in [literature](literature.md) | Applies the known harmonic representation to this explicitly normalized family and the actual deposited designs; no harmonic-priority claim. |
| Specific measured joint result | The two papers analyze dimension tests; neither cited analysis supplies this shared-theta profile plus physical ordinary joint certificate | Connects two different measured geometries, shows removal effects in a declared nuisance family, and brackets the cost of an ordinary explanation with dependence sensitivity. Literature priority is not established. |

The greatest-value new measurement is **preparation-resolved leakage population and
leakage-level response, interleaved with randomized full tables and time-tagged
repetitions**. This addresses the surviving explicit construction and drift. Tighter
assignment fidelity alone or another uncalibrated determinant does not remove it.
The result is intentionally conditional and finite-design: no universal impossibility
or established new-physics conclusion follows. The original optical feasibility
and its limitations remain as a secondary result. An elementary-only Lean companion
would not add scientific value here.
