# Classical-first fits pass the complete fixed region

**The fixed 95% joint region contains explicitly physical classical models
under both conditional IBM record interpretations. Separately constructed
quantum-memory counterparts also pass every component of that same region.**
The new search removes the earlier fixed quantum-counterpart restriction and
permits dissipative dynamics. This supersedes the failure of the four earlier
candidate points; it does not invalidate those exclusions.

This is compatibility with the declared, conservative confidence region.
It is not proof of the explanations' correctness, a quantum-memory detection,
or a guarantee of survival under every sharper statistical analysis. The
independent-row and fixed-exposure assumptions remain unverified, and neither
historical record interpretation has been established from original jobs.
No claim of complementary measured-source gain is made.

## Architecture and the classical families actually searched

For delay j, the physical classical model is

\[
 p_j(b,c|a,y,z)=\sum_e\operatorname{Tr}\bigl[
 E_{c|z}\,C_{e,j}\bigl(B_{b|y}(A_{e,j}(\rho_a))\bigr)\bigr].
\]

The A_e,j are CP maps with TP sum. Their index e is an ordinary classical
record retained until the CPTP return map C_e,j. It does not receive the
preparation label or future setting. Both observed flags b and c are retained.
The system is a qubit; no leakage or instrument/environment interaction is
needed for these admitted witnesses. Allowing those imperfections would not
remove them.

| Quantity | Sharing and status |
|---|---|
| Six preparation states and three binary final POVMs | Shared across all delays and settings; physical nuisance choices, not independently measured calibration. |
| Eighteen two-flag middle instruments B | Shared across delays, input preparations and final axes. The first POVM effect for fixed first axis and flag is exactly shared across re-preparation labels. |
| Pre/return dynamics | Separate physical CP/CPTP maps at each delay. No stationary generator, common unitary, or bounded-angle law is assumed. |
| Literal-row witness | One record value: memoryless general CPTP pre- and post-channels. Shared SPAM starts from the previously fitted physical choice. |
| Regrouped witness | Thirteen record values. Conditional return channels are a fixed physical bank at each delay: one fitted dissipative channel, six small unitary perturbations of it, and six Pauli-state reset channels. The pre-instrument is a freely fitted 13-outcome CP instrument; shared B and SPAM can be updated. |
| Larger null | Arbitrary classical memory, more record values, b-dependent controls and more general conditional returns remain allowed. These searches are inner families, not globally exhaustive null optimizations. |

The unrestricted instrument constraint is J(B_b|y)>=0 with TP sum and shared
first effects. **There is no constraint of the form J−wqJ_Id>=0 during the
classical search.** The earlier `fit_delay_candidates.py` retains its historical
fixed-w search for reproducibility; the current discovery entry point is
`fit_classical_models.py`.

Alternating convex blocks fit full multinomial likelihood. General CPTP
channels replace bounded rotations; optional shared states and effects range
over their physical convex sets. A later conditional-record search gives the
regrouped witness. This is still local discovery. Its success supplies a
physical point; its failures would not prove exclusion of a family.

Rounded matrices are repaired to exact TP/shared-effect identities and given
small, explicitly recorded depolarizing repairs where needed for exact strict
positivity. This is a post-fit physicalization step, not a quantum-memory
constraint. Every final probability is checked after repair. Optimizer status,
including `optimal_inaccurate`, is never a certificate.

## Complete statistical acceptance

The [statistical construction](joint-statistics.md) is unchanged: 0.025 error
probability for the aggregate Cantelli region and 0.025 for a simultaneous
Clopper–Pearson cell box. Each cell tail receives 1/(80×11664)=1/933120.
Independent fixed-exposure multinomial rows are required for the aggregate
component. Within-row flag correlations are included in the four-outcome
multinomial calculation.

| Mapping / model | Pearson X | Cantelli bound (approx.) | Complete fixed region |
|---|---:|---:|---|
| Regrouped, classical-first | 9411.050346 | 0.03849769 | Inside; all 11,664 cells pass |
| Regrouped, post-fit pair | 9449.896172 | 0.03449065 | Inside; all 11,664 cells pass |
| Literal rows, classical-first | 8759.181634 | 0.99294743 | Inside; all 11,664 cells pass |
| Literal rows, post-fit pair | 9069.988454 | 0.14506974 | Inside; all 11,664 cells pass |

The aggregate cutoff is 0.025; verified directed lower bounds exceed it.
The regrouped pair remains relatively close to this conservative boundary.

Every row in this table passes **all 11,664 cell constraints**, with zero
excluded or unresolved cells. Acceptance uses directed lower bounds on both
binomial tails, not inverse-CDF rounding. The checker now evaluates both
components even when a point is already excluded by one. A point can be called
inside the joint region only when both independently return `inside`.
Tests include an example that passes the aggregate component and fails the
cell component, as well as exact enumeration of binomial tails.

The two fitted control/process choices are alternatives for two conditional
mappings; no selection of a convenient mapping establishes a violation.
The opposite first-bit convention is the corresponding re-preparation-label
permutation. Fitting and post-fit model selection do not change the confidence
region. They do not turn its membership into a composite-model p-value or
provide a guarantee for an unrestricted retrospective search over tests.

## Quantum counterparts are a separate calculation

Only after an unrestricted classical point passed both components did we
explore a physical instrument direction

\[
 \widetilde B_b=(1-t)B_b+tD_b,\qquad D_b(X)=\operatorname{Tr}(X)I/4.
\]

The selected t values are 0.0003 (regrouped) and 0.001 (literal rows). These are
post-fit model choices, not measured control errors. The altered classical
points were independently rechecked against the entire unchanged region.
Then use K_b=(Btilde_b−w q_b Id)/(1−w), with q_b=Tr[Btilde_b(I/2)].
Exact rational LDL checks verify each K block and TP sum. The selected w is
99% of the exact rank-one CP endpoint for that chosen instrument. It is not a
globally optimized bound on memory permitted by the observations.

For each classical record e, a quantum branch swaps A_e,j(rho) into an
environmental qubit, applies K to I/2, swaps back and applies C_e,j. The other
branch transmits through K normally. Summing the branches implements Btilde
exactly, for every input operator and retained outcome, including reference
correlations. Thus both physical models have exactly the same full probability
table and statistical decisions. Their instruments differ.

Dissipative wrapping does **not** automatically preserve the quantum resource.
We therefore construct and test the full process directly. In wire order
A,B,O,C, put

\[
 W_{\rm cl}=\sum_e J(A_e)_{AB}\otimes J(C_e)_{OC},\qquad
 W_{\rm sw}=J\!\left(\sum_e C_e\circ A_e\right)_{AC}
                   \otimes(I_B/2)\otimes I_O.
\]

The quantum process is (1−w)W_cl+wW_sw, with trace 4. Classical W_cl is
explicitly temporally separable. A numerically discovered test vector is
rounded to integer real/imaginary coordinates; its expectation against the
normalized process's partial transpose is evaluated as an **exact rational
number**. Negative expectation, not a floating eigenvalue, certifies the
constructed quantum process. No equivalence of PPT and classical memory is
assumed. The certificate concerns at least the named delay, not necessarily
all nine processes.

| Mapping | Memory weight w | NPT delay | Exact negative expectation (decimal display) | Probe TV gap | Shots per candidate* |
|---|---:|---|---:|---:|---:|
| Regrouped | 0.000167241111 | 28.444,24.889 | -9.02756131e-06 | 0.000164465748 | 324,006,768 |
| Literal rows | 0.000558450137 | 28.444,21.333 | -6.88915267e-05 | 0.0005541207 | 28,542,781 |

*Conservative shot scale for two simultaneous Hoeffding radii at least
95% confidence to sum to less than the predicted binary-event gap.
This assumes ideal external probes and is not an experimental power analysis.

The [rational coefficients](results/ibm-classical-candidates.json) and
[complete audit](results/ibm-classical-audit.json) contain the exact endpoints,
vectors, negative expectations and each statistical component. Independent
tests compare all predicted probabilities with direct dissipative channel
composition, check the Kraus/SWAP identity, and reject invalid process records.

## Remedy, acquisition assumptions and remaining complementary-data goal

An isolated-instrument probe distinguishes each pair, assuming independently
anchored external preparations/readout and prevention of environmental return.
The predicted TV gaps and conservative shot-precision estimates are saved in
the audit. They are small and require substantial precision; a practical
calibration remedy has **not** been established. Neither instrument's distance
from the real device is known. No proposed probe is presented as measured data.

The acquisition check examined all 13 commits reachable from pinned NMN main
`154235f8bbf5e70eb71c325370a67b1894490452`, including five distinct notebook
contents. The notebook versions first encountered at `154235f`, `0336da1`,
`5545a5a`, `bd4bf0a`, and `f5c19b4` contain no recovered acquisition call,
backend identifier, job ID, retrieval call, or acquisition timestamp. The
current public refs inspected identify that main branch and no tags.
The `n=8000 # num shots` line is in synthetic bootstrap generation, beside
`g=1000`; it is not a saved hardware job configuration. A matching character
sequence in an encoded plot is not IBM job metadata. This bounded source check
therefore **does not verify** independent rows, predetermined exposures,
interleaving, drift control, or the historical regrouping.

The paper's assignment point estimates do not characterize the whole flagged
instrument and lack deposited calibration trials and epoch linkage. The q0
IQ/QPT archive belongs to another acquisition and cannot calibrate these IBM
operations. No independently measured, matched instrument calibration was
recovered. The [acquisition specification](missing-acquisition.md) remains
necessary; contacting authors or accessing private hardware jobs was not
undertaken. Missing provenance is an empirical limitation, not a mathematical
no-go theorem.

**Decision: retain the conditional full-region ambiguity result; continue the
acquisition/complementarity objective.** A positive quantum-memory exclusion
from this fixed region is unavailable because explicit classical points lie
inside it. The new question answered is whether the fixed region contains explicitly
physical classical explanations. The unrestricted dissipative/record models
supply such explanations. This is a
stronger conditional result than the earlier rejected examples, but it still
does not establish a consequential combined-measurement resource bound or a
practical distinguishing calibration. No novelty is claimed for instrument
ambiguity, classical-record models, the SWAP identity or the statistics.

## Reproduction

```bash
python studies/calibrated-quantum-memory/classical_audit.py --source-dir /tmp/nmn --check
python studies/calibrated-quantum-memory/joint_statistics.py --source-dir /tmp/nmn --check
python -m unittest discover -s tests -p 'test_calibrated_quantum_memory.py' -v
```

Optional discovery starts from the saved, data-fitted classical seed and
imposes ordinary CP/TP constraints. It does not require any quantum certificate:

```bash
python -m pip install -r studies/calibrated-quantum-memory/requirements-fit.txt
python studies/calibrated-quantum-memory/fit_classical_models.py --source-dir /tmp/nmn --mapping literal_rows --rounds 1 --output /tmp/classical-candidate.json
python studies/calibrated-quantum-memory/classical_audit.py --source-dir /tmp/nmn --candidate-file /tmp/classical-candidate.json --output /tmp/classical-candidate-audit.json
```

For the 13-record family, `--fixed-return --fit-spam` retains the chosen return
bank and fits shared SPAM. Without `--fixed-return`, each conditional return
is also optimized over CPTP channels. No solver trajectory is used as proof;
the committed rational points are the reproducible certificate objects.
Raw sources remain external. The PR remains draft, with no merge or author
contact performed.
