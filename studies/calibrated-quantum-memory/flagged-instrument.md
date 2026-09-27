# Classical null with an observable instrument register

This is an architecture-level extension of [the single-outcome argument](theory.md),
not a certificate evaluated on the deposited experiments and not a novelty
claim. It specifies the model that additional measurements would test.

## Wires, flags, and allowed controls

For fixed physical intervention setting y, let B_(r|y) be CP maps from the
accessible system S to S', with their sum trace preserving. r contains **every
classical output available to later dynamics or controls**, not just the bit
saved in a spreadsheet. Package the instrument as

\[
\widehat B_y(X)=\sum_r |r\rangle\langle r|_R\otimes B_{r|y}(X).
\]

The null allows arbitrary classical environmental memory l at the cut,
arbitrary quantum system propagation on either side, and outcome-dependent
feed-forward:

\[
p(r,c|a,y,z)=\sum_l\operatorname{Tr}\!\left[
 E_{c|z}\,C_{l,r,y}\!\left(B_{r|y}(A_l(\rho_a))\right)\right].
\]

Here {A_l} is any CP instrument with TP sum; each C_(l,r,y) is CPTP.
There is no bound on the number of classical states. A_l cannot depend on the
secret preparation label a or a later setting; C may use l,r,y but not a.
Final z selects an anchored POVM. Common causes independent of a can be
included in l. Allowing C to receive r is conservative for a device with no
physical feed-forward. Outcome-conditioned data selection is not permitted:
all branches enter the likelihood with their original denominators.

For NMN, r includes the first outcome, y the **fixed rotation** setting, and
actual re-preparation labels are functions of r,y. The reconstructed output
label must not be silently reinterpreted as a deterministically chosen input.
For an active reset, r must include any reset attempt/acceptance history that
remains accessible. A discarded record can be traced out only after physical
erasure/isolation is justified; merely omitting it from the archive is not
erasure. A hidden quantum output or leakage level requires a larger system
model. Direct instrument–environment coupling violates the factorization
used here and needs its own justified error bound or explicit model.

## Globally valid sufficient inequality for a calibrated breaking instrument

Let an ideal reference K be any **flagged entanglement-breaking channel**.
For example, a fixed measure/reprepare instrument has

\[
 K_r(X)=\operatorname{Tr}(F_r X)\sigma_r,\quad
 F_r\succeq0,\quad\sum_rF_r=I,
\]

with fixed normalized output states sigma_r. They need not be maximally mixed
or independent of r. More general EB blocks with internal classical indices
are allowed. Define epsilon = (1/2)||hat B - hat K||_diamond on the full S'R
output, including flags available to downstream controls.

Replacing B by K in the null gives

\[
\Lambda_K(X)=\sum_{l,r}
 \operatorname{Tr}\!\left[A_l^*(F_r)X\right]C_{l,r}(\sigma_r).
\]

The effects A_l^*(F_r) are positive and sum to I; the output states are
normalized. Thus Lambda_K is EB even with arbitrary l and feed-forward.
For qubit input/output, its six-state average recovery fidelity is at most
2/3, by the benchmark proved in the existing theory. Preprocessing, keeping l,
and postprocessing controlled by l,r are CPTP maps, so diamond-norm
contractivity gives

\[
 F_6(\Lambda_B)\le\min\{1,2/3+\epsilon\}.
\]

The same fixed intervention and subsequent recovery must be used for all six
input states. Input-dependent recovery would reveal a and invalidate the
benchmark. Closeness of the **averaged** channel to an averaged EB channel
does not suffice: the flagged Pauli example in the existing tests violates
that substitution. The claim is valid for a qubit and the declared cut, not
arbitrary leakage or spillage. Unknown sigma_r must be included in a jointly
valid region or replaced by a preregistered reference; estimating it and then
treating it as exact would be unjustified.

## A conservative count-level route once the required tables exist

With trusted six-state preparations and terminal Pauli readout, record the
joint probabilities p(r,c|s,k) of each intermediate flag and final bit, at
every calibration input s and output axis k. Do **not** normalize conditional
on r. Their sums/differences determine

\[
 B_r(\rho_s)=\tfrac12\left[p(r|s)I+
       \sum_{k=x,y,z}\{p(r,0|s,k)-p(r,1|s,k)\}\sigma_k\right].
\]

The p(r|s) are shared across k. For the feasible instrument region impose
J(B_r)>=0 and sum_r Tr_output J(B_r)=I, with one Choi convention throughout.
For each retained probability cell use a simultaneous Hoeffding interval
with radius sqrt(log(2M/alpha)/(2N)), or a tighter interval with proven
simultaneous coverage. Independent trials **within** each table are required;
a union bound needs no independence between tables. Include the full dynamics
cells and six success events in M, as specified in the acquisition schema.
SPAM anchors, drift, missing trials and classifier selection require additional
constraints/budgets rather than an inverse-confusion correction treated as exact.

One fully conservative norm bound is available without a local optimizer.
Writing sigma_0=I and Delta_r=B_r-K_r, the Pauli expansion gives

\[
\epsilon\le \min\!\left\{1,\frac14
  \sum_{r,\mu=0}^3\|\Delta_r(\sigma_\mu)\|_1\right\}.
\]

Indeed each map X -> Tr(sigma_mu X) Delta_r(sigma_mu)/2 has diamond norm
at most ||Delta_r(sigma_mu)||_1/2, and classical output blocks add in trace
norm. Form B_r(I) from a pair of opposite probe states and B_r(sigma_mu)
from their difference. Taking interval upper bounds on these trace norms
(e.g. bounding every scalar Pauli coefficient and using the triangle
inequality) yields a globally valid epsilon_U. A tighter channel-constrained
upper bound needs a certified **global** bound; a local optimum is insufficient.
No numerical flagged bound has been evaluated here.

On the simultaneous event, the dynamics lower bound F_L requires
`epsilon_required >= max(0, F_L - 2/3)` for every null model. If calibration
transfers with system-only half-diamond error at most d and anchored probability errors
at most eta, the exclusion condition is

\[
 F_L>2/3+\epsilon_U+d+\eta.
\]

Direct environment coupling is not covered by a system-only d: it needs a
justified bound on the full coupled operation or a different model.

Keep the same original confidence budgets when either source is deleted.
A positive inequality alone is not the requested scientific result: also
exhibit a physical quantum process and instrument satisfying **every** joint
constraint and complete, physical source-only models. The available records
have not supplied the premises for these empirical tasks. The prospective
single-outcome synthetic witnesses do not certify this measured architecture.

## Relation to closest work

[Giarmatzi–Costa](https://doi.org/10.22331/q-2021-04-26-440) already treats
quantum environmental-memory witnesses; [Taranto et al.](https://arxiv.org/html/2307.11905v2)
provides the operational classical-memory hierarchy (classical memory is not
equated with separability or PPT). [NMN tomography](https://arxiv.org/html/2308.00750v3)
already obtains complete measure/reprepare statistics through retrospective
labeling. [White et al.](https://doi.org/10.1103/PhysRevX.15.021047)
already includes noisy control, self-consistency, gauge and spillage: the
system-only factorization here is a stricter premise, not an improvement on
that characterization framework. [Rosset et al.](https://doi.org/10.1103/PhysRevX.8.021033)
already develops quantum-storage verification. The flagged inequality follows
from those established EB and channel-norm ideas; no novel theorem is claimed.
The full [prior-art matrix](prior-art.md) retains the distinct storage and
recent device-independent targets. A publication-level novelty search remains
incomplete and cannot be replaced by the limited source audit.
