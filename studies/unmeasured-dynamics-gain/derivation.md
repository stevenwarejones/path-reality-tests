# From count constraints to an unmeasured prediction

## Declared class and statistical region

A single qubit is reset to a shared density operator `rho`. Three arbitrary stationary CPTP maps `Gi,Gx,Gy` act in written order; the terminal binary POVM is `{E,I-E}`, shared by all words. There are no intermediate measurements in these records. Thus

\[
p(w)=\operatorname{Tr}[E\,G_{w_L}\cdots G_{w_1}(\rho)].
\]

Preparation, readout and all channel entries are nuisance parameters. No unitality, gate ideality, ideal inverse, computational readout contrast, or RB depolarizing law is assumed. The theorem ranges over the whole physical class, including singular boundary states and channels. Numerical searches use the inherited 52-factor-coordinate box with three redundant process scales fixed, but **that box is not used to certify coverage of the physical class**. The analytic proof does not rely on its completeness.

We reuse `gain_core.region` without modification: simultaneous two-sided Clopper–Pearson cell intervals for 6,661 rows consume error probability 0.025; simultaneous Hoeffding intervals for 12 predetermined, shot-weighted groups consume 0.025. Individual cell tail budget is `1/(80*6661)`; group one-sided budget is `1/960`. Overlapping groups need not be independent of each other; independence of individual Bernoulli shots under the fixed model is the relevant assumption. Endpoints are conservatively rounded with exact rational checks inherited from PR17. Counts, denominators, row identity and order are unchanged.

A source deletion uses exactly the same endpoints for the retained cells and groups. The region was fixed before this new target screen, but the analysis is retrospective. On the one simultaneous coverage event, all inequalities below hold for every target to which their hypotheses apply. Data-dependent choice of a target cannot invalidate that implication. It does not establish held-out predictive performance, remove discovery bias, or establish novelty. Unrestricted within-row dependence removes the asserted coverage.

## Qubit composition lemma

For a prefix `w` and suffix `f`, write `sigma=G_w(rho)`, and

\[
p_0=\operatorname{Tr}E\rho,\quad p_w=\operatorname{Tr}E\sigma,
\quad p_f=\operatorname{Tr}E G_f(\rho),\quad
p_{wf}=\operatorname{Tr}E G_f(\sigma).
\]

Let `l <= u` be the eigenvalues of `E`. If `u=l`, all four probabilities equal that constant and the assertion follows. Otherwise write `E=l I+(u-l)P`, where `P` is a rank-one projector, and set `q0=(p0-l)/(u-l)`, `qw=(pw-l)/(u-l)`.

For two qubit states with those fixed projector probabilities, the maximum trace distance is

\[
D(\rho,\sigma)\leq
\sqrt{q_0(1-q_w)}+\sqrt{q_w(1-q_0)}.
\]

To prove it, choose the Bloch z axis along `P`. The z components are `2q0-1` and `2qw-1`; each transverse length is at most `2 sqrt(q(1-q))`. Opposite maximal transverse components maximize the Euclidean distance. Substitution into `D=|r-s|/2` gives the displayed expression. Pure states attain this geometry. This uses dimension two.

The dual CPTP map is positive and unital, so `l I <= G_f*(E) <= u I`. For a traceless Hermitian difference of states,

\[
|p_{wf}-p_f|\leq (u-l)D(\rho,\sigma)
\leq \sqrt{(p_0-l)(u-p_w)}+\sqrt{(p_w-l)(u-p_0)}
\leq h(p_0,p_w),
\]

where

\[
h(a,b)=\sqrt{a(1-b)}+\sqrt{b(1-a)}.
\]

The last inequality follows termwise from `0<=l<=p0,pw<=u<=1`. If `a+b<=1`, `h` is nondecreasing in both arguments throughout `[0,a] x [0,b]`; equivalently write it as `sin(arcsin(sqrt(a))+arcsin(sqrt(b)))` on that domain. Therefore `p0<=a`, `pw<=b`, `pf>=c` imply

\[
p_{wf}\geq \max(0,c-h(a,b)).
\]

Every full physical model in the joint confidence region maps into these three scalar interval constraints, so this is an **outer** bound. No scalar-feasible point is promoted to a physical witness.

## Exact application and numerical allowances

Take `w=GiGi`, `f=GyGy`. The pinned rows are:

| Global row | Source / literal word | Plus / total | Endpoint used |
| --- | --- | --- | --- |
| 0 | GST / empty | 0 / 50 | upper 0.2318121895 |
| 10 | GST / `GyGy` | 50 / 50 | lower 0.7681878105 |
| 4662 | RB / `GiGi` | 0 / 285 | upper 0.0452128544 |

The new target hash is `dcc185c9119d8ed71fca5a402c4353202bad7f19493372053cf68bf8c5186350` under `sha256(bytes([0,0,2,2]))`. It is absent from every recorded expanded word, including duplicates and compressed expressions.

Let the endpoint integers on the inherited `D=10^10` grid be `A,B,C`. Integer square roots give

\[
L=\frac{C-\lceil\sqrt{A(D-B)}\rceil-\lceil\sqrt{B(D-A)}\rceil}{D}
=0.1113640705.
\]

This is downward conservative with **no floating-point square root in the certificate**. It is not claimed sharp. The bound involves probability predictions, not observed leakage populations or a probability-law deformation.

Witness factors are the exact dyadic real numbers stored in JSON. Reused directed-interval calculations certify their normalized Kraus construction and bound the difference between exact physical predictions and the Bloch propagation kernel by at most `1e-9` through length 8,198. The four-gate target lies within that envelope. Saved probabilities are used for exact rational retained-constraint arithmetic. Full-source reproduction recomputes them and adds any artifact/kernel discrepancy to the kernel certificate before accepting the common envelope. Independent complex density-superoperator propagation provides a separate implementation check; it is not substituted for the rigorous interval envelope.

Outward `1e-8` target quantization, plus `1e-9` on either side, stabilizes the report across numerical-library implementations. Separation uses the upper endpoint of the countermodel and the downward analytic lower bound. Every retained inequality has strictly positive margin after the envelope. Local search stopping flags are not certificates of minima.

## Aggregate and composition controls

Independent variables for all recorded word probabilities form an outer relaxation of the physical model. At the joint feasible point the measured variables obey every count constraint. The unrecorded variable can be independently set to zero or one, proving its full interval is [0,1] in that relaxation. This demonstrates the need for composition, not a physical source-only construction.

For the aggregate-only physical countermodel, set `Gi=I`, `Gx=Rx(pi+delta)`, `Gy=Rx(pi-delta)`, and pure +Z preparation. For a word with `nx,ny` occurrences, its ideal-projective plus probability is

\[
B_w=\sin^2\!\left(((n_x+n_y)\pi+(n_x-n_y)\delta)/2\right).
\]

A diagonal effect produces `p_w=l+(u-l)B_w`. At `delta=0.005189999999999999`, solve two linear equations for `l,u` to match the two overall observed mean centers. The result is `l≈0.05941543077`, `u≈0.99996377502`, both physical. The saved normalized Kraus factorization includes the inherited tiny positive regularizer, so **the saved factors, not an ideal-unitary approximation, define the certified model**. Its complete word probabilities are recomputed and the two retained overall-mean constraints checked. The construction is well inside the aggregate region but does not satisfy all detailed records.

## Assumption enlargements

A qutrit cyclic shift `V|j>=|j+1 mod 3>` with `rho=|0><0|`, `E=|2><2|`, `Gi=V^2`, `Gy=V` has `p0=pII=0`, `pYY=1` and `pIIYY=0`. It violates the qubit lemma while remaining stationary and CPTP. It only refutes dimension-independent use of the lemma; it is not a full-data leakage explanation.

For an explicitly hypothetical context model, suppose each realized gate occurrence is within half-diamond distance `delta` of a fixed stationary qubit gate, with common exact SPAM. Circuit-dependent operations are fixed across shots, preserving Bernoulli sampling within each circuit. A length-L word changes probability by at most `L delta` by telescoping and trace-distance contractivity. The reference probabilities satisfy `pII<=b+2delta`, `pYY>=c-2delta`; the actual target differs from its reference by at most `4delta`. Thus, where `a+b+2delta<=1`,

\[
p_{IIYY}^{\rm actual}\geq\max\{0,c-6\delta-h(a,b+2\delta)\}.
\]

The grid values are sensitivity hypotheses, not measured drift or memory bounds. Random correlated drift across shots requires a different confidence construction. Finite classical memory and unrestricted qutrit channels are outside the primary theorem; no resource ordering between them is asserted.

## Exact all-word equivalence and the resource task

In trace/Bloch coordinates, `S_s=diag(1,s,s,s)` preserves the trace coordinate. For any word,

\[
(E S_s^{-1})(S_s G_{w_L}S_s^{-1})\cdots
(S_s G_{w_1}S_s^{-1})(S_s\rho)=E G_w\rho.
\]

This is exact cancellation for all finite words, not a small-error or finite-data approximation. `S_s` need not itself be a quantum channel. Instead we independently verify physicality of **each transformed native channel**, state and effect at `s=1` and `s=1009/1000`.

`resource_certificate.py` reconstructs exact normalized Kraus maps with 60-digit directed interval arithmetic, forms their Bloch matrices, and applies the rational similarity. Choi matrices use `J=sum_j sigma_j^T tensor G(sigma_j)/2` (input subsystem first); their trace is two. The trace row fixes TP exactly, and interval LDL verifies `J-10^-6 I` positive definite. Preparation length is below one; both effect eigenvalues are in [0,1]. Hence both transformed models are physical and inherit every joint count constraint by the exact identity above.

For the resource target, repeat `Gy` **10,749** times, with interval exponentiation by squaring. This computation has its own interval enclosure and does **not** reuse an 8,198-gate scalar error budget. The trusted Bell-input output is `J/2`. At `s=1`, interval LDL verifies `(J/2)^T_output -3e-6 I` positive definite. At `s=1009/1000`, a stored fixed vector gives an interval Rayleigh quotient strictly below `-3e-6`. For two qubits, PPT is equivalent to separability; a channel is entanglement breaking exactly when its Choi state is separable. These established results yield the two different classifications.

The fixed trusted Bell input is outside the original preparation alphabet and is **not** transformed along with `S_s`. That is why the enlarged operational task can distinguish the two realizations even though all original terminal-word probabilities coincide. If the ancillary preparation and measurement were also unknown and transformed, this argument would not certify an observed resource.

For a trusted Z measurement of the reset, the plus-probability difference is `(s-1)c/2`, interval-certified between 0.0044564 and 0.0044565. This identifies `s` on this specific orbit when the measurement is externally fixed. It is not a general informational-completeness claim. No unique resource identification, general no-gain theorem, or robust memory advantage follows from the example.
