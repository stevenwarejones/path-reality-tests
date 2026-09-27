# IBM: a physical shared-instrument ambiguity, with a statistical limit

**Decision: archive the search for a positive memory claim within the current interval region.**
A focused numerical search now supplies explicit classical and quantum-memory
models satisfying every selected IBM count constraint, with the same shared
instrument structure across all nine delay pairs. An exact operator identity
makes their complete flagged probabilities identical. This is an empirical
confidence-region compatibility result and a defined-family identifiability
limit, not an empirical detection of quantum memory or a general no-go theorem.

**The statistical limitation is material:** the saved points have substantial
aggregate likelihood residuals. Passing the simultaneous cell intervals does
not certify an adequate joint likelihood fit. A stronger confidence region
could reject these particular points, so this update does not claim that all
statistically efficient analyses of NMN must remain ambiguous. The exact
physical equivalence itself holds without sampling noise.

## Selected records, sharing and assumptions

We retain **all 9 × 324 × 4 = 11,664 IBM cells** in the pinned
`NMN_tomog_rerun.json`; no outliers, outcomes, axes or delays are removed.
UQ and the separate q0 apparatus are not calibration sources for IBM.
The [existing semantics audit](acquisition-semantics.md) remains authoritative.
No original job-to-archive transformation has been newly verified.

| Item | Construction and status |
|---|---|
| Preparations and final effects | Six ideal Pauli states and three Pauli measurements form one admitted physical SPAM choice. This is a witness choice, not a calibration inferred from labels. Enlarging a nuisance region that contains this choice retains the witnesses. |
| Middle settings | y=(first axis m, fixed re-preparation rotation label p), 18 settings; two **retained** outcomes b. Each setting has two CP qubit maps with TP sum. |
| Shared operations | One instrument per y is used across all inputs, final axes and nine delays. The first POVM effect for fixed m,b is exactly shared across all p. No dependency on a later final axis or secret preparation label is allowed. |
| Processes at different delays | The null permits separate arbitrary processes. The exhibited classical point happens to choose identity propagation at every delay. Equality is a property of this one witness, not an imposed dynamical law or a conclusion about the apparatus. |
| Flags and feed-forward | Both outcomes b are retained in the probability tables. The allowed null includes arbitrary classical environmental records and b-dependent controls. Our classical witness needs neither environmental memory nor feed-forward. No postselection discards b. |
| System boundary | Both constructions use a qubit system, no leakage, and system-only instruments. The quantum model additionally uses an environmental qubit and a classical branch selector. No preparation label leaks to them. Allowing more imperfections cannot remove these admitted witnesses. |
| Unmeasured control assumptions | Preparations/readout, all instrument CP maps, and the process split are candidate models. No independent IBM calibration-count region certifies those controls. Distances between the two exhibited instruments are not distances from the real apparatus. |

We verify two interpretations separately:

* **Inferred first-bit-1 regrouping:** every physical row has 8000 trials.
* **Literal rows:** row totals are treated as predetermined, outcome-independent
  trial counts. This conditional interpretation uses each actual row total;
  it does not claim that retrospectively assembled rows really are multinomial.

The first-bit-0 regrouping is obtained from the first construction by swapping
p with its opposite. It has identical constraints, probabilities and bounds
under that relabeling. Second-character regroupings conflict with the
notebook's declared analysis bit order and the fixed-shot inference; they
remain acquisition alternatives in the prior audit, not verified physical
sampling models. No exclusion is selected from a convenient mapping. Because
the uncertainty union already contains a physical classical point, a uniform
exclusion over the larger mapping union is unavailable for this region.

## Preliminary calculation before the certificate

A 66-parameter physical model with per-delay unitary propagators and shared
assignment errors fits the Hoeffding box (largest residual 0.027224 versus
radius 0.028563), but not the tighter exact-binomial box. Adding preparation
and output-state errors did not find a point in that smaller restricted
family. These failed searches are **computational uncertainty**, not exclusions.

We then made the shared instruments general CP maps and solved a small convex
feasibility problem. Each flag map has a 4×4 Choi matrix J_(b|y); impose
J>=0, sum_b Tr_out J=I, and common first effects across p. The known probe
probability is Tr[J_(b|y)(rho_a^T ⊗ E_(c|z))]. The optimizer only finds a
candidate. The committed coefficients are rationalized and repaired to exact
trace preservation/shared effects; an independent verifier decides whether
that candidate is usable. No PPT relaxation point is called classical.

| Conditional interpretation | Minimum cell-interval slack | Quantum branch weight w | Exact NPT test expectation |
|---|---:|---:|---:|
| Inferred regrouping | 0.000373783224 | 0.003 | −3/8000 |
| Literal rows | 0.000374024277 | 0.003 | −3/8000 |

All nine delay pairs are checked, not pooled into a single multinomial trial.
The common candidate probability must satisfy every run's interval separately.
The decimal slack is descriptive; acceptance does not depend on its precision.

## Exact physical equivalence: instrument transmission or environmental memory

Let K_(b|y) be any qubit instrument and let tau=I/2. Write
q_(b|y)=Tr[K_(b|y)(tau)]. For 0<w<1 define

\[
 B_{b|y}(X)=(1-w)K_{b|y}(X)+w q_{b|y} X. \tag{1}
\]

B is a valid flagged instrument: its blocks are CP and their sum is TP.
Consider two full physical models:

1. **Classical/no environmental memory:** identity evolution before and after
   the middle instrument B. A fresh local ancilla can implement B and is
   discarded within the instrument. No environment carries a quantum state
   across the cut.
2. **Quantum environmental memory:** use instrument K. An independent branch
   selector chooses, with probability w, to SWAP the incoming state into an
   environmental qubit initially in tau, apply K to the resulting system
   state tau, then SWAP back. On the other branch, propagate the system
   directly with no SWAP. The selector remains available to control the
   second SWAP. Trace out the environment after that operation.

In the quantum branch, outcome b occurs with probability q_(b|y), and the
returned system state is X. Thus both constructions implement exactly the
same subnormalized output map (1), for **every input operator, setting and
flag**, including inputs correlated with a reference. All final POVMs,
coarse-grainings and retained outcome tables agree. Arbitrary fixed pre/post
unitaries could be included in both constructions; they are unnecessary for
the committed witnesses.

This is a physical ambiguity, not an unphysical similarity transform. Given
the fitted classical instrument B, our quantum instrument is explicitly

\[
 K_{b|y}=\frac{B_{b|y}-w q_{b|y}\operatorname{Id}}{1-w},\qquad
 q_{b|y}=\operatorname{Tr}[B_{b|y}(I/2)]. \tag{2}
\]

The inverse is not automatically CP; **the verifier proves CP for every
retained block**. Both models share their own instruments across all delays
and have identical observed probabilities. Their instruments differ: knowing
only the multi-time probabilities cannot determine which instrument was used.

The underlying quantum process is genuinely nonclassical as a process tensor.
In wire order A,B,O,C (input, pre-instrument system, post-instrument system,
final output), its unnormalized Choi operator is

\[
 W_w=(1-w)\Phi_{AB}\otimes\Phi_{OC}
       +w\Phi_{AC}\otimes(I_B/2)\otimes I_O,
\]

where Phi=|00+11><00+11| and Tr W=4. The normalized test vector
v=(|0001>−|1000>)/sqrt(2) gives

\[
 \langle v|(W_w/4)^{T_{AB}}|v\rangle=-w/8<0.
\]

Classical-memory combs are temporally separable and hence PPT; this negative
expectation certifies the resource of this **constructed process**, with its
specified slot architecture. It does not certify the actual device's memory,
nor exclude the larger allowed null on these data: model 1 already fits them.
No equivalence of PPT with classical memory is used.

## Certified boundary and which nuisance direction matters

For a fixed positive definite fitted Choi block J_b, let v_I=vec(I) and
q_b=Tr J_b/2. The rank-one update criterion gives the exact endpoint

\[
 J_b-wq_b v_Iv_I^\dagger\succeq0
 \quad\Longleftrightarrow\quad
 w\le [q_b v_I^\dagger J_b^{-1}v_I]^{-1}.
\]

Taking the minimum across all 36 blocks gives 0.003000204821 for the regrouped
point and 0.003000203895 for the literal point; the stored fractions are exact.
Every w from zero to this endpoint (and below 1) gives an observationally
identical physical family; strictly interior w has the exact PD certificate.
These are **endpoints for the fixed exhibited instruments**, not globally
optimized upper bounds on memory allowed by the data. We do not interpret a
failed larger-w optimization as a resource bound.

The relevant nuisance direction is coherent identity transmission with a
flag drawn from q_b. Its full flagged half-diamond distance satisfies

\[
 \tfrac12\|\widehat B-\widehat K\|_\diamond
 =w\tfrac12\|\widehat F-\widehat K\|_\diamond\le w,
 \qquad F_b(X)=q_bX.
\]

This bound covers arbitrary reference correlations and every state the null
could send to the instrument. It is not inferred from six isolated probes.
The selected isolated-probe TV differences below reach 98.24% and 99.21% of
that ceiling: for this measured-fit direction, the generic norm allowance is
already close to tight. No witness on the existing full tables can distinguish
this pair; their probabilities agree exactly. A different witness-sensitive
bound cannot fix this particular ambiguity unless its calibration constraints
exclude at least one instrument.

## A feasible discriminating measurement and source deletion

Isolate the middle instrument on the same device/epoch, retaining its outcome
and preventing return through the external environment. Use the external
state preparation and final readout anchors declared above, or account for
their independent uncertainty.

| Interpretation | Input / first axis / rotation label / final axis | Event | Difference in probability |
|---|---|---|---:|
| Regrouped | Z+ / z / Z− / z | final c=0, sum over b | 0.002947284930 |
| Literal | Z− / y / Z+ / z | final c=1, sum over b | 0.002976431573 |

These are predictions of the fitted physical models, **not measured new data**.
The exact event probabilities differ even though every original multi-time
probability agrees. About 1.01 million or 0.99 million shots per candidate
would make two Hoeffding 95% radii sum to less than the corresponding predicted
gap. This is a conservative precision calculation, not a power guarantee;
actual separation needs additional shot margin and SPAM/drift error below the
gap. The earlier acquisition specification remains the larger full-class test.

As concrete internal source deletions, retain only m=z settings (3888 cells),
or only m=x,y (7776 cells). **The same two complete physical models** satisfy
every retained constraint under either deletion, with the original 11664-cell
confidence budget. More generally, any subset of these same single-slot
operations preserves the exact equality (1). Calling a subset “calibration”
cannot resolve the ambiguity within this model family. This is not a claim
that every conceivable operation or independently anchored calibration is
uninformative.

## Statistical coverage, calibration consistency and limitations

At alpha=0.05 use two-sided Clopper–Pearson intervals for each retained
binomial cell, assigning each tail alpha/(2M), M=11664. The union bound gives
simultaneous coverage at least 95% if trials are independent with stable
probabilities within each physical row and the chosen mapping/denominators
are correct. Cross-row independence is unnecessary. Each mapping has its
own conditional region; their union covers the correct interpretation without
selecting a mapping to obtain a violation. No alpha is reassigned on deletion.

The verifier uses exact rational LDL positivity, exact TP/shared-effect
identities, and **downward-rounded partial binomial-tail sums compared with an
upward-rounded threshold**. It does not trust the SDP's status or floating-point
inverse-CDF values. Saved decimal CP slacks are display-only.

Two limitations prevent stronger claims:

* **Aggregate fit:** the multinomial deviance from saturated per-row
  probabilities is 17570.91 (regrouped) or 16409.19 (literal), versus 8748
  saturated probability coordinates. Even unrestricted probabilities shared
  across delays give substantial residuals (saved in the result). These are
  descriptive diagnostics, not calibrated composite-null p-values. They make
  clear that interval-box compatibility is weaker than a satisfactory aggregate
  fit. A stronger region and more general delay-dependent classical models
  remain an unresolved calculation; failed restricted searches do not exclude
  them. We do not claim universal finite-data indistinguishability of NMN.
* **Published assignment points:** if Table 4's IBM values 0.976/0.973 are
  imposed as an exact, shared forward misclassification channel, every
  first-bit-0 joint event has probability <=0.976. The regrouped row
  `zp,z,zp,z`, delay `28.444,21.333`, has 7948/8000 for 00. The saved directed-
  rounding binomial-tail bound rejects that fixed ceiling within the same
  simultaneous budget. This would be a calibration/record-consistency failure
  for quantum and classical models alike. No independent calibration trials,
  uncertainties or epoch linkage justify imposing those point values as exact.
  We have not fitted them as additional count records or silently treated them
  as trusted error bars.

## Prior work and scientific increment

Self-consistent instruments, SPAM gauge freedom and uncertain process
reconstruction are established: [Li et al., QST 9 (2024)](https://doi.org/10.1088/2058-9565/ad3d80),
[White et al., PRX 15 (2025)](https://doi.org/10.1103/PhysRevX.15.021047), and
[Taranto et al., Quantum 8 (2024)](https://doi.org/10.22331/q-2024-05-02-1328)
are direct comparisons. The temporal entanglement witness follows
[Giarmatzi–Costa (2021)](https://doi.org/10.22331/q-2021-04-26-440).
The original [NMN paper](https://arxiv.org/html/2308.00750v3) uses a different
trusted-operation reconstruction and bootstrap; our conditional box analysis
does not invalidate its reported results.

The increment is a **full selected-record physical certificate**, including
shared first effects, two mapping interpretations, an exact flagged
instrument/environment exchange, a fixed-model robustness endpoint and a
concrete separating probe. No novelty is claimed for gauge ambiguity or the
SWAP mechanism. No synthetic count family was added. Whether this constitutes
a consequential empirical ambiguity under a stronger statistical analysis
remains open; the original positive research mission is not complete.

## Reproduction

```bash
python -m pip install -r studies/calibrated-quantum-memory/requirements.txt
python studies/calibrated-quantum-memory/audit.py --source-dir /tmp/nmn --download --check
python studies/calibrated-quantum-memory/identifiability.py --source-dir /tmp/nmn --check
python -m unittest discover -s tests -p 'test_calibrated_quantum_memory.py' -v
```

The verifier uses the pinned external IBM counts and
[rational model coefficients](results/ibm-instrument-witnesses.json).
The [result](results/ibm-identifiability.json) records each delay, all deletion
checks, physicality, exact endpoints, aggregate diagnostics and the probe.
Optional discovery is separate from acceptance:

```bash
python -m pip install -r studies/calibrated-quantum-memory/requirements-fit.txt
python studies/calibrated-quantum-memory/fit_identifiability.py --source-dir /tmp/nmn --output /tmp/candidate.json
python studies/calibrated-quantum-memory/identifiability.py --source-dir /tmp/nmn --witness-file /tmp/candidate.json --output /tmp/candidate-result.json
```

A new optimizer run need not reproduce identical coefficients. It must pass
all independent checks. No author contact, new hardware acquisition or merge
was performed.
