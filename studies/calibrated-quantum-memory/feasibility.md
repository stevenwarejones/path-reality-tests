# Feasibility gate: calibrated environmental memory

Decision (2026-09-27, third calculation): **archive the search for a positive exclusion within
the current IBM interval region; retain the draft as a scoped physical
identifiability result. The empirical mission remains incomplete.**
The [new calculation](ibm-identifiability.md) verifies full classical and
quantum models for every selected IBM cell under two conditional mappings.
The [joint follow-up](joint-statistics.md) now excludes those points and two
improved delay-dependent pairs using a valid finite-sample aggregate region.
The full delay-dependent classical class remains unresolved.
No independent matched causal-break calibration has been verified. This is
not a general impossibility theorem or a refutation of the original experiment.
The prior acquisition search below is retained as history, not expanded into
another terminal-QPT inventory.

| Route | Decision and reason |
|---|---|
| Positive claim from this simultaneous cell box | Archive: explicit physical classical points satisfy every constraint. |
| Shared-instrument identifiability | Completed at the declared-region level: exact full-record pair, fixed-point boundary, source deletions and separating probe. No aggregate fit or broad novelty claim. |
| Stronger statistical NMN analysis | Continue: a valid joint region excludes four exhibited points; improved delay-dependent unitary/shared-SPAM pairs are independently physical but still rejected. The larger classical class remains computationally unresolved. |
| New empirical architecture | Reopen only with actual intermediate intervention records and usable control constraints; adjacent terminal calibration/QPT is insufficient. No better qualifying dataset was verified in this update. |

## Target and proposed measurement combination

The target is quantum information carried by the environment across a system
intervention. The null allows arbitrary classical environmental memory and
quantum system propagation. Non-Markovianity or a non-entanglement-breaking
storage channel alone does not answer this question.

| Measurement family | Records and exposures | Apparatus, times and intervention | Role and matching status |
|---|---|---|---|
| NMN multi-time outcomes | `NMN_tomog_rerun.json`: 9 × 324 four-cell rows; opposite-repreparation regrouping gives 8,000 per physical circuit | IBM Perth; delay keys 21.333, 24.889, 28.444 ns; preparation, intermediate measurement/rotation, final measurement; no acquisition timestamps in JSON | Dynamics; usable with explicitly inferred regrouping semantics, not a calibration source |
| NMN UQ outcomes | `NMN_lab_rslts.json`: 324 four-cell rows; same regrouping gives 9,216 per circuit | University of Queensland; 97,97 ns; different apparatus from IBM | Independent dynamics experiment; never pooled with IBM |
| Matched readout information | Original paper Appendix Table 4: six assignment probabilities; calibration counts, denominators, drift records and timestamps not deposited in NMN repository | UQ: first/second measurement separately; IBM: common reported pair | Point estimates in published bootstrap; cannot form an independent finite-data region or bound the whole intervention |
| Matched instrument characterization | Desired: all 18 six-state-input/three-Pauli-output tables, denominators, timing and shared operation identifier | Must characterize the same system-only intervention used in dynamics; preparation and readout anchors also needed | Missing for the prospective certificate; cannot substitute another laboratory's GST |

The first two rows are **alternatives/replications**, not a complementary core
pair. No actual two-source memory certificate is claimed. Derived process
matrices and bootstrap outputs are not additional measured datasets.

## Mathematical bridge and feasibility calculation

For the intended general experiment, shared instruments \(\mathcal I_{b|y}\)
and a process \(W_j\) generate probabilities through the Choi link product
\(p_j(b,c|a,y,z)=W_j\star\rho_a\star J(\mathcal I_{b|y})\star E_{c|z}\),
with the link product contracting matched input/output wires and its partial
transposes fixed by a single Choi convention. The unrestricted nonlinear problem
is not solved here. A shared-instrument physical subfamily is now implemented
and independently certified in the IBM report; identity delay processes are
a witness choice, not a restriction of the full null. Calibration probabilities constrain these same maps through
trusted input states and output effects. Different delays have separate
processes; no stationary generator is imposed. Device identity and matching
times are required evidence, not inferred from similar fidelities.

Before fitting any such model, the first recorded outcome must not depend on
later controls if the earlier operations and acquisition distribution are
shared. The NMN regrouping recovers identical circuit totals in all 3,240
circuits. Under this mapping and independent trials within each circuit, a
simultaneous Hoeffding region gives a **5.7844 percentage-point minimum
worst-case marginal discrepancy** for UQ from any shared past-marginal model.
IBM's nine settings give zero lower bounds with this conservative diagnostic.
Zero alone does not establish compatibility; the later IBM analysis
constructs physical points in a declared exact-binomial cell region. This is a process/instrument stability
and record-semantics gate, not a quantum-memory witness. The original job
records and regrouping implementation are not present, so the interpretation
remains conditional. See [the report](README.md).

Two mechanisms a terminal recovery experiment cannot distinguish are a qubit
preserved in the system through a failed reset, and a qubit stored in the
environment through a working reset. Matched intervention tomography separates
them. [Theory](theory.md) proves a full-class bound for a narrower single-outcome
reset protocol: \(F_6\le 2/3+\epsilon\), where \(\epsilon\) is half the diamond
distance from a completely depolarizing reset. Eighteen calibration tables
provide a conservative global upper bound on epsilon without a nonconvex fit.
At 8,000 shots per setting this can resolve an ideal SWAP-memory example.
This is **synthetic feasibility**, not a calculation of experimental power for
NMN, whose operations differ and whose calibration is missing.

## Alternative acquisitions checked before narrowing

1. **44hsiang/quantum_memory_data**, commit
   `25459810b88553dc5bc8641956e9a72e175008af`: follow-up now reconstructs
   measured q0/q2 IQ calibration node 6047 and q0 terminal QPT nodes 6048/6099,
   including HDF5 counts, raw Bloch vectors, parent IDs, timestamps and stored
   readout-operation equality. This is a real adjacent readout/QPT match;
   it does not characterize a middle instrument. Scheduling fields differ,
   classifier validation is unverified, and no numerical transfer bound is
   supplied. See the [exact paths and source-matching table](acquisition-semantics.md).
   Derived GST/QPT/ellipsoid estimates and simulations are not independent
   measured sources. The [APS 2026 abstract](https://meetings-archive.aps.org/smt/2026/mar-g07/14/)
   matches the storage/geometric objective; an exact journal-publication
   association remains unverified. Source files stay external; no explicit
   repository license was inferred.
2. **White et al., PRX 15, 021047 (2025)** and its linked
   [public code](https://github.com/gwhitequantum/process-tensor-network-tomography/tree/e3f01b011f28c1497f7ebb3dc4f98a782a75095f):
   follow-up safely decodes all 20 SU4 training/validation float lists and
   inventories representative DD/optimizer pickle opcodes. The lists contain
   terminal probabilities, not original counts, setting/outcome keys or
   matched intermediate instrument tomography. Two training lists have 11999
   entries; supplied preprocessing can omit absent outcomes, so reshaping and
   zero imputation without keys are unjustified. DD object files were not
   executed or fully decoded, and this is not an exhaustive archive search.
   Exact paths, versions and limitations appear in the
   [focused schema audit](acquisition-semantics.md). Existing self-consistent
   fits do not supply the desired global classical-memory exclusion.
3. **Aloy et al., Nature Communications 17, 2474 (2026)**:
   [primary paper](https://www.nature.com/articles/s41467-026-69030-x) describes
   terminal prepare-and-measure delay tables and readout calibration (2,000
   shots per prepared computational state). Its
   [Zenodo archive](https://doi.org/10.5281/zenodo.17979407) was identified but
   not downloaded. The stated protocol lacks a mid-evolution causal break and
   therefore does not implement this certificate. No claim is made about every
   file in that archive. It is an independent comparison, not borrowed calibration.

## Bounded intervention-architecture check

A short fallback search considered actual intermediate interventions, rather
than more terminal QPT. The ranking below is for further acquisition checking,
not a claim that any candidate already supplies a usable joint dataset.

| Rank / architecture | Multi-time intervention and records | Flagged instrument / matched controls | Boundary and complementarity | Decision |
|---|---|---|---|---|
| 1. Guo et al., [photonic common-cause process](https://arxiv.org/abs/2003.14045) | Experimental instrument-specific memory and process recovery; intermediate operation choice is part of the reported experiment. Public raw outcome tables were not verified in this bounded check. | An independent finite-count region for outcome probabilities and postmeasurement states, matched to SPAM and epoch, was not verified. | A distinct photonic architecture; require explicit surviving optical modes, heralds and registers. Which measured family breaks instrument/process ambiguity remains unresolved. | Most relevant next acquisition check; do not promote the paper's reconstructed process to raw calibration. |
| 2. Gong et al., [five-qubit temporal correlations](https://doi.org/10.1038/s41534-020-00329-3) | Intermediate reset/probe construction and final postselected tomography; paper links [Figshare 12715493](https://figshare.com/articles/dataset/12715493). Linked data were not successfully retrieved. | No matched independent flagged-reset characterization was verified. | Probe ancillas, postselection and controller records must be explicit; final probe correlations do not themselves certify the reset's error. | Potentially useful intervention architecture, blocked at record verification. |
| 3. Gosling et al., [fluxonium/TLS quantum-jump monitoring](https://arxiv.org/html/2603.11889v1) | Stroboscopic readout every 6 microseconds; [Zenodo 21908500](https://zenodo.org/records/21908500) lists a trace archive, analysis notebook and Solomon-equation code. The 213.8 MB trace archive was not downloaded. | Readout time series and IQ plot are reported; no complete instrument-state/SPAM confidence region was established here. | The paper explicitly finds a semiclassical population-rate model sufficient; no qubit manipulations provide independent coherence probes. | Archive for this quantum-environmental-memory target. It is useful evidence of memory, not the needed distinction from classical environmental records. |

No author was contacted and no cross-device calibration was borrowed. This
check does not establish that the first two archives lack the required files;
it records the limits of inspection. Neither becomes a substitute empirical
result for the IBM calculation. A pivot should proceed only after the required
records are actually recovered and a preliminary full-record feasibility test
shows that the new operations constrain the relevant ambiguity.

## What would reopen the empirical gate?

The [exact acquisition specification](missing-acquisition.md) now freezes a
q0 cut, 18 flagged calibration tables and 18 full dynamics tables per matched
block/delay pair, with all outcomes, SPAM anchors and transfer controls.
The [flagged null](flagged-instrument.md) incorporates every available record
and feed-forward. Neither has been empirically certified from this subset.


For NMN: original physical-circuit records or documented invertible regrouping;
timestamps/interleaving to investigate the UQ past-marginal variation; matched
instrument/readout/preparation calibration counts; and a justified bound on
leakage and intervention–environment coupling. Then quantify whether the
classical-memory-required instrument error exceeds the calibrated allowance.
The UQ 0.057844 bound concerns context dependence, **not** the reset epsilon in
the prospective theorem; the two numbers must not be compared as if identical.

There is no empirical null exclusion, no verified empirical combination gain,
and no claim of novelty for the prospective bound. Missing files are not an
identifiability theorem. The earlier restricted ambiguity in [theory](theory.md) is synthetic.
The later [IBM ambiguity](ibm-identifiability.md) is fitted to all selected
measured counts, but only certifies membership in its declared interval region;
its substantial aggregate residuals prevent a stronger empirical claim.
