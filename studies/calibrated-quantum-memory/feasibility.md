# Feasibility gate: calibrated environmental memory

Decision (2026-09-27): **narrow to an acquisition audit and a prospective,
calibration-robust certificate. The empirical mission remains incomplete.**
Repository base inspected: `db1a02dcab5cc12107cc86fdb371dfba95353cc2`.
No independent matched causal-break calibration has been verified in the
sources inspected below. This is a bounded search result, not an impossibility
theorem or a refutation of the original experiments.

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
transposes fixed by a single Choi convention. This general nonlinear problem
is not implemented here. Calibration probabilities constrain these same maps through
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
Zero does not establish compatibility. This is a process/instrument stability
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
   `25459810b88553dc5bc8641956e9a72e175008af`: downloaded the archive and inspected
   the README, figure scripts, simulation code and representative `node.json`
   acquisition records. The 152 ns family contains 245 node records (QPT,
   ellipsoid and interaction-control acquisitions). For example node 6099
   specifies q0, 10,000 runs, a 102 ns interaction and a 2026-02-21 timestamp.
   The GST/ellipsoid family has 189 node records; iSWAP T1/T2 has two dated
   2026-03-02. These are not automatically contemporaneous causal-break
   calibrations. The [APS 2026 abstract](https://meetings-archive.aps.org/smt/2026/mar-g07/14/)
   by Shih-Hsiang Chao et al. matches the geometric/storage objective, but the
   README supplies no paper DOI; an exact journal-publication association was
   not established. Do not assert that the conference abstract proves that
   linkage. Inspected protocols concern terminal channel characterization;
   no matched intermediate instrument calibration was verified. No HDF5 count
   inference or environmental-memory conclusion is made. No explicit repository
   license was found in the inspected top-level files; source data stay external.
2. **White et al., PRX 15, 021047 (2025)** and its linked
   [public code](https://github.com/gwhitequantum/process-tensor-network-tomography/tree/e3f01b011f28c1497f7ebb3dc4f98a782a75095f):
   downloaded and inventoried the archive, read the paper's self-consistency,
   gauge and spillage treatment, and inspected file structure. It includes
   fitted process tensors, validation statistics, circuit job files, RB and
   decoupling families. These are promising for a dedicated adaptation, but
   no independent matched reset-calibration dataset was verified here. This
   source is **not exhausted**; a complete job-schema and instrument analysis
   remains open. Existing self-consistent fits do not supply a global composite
   classical-memory exclusion. Pickles were not executed.
3. **Aloy et al., Nature Communications 17, 2474 (2026)**:
   [primary paper](https://www.nature.com/articles/s41467-026-69030-x) describes
   terminal prepare-and-measure delay tables and readout calibration (2,000
   shots per prepared computational state). Its
   [Zenodo archive](https://doi.org/10.5281/zenodo.17979407) was identified but
   not downloaded. The stated protocol lacks a mid-evolution causal break and
   therefore does not implement this certificate. No claim is made about every
   file in that archive. It is an independent comparison, not borrowed calibration.

## What would reopen the empirical gate?

For NMN: original physical-circuit records or documented invertible regrouping;
timestamps/interleaving to investigate the UQ past-marginal variation; matched
instrument/readout/preparation calibration counts; and a justified bound on
leakage and intervention–environment coupling. Then quantify whether the
classical-memory-required instrument error exceeds the calibrated allowance.
The UQ 0.057844 bound concerns context dependence, **not** the reset epsilon in
the prospective theorem; the two numbers must not be compared as if identical.

There is no empirical null exclusion, no verified empirical combination gain,
and no claim of novelty for the prospective bound. Missing files are not an
identifiability theorem. The exact restricted ambiguity in [theory](theory.md)
is a constructive mathematical example and is explicitly not fitted to NMN.
