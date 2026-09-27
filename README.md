# Quantum Data Identifiability

**What can quantum experimental data distinguish?**

Reproducible analyses of competing physical models, calibration ambiguities, and
constraints from combined datasets, plus experiment designs for unresolved
questions. Studies include [Wen et al. 2026 photon propagators](studies/wen-2026-propagator/README.md),
[NIST Bell test records](studies/nist-bell-causal-audit/README.md), and
[gate characterization (GST/RB)](studies/shared-quantum-dynamics/README.md).

* **Wen figure reproduction:** reproduces the deposited Figure 3 and 4 summaries.
  Comparing the coherent reconstruction with the deposited theory gives **5.57%**
  symmetric relative error; the paper's **4.45%** remains unreproduced under the
  tested conventions. A different comparison, endpoint data versus coherent
  reconstruction with each normalized to sum to one, gives **4.49%**.
  [Methods and unresolved inputs](studies/wen-2026-propagator/method-comparison.md).
* **NIST raw record reconstruction:** independently decodes **107,109,596 paired
  records** and reproduces the archived settings and click words. Preserving all
  distinct eligible pulse bits instead changes **105 Alice and 94 Bob words**,
  adding 31 and 27 clicks in the primary window. The effect on the published Bell
  analysis has not been assessed.
  [Reconstruction details](studies/nist-bell-causal-audit/reconstruction.md).

Each study documents its data, assumptions, reproducible calculations, and limits
on what the observations identify. The tables separate analyses of existing
records from prospective designs and their synthetic checks.

The companion [ontology-separation](https://github.com/stevenwarejones/ontology-separation)
provides the formal Lean framework for proving predictions from physical
assumptions. This repository focuses on experimental records, numerical inference,
and measurement design. See the companion's
[path interference case study](https://github.com/stevenwarejones/ontology-separation/blob/ab1eeb81119ff3fdcb46a771d9bcdf5e39ea3d29/docs/PATH_INTERFERENCE_CASE_STUDY.md)
for the distinction between equivalent descriptions and observable model separation.

## Analyses of existing data

| Study | Question | Status |
|---|---|---|
| [Collective interference](studies/collective-interference-identifiability/README.md) | Do singleton calibration and three-atom records exclude mixtures of pair-sized groups? | Shared-drift ambiguity, exact pair-control evasion and a prospective class-wide triple control; apparatus-supported dynamics remain open |
| [Wen 2026 propagator](studies/wen-2026-propagator/README.md) | What can reconstructed propagators and endpoint statistics tell us about path descriptions? | Public-data audit available; full reconstruction needs additional records |
| [NIST Bell-record causal audit](studies/nist-bell-causal-audit/README.md) | What conditional influence bounds do public trial records support? | Full raw reconstruction, measured marginals and explicit calibration limits |
| [Born-rule identifiability](studies/born-rule-identifiability/README.md) | Can additional datasets separate probability changes from control and detector errors? | Measured scan–Viviani combination, ordinary qutrit compatibility certificates and conditional distance bounds; no identified probability-rule violation |
| [GST/RB discrimination resource](studies/unmeasured-dynamics-gain/README.md) | Can both sources restrict an operational discrimination resource? | Certified two-sided resource gain, largely detector contrast; fixed-input prediction remains one-sided |
| [Detector-normalized dynamics](studies/detector-normalized-dynamics/README.md) | Does the gain survive detector normalization, aggregate adequacy and stability checks? | Exact detector-retuning obstruction, finite-sample rejection of exhibited tables, and matched-probe stability limits; dynamics gain remains open |
| [Certified GST/RB combination gain](studies/gst-rb-combination-gain/README.md) | Do both families strictly restrict one observable over all stationary qubit CPTP models? | Conditional full-class bound and two physical source-removal witnesses; recorded-context contrast, not unseen-sequence prediction |
| [Correlated qubit mechanisms](studies/correlated-qubit-mechanisms/README.md) | Can complementary controls predict correlated readout response under a different intervention? | Matched-control forcing model tested on shock records; both extensions fail transfer, and adequate physical countermodels remain open |
| [Shared quantum dynamics](studies/shared-quantum-dynamics/README.md) | What do matched GST/RB sequences identify about leakage and finite memory? | Shared CPTP fits, complementary measured constraints and an exact restricted resource obstruction; long-sequence prediction limits retained |
| [Calibrated environmental memory](studies/calibrated-quantum-memory/README.md) | Can matched intervention calibration distinguish environmental quantum memory from system survival? | Cell-region ambiguity; a joint test rejects exhibited points; full classical class unresolved |
| [Collapse compatibility](studies/collapse-compatibility/README.md) | What can interference, mechanical, space and radiation data jointly determine about collapse noise? | Reproduced sources; two conditional feasibility studies and explicit inference limits |
| [Moving-noise feasibility boundary](studies/moving-noise-boundary/README.md) | Can complementary archives constrain one moving dephasing field? | Tested quantum/measurement bounds and public CAL image audit; empirical complementarity gate blocked |
| [Timing structure and identifiability](studies/timing-structure/README.md) | What can cancellation and temporal diagnostics reveal in existing NIST records? | Conditional archive analysis and exact invisibility witnesses; no detected effect and limited sensitivity |
| [Fine timing information](studies/fine-timing-structure/README.md) | What do pulse identity and finer timing add beyond four categories? | Conditional archive bounds and paired recovery; no measured departure |
| [Timing sensitivity and actual settings](studies/synthetic-timing-shifts/README.md) | Could sparse records reveal timing shifts or effects that reverse sign? | Conditional bounds using actual NIST settings, plus separately labeled synthetic and empirical background injections |

The timing sensitivity study includes both analyses of actual settings and
injected effects. Its recovery simulations are sensitivity checks, not observed
experimental departures.

## Prospective designs

These studies specify measurements and calibration requirements, with analytical
bounds and synthetic checks. They do not report experimental exclusions.

| Study | Question | Status |
|---|---|---|
| [Path contextuality](studies/path-contextuality/README.md) | Can a finite weak path probe reject explicit noncontextual representations? | Prospective conditional design; exact instrument and calibrated-trial budget |
| [Sharp compatibility boundaries](studies/path-compatibility/README.md) | Which complete path probe tables admit the specified models? | Exact bounds and attaining models supporting the contextuality design; no measured violation |
| [Phase-intervention design](studies/phase-intervention-design/README.md) | Which calibrated dephased models would a four-phase experiment exclude? | Proposed; synthetic checks and conditional power |
| [Spacetime causal influence](studies/spacetime-causal-influence/README.md) | Can a remote intervention change a complete local record outside its future light cone? | Conditional protocol, formal companion and synthetic sensitivity study |
| [Continuum and finite models](studies/continuum-finite-models/README.md) | Which specified finite dynamics differ from a continuum ring? | Synthetic conditional design with formal companion results; apparatus calibration pending |

## Data and reproduction

The Wen study includes the pinned CC0 Dryad originals for reproducible offline
checks, with the official source, version and SHA-256 manifest retained. Its
[completion report](studies/wen-2026-propagator/completion-report.md) links the
reproduced results and remaining input requirements. Paper PDFs and complete
workbook extracts are not committed. The NIST audit keeps its source inputs
external and pinned by hash; committed tables are derived analysis artifacts.
New source datasets or event extracts require explicit authorization and review
of their source and license before inclusion.

Repository code is distributed under the [Apache 2.0 license](LICENSE); external
data and papers retain their stated licenses. Changes are reviewed through pull
requests before merging.
