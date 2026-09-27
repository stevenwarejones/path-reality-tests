# Continuum dynamics and finite alternatives

**Synthetic prospective study; no experimental
exclusion.** A finite spectral model reproduces restricted continuum experiments
exactly. A specified nearest-neighbor spatial lattice instead changes relative
mode phases, with a detectable difference only under adequate access and calibration.

The [model](model.md) gives the circle dynamics, spectral counterexample,
approximation/tail derivation, finite-resource quantifiers and exact blind spots.
The executable [decision rule](decision.py) converts complete calibration/main
counts into rejected candidate site counts with coordinate witnesses. The
[shared-nuisance refinement](joint.py) additionally enforces the common scale
and phase across momenta, with six target-gap refinements and [certified sensitivity bounds](results/shared-sensitivity.json).
The [protocol](protocol.md) specifies a common 24-setting experiment, complete
outcomes, nuisance boxes, simultaneous rejection and conditional acquisition cost.
The [literature review](literature.md) assesses primary sources and dataset fit.
No source data are used in the generated [sensitivity map](results/sensitivity.json).

![Synthetic robust coordinate separation](results/sensitivity.svg)

Run from any working directory:

```sh
python /path/to/repository/studies/continuum-finite-models/design.py --check
python /path/to/repository/studies/continuum-finite-models/joint_design.py --check
```

From the repository root, `python -m unittest discover -s tests -v` includes
independent matrix diagonalization/evolution, complete POVM tables, convergence,
stable tiny-tail control, calibration boundaries, strict count-based decisions,
complete failures, heterogeneous joint distributions and sampling
checks. The suite contains 129 tests, including 30 continuum tests; the new
connection checks use matrix exponentials, nonunit hbar and arbitrary states
outside the accessible measurement subspace. Running `design.py` without `--check` regenerates JSON and SVG. Source
hashes identify the numerical code and its model/protocol assumptions; the Git
commit identifies the complete revision. Floating values compare with declared
tolerances; integer counts and categorical decisions compare exactly.

The map covers 660 synthetic parameter combinations. A zero **certified** gap
is inconclusive unless an overlap construction applies. It is not labeled as
an optimized minimum over jointly shared nuisance parameters. The separate
shared-scale study recovers certified gaps for some of these inconclusive rows.
Calibration and main trials are reported separately. The total eligible
count includes 24 main and 3 calibration strata. Total source attempts remain
unknown because preparation efficiency and physical scale/phase certification
are missing. Prospective unconditional power is at least 87.50%, conditional
on the declared external certificates, at type-I risk at most 5.00%.

## Completion/disposition register

The machine-readable [register](questions.json) records dependencies and next
steps. A `proved` status refers to the stated formal result; apparatus
certificates are recorded separately as external requirements. Finite-menu
numerical certificates remain computational evidence. None of the numerical
corroboration is a substitute for a Lean proof.

| ID | Disposition |
|---|---|
| CF01 | Common operational interface and finite objects specified |
| CF02 | Infinite spectral evolution and full chain compiled and audited |
| CF03 | Full Fourier inversion, physical Schrödinger dynamics and uniqueness, alias-free band transport; independent matrix checks |
| CF04 | Spectral finite counterexample checked with complete POVMs; formal finite embedding compiled |
| CF05 | Global 1/24 approximation, normalized tail and Born bounds implemented formally |
| CF06 | Heterogeneous-menu product TV and finite-resource test errors, including unequal/zero counts and a common convergence threshold |
| CF07 | Complete noisy readout in infinite and site spaces, phase robustness, wrapping and degenerate controls; apparatus certificates external |
| CF08 | Exact one-mode overlap and two-mode ratio obstruction proved; finite-menu interval certificates checked separately |
| CF09 | Executable count-based decision and conditional eligible budgets complete; source and physical calibration cost blocked |
| CF10 | Scoped primary-source comparison completed; unavailable version details not used |
| CF11 | No inspected compatible dataset; prospective outcome selected |
| CF12 | Literal path occupancy outside scope; no conclusion relies on it |

## Formal evidence

The [formal source, axiom audit and theorem report](https://github.com/stevenwarejones/ontology-separation/tree/f9f7ae2905827117318149cd51aea27867be267a)
cover physical site dynamics, full-space noisy measurements, phase/scale controls
and heterogeneous independent menus. The audit covers all 201 continuum roots;
the report contains 164 exported theorems. The only axioms used are `propext`,
`Classical.choice` and `Quot.sound`. The model and protocol identify the external
physical assumptions and the separately derived sharper truncation bound.
