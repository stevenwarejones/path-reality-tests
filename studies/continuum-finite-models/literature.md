# Literature and evidence assessment

Review date: 2026-09-26 UTC. This is a scoped primary-source review, not a claim
of exhaustive novelty clearance. The contribution sought is formal verification,
explicit access/calibration assumptions and numerical validation of known physics.

## Search record

Questions: time slicing versus spatial approximation; spectral versus local
finite dynamics; observable error; dispersion-sensitive apparatus; usable data.
Searches included “Trotter 1959 product semi groups operators strong convergence”,
“finite difference Schrödinger cosine dispersion”, “exact discretization
Schrödinger equation”, “Brun Mlodinow 1802.03911”, “ring guided atom interferometer
data”, and “atom interferometer dispersion calibration”. Local mathlib source
inspection covered Analysis/Fourier/AddCircle, Analysis/Fourier/ZMod,
Analysis/Complex/Trigonometric and Trigonometric/Bounds. The last includes
`abs_cos_sub_cos_le`; `Real.cos_bound` supplies a conservative 5/96 fourth-order
bound on |x|≤1, not the desired global sharp 1/24 remainder.

## Sources read and implications

| Primary source / version | Derivation inspected | Finding and implication |
|---|---|---|
| [Guth, MIT 8.323 notes, 2008](https://web.mit.edu/8.323/spring08/notes/ft1ln05-08-2up.pdf) | Eqs. 5.3–5.10 and imaginary-time regulator discussion | Finite time slices retain integrations over continuous positions; real-time kernels are oscillatory. This does not establish arbitrary spatial convergence. |
| [Tarasov, PLA 380 (2016), 68–75, author copy](https://theory.sinp.msu.ru/~tarasov/PDF/PLA2016.pdf) | Eq. 2, Fourier construction, discussion of long-range differences | An exact discrete representation can use long-range couplings. Its infinite lattice differs from the finite accessible-band counterexample, but blocks conflation of discreteness with nearest-neighbor dispersion. |
| [Brun–Mlodinow, arXiv:1802.03911v1](https://arxiv.org/html/1802.03911v1), [latest author PDF](https://arxiv.org/pdf/1802.03911), [PRD 99, 015012 (2019)](https://doi.org/10.1103/PhysRevD.99.015012) | Eqs. 17–23: spin/orientation-dependent correction and interferometer phase | Their anisotropic Dirac quantum walk is a different physical model. Neither their proposed sensitivity nor Lorentz-violation limits are bounds on this ring's isotropic a²k⁴ correction. The published title says “Detecting … via”; the arXiv title says “Detection … by”. No complete published-version equation comparison was available. |
| [Watrous, 2018 author pre-publication text](https://cs.uwaterloo.ca/~watrous/TQI/TQI.pdf) | Theorem 3.4 and Proposition 3.5 | Half trace norm is the state-discrimination convention; a complete measurement cannot increase distinguishability. The optimum is not a calibrated detector implementation. |
| [Gautier et al., Science Advances 8 (2022), eabn8009](https://pmc.ncbi.nlm.nih.gov/articles/PMC9187224/) | Eqs. 2–5, Table 1, preparation/detection methods, data availability | A closest inspected precision matter-wave experiment uses Raman pulses, free-fall trajectories and a rotation-dependent phase. It does not measure free-ring mode dispersion. Its public data availability points to Zenodo, but its quoted Sagnac precision cannot be transferred to this model. |

Trotter's [1959 original](https://doi.org/10.1090/S0002-9939-1959-0108732-6)
was sought; the publisher PDF could not be retrieved. Its exact operator-domain
hypotheses have therefore **not** been verified here. No theorem depends on it:
the selected free model uses exact spectral evolution. General product-formula
convergence, Euclidean Feynman–Kac measures and real-time path measures are outside
scope. The Fourier L² identification remains an external analytic bridge.

Following the Brun–Mlodinow companion-paper lead located
[arXiv:1802.03910](https://arxiv.org/abs/1802.03910), on quantum-walk dynamics.
It is a lead, not a fully read derivation here. Searches for related later work
located guided/multi-loop Sagnac proposals; they did not establish a dataset
with this study's common ring preparation/control/readout calibration. The forward-citation search additionally located Chou (below), which directly
addresses the finite spectral Hamiltonian and strengthens the closest-prior-art
comparison. No exhaustive priority claim is made.

## Dataset disposition

- Wen repository originals are processed optical propagator summaries. Their
  sampling pixels are detector/instrument coordinates, not fundamental sites.
  The massive free-ring Hamiltonian is inapplicable; no photon mass is inserted.
- NIST records concern Bell outcomes and setting assignments. They do not identify
  the free-ring kinetic spectrum or its phase/time/momentum nuisance parameters.
- Gautier et al. explicitly provide [data/code at Zenodo 6372385](https://doi.org/10.5281/zenodo.6372385).
  The landing page was not retrievable in this environment. The paper alone
  already identifies incompatible dynamics (gravity, Raman pulses and rotations),
  so no exclusion is extracted. A new forward model would be required before
  examining its residuals as a lattice test.

No inspected dataset supports this proposed exclusion. This is a scoped finding,
not a claim that no suitable dataset exists anywhere. No new source files, PDFs
or event extracts have been imported. The selected deliverable is prospective.
The missing records are: eligible-trial ledger including failures and heralds;
mode preparation/tail and contamination bounds; recombiner phase characterization;
time, circumference, mass and rotation calibration; trap/interactions error
budget; calibration transport across all strata; source preparation efficiency
and its confidence statement. Their absence blocks any apparatus exclusion and
any unconditional source-attempt budget.

## Final model-specific comparison

The later [Chou, arXiv:2008.03698v1 (2020)](https://arxiv.org/pdf/2008.03698)
was read through its finite-grid definitions and Fourier Hamiltonian derivation
(Eqs. 7, 11–16). It explicitly starts with quadratic momentum energies and a
normalized finite DFT, contrasting this with the central-difference dynamics.
This is closer prior art for the spectral counterexample than an arbitrary
finite-dimensional witness. The present study adds common detector/access
conditions and conservative calibration/sampling analysis; the underlying
dispersion distinction is established. Different Lorentz-violation laws and
unrestricted dimension witnesses cannot supply a bound on the selected ring.

Final searches on 2026-09-26 UTC included the exact Chou/Tarasov titles,
“1802.03911 dispersion interferometry”, and “Brun Mlodinow 015012 correction”.
No inspected correction changes the chosen nonrelativistic ring calculation.
Published-version details unavailable above remain unavailable; no transferred
estimate depends on them. The outstanding formal mathematics and apparatus
certificates are separately recorded, not treated as literature conclusions.


A final refresh also found Pachón and Gómez,
[arXiv:2604.20776v2](https://arxiv.org/html/2604.20776v2)
(28 June 2026 version). The primary HTML introduction and construction overview
were inspected on 2026-09-26 UTC. They formulate finite-dimensional Hamiltonian
path sums in discrete phase space for odd prime dimensions, with qutrit examples.
This reinforces the distinction between a finite path representation and a
nearest-neighbor spatial-lattice hypothesis. It does not provide a free-ring
apparatus dataset or a calibration certificate for this protocol. No numerical
estimate is transferred from that paper; Chou remains the closer comparison for
the specific quadratic-spectrum versus finite-difference Hamiltonians here.

The completion check revisited Chou v1, Eqs. 7 and 11–19, Tarasov's published
2016 Eqs. 1–4, Watrous's author text, Theorem 3.4 and Proposition 3.5, and
Brun–Mlodinow v1's Hamiltonian/dispersion construction. The finite Fourier
propagator here follows the selected central-difference Hamiltonian; it is not
Tarasov's infinite long-range exact discretization. Measurement probabilities
use complete outcome tables and TV is half L1. Product bounds concern independent
fixed acquisition, not independent bins within one multinomial observation.
The numerical implementation stores energies and evolves with exp(-i*t*E/hbar);
the formal spectral functions store frequencies E/hbar. The new derivative
check uses nonunit hbar to test this convention explicitly. Neither the shared
energy-ratio obstruction nor the different quantum-walk dispersion supplies a
finite-menu statistical certificate or apparatus sensitivity for this ring.
