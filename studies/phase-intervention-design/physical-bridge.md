# Connecting the phase-intervention apparatus to its model

This is a general design analysis. No apparatus has yet been selected or calibrated. None of the following conditions is assumed to have been experimentally established.

| Model object | Required physical meaning | Evidence needed |
|---|---|---|
| State ρ | Prepared state at the selected intervention plane | Preparation stability, spatial/polarization modes and contamination bounds |
| P and Q=I−P | Orthogonal selected-region and complementary subspaces | Region edges, finite resolution, mode overlap and phase-mask response |
| Vθ=Q+exp(iθ)P | Phase change with otherwise fixed apparatus | Phase calibration and bounds on setting-dependent transmission, diffraction or polarization changes |
| Effect 0≤E≤I | Fixed downstream outcome, including efficiency | Detector response and complete herald-denominated outcomes |
| D(ρ)=PρP+QρQ | The precisely specified dephased comparison class | Model assumptions plus characterization of any randomized-phase control |
| Intermediate projections | Actual filters or justified truncation of propagation | Intermediate leakage bounds and separate quadrature error control |

The decision rule in [follow-up-protocol.md](follow-up-protocol.md) is conditional on this correspondence. Successful algebra checks do not establish these physical premises.

## Finite-window error: what can be bounded

Let the input be normalized, the propagation steps Uj unitary, and Pj orthogonal projections at intermediate planes. Define

\[
A=U_M\cdots U_1,\qquad B=U_M P_{M-1}U_{M-1}\cdots P_1U_1.
\]

The truncated construction B describes successful passage through all those projections. It is generally different from A. For ideal, unfiltered states ψj=Uj⋯U1ψ, set ℓj=‖(I−Pj)ψj‖². A telescoping expansion, with contractions after each omitted component, gives

\[
\|(A-B)\psi\|\leq\sum_{j=1}^{M-1}\sqrt{\ell_j}=\varepsilon.
\]

For any fixed detection effect 0≤E≤I, write u=Aψ and v=Bψ. Expanding ⟨u,Eu⟩−⟨v,Ev⟩ and using ‖u‖,‖v‖≤1 gives

\[
|p_A(E)-p_B(E)|\leq\min(1,2\varepsilon).
\]

These are unconditioned probabilities. Renormalizing B by its success probability changes the comparison; the same bound cannot simply be carried over. The bound may be loose, but it states exactly which missing measurements would make a truncation claim defensible. It also requires genuine projections and controlled propagation; a finite-node quadrature approximation needs its own discretization bound.

A two-mode counterexample is enough to expose the issue. Let U1=U2=H (the Hadamard matrix), ψ=|0⟩ and P1=|0⟩⟨0|. Unrestricted propagation returns |0⟩ with probability 1. With the intermediate projection, the final |0⟩ probability is 1/4, total survival is 1/2, and the conditional |0⟩ probability is 1/2. Even zero population outside the final output window does not establish that intermediate omissions were harmless: amplitudes can leave and return. This calculation is checked in `intervention_checks.py`; it is a synthetic model example, not an apparatus measurement.

## What acquisition can resolve

For each intermediate plane, measure or bound population outside the selected window using the same initial preparation and unrestricted earlier propagation. Record sufficient spatial range and detector calibration to bound the missing tails. Separately document real apertures, mode profiles and integration weights. Recompute both projected and unrestricted predictions, preserving absolute survival as well as conditional shape. A numerical window scan without tail/calibration control is a sensitivity analysis, not proof of convergence.
