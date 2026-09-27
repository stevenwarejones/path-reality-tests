# Claim-by-claim prior-art audit

| Claim / closest primary source | Precise increment and limitation |
|---|---|
| [Blume-Kohout et al. 2017](https://www.nature.com/articles/ncomms14485): interleaved GST/RB, Markovian discrepancy, alternating two-state construction | This is the closest experimental baseline. Their alternating-error explanation is explicitly reproduced as a model class, not claimed as new. We use held-out length families and retain full per-sequence predictions instead of comparing only decay rates. Our restricted depolarizing gate family is less general than their GST estimation. |
| [Nielsen et al., GST review, Quantum 5, 557](https://quantum-journal.org/papers/q-2021-10-05-557/): self-consistency, gauge, predictive characterization | Common-gate propagation and local Jacobian diagnostics are standard. Our operational comparisons use probabilities; fitted coordinate differences are not device differences. No Born-based gate estimate is offered as independent proof of Born's rule. |
| [Wood and Gambetta, 1704.03081](https://arxiv.org/abs/1704.03081): leakage and seepage characterization | Loss/return recurrence and the importance of computational-subspace measurement are established. The new artifact specializes them to an exact all-word observable-coordinate fiber, supplies its conditional physical interval, and verifies equivalent fits on both selected acquisition families. This is not a new definition of leakage or a general leakage-certification theorem. |
| [Wallman, Barnhill, Emerson, 1412.4126](https://arxiv.org/abs/1412.4126): leakage benchmarking with extra control assumptions | Ordinary binary terminal records are not relabeled as their leakage measurement protocol. The minimal monitor proposed here is conditional on calibrated contrast and a fixed computational subspace. |
| [Chen and Baldwin, 2502.00154v2](https://arxiv.org/html/2502.00154v2): leakage readout ambiguity and reanalysis of Quantinuum records | Particularly close to the rejected Quantinuum candidate: readout can distort inferred fidelity; population/postselection protocols need explicit assumptions. Our selected result concerns exact word probabilities and GST plus RB, not a claim that readout confounding was previously unknown. |
| [Białecki et al., 2308.11246](https://arxiv.org/abs/2308.11246): repeated-operation dimension witnesses | Repetition assumptions motivate the move beyond independent table states/effects. We do not transfer the Belem/Lima/Quito records to Nairobi or claim a new dimension witness. |
| [Strikis, Datta, Knee, 1811.05220](https://arxiv.org/abs/1811.05220): leakage and dimension certification | A table-level leakage alternative need not be one shared gate set. This study supplies shared CPTP constructions, but does not solve all gates for the Nairobi tables. |
| [Gullans et al., PRX Quantum 5, 010306](https://doi.org/10.1103/PRXQuantum.5.010306): characterization with time-correlated noise | Finite-memory process descriptions and the distinction from Markovian gates are established. Our dormant flag and clock pair are explicit restricted realizations with an exact missing intervention, not a new general hidden-state realization theorem. |
| [Blume-Kohout et al., 2012.12231](https://arxiv.org/abs/2012.12231): wildcard error budgets | Context-cost interpretation must distinguish a physical mechanism from unexplained model discrepancy. Our shared-word confidence certificate is conservative and returns zero; no fit residual becomes a necessary leakage budget. |
| [Hyyppä et al., PRX Quantum 5, 030353](https://doi.org/10.1103/PRXQuantum.5.030353): leakage-control experiment | Its measured leakage arrays are useful controls for methods, but its pulse calibration is not a constraint on this ion or Nairobi. |

The meaningful increment is a bounded resource nonidentification statement tied to
actual shared-operation records, with a sharp conditional fiber and explicit
source-removal witnesses. A literature sweep did not establish priority for this
specific specialization. Publication-level importance is **not established**.
The known ingredients and existing analysis are sufficient reason not to open a
Lean companion merely for elementary recurrence algebra. A formal companion would
make sense only if a materially stronger general theorem is developed.

## Escalation comparison

| Escalated claim | Closest result | Increment actually supported |
|---|---|---|
| General stationary CPTP baseline | Blume-Kohout et al. 2017 and the deposited GST analysis; Nielsen et al. GST review | Closes our implementation's restricted-channel gap under the retained held-out split. General gate-set fitting and the GST/RB pairing are not new. Local failures do not independently reproduce a global rejection. |
| All-word gate-dependent, nonunital leakage orbit | GST similarity/gauge freedom; Wood–Gambetta leakage/seepage framework; Chen–Baldwin 2025 readout dependence | Supplies an explicit physical similarity, transformed return-state polarization, and its positivity conditions. It proves equality of arbitrary word probabilities, not only averaged RB decays. Neither generic readout confounding nor the use of a GST similarity is claimed as new. Priority for this specific orbit has not been established. |
| Five doubled-word outer bounds | [Gutoski–Wu, section 2.3](https://arxiv.org/abs/1011.2787): Bures-angle triangle inequality and channel contractivity | Applies standard geometry to unknown qubit SPAM and the deposited GST/RB counts, with exact outward count-tail arithmetic. This is a full-class certificate; it is not a new general inequality priority claim or proof of strict complete-region contraction. |
| Singular observable-basis interval obstruction | GST fiducial/Gram-matrix inversion and gauge-free realization | Exhibits opposite exact determinant signs inside this acquisition's entrywise confidence box. Only that interval-inversion relaxation fails; CPTP-aware certification may be stronger. |

The closest leakage paper explicitly changes measurement/control assumptions to
recover interpretable RB quantities. Our complementary negative statement asks
what remains invariant if those calibration anchors are absent. The added
nonunital result makes return-state polarization part of that ambiguity. It still
assumes incoherent block dynamics and one leakage level; coherent leakage and
independently fixed instruments are not covered.

No publication-level novelty claim follows from these additions. A companion
formalization is not opened: the proof is a small similarity calculation plus
standard contractivity, and independent physical/math tests are more useful at
this stage than an elementary formal PR. A stronger all-source feasible-region
certificate would be a materially different result.

The final escalation additionally uses a **nonminimal physical realization**,
rather than an invertible similarity: QG_g=Ψ_gQ collapses a qubit-plus-leakage
model onto an arbitrary interior qubit model. The sharp CP budget is
λ_g≤2λ_min(J_g) for the fixed representation and maximally mixed replacement.
Its return-state condition allows nonunital channels. This is attached to the
compatible all-length measured prediction vector. It therefore fixes the earlier
weakness that our leakage examples belonged only to a failed held-out fit.

The nearest conceptual prior art remains GST gauge/realization nonuniqueness,
dimension-witness limitations, and leakage/readout calibration dependence. We do
not claim to have discovered hidden realizations or invisible dimensions. The
precise implemented increment is a positive-transfer CPTP lift with a conditional
CP budget, explicit return polarization, finite two-state-memory realization,
and a calibrated monitor that separates it on this acquisition's operations.
The exact readout condition z=Tr(Eτ) is substantive; a calibrated leakage response
can rule the lift out. Neither its priority nor a major-advance claim is established.

The original paper's non-Markovianity analysis uses likelihood-based aggregate
badness of fit. Passing our conservative simultaneous cell/mean region does not
reverse that analysis, and it is not equivalent to passing its likelihood test.
The new all-length example only certifies compatibility with our explicitly
stated region; frozen prediction still fails. This difference in inferential
standard must accompany any comparison of the two analyses.
