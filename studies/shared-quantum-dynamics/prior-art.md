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
