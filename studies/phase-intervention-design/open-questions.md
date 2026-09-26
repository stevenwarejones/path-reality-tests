# Open questions — phase-intervention design

No apparatus has been calibrated and no empirical result is claimed. The
two-mode implementation in the protocol is a candidate, not a demonstrated setup.

| ID | Question | Current evidence and missing inputs | Closure condition |
|---|---|---|---|
| PH01 | Can the selected apparatus implement the intended regional phase map? | Ideal map is specified; compare actual phase-mask resolution, spatial/polarization coupling and settling against it. | Quantified implementation error or documented failed feasibility. |
| PH02 | Can preparation, loss and detector setting dependence be bounded? | Full trial outcomes and nuisance allowance B are specified; collect independent calibration and contamination bounds. | Defensible nuisance budget, including calibration uncertainty, or explicit inability to reject the null. |
| PH03 | Does the spatial representation approximate the physical propagation? | A projection error bound and leave-and-return counterexample are derived; obtain intermediate leakage and quadrature controls. | Stated domain/mode model with quantitative approximation error. |
| PH04 | Is the predicted contrast detectable at a justified sample size? | Hoeffding and CP power certificates, including occupation calibration, are in design.py; achieved signal, nuisance and independence need pilot evidence. | Frozen selection, multiplicity plan, nuisance margin, sample size and stopping rule. |
| PH05 | Can an extension exclude a stronger noncontextual class? | A relevant finite-pointer theorem is identified; audit its proof and exact instrument construction, then design equivalence/disturbance controls and finite-statistics treatment. | Applicable theorem plus justified controls and data analysis, or a clear negative feasibility result. No claim is inherited from ordinary phase contrast. |
| PH06 | What would this study add to established interference experiments? | Closed at literature scope: Hardy’s interferometer argument and Grangier–Roger–Aspect’s triggered experiment are close precedents. This is formal/statistical re-analysis of established physics; the design stops at a conditional template. | Explicit comparison identifying the contribution, even if it is reproducibility rather than novel physics. |
| PH07 | Is the occupation readout related to the phase-run ontology by the declared premise? | The protocol separates a fifth region-measurement context and its interval from the ontic identification μ(P)=w. Missing-region detection and context disturbance remain physical uncertainties. | Calibrated complete region outcomes and an explicit defended context-identification premise, or a test conditional on that unverified premise. |
| PH08 | How much setting information reaches Q responses? | Leakage ℓ contributes at most (1−w)ℓ, conservatively included as ℓ in B. No apparatus leakage measurements exist. | Quantified leakage and nuisance with coverage, or no exclusion. |
