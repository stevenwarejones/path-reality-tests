# Open questions — phase-intervention design

No apparatus has been calibrated and no empirical result is claimed. The
two-mode implementation in the protocol is a candidate, not a demonstrated setup.

| ID | Question | Current evidence and next action | Closure condition |
|---|---|---|---|
| PH01 | Can the selected apparatus implement the intended regional phase map? | Ideal map is specified; compare actual phase-mask resolution, spatial/polarization coupling and settling against it. | Quantified implementation error or documented failed feasibility. |
| PH02 | Can preparation, loss and detector setting dependence be bounded? | Full trial outcomes and nuisance allowance B are specified; collect independent calibration and contamination bounds. | Defensible nuisance budget, including calibration uncertainty, or explicit inability to reject the null. |
| PH03 | Does the spatial representation approximate the physical propagation? | A projection error bound and leave-and-return counterexample are derived; obtain intermediate leakage and quadrature controls. | Stated domain/mode model with quantitative approximation error. |
| PH04 | Is the predicted contrast detectable at a justified sample size? | Synthetic algebra and a conservative conditional power calculation are available in design.py; actual signal and nuisance bounds require independent pilot data. | Frozen selection, multiplicity plan, nuisance margin, sample size and stopping rule. |
| PH05 | Can an extension exclude a stronger noncontextual class? | A relevant finite-pointer theorem is identified; audit its proof and exact instrument construction, then design equivalence/disturbance controls and finite-statistics treatment. | Applicable theorem plus justified controls and data analysis, or a clear negative feasibility result. No claim is inherited from ordinary phase contrast. |
| PH06 | What would this study add to established interference experiments? | The phase-unitary/POVM model is established; Biswas et al. and Hacker et al. provide close theoretical and apparatus references. The present contribution is formal verification and reproducibility, with no experimental novelty claim. | Explicit comparison identifying the contribution, even if it is reproducibility rather than novel physics. |
