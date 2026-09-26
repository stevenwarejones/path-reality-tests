# What the public-data audit establishes

Reviewed 2026-09-26 UTC. The deposited Figure 3/4 summaries are reproducible under
the documented mappings. The independent endpoint series E is closer to coherent
reconstruction Q than to the deposited incoherent comparator C in the checked
descriptive metrics. This supports that limited comparison. It neither certifies
individual photon trajectories nor supplies a significance test from missing
repeat-level records.

The public-data analysis and ideal formal comparison have been carried out.
End-to-end reconstruction from photon acquisitions, a retrospective uncertainty
audit, and an experimental implementation of the follow-up remain incomplete.
The distinctions below make those limits explicit.

## Requirement-to-evidence checklist

| Requirement | Evidence | Disposition |
|---|---|---|
| Identify paper, supplement, data version, license and original bytes | [Provenance](provenance.json), [manifest](manifest.json), [download instructions](raw/README.md), all 19 digest checks | Complete for the pinned version; new releases need a separate audit |
| Inspect every sheet and distinguish signals, means and uncertainties | [Dictionary](data-dictionary.md), [column map](column-map.json), replication inventory | All 18 workbooks / 19 sheets mapped; processed values are not called raw trials |
| Reproduce central probability and phase/action plots | [Figures, residuals and summary](results/replication/summary.json), [validation](results/validation.md) | Dedicated Figure 3/4 deposit replotted in original units; original figure-generation pipeline unavailable |
| Quantify and investigate numerical differences | [Method comparison](method-comparison.md), interpolation/normalization/rounding checks | Conditional results documented; the actual implementation of the reported 4.45% remains unidentified |
| Reconstruct the inferential chain and inspect the physical bridge | [Study explanation](README.md#from-measurements-to-path-claims), [physical correspondence](physical-bridge.md) | Each acquisition/reconstruction step separated; calibration and intermediate truncation unverified |
| Compare coherent, incoherent and broader contextual models | Model table below; pinned formal companion | Exact ideal exclusion and an explicit matching response; no statistical exclusion of an empirical trajectory class |
| Assess retrocausal and superluminal interpretations | Ordered paraxial propagation; no causal intervention in the deposit | This is not an identifiable causal-response test; no arbitrary retrocausal rival introduced |
| Propagate justified statistical uncertainty | Inventory of supplied SE/SD; [missing-record specification](required-records.md) | Blocked by joint repeat/calibration records; no manufactured counts, covariance or p-values |
| Explain why summaries cannot replace the missing inputs | [Exact counterexamples](identifiability_checks.py), [checked output](results/identifiability.json) | General non-identifiability established by explicit different joint records with the same marginal summaries |
| Supply an exact normalized formal example and correspondence | [Merged Lean case study](https://github.com/stevenwarejones/ontology-separation/blob/ab1eeb81119ff3fdcb46a771d9bcdf5e39ea3d29/docs/PATH_INTERFERENCE_CASE_STUDY.md); general phase extension linked from the follow-up validation | Formal models checked separately from apparatus claims |
| Select the strongest immediate follow-up | [Independent phase-intervention study](../phase-intervention-design/README.md) | Full proposed protocol, complete outcomes, nuisance rule and conditional power; no data acquired |
| Review supportive and adverse literature at multiple stages | [Search log](literature-search-log.md), [comparison](literature-comparison.md), follow-up review | Bounded primary-source review; no exhaustive priority or absence claim |
| Account for every open question | [Question register](open-questions.md), grouped disposition below | Every current entry retained; external dependencies are not marked scientifically resolved |

## Actual model comparison

| Defined model | Prediction / construction | What can be concluded |
|---|---|---|
| Coherent finite amplitude composition | Multiply the specified K entries along each history, add amplitudes, then apply the stated detection rule | Deposited Q compared with independently acquired E; raw-factor reconstruction and calibration still missing |
| Incoherent sum for the same declared histories | Delete cross terms and sum squared magnitudes, with the documented response and normalization | Deposited C is descriptively farther from E in the checked metrics. Exact ideal class exclusion lives in Lean; empirical statistical exclusion is not established |
| Transfer-matrix representation of the same finite model | The matrix-product entry equals the corresponding finite intermediate-index sum | Same complete prediction family under the same boundaries and detector rule; this representation cannot be distinguished merely by redescribing the sum |
| Unrestricted stochastic setting-response model | A one-state classical model draws outcome o with probability p(o|s) for each setting s | Exactly matches any normalized finite table. It is an explicit response construction, not a dynamical trajectory theory or a proof about Bohmian mechanics |
| A retrocausal alternative | No specified rival and discriminating intervention supplied by these data | No model-class exclusion or causal direction can be inferred by inventing a comparator after the fact |

An arbitrary-unit intensity curve is not itself a normalized behavior. Applying
the response construction to a unit-sum curve reproduces that conditional shape
only; it does not identify absolute success probability. The follow-up instead
uses selected-bin, other-detection and failure outcomes per eligible herald.

## Why more processing cannot supply the missing records

`identifiability_checks.py` uses exact rational arithmetic and no experimental
inputs. Its examples are counterexamples to general reconstruction claims:

1. Two shared propagator factors each have mean 1 and variance 1/100 in both
   records. Pairing their fluctuations in the same or opposite direction gives
   product means 101/100 or 99/100. Marginal error bars do not determine a product's
   expectation or its covariance with other reconstructed paths.
2. Amplitude magnitudes {1,2,3} and phases {0,π/2,π} are unchanged as marginal
   lists. Different magnitude–phase pairing gives coherent intensities 8 and 10.
   Separate group summaries do not recover the joint path-amplitude vector.
3. The same conditional detected shape (1/4,3/4) occurs at efficiencies 1/5 and
   4/5, with different no-detection probabilities. Shape normalization loses the
   information needed for an absolute-loss comparison.

These examples establish non-injectivity of those summary maps. They do not show
that either constructed record was used in the paper, nor quantify the missing
records' actual effect. The empirical magnitude of that effect remains unknown.

## Disposition of the question register

| Questions | Status and consequence |
|---|---|
| Q01–Q03, Q05 | Acquisition/inventory and method-level reconstruction are settled at the recorded scope. Shared propagator factors are not independent path trials. |
| Q04, Q08, Q14 | Repeat-level reconstruction, absolute calibration and the original noise pipeline require additional records. The exact counterexamples explain why processed means cannot replace them. |
| Q06–Q07, Q09 | The 4.45% implementation, two reported fidelity values and printed normalization convention need the actual evaluated arrays and code. Checked alternative conventions have not identified the implementation. Findings remain conditional. |
| Q10–Q11, Q17 | Bin centers/weights, phase/action grouping and finite-aperture correspondence need instrument and grouping metadata. A plausible half-bin convention and a projection counterexample narrow the issue without resolving the apparatus. |
| Q12–Q13, Q16 | Operational equivalence, lack of a causal intervention, and the scope of the formal certificate have precise answers. They limit the ontological conclusions independently of the numerical discrepancies. |
| Q15 | Focused version/code/correction searches found no usable additional release. Coverage limits remain explicit; a search cannot prove absence. |
| Q18 | Repository verification is executable in CI; its submitted-commit status is reported by Actions. Software correctness does not establish experimental calibration. |
| Q19–Q21 | Histogram-tail and figure-column/bin-label ambiguities are recorded with exact affected cells/series. They do not change the dedicated Figure 3/4 replot. Generation code is needed to identify intent. |

The concrete remaining inputs are collected in [required-records.md](required-records.md).
New public data or code can be checked against those requirements.
