# Collective interference: fixed-channel exclusion and a drift ambiguity

The [shared-dynamics and control-design extension](control-dynamics-report.md) proves that the proposed pair statistic has **zero worst-case separation**: a physical cluster mixture matches it exactly while passing the expanded measured region. A matched labelled-triple control gives a sufficient class-wide test with conditional source-removal witnesses. A joint July 7 Hamiltonian diagnostic fails calibration constraints; apparatus-supported robustness remains open.


**The selected records exclude the declared fixed-channel cluster null, but
admit an explicit shared fluctuating cluster model in the expanded confidence
region.** Both preparation families use the same two-channel law. This is a
physical countermodel for the broad passive-channel class; its large
fluctuation is not established as plausible for the apparatus's lattice
Hamiltonian. The apparatus-robust research goal remains **incomplete**.

| Current robustness result | Exact-check conclusion |
|---|---:|
| Shared-drift budget excluded through (sum of per-input TV deviations) | 0.0397 |
| Physical shared-drift model feasible at | 0.89852380063206 |
| Mean singleton, all parity-pattern and 12 row-hit constraints | All pass |
| Old focusing model rejected by added row-hit constraints | 8 rows |
| Independently calibrated settings with physical joint cluster models | 16 / 16 |
| Proposed matched distinguishable-pair event: fixed quantum / drifting cluster | 0.143236 / 0.200070 |

The critical drift budget is bracketed, not solved. The feasible endpoint is
about 30% TV per input; it does not show that tiny errors explain the data.
A separate event-specific envelope excludes fixed preparation changes through
sum TV .00554, about 17.7 times the old generic sufficient allowance.

Read the [robustness decision table and derivation](robustness-report.md) for
the nuisance coordinates, shared acquisition law, complete setting checks,
confidence allocations, limitations and proposed separating control. The
selected-data and 16-setting analyses have separate declared confidence
families; they are not combined into an archive-wide 95% claim.

## Preserved fixed-channel source-combination result

The two measurement families come from Young et al.,
[An atomic boson sampler](https://arxiv.org/html/2307.06936v2),
[Zenodo 10453016](https://zenodo.org/records/10453016): 931 unique singleton
shots and 2,999 three-atom shots, with 93 collisionless full-row events.
Three supplied singleton arrays are translated copies of that one record.

## Fixed-channel results

The final declared simultaneous region includes singleton row and conditional
cell probabilities, loss/crop tails, all 497,785 possible parity patterns with
at most three occupied sites in the selected crop, and the bunching event.
Of those patterns, 2,539 were observed. Exact marginal binomial intervals and
union bounds give nominal 95% familywise coverage across the declared n=2–5
family, conditional on IID sampling and the stated physical bridge. This is
exploratory archive reanalysis, not preregistered confirmation or coverage over
an arbitrary search of other settings and archives.

| Final expanded-region quantity | Exact-check result |
|---|---:|
| Bunching lower endpoint | 0.023025033926510948 |
| Distinguishable exact-event upper bound | 0.0112 |
| Size-two cluster ceiling | 0.0224 |
| Separation margin | **0.000625033926510948** |
| Sufficient sum of full-cell transfer-TV allowances, other effects zero | **less than 0.000312516963255474** |
| Conditional non-cluster mixture-weight lower bound | 0.0139516501453 |
| Joint quantum model bunching probability | 0.0259974396797 |
| Higher-source-only cluster model bunching probability | 0.0289752778116 |

| Retained sources | Conclusion for the same declared region and physical assumptions |
|---|---|
| Singleton + complete selected three-atom parity record | C2 excluded; a physical identical-boson model satisfies every retained constraint. |
| Singleton only | The preserved fully distinguishable physical witness reproduces all empirical singleton cell frequencies and belongs to the full singleton region. |
| Complete selected three-atom record only | A physical fixed-channel pair-plus-singleton model satisfies every retained pattern interval and bunching constraint. Exact translation of its unconstrained channel is still enforced. |
| Additional controls | No extra dataset is asserted indispensable or used to shrink the selected channel region. |

Thus both selected sources are necessary for this preserved **region-based**
exclusion. This deletion claim concerns the original full-pattern region,
before adding row-hit constraints or enlarging the channel class.
The very sparse full-pattern confidence region is broad: compatibility is not
a likelihood fit, proof that every other statistic is inconclusive, or a fit
to the entire multi-setting archive or a specific lattice Hamiltonian.

## What improved and what did not

| Analysis design | Distinguishable ceiling | Cluster margin |
|---|---:|---:|
| Preserved original row-only certificate | 0.01166 | 0.00032675 |
| Row-only relaxation using the revised singleton region | 0.0119 | −0.00015325 |
| 2D exact-event certificate with scalar bunching region | 0.0112 | 0.00124675 |
| Final 2D certificate with complete parity region | 0.0112 | 0.0006250339265 |

These are alternative confidence designs for comparison, not four additional
simultaneous 95% discoveries. Additional cell categories widen the row
intervals enough that the revised row-only certificate does not exclude C2.
Retaining distinct x positions restores the exclusion. Including all higher
patterns then spends part of the bunching budget. Neither upper bound is
claimed globally optimal. [revision.md](revision.md) gives the proof,
allocation, runtime, source-sharing table and discovery history.

Independent input-site records exist at 2.428571 ms on July 7 and at 4.65 ms
on July 2, 2022. The selected 2.45 ms record is from July 2. These controls do
not bound the selected channel without an additional justified dynamical and
date-transfer model. No measured bound below the sufficient TV allowance has
been established; no arbitrary interaction or detection error is inserted.
A shared passive-channel drift model now fits the expanded selected region,
and a proposed matched pair control separates it from the exhibited quantum
model. Its restriction to an apparatus-specific Hamiltonian remains open. The physical null is not every definition of absent
three-particle interference, entanglement depth or a universal resource.

## Reproduce

Fast offline exact checks use Python's standard library:

```sh
python studies/collective-interference-identifiability/certificate.py --check
python studies/collective-interference-identifiability/witness.py --check
python studies/collective-interference-identifiability/cell_certificate.py --check
python studies/collective-interference-identifiability/full_models.py --check
python studies/collective-interference-identifiability/robustness.py --check
python -m unittest discover -s tests -p 'test_collective_interference.py' -v
```

Reconstruct the pinned originals and the complete derived histogram:

```sh
python -m pip install -r studies/collective-interference-identifiability/requirements.txt
python studies/collective-interference-identifiability/sources.py --cache /tmp/collective-raw --download --check
python studies/collective-interference-identifiability/build_full_models.py --cache /tmp/collective-raw --check-source
python studies/collective-interference-identifiability/control_audit.py --cache /tmp/collective-raw --check
python studies/collective-interference-identifiability/build_robustness.py --cache /tmp/collective-raw --check-source
```

The 59 MB archive and extracted NC files remain outside git. The deposited
parity histogram contains derived aggregate counts, not shot chronology.
PR CI now runs both exact checks and full source reconstruction. Solver-based
regeneration uses `build_cell_certificate.py` and `build_full_models.py --cache
/tmp/collective-raw`; optimizer versions can change proposed witnesses. The
independent exact check decides validity, not byte equality of new proposals.

Further reading: [physical class and preserved original proof](theory.md),
[revision and full-region derivation](revision.md),
[source audit and primary prior art](source-audit.md),
[claim status](claim-status.md), [source manifest](sources.json).
