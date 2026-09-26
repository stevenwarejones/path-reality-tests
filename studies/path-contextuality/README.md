# A finite path-projector contextuality experiment

This is a **prospective, conditional design**, with exact arithmetic and synthetic
power calculations. It is not an executed experiment. The formal companion is
[ontology-separation PR #87](https://github.com/stevenwarejones/ontology-separation/pull/87).
The closest theorem is Kunjwal–Lostaglio–Pusey (2019), Theorem 3; the contribution
is its finite formal implementation and a reproducible costed example.

The test rejects a common-preparation stochastic model with a bounded negative
probe response and an identity-plus-disturbance transition representation. Those
premises arise from specified measurement/transformation equivalences together
with noncontextuality. Definite path occupancy is unnecessary for the bound.
Thus it also excludes the subclass assigning definite occupancy, but says
nothing against all trajectory theories. Ordinary phase fringes do not establish
these equivalences.

- [Derivation and model comparison](derivation.md)
- [Complete protocol, calibration and statistical assumptions](protocol.md)
- [Prior-work comparison and annotated bibliography](literature.md)
- [Search log](literature-search-log.md) and [structured evidence](evidence.json)
- [Question register](open-questions.md) and [requirement-to-evidence report](completion-report.md)
- [Reproducible calculation](design.py) and [machine results](results.json)

## Reference calculation

Use path basis (Q,P), K₋=diag(4/5,3/5), K₊=diag(3/5,4/5), source (4,3)/5,
and final success (4,−3)/5. The witness uses the joint negative-pointer/success
probability a, bypass success f, pointer cap q and disturbance weight d:
W=a−qf−d(1−f). Exact values are a=1369/15625, f=49/625, q=16/25,
d=1/50, W=297/15625. All failed postselections are retained.

| Final efficiency | Probe efficiency | Final bit flip | Witness gap | Total eligible trials, including calibration audit |
|---:|---:|---:|---:|---:|
| 100.00% | 100.00% | 0.00% | 0.019008 | 4,351,904 |
| 90.00% | 95.00% | 0.00% | 0.0134224 | 6,411,796 |
| 80.00% | 90.00% | 1.00% | 0.007081216 | 17,065,004 |
| 50.00% | 80.00% | 2.00% | −0.0054768 | No violation |

These are conservative fixed-count Hoeffding budgets at conditional false-rejection
risk 0.005 and power at least 90.00%, plus 44 audit contexts at simultaneous
coordinate radius 0.01 with error probability 0.005. Audit coverage is **not**
evidence that exact operational equivalences or ontic closeness hold. Calibration
control characterization is a separate prerequisite. In particular, the d-control
requires a characterized complete X readout; unknown readout losses cannot be
called d. No run duration or emitted-photon budget follows without source,
herald and calibration-device measurements.

The declared rational grid has 1,215 candidates. With final efficiency 90.00%,
probe efficiency 95.00% and final flip 1.00%, its best cost is 3,744,976 eligible
trials at (t,u,v)=(3/10,1/4,1/4). This is a finite-grid result, not a global optimum.
The source and postselection parameters as well as strength were searched;
selection precedes any confirmatory observations.

## Reproduce

Use Python 3.12 and an isolated environment:

```sh
python -m venv .venv-contextuality
.venv-contextuality/bin/pip install -r studies/path-contextuality/requirements.txt
.venv-contextuality/bin/python studies/path-contextuality/design.py --check
.venv-contextuality/bin/python -m unittest discover -s tests -p test_path_contextuality.py -v
```

The exact Fraction evaluator is checked against independent complex matrix
multiplication, including imaginary off-diagonal inputs. Seed 260926 is fixed.
The Monte Carlo result is accompanied by its own binomial interval and does
not replace the analytic power certificate. Snapshot comparison keeps counts,
keys, decisions and strings exact, with explicit tolerance only for floats.

The extension now includes `dual_decision` and `dual_certified_budget` for both
success-cell witnesses; see the joint-decision section of [protocol.md](protocol.md).
The original grid and random-allocation figures remain single-witness baselines.

The loss comparison in [the protocol](../path-contextuality/protocol.md) gives
joint-rule budgets of 4.43M, 6.22M, 10.23M and 56.41M eligible trials. The last
three improve on the negative witness under those specified noise models.
The q calibration must attain the actual negative effect's largest eigenvalue;
tomography must confirm the eigenbasis and ordering, including uncertainty.
