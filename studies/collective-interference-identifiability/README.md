# What does calibrated three-atom bunching exclude?

**The public atomic records conditionally exclude mixtures of distinguishable
groups containing at most two atoms. The exclusion is too fragile to claim
that the released calibration independently establishes the required apparatus
assumptions.** This is a verified conditional result, not completion of the
broader apparatus-robust research goal.

The combination is **two measurement families in one apparatus**: a singleton
propagation record and three-atom parity images from Young et al.,
[An atomic boson sampler](https://arxiv.org/html/2307.06936v2),
[Zenodo 10453016](https://zenodo.org/records/10453016).
The three singleton arrays are shifted copies of the same 931 shots. They do
not independently establish translation invariance or provide 2,793 trials.

| Quantity | Verified result |
|---|---:|
| Three-atom prepared shots / bunching events | 2,999 / 93 |
| Three-atom shots with three detected occupied sites | 2,329 |
| Unique singleton calibration shots | 931 |
| Simultaneous lower bound on bunching per prepared shot | 0.02364675 |
| Certified upper bound on distinguishable row coincidence | 0.01166 |
| Ceiling for the size-two cluster class | 0.02332 |
| Separation in absolute probability | 0.00032675 |
| Sum of allowable singleton-transfer TV errors, with other errors zero | **less than 0.000163375** |
| Lower bound on the defined non-cluster mixture weight, under exact sharing | 0.007005789 |

The statistical construction allocates a total error budget of 0.05 across
four inspected particle numbers (2–5), using exact binomial tails and a union
bound. Coverage requires IID shots within each genuine record and the stated
measurement model. It does not require independent copies of calibration data.
The result is exploratory archive reanalysis, not preregistered confirmation.

The exact certificate covers **every probability vector in the declared
confidence region**, with unobserved crop-tail mass retained. It is not a
finite model bank, an optimizer's failure to find a counterexample, or a
bootstrap interval. Four boxes and rational dual bounds suffice; the checker
uses only Python's standard library and exact integer/fraction arithmetic.

## What the combination adds, and what remains open

Singleton images alone cannot constrain the internal distinguishability
resource: an explicit fully distinguishable physical model reproduces the
entire empirical singleton image distribution. Its bunching probability is
0.00462813088179467. Using identical internal states with the same propagation
gives a physical example with bunching probability in
[0.0252913383446, 0.0252913383815], inside the declared joint summary region.
Thus the joint summary region is nonempty, while the size-two-cluster part is
empty under exact sharing.

The bunching statistic alone also admits fully distinguishable explanations
with different propagation. **We have not constructed a model satisfying every
constraint of the complete higher-particle image distribution after dropping
calibration.** Consequently this is a certified joint summary exclusion plus
a full-singleton-source witness, not a proof that both complete sources are
indispensable. Nor is the compatible quantum example a fit to every
higher-order output count. These distinctions prevent a coarse ablation from
being mistaken for complete-source separation.

The numerical separation allows only 327 extra bunching events per million
prepared shots from all unmodelled effects combined. Independent calibration
of spatial transfer, drift, interactions, false positives and preparation
errors to that level has **not** been established here. Exact translation is a
particularly strong assumption: the archive's copied calibration arrays do
not test it. No broad empirical certification, violation of quantum mechanics,
computational advantage, or consequential novelty is claimed.

## Read and reproduce

* [Feasibility gate](gate.md): archive choice and exact measurement families.
* [Theory and assumptions](theory.md): physical class, parity channel, confidence proof and robustness budget.
* [Source audit and prior art](source-audit.md): dependency/cropping findings and comparison with existing witnesses.
* [Claim status](claim-status.md): measured, proved, conditional and open claims.
* [Source manifest](sources.json): pinned URLs, hashes and licensing notes.

Fast offline checks (Python 3.12; no third-party packages):

```sh
python studies/collective-interference-identifiability/certificate.py --check
python studies/collective-interference-identifiability/witness.py --check
python -m unittest discover -s tests -p 'test_collective_interference.py' -v
```

Full source reconstruction (59 MB download, approximately 630 MB extracted;
raw data remain outside git):

```sh
python -m pip install -r studies/collective-interference-identifiability/requirements.txt
python studies/collective-interference-identifiability/sources.py --cache /tmp/collective-raw --download --check
```

To regenerate certificates after an intentional analysis change:

```sh
python studies/collective-interference-identifiability/build_certificate.py
python studies/collective-interference-identifiability/build_witness.py --cache /tmp/collective-raw
python studies/collective-interference-identifiability/witness.py
python studies/collective-interference-identifiability/diagnostics.py
```

Solver versions can change the witness/dual parameters; trust the independent
exact check, not byte equality of newly optimized parameters. Source counts
and checked committed outputs must reproduce exactly. The dedicated workflow
runs fast checks on PRs; its manually requested full-source job downloads and
reconstructs the original records.
