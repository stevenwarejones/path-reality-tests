# Spacetime causal influence: a conditional optical protocol

Can a randomized local intervention change a complete receiver record outside
its future light cone, or a record fixed before the choice? This study supplies
an operational null, countermodels, a prospective protocol and reproducible
sensitivity calculations. It contains **no apparatus observations**. Standard
local quantum operations predict zero spacelike influence; the tunable alternatives
below are probability tables for sensitivity analysis, not proposed field theories.

The recommended implementation compares fresh amplitude-modulator settings with
an independently latched receiver bit across a surveyed free-space separation.
The [protocol](protocol.md) compares this with a later-choice/earlier-record test.
The optical version has an informative positive control after the causal arrival
time; the earlier-record version remains a useful independence/storage control
but offers no stronger interpretation and no analogous backward-time positive control.

## Exact question and rejected conjunction

For setting x and complete local outcome b,

\[
\Delta=\frac12\sum_b|P(b\mid do(x=1))-P(b\mid do(x=0))|.
\]

The binary record is one eligible click versus every other detection category;
no-click and invalid-detection trials remain in the denominator. Missing logs
require a declared recording-failure treatment, not silent deletion. For this
record Δ=|p₁−p₀|. Coarsening can hide changes among other raw categories.

The finite causal class has setting-independent preparation μ(λ), normalized
sender operations Kₓ(a|λ) and a receiver R(b|λ) with no setting input. Summing
all a in μKₓR gives the same receiver marginal for both settings. A positive
calibrated result rejects the conjunction of that causal response structure,
intervention independence, complete sampling and the stated calibration bounds.
It does not uniquely identify which premise failed. Correlated shared causes,
local quantum operations and operationally nonsignaling hidden-path or retrocausal
descriptions remain compatible with a null result. A finite null result bounds
influence for the stated binary record and preparation, rather than proving zero.

The [physical bridge](physical-bridge.md) separates the finite proof, conditional
quantum-operation identity and real-device assumptions. It derives a conservative
probability allowance B from calibrated mismatch couplings. Timing errors in
seconds and optical tail amplitudes cannot simply be added to a probability gap.
The complete [question register](open-questions.md) gives dispositions ST01–ST12.

## Sensitivity and feasibility

The primary [statistical derivation](statistics.md) uses a predeclared fixed
horizon, fresh fair settings and a bounded randomized score. Device memory is
allowed. Its estimand is the average signed conditional effect; interpreting its
absolute value as one binary TV requires a constant effect. A separate exact
binomial comparison declares the stronger IID assumptions, including conditional
coverage for random per-setting sample sizes.

The following **hypothetical observed gaps** include losses. At α=0.01 and
miss probability β=0.10, the conservative inequality guarantees at least 90.00%
power under a signed effect of at least d on every history. It assumes conditional
setting bias ζ=0. Nonzero ζ needs the additional correction in the statistics guide.

| Observed gap d | Nuisance allowance B | Sufficient total trials | Time at hypothetical 10 kHz |
|---|---|---:|---:|
| 0.10 | 0.002 | 3,038 | 0.304 s |
| 0.03 | 0.002 | 37,211 | 3.721 s |
| 0.01 | 0.002 | 455,830 | 45.583 s |
| 0.003 | 0.002 | 29,173,105 | 2,917.311 s |
| 0.01 | 0.01 | No uniform separation | — |

These times exclude calibration, controls, resets and logging overhead. More
shots cannot uniformly separate an alternative already admitted by the nuisance
null. If d instead describes a clean effect, worst-case nuisance cancellation
must also be subtracted before applying this table.

A synthetic null example with 100,000 trials and 50,000 setting/outcome matches
gives a 99.00% absolute observed-effect upper bound of 1.03%. With B=0.002 the
clean-effect upper bound is 1.23%, before adding calibration failure risk. This
is a calculation, not a measured constraint.

The illustrative geometry uses 30 m separation, 0.5 m support radii, 0.1 m survey
allowance, complete setting support [−5,20] ns and receiver support [35,65] ns.
The strict vacuum-light-cone margin is 7.91452794 m (26.40002351 ns). These are
requirements to demonstrate, not claimed instrument specifications. Achieved
feasibility is externally blocked by missing surveys, timing-tail data, randomizer
predictability bounds, leakage controls and record-integrity measurements.
The prospective resource study is complete conditional on those inputs.

## Evidence and contribution

| Kind | Result | Evidence / limit |
|---|---|---|
| Established theory | Local normalized operations preserve a remote marginal under the stated causal premises | [Primary evidence and closest-work comparison](evidence-table.md); not a new no-communication theorem |
| Formal verification | Finite marginal identity, normalized g-family, exclusion, toy access equivalence, countermodels, quantum algebra, interval and geometry implications | [Lean guide](https://github.com/stevenwarejones/ontology-separation/blob/study/spacetime-influence/docs/SPACETIME_INFLUENCE.md); no continuum or concentration theorem is formalized |
| Mathematical analysis | Fixed-horizon memory interval, IID difference interval, nuisance transport and power inequality | [statistics.md](statistics.md), [physical-bridge.md](physical-bridge.md) |
| Numerical evidence | Seeded synthetic repetitions, numerical binomial powers and strict sample counts | [config.json](config.json), [results.json](results.json), [design.py](design.py) |
| Proposed measurement | Complete locally latched binary outcomes, fixed windows, in-cone positive controls and negative controls | [protocol.md](protocol.md); unexecuted |
| Open physical bridge | Actual support, predictability, leakage and calibration transport | ST04, ST06, ST07 in the [register](open-questions.md); no empirical exclusion |

The closest reviewed statistical proposal is Albanese's 2026 binary-channel
criterion. Detector causality, optical precursors and loophole-controlled Bell
experiments supply distinct precedents; none is treated as this study's data.
The contribution is formal verification and a reproducible conditional protocol
analysis. No priority, improved measured bound or new-physics claim is made.
The [search log](literature-search-log.md) records coverage through 2026-09-26 UTC
and the remaining source-access limits. The [completion report](completion-report.md)
maps requirements to exact evidence.

## Reproduction

From the repository root, in an isolated Python 3.12 environment:

```sh
python -m pip install -r studies/spacetime-causal-influence/requirements.txt
python studies/spacetime-causal-influence/design.py --check
python -m unittest discover -s tests -p 'test_spacetime_influence.py' -v
```

`design.py` without `--check` regenerates the structured snapshot. PCG64 seed
20260926 and all synthetic inputs are recorded. Float comparisons allow only
1e-10 relative / 1e-12 absolute rounding differences; discrete counts, decisions,
keys and labels must match exactly, with mismatch paths reported. Tests include
independent binomial-tail coverage calculations, exhaustive memory-null trials,
strict boundaries, selection/common-cause examples and geometry failures.
No external data, source PDFs or apparatus pricing are included.
