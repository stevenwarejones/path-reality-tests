# Synthetic timing-shift sensitivity

Can a sparse local record reveal an artificial remote-setting-dependent arrival
shift even when total clicks do not change? Can a frozen learned score recover an
effect whose sign reverses between observable regimes?

The original Gaussian study answers those questions **only for simulated records**.
The [generated report](report.md) contains detection frequencies, null controls,
and Monte Carlo uncertainty. Its `design.py` reads no experimental event files.
All experimental, causal and spacelike claim flags remain false.

The follow-up [empirical-background study](empirical-protocol.md) now injects
shifts into copied Bob timing records under fresh artificial assignment labels.
Its [separate report](empirical-report.md) measures recovery on the actual local
timing background, without testing the actual remote labels. The original
Gaussian simulations below remain separate and unchanged.

## Connection to the available data

The [NIST audit](../nist-bell-causal-audit/README.md) reconstructs 107,109,596
paired trials; its [report](../nist-bell-causal-audit/results/report.md) lists
primary arm click rates of roughly 130–429 per million. This simulation uses
200 clicks per million, 10 million trials for a smaller exposure, and 53,554,798
trials for an illustrative fixed-receiver-setting exposure of half that archive.
The half-size is a planning scale, not a measured receiver-setting count.

The [raw reconstruction description](../nist-bell-causal-audit/reconstruction.md)
gives a timetagger bin of 78.125 ps. We use that unit to label injected shifts.
The Gaussian width of two bins (156.25 ps), the Gaussian distribution itself,
and the persistent two-regime structure are **model choices, not fits or
calibrations to the raw timestamps**. Sensitivity in ps is conditional on them.
The sparse primary-window click rate is used as an illustrative record rate;
this is not a model of all raw detector events or of NIST's window acceptance.

## Generative model

There is one receiver with a fixed local setting and a synthetic remote bit X.
Every trial has an independent fair X, followed by a Bernoulli click and, if it
clicks, an independent Gaussian arrival time. A no-click is a recorded outcome.
There is at most one click per trial. Timing histograms have one-tag-bin spacing,
with both overflow tails kept. Shifts never remove or add events through a
timing-window cut. The model excludes detector dead time, multiclicks, timestamp
excursions, setting feedback, missing records and cross-trial device memory.

Thirty-two chronological blocks have an observable pretrial regime `block % 2`.
The first eight train the learned detectors; the remaining twenty-four test them.
Both regimes occur equally often in each segment. Block sizes differ by at most
one trial where the total is not divisible by 32. The regime label is provided
by the design; discovering a useful regime in real data is an additional task.

- **Fixed timing shifts:** arm means are -delta/2 and +delta/2. Delta is 0.125,
  0.25, 0.5 or 1 tag bin: approximately 9.77, 19.53, 39.06 or 78.13 ps.
- **Reversing shifts:** the same arm-mean difference alternates sign between
  regimes. Equal mixing cancels the pooled timing contrast in expectation.
- **Rate controls:** p0 = p(1-g/2), p1 = p(1+g/2), with g = 0.5. The reversing
  version alternates that gap, retaining the same mean click rate.
- **Null with drift:** the common baseline click probability ramps from 0.5p to
  1.5p and the common timing mean from -3 to +3 bins. Neither depends on X.
- **Transfer failure:** both training regimes have the same positive shift;
  test regimes have opposite shifts. The feature learned in training no longer
  captures the test dependence. This is an intentional low-power counterexample.

Sampling binomial setting counts, binomial clicks, and multinomial timing bins
is exactly equivalent to compressing the independent trial model into counts.
It does not replace 53 million Bernoulli trials with 53 million detected photons,
nor use a Poisson approximation. Trials and no-click counts remain recoverable
for every block and setting. No real event excerpts are committed.

## Detectors and finite-sample calibration

All four detectors use the **same held-out test segment**:

| Detector | Frozen receiver-local feature f |
|---|---|
| Click count | +1 for a click, 0 for a no-click |
| Pooled timing | -1 for arrival before zero, +1 otherwise, 0 for a no-click |
| Trained timing | Pooled timing multiplied by a separately learned sign for each regime |
| Trained histogram | Separately learned -1/0/+1 weight for each regime and timing bin |

Training uses differences of counts divided by the number of assignment-arm
trials. The histogram learns the sign of each bin's difference. The timing
detector learns the sign of the weighted early/late difference. Ties get zero
weight. Neither learner receives the injected shift, its direction, or any test
outcomes. After training there is no feature selection, retuning, early stopping
or threshold tuning on the test segment. All regimes and bins were specified
in the simulation design. These are simple trained scores, not a claim to have
optimized an arbitrary-record classifier.

For each test trial, the signed score is (2X-1)f(record, regime). Under the
synthetic null, the held-out setting bits are independent fair coins independent
of the entire held-out local record and training data. Conditional on those
records and frozen features, each nonzero score is an independent fair sign.
If m scores are nonzero, their number of positive signs is Binomial(m, 1/2).
We compute the exact two-sided binomial p-value. With m = 0 it is 1.

The four tests each reject at 0.01/4. By the union bound, the chance of any
rejection is at most 0.01 under this null, without assuming independence between
the tests. This conditional argument allows arbitrary *fixed, setting-independent*
clock/rate drift. It does **not** establish validity for arbitrary device memory,
assignment bias, or feedback from settings into subsequent local records. It is
not the history-robust confidence sequence used by the NIST causal audit.

No-click trials have zero score but are kept in the trial model and training
denominators. Conditioning on the number of nonzero scores here is justified
by the synthetic null's assignment independence; it is not a coincidence cut.
The record feature never contains the sender's outcome. These tests detect
dependence; they do not measure a causal effect or give an upper confidence bound
on total variation. A test that fails to reject supplies no general exclusion.

## Reproduction and reporting

```sh
python -m pip install -r studies/synthetic-timing-shifts/requirements.txt
python studies/synthetic-timing-shifts/design.py
python studies/synthetic-timing-shifts/design.py --check
python -m unittest discover -s tests -p 'test_synthetic_timing_shifts.py' -v
```

The fixed seed and per-scenario seeds reproduce 400 independent simulations for
each of 26 scenarios. JSON includes integer rejection counts and pointwise 95%
Clopper-Pearson intervals for every detector and the any-test decision. A zero
rejection count does not imply a zero false-alarm probability. The scenario grid
is a power study, not 26 opportunities to claim an experimental discovery.
Pinned NumPy/SciPy versions are required for byte-for-byte numerical artifacts.

## Findings in the committed run

At the larger exposure, the trained timing score detects the 0.25-bin reversing
shift (19.53125 ps under the assumed unit) in 352/400 runs, or 88%; pooled timing
rejects in 2/400 and counts in 3/400. At 0.5 bin, trained timing detects all
400 injections. These rates are conditional on the assumed Gaussian width and
the supplied regime label; 400/400 is not a guarantee of perfect future power.

The more flexible histogram learner detects the 0.25-bin reversing shift in only
106/400 runs (26.5%). Its additional learned weights cost power with sparse
training clicks, although it detects all 400 larger-exposure reversing rate-gap
injections. A simple physically motivated feature can outperform a flexible one.

When the training-to-test regime relationship breaks, trained timing detects
0/400 larger-exposure injections despite the large one-bin shift. The histogram
detects 13/400. Thus this study demonstrates a useful *conditional* detector and
an explicit blind spot, not comprehensive detection of all sign-changing effects.

Across the four null controls, any-test rejections are 5, 3, 1 and 1 out of 400.
The 5/400 estimate (1.25%) exceeds the nominal 1% numerically but its 95% Monte
Carlo interval includes 1%; no seeds or thresholds were changed to suppress it.
Finite-sample calibration follows from the fair-sign argument above, with unit
tests enumerating small null distributions. These small simulation samples
alone cannot validate a 1% tail probability precisely.

## What would be needed next

Before a real-record analysis: justify the setting-assignment model; replace this
single-click timing law with calibrated local records; preserve no-clicks and
multiclicks; account for timestamp failures and missing records; and bound timing
alignment and ordinary leakage. Choose features and the test population before
opening a fresh held-out segment or independent run. The existing archive has
already been explored, so simply naming a segment "held out" cannot undo prior
inspection. A split alone does not provide independence under device memory.

The useful finding here is whether particular injected effects can survive
sparse counts and training costs. Nothing here establishes that the real
experiment contains such effects, that the chosen regimes exist there, or that
these shifts would be spacelike influences.
