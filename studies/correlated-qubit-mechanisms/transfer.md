# Intervention transfer is now the primary target

**The first explicit forcing-to-response transfer fails.** A model fitted to
matched pulse-tube-on/off records predicts far too many observed flags under the
controlled-shock reference. A saturating extension is exactly indistinguishable
on all pulse-tube calibration inputs, yet gives a very different shock prediction;
it also fails. These are neither successful physical predictions nor adequate
physical countermodels establishing non-identifiability.

This is a test of a specific transfer assumption, beyond improving the earlier
0.39% predictor gain. The remaining goal is an operational intervention response
that complementary records jointly determine, or empirically adequate physical
models that agree on retained observations but disagree on that target.

## Which acquisitions are actually matched?

The previous common-flag audit established a shared observed coordinate. It did
not establish a shared forcing trace or a physical observation law.

| Records | Matching supported by deposited notebook and arrays | Remaining limitation |
|---|---|---|
| Fig. 2 run 562 | 8,192 qubit samples, 1 ms apart, and a same-number acceleration/trigger file; notebook identifies PT on | One short acquisition; no full preparation/readout instrument log |
| Fig. 2 run 567 | Same schema and same-number acceleration/trigger file; PT off | One short acquisition; no randomized intervention sequence |
| Fig. 2 runs 774/782 | Longer qubit records used by the preceding audit | The notebook pairs their spectra with acceleration runs **483/480**, not same-number forcing records; these cannot be silently joined sample by sample |
| Fig. 6 runs 832/833 | Run 832 contains acceleration and a trigger; run 833 contains 50 paired qubit acquisitions with a common 1 ms clock | **Separate run numbers.** The notebook uses one forcing reference and folds qubit records at the programmed period; it does not supply a forcing trace for each qubit acquisition |

Both short PT qubit traces start at zero; their accelerometer clocks span almost
10 s, sampled every 0.1 ms, with the recorded trigger at 0.5 s. We align qubit
sample `t` with accelerometer time `t+0.5 s`. The shock reference samples every
0.05 ms and has its trigger at 0.5151 s. The notebook specifies a 1.5 s period and
folding origin 1.7 s on that reference clock, corresponding to 1.1849 s on the
qubit clock (nearest deposited sample: 1.185 s).

This is a **protocol-relative** time coordinate, not a calibrated per-acquisition
shock-onset time. A separate reference may describe repeatable forcing, but its
repeatability across the 50 qubit acquisitions is not measured by these files.

The common flag remains `u>0.5` in the G/E-referenced first projected quadrature.
Balanced G/E calibration and PCA centering place its boundary at zero. The
notebook's positive normalization constants do not change the flag. Identical
reference arrays establish this coordinate, not equal state preparation,
classification errors, pulse duration/amplitude, backaction, or intrinsic rates.
Saved projected traces and clocks do not reconstruct those instrument settings.
No causal classifier latency is established.

## Fixed calibration and transfer test

The [protocol](transfer-protocol.json) was recorded before computing the initial
fit and timing outputs. Earlier aggregate intervention counts were already known;
this is retrospective and not a newly blind experiment.

After the initial uncentered test failed, an input audit found pre-trigger voltage
means of approximately −0.00949 V (PT on), −0.01596 V (PT off), and **−0.60368 V**
(shock reference). Squaring an unremoved record offset changes the forcing proxy.
We therefore repeated the entire fixed fit after subtracting each record's
pre-trigger mean, which is available before its qubit samples. This is a
**subsequent baseline diagnostic**, not a prespecified correction. The main tables
below use it; both complete outputs are retained. Baseline estimation uncertainty
is not calibrated.

The initial uncentered model predicted 740,937 paired shock flags, and its capped
extension predicted 1,277. After baseline removal the predictions become 648,620
and 1,073. Removing the offset matters, but does not rescue either transfer.

- Fit only samples 0–4095 of each matched PT record; retain samples 4096–8191 for
  calibration diagnostics. No shock qubit outcome enters parameter fitting.
- Convert deposited acceleration voltage using the notebook's factor of
  5 m/s² per volt. Apply a causal 1 ms amplitude low-pass, square it, and apply a
  causal exponential reservoir filter. This is an **acceleration-derived proxy**,
  not measured absorbed power or phonon energy at the chip.
- Select one common reservoir time from 1, 3, 10, 30, 100, 300 and 1000 ms using
  training likelihood. The selected time is 100 ms. Rescale the proxy using only
  training forcing values.
- Model each binary flag with nonnegative continuous-time excitation and recovery
  rates, `u_q=b_q+k_q E(t)` and `d_q`, held constant within each 1 ms transition.
  Use the exact two-state transition matrix. Parameters differ between qubits;
  onsets and flagged durations need not coincide or have identical marks.
- Propagate forward from each acquisition's observed initial pair. Subsequent
  target flags are not fed back into the forecast. The separate conditional
  transition likelihood diagnostics do use observed previous flags and are
  labeled accordingly.

The fitted triples `(b,k,d)` are approximately `(2.347,59.653,2431.324)` and
`(4.722,78.999,2733.250)` in the model's rate units. They are parameters for the
observed flag process, not calibrated intrinsic qubit rates. Four deterministic
initializations per grid point converged without contacting the numerical search
bounds; this is not a global optimum certificate.

If the binary state were a physical two-level population, these nonnegative
rates admit ordinary excitation and relaxation Lindblad jump operators. This
supplies coarse dynamical realizability. It does **not** establish that the saved
threshold flag is that physical state, that its observation response transfers,
or that the proxy determines those rates through known material physics.

## An exact calibration-support ambiguity

Define `M` as the maximum filtered proxy over both complete PT control records.
Compare the fitted linear response with the extension

`u_q=b_q+k_q min(E,M)`.

The parameters, initial conditions, memory filter and recovery rates are
unchanged. These two models have **identical transition kernels and complete
conditional path laws on every retained PT input**, because `E<=M` there. This is
not agreement of only one aggregate count. Above that measured range one response
continues to rise and the other saturates. The saturation point is a declared
extension, not a measured property of the device.

Here the calibration inputs are the matched 562/567 traces. Neither model has
been fitted to every other published PT record or to the full prepared-state
reference distributions.

After pre-trigger centering, `M=1.029518`; the shock-reference maximum is
6230.408, and **74.86% of shock transitions exceed the PT maximum**. The initial
uncentered diagnostic placed every shock transition above the PT range. The
corrected input still demands substantial extrapolation. A response law beyond
that range is not identified just from the observed PT forcing values.

Exact calibration equivalence is insufficient for the requested negative result:
a base model must also be empirically adequate and physically connected to the
measurement. The next two tables expose why these particular fits do not supply
that result.

## Calibration validation and exposure-normalized transfer

Each short-record validation segment supplies 4,095 transitions after its measured
initial state. The two extensions have identical predictions here.

| Validation record | A flags: observed / predicted | B flags: observed / predicted | Paired flags: observed / predicted |
|---|---:|---:|---:|
| PT on | 55 / 40.63 | 60 / 50.12 | 2 / 0.596 |
| PT off | 6 / 4.05 | 8 / 7.19 | 0 / 0.007 |

This short validation does not establish aggregate adequacy. The model also
predicts more multi-sample flagged runs than occur. No independence-calibrated
confidence region is attached to these comparisons.

The shock target has **1,638,350** post-initial samples across 50 acquisitions.
All rates below use that same exposure.

| Shock target | Observed | Linear transfer | Saturating extension |
|---|---:|---:|---:|
| A flags | 4,216 | 763,621.61 | 37,422.06 |
| B flags | 7,594 | 784,531.45 | 44,812.99 |
| Paired flags | 32 | 648,619.92 | 1,073.16 |
| Paired flags per million samples | **19.532** | **395,898.26** | **655.02** |

The exact-calibration-equivalent extensions differ by about 604-fold in predicted
paired events, but **both miss the measured intervention response badly**. They
cannot be described as empirically adequate competing explanations. Nor does
this failure reject all mechanical mechanisms: proxy choice, spectral coupling,
reference-run mismatch, saturation, observation response and dynamics could all
change the transfer.

## Timing and recovery targets, beyond the total

The protocol specifies six phase intervals, cross-qubit onset lags from −10 to
+10 ms, and complete flagged-run durations of 1, 2, 3–5 or at least 6 samples.
The [full output](results/transfer.json) retains exposure-normalized lag rates,
duration fractions, boundary-censored counts and every acquisition. A positive
lag means B starts later than A. A flagged run is not equated with a microscopic
excited-state lifetime.

The phase fold retains 21 complete cycles per acquisition, totaling 1,575,000
samples and 31 paired flags. The all-record table above includes the additional
unfolded samples and their one further paired flag.

| Protocol phase (s) | Exposure | Paired flags | Observed per million |
|---|---:|---:|---:|
| 0–0.05 | 52,500 | 0 | 0 |
| 0.05–0.10 | 52,500 | 0 | 0 |
| 0.10–0.20 | 105,000 | 10 | 95.24 |
| 0.20–0.40 | 210,000 | 11 | 52.38 |
| 0.40–0.80 | 420,000 | 7 | 16.67 |
| 0.80–1.50 | 735,000 | 3 | 4.08 |

The centered capped model predicts approximately 780 paired flags per million
throughout 0.10–0.80 s, and 598 per million in the late interval. It does not
reproduce the observed response shape. The linear model predicts
long, highly occupied intervals. Actual complete runs are much shorter: A has
4,110 one-sample runs, 50 two-sample runs and two 3–5-sample runs; B has 7,316,
133 and four respectively. Neither has a complete run lasting six samples. The observed one-sample fractions are 98.75% and 98.16%; the linear model
predicts 48.19% and 49.41%, and the capped model 89.52% and 91.31%. Both
models predict excessive persistence even after normalizing out the total count. Zero-lag onset pairs number 29, against
approximately 19,312 and 876 predicted by the linear and capped models. The
remaining lag bins are reported rather than assigning that mismatch to a cause.

These diagnostics prevent calling either model successful merely by adjusting
a coincidence normalization. They do not provide a mechanism separation with
calibrated uncertainty or a full-record physical exclusion.

## What each source contributes—and does not resolve

| Information | What it supplies | What remains unresolved alone |
|---|---|---|
| Prepared-state references | Common observed flag coordinate | Intrinsic state, preparation error and backaction separation |
| Qubit histories | Conditional observed transition and duration information | Dependence on an intervention's forcing and its extrapolation |
| Matched PT acceleration and trigger | Time-associated forcing proxy for calibration | Chip-level absorbed spectrum and high-forcing response law |
| Shock forcing reference | A conditional input for a proposed intervention forecast | Actual input to each separate qubit acquisition and repeatability |
| Combined retained inputs | A concrete conditional forecast and its failure; exact PT-support equivalence of two extensions | An adequate transferable physical likelihood or a demonstrated combination payoff |

Radiation and injection records still lack a justified bridge into this model.
They are not silently counted as constraints on its prediction.

The earlier shared-latent mark model remains an observation-level benchmark. Its
identical-mark joint-event competitor does not represent general common physical
disturbances. A later physical comparison must permit asymmetric excitation,
relative delay distributions and recovery, and give its shared environment
temporal dynamics. The present driven Markov attempt supplies explicit temporal
and asymmetric local dynamics, but **does not complete that competing-mechanism
comparison**. No physical winner is inferred from either analysis.

Controller comparisons are secondary and currently unsupported. The earlier
unequal-coverage selections remain descriptive. Any operational test must use a
common coverage or cost criterion, a logged causal phase estimate, and measured
classification/actuation latency; zero retrospective events is not a risk bound.

**Next decision:** constrain the forcing/observation bridge before adding more
fit flexibility. Recover the relevant preparation/readout settings and whether
the shock reference applies to each retained run; test frequency-sensitive
response only with parameters constrained by calibration records. A successful
transfer must predict phase response, onset lags and recovery together. If that
bridge remains unidentified, the next countermodels must fit the retained
observations adequately—not merely share the same inadequate calibration law.
