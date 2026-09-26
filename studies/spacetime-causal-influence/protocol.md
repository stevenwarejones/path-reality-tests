# Prospective optical intervention protocol

This is a protocol for discussion with an optics collaborator. It is not a
registered experiment, a claim of apparatus access, or an observation of a
superluminal effect. Every numerical specification below is hypothetical.

## Choice of protocol

| Feature | Separated optical source/modulator/receiver | Earlier receiver record, later choice |
|---|---|---|
| Intervention | Fresh local bit selects a short amplitude-modulation pulse at A | Fresh bit actuates the same modulator only after B's record is secured |
| Primary record | Any predeclared B click in a fixed early window, with all other cases mapped to zero | Previously latched B bit; no subsequent relabeling |
| Causal null | Local nonselective operation outside B's causal past leaves its marginal invariant | Future choice neither changes the earlier response nor shares/predicts its preparation |
| Closest reviewed work | Optical precursor observations; fast randomized spacelike optical Bell geometry; Fermi detector response | Delayed-choice causal models and no-backward-signaling marginal logic; ordinary quantum eraser coincidences measure another quantity |
| Measurement-to-model premises | Complete local records, fresh assignment, full support certificate and nuisance coupling bounds | The same assignment/record premises plus earlier record fixation and no B→RNG data path |
| Sensitivity | Early interval can bound binary influence; same hardware gives a delayed in-cone positive control | IID association intervals have the same count scaling; backward-effect identification needs an additional causal model; no causal positive control on the fixed past bit |
| Practical limitation | Separation, synchronization, support tails, leakage and detector recovery need measured bounds | Record attestation does not establish hidden setting independence; delayed optics does not predict a nonzero past marginal |

Develop the separated optical protocol. The reason is measurable diagnostic
access: a later window can test the actual control-to-detector response on the
same link. The later-choice protocol is retained as a legitimate null-control
variant, with exact record-order conditions and countermodels in the
[physical bridge](physical-bridge.md). It supplies no stronger ontology
exclusion without additional assumptions, and is not a proposed way to observe
retrocausal trajectories. Neither route has a standard-QM positive signal
outside the relevant causal cone. The normalized g-family specifies statistical
targets only, not a new dynamical theory.

## Optical layout and required equipment

Use a stable optical source feeding an amplitude modulator at A, an attenuated
free-space transmission path, and a fixed detector/readout at B. A fresh
hardware random choice selects one of two specified voltage pulses, with a
common baseline restored by the end of the intervention support.
Keep the preparation entering the modulator fixed. A phase-only element with
an intensity-insensitive receiver is a poor positive-control choice; amplitude
modulation gives a directly measurable delayed response. A weak coherent source
is sufficient for this causal engineering test. It does not test uniquely
single-photon or entanglement physics. A heralded source is an optional variant
only if eligibility is fixed before X and independent of it.

Required components: source and calibrated attenuation, fast amplitude modulator
and driver, characterized fresh randomness source, surveyed free-space path,
fixed receiver detector, local time-tagger/FPGA latch, independent local storage,
clock synchronization with bounded errors, pulse/RF probes, shielding and
dummy-load controls. Detector bandwidth, efficiency, dark counts, dead time and
afterpulsing must be measured. No manufacturer, price, laboratory availability
or achieved performance is asserted.

The source can illuminate B before the current intervention. Such detections
are allowed background; dependence on the new setting is the target. Include
all earliest setting leakage and the entire controlled drive pulse in A's
spacetime support. Residual setting-dependent drive after the nominal pulse
must be bounded or included by extending that support. Later storage of X is
needed for reconciliation; it does not retroactively extend the intervention's
earliest emission time. Any storage/communication activity that could reach B
before its latch is part of the leakage analysis. Clock/trigger signals
are distributed independently of X. Disable feedback from B into the RNG.
Do not communicate the current setting to the B acquisition computer during
the primary support interval. Joint analysis of separately secured records
happens later; B's bit was locally available when fixed.

## Hypothetical timing and distance budget

Times refer to a shared calibrated clock relative to nominal fresh choice.
Each interval includes its proposed uncertainty allowance; the entire interval
must be calibrated, not inferred from nominal component speed.

| Item | Hypothetical support / target | Needed evidence |
|---|---|---|
| Fixed input preparation and clock | Established before −100 ns; no present X dependence | Source stability and isolation; prior settings included as history |
| Setting-bearing RNG, driver and modulation | −5 to +20 ns; includes earliest possible leakage | Bit-generation latency, predictability bound, driver traces and settling |
| B detection gate | +40 to +60 ns nominal | Time tags tied to detector interaction, window calibration |
| Complete detector-to-local-latch support | +35 to +65 ns including errors | Worst-case latency and tail probability; failed latches accounted for |
| Center separation | 30.0 m | Survey |
| A/B support radii | 0.5 m each | Full cable/electronics/interaction envelope |
| Distance uncertainty | 0.1 m | Calibration coverage |
| Primary trial rate | 10,000 Hz (100 μs spacing) | Recovery/afterpulsing and storage measurements |
| Positive-control gate | +105 to +130 ns in separate predeclared control runs | Overlap with the delayed modulation pulse; verified latch timing |

The early support has a 26.40 ns conservative spacelike margin under these
inputs. At 3 m the same setup fails the geometry condition, so tabletop work
is a preliminary electronics/control demonstration. At roughly 100 ns free-space transit, a modulation pulse during 0–20 ns
reaches the receiver around 100–120 ns; the positive-control gate must overlap
that transient. Positive controls use separate predeclared runs to avoid an
early detection's dead time suppressing the delayed response. This window is
inside the future cone and is a functionality check. It is not included in
primary counts.

The trial spacing is not evidence of independence. It allows time for a pilot
to characterize recovery. The primary memory-valid score remains conditional
on fresh present randomization; an afterpulse that depends only on earlier
history does not violate its null. Present-setting electronics cross-talk does.

## Trial definition and immutable records

Choose N clock-triggered trial IDs in advance. Each initiated trial records X,
the complete local raw category and the predeclared binary mapping at B:

| Raw category | Binary B |
|---|---:|
| Exactly one eligible selected-detector click in the primary gate | 1 |
| No click | 0 |
| Multiple clicks, out-of-window events, or flagged invalid detection | 0 |
| Missing local bit or irrecoverable record | No silent deletion: terminate the claim or use a predeclared worst-case corruption bound |

Store raw categories, timestamps, trial ID, source/control status and diagnostics.
The binary record is complete for every properly recorded initiated trial;
multiple/invalid outcomes are not dropped from denominators. Instrument failure
is not harmless just because its category is mapped to zero: setting-dependent
failure can produce a gap and must enter the response/null or nuisance model.

Prefer no heralding for this initial protocol. In a heralded implementation,
eligibility must be fixed by a source event before current X, with a common
preparation distribution after heralding. A future remote coincidence is never
an eligibility criterion. Retain all eligible heralds including missing photons.
No subtraction of fitted backgrounds or fractional pseudo-counts is permitted
in the primary binary analysis.

Secure each local latch and storage log before cross-station data reconciliation.
Digitally signed or chained records detect later edits only under their stated
trust model. Their timestamps and acquisition electronics still require physical
calibration. An unbounded corruption rate prevents a spacetime interpretation.

## Controls and calibration

Before freezing the design, use independent pilot trials to check the delayed
positive response and determine a defensible effect/nuisance range. Do not use
pilot outcomes in the confirmation sample. Freeze source intensity, modulator
voltages, window edges, category rules, N, α, and the setting-predictability
allowance ζ before acquisition.

Run source-off, blocked-beam, dummy-driver-load, disconnected-optical-path and
fixed-modulator controls, retaining every trial. Compare randomized electrical
drive versus actual optical switching to isolate RF/power-supply leakage.
Change cable routing and shielding in calibration, then freeze the primary
configuration. Check detector recovery and correlations with previous settings,
not only the current one. These controls diagnose mechanisms; transporting a
bound from a control to the live optical system remains an explicit premise.

Measure timing error tails and event mass near both fixed gate edges; this
converts window uncertainty to outcome mismatch. A Gaussian fit or nominal
bandwidth alone is inadequate. Allocate calibration failure probability across
every reported bound, including shared calibrations. A zero-failure IID
calibration with m trials gives one-sided CP upper failure rate
1−αcal^(1/m); a small bound therefore costs real calibration trials. Neither
zero observed failures nor a manufacturer jitter number supplies an exact zero.

Analyze one primary early gate with `score_interval` and `classify`. The
[statistics derivation](statistics.md) specifies memory, bias, calibration and
power assumptions. Report an upper bound even if no rejection occurs. Optional
IID two-sample intervals are a labeled secondary comparison. If additional
confirmatory windows are selected, freeze them and split α before collecting
confirmation data. Raw unsmoothed local records define the primary bit; centered
filters can import later information into an apparently earlier window.

## Decision and feasibility

If geometry fails, report an ordinary causal-channel/control experiment. If
nuisance or predictability cannot be bounded below the target gap, stop the
exclusion claim: more shots cannot resolve a model-class overlap. If the
primary interval strictly exceeds the justified nuisance ceiling, first test
timing, selection, source dependence, record integrity and ordinary leakage
using the frozen controls; call it a candidate operational anomaly pending
independent reproduction. No single positive table identifies a consistent
relativistic field theory or proves retrocausality.

The software establishes a justified conditional resource study. It does not
establish achieved optical sensitivity. The external inputs still missing are
the surveyed component supports, local latch timing, conditional RNG model,
tail/response calibration with coverage, control-to-primary transport bounds,
detector recovery records and the available complete-trial rate. A collaborator
can use these concrete requirements to accept, enlarge or reject the proposed
30 m envelope before any physics claim is made.

The primary memory analysis needs the calibration envelope to hold conditionally
on device history, or to bound its realized average with stated coverage. The
IID zero-failure calibration calculation above does not establish that transport
by itself. A characterized device envelope or a separate sequential calibration
argument is required; otherwise restrict the entire inference to a justified
IID regime. This missing apparatus input is included in ST06/ST07.
