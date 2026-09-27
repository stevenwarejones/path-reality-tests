# Matched gamma follow-up: a quantified acceptance requirement

**The requested research outcome remains incomplete.** The matched gamma records
support a useful conditional sensitivity calculation, but not a full-record
physical ambiguity or certified susceptibility separation. The missing pieces
are now more specific than “no counts”: the publication supplies the counts;
incident-event acceptance and the observation kernel remain uncalibrated.

## Corrections to the first audit

[Larson et al., Appendix A.3 and Table 2](https://arxiv.org/html/2503.07354v1#Sx4.T2)
state that both chips used the same cold finger and cooldown, and provide all
12 footprint numerators and accepted-window denominators. The acquisition
exposures are 11.78 h (bare) and 5.35 h (1 µm Cu). The earlier gate incorrectly
said these denominators were not supplied. They are absent from the figure CSVs,
but present in the paper. The acquisition now extracts and hash-checks that
HTML table, parses every count/error field, and reconstructs all 12 published
contrasts to their rounding precision.

| Qubit | Bare: switches / accepted windows | Cu: switches / accepted windows |
|---|---:|---:|
| Q1 | 1351 / 2889 | 123 / 1149 |
| Q2 | 1306 / 2621 | 219 / 1398 |
| Q3 | 1342 / 3150 | 248 / 978 |
| Q4 | 1053 / 2593 | 160 / 1208 |
| Q5 | 904 / 2344 | 718 / 1445 |
| Q6 | 1433 / 3014 | 145 / 1238 |

The charge sensors are bare Q2 and Cu Q5. These are **different accepted event
sets**, because an event is excluded for a qubit if any part of its coincidence
window is masked. Dividing the columns by their maximum denominator would not
recover retention: that maximum is neither the total trigger count nor the
number of incident impacts. Nor do the twelve marginals supply their eventwise
intersections or a multinomial likelihood.

There is a second correction. Our previous result JSON applied a net-parity
moment bound directly to the footprint contrasts. The paper describes counting
whether a switch occurs anywhere in an accepted window. That indicator is not
generally endpoint parity. We removed those empirical any-tunneling bounds;
the theorem in [model.md](model.md) remains valid for its specified observation.
This is a correction to our application, not a demonstrated error in the paper's
processing: the released excerpts cannot establish how multiple nearby events
were resolved in the complete acquisitions.

## What the temporal records actually constrain

The complete gamma ZIP has 57 data files. It supplies figure traces, rate/error
curves, spectra, simulations and conventional T1 summaries, not the complete
footprint event catalog, masks, raw shots, extraction code or run log.

- Each Fig. 4 trace contains 21,001 displayed samples spanning 15.452388 s,
  separated by 0.735828 ms. The displayed coordinate omits charge-bias reset
  interruptions; it is not a continuous acquisition wall clock.
- The footprint processing uses a 100-point moving average and coincidence
  window (nominally 73.5828 ms), plus a 200-point charge-step filter. Charge bias
  is reset every 5000 shots. Readout is staggered between qubits, with 50–100 µs
  waits. These delays matter relative to the cited approximately 100 µs recovery.
- A diagnostic partition into 100-interval blocks gives 188/197 valid blocks for
  the two bare streams and 210/144 for Cu. None of these valid blocks contains
  two digital switches. Thus OR and endpoint parity agree in this small
  partition. This is **not** a bound on their disagreement in the full data,
  a reconstruction of Table 2, or a detection-efficiency calibration. Midlevel
  digital values are treated as masked, and the partition is explicitly ours.
- The dose-sweep parity processing instead uses 40-point windows (with a stated
  exception), dose-dependent repetition periods, and 53–93 ms coincidence windows
  with source versus 400–413 ms without. These are different dose populations and
  filtering settings, not several known integration windows applied to one
  stationary latent rate distribution. Repeated-dose rows have no deposited run
  chronology that would justify treating them as independent replicates.
- Conventional T1 versus distance uses separate repeated delay sweeps. It is not
  charge-triggered recovery, and cannot supply a second Laplace-transform value
  for the same selected impacts. The approximately 100 µs recovery quoted in
  the discussion comes from earlier injection work, not these T1 summaries.
- The four-way charge calibration in Fig. 12 was acquired months after the
  earlier charge sweeps. Same cooldown alone does not validate a shared
  charge-acceptance function across these protocols or epochs.

Consequently the proposed multiple-window moment inference cannot be instantiated
from these records without additional assumptions. We do not optimize a moment
problem with falsely matched windows. The 10 µm injection archive and mechanical
apparatus supply no calibration of this gamma observation kernel; their previous
reconstructions are unchanged.

## Conditional background-timing sensitivity

To quantify the consequence of the observation distinction, temporarily assume
faithful digitization and one signal-bearing sampling interval j. Let S be its
odd signal switch, with probability q≤1/2, independent of a background switch
vector B over the coincidence window. Background entries can be dependent. Put
b=P(B≠0), r=P(B=e_j), and let Y indicate any switch in B XOR S e_j. Then exactly

\[
 p=P(Y=1)=b+q(1-b-r),\qquad 0\le r\le b,
\]

because toggling j creates a detection only from B=0 and erases one only from
B=e_j. For c=2q this implies

\[
 \frac{2(p-b)}{1-b}\le c\le
 \min\!\left(1,\frac{2(p-b)}{1-2b}\right).
\]

The upper endpoint puts all background switches in interval j; the lower puts
them outside j. Both are nonnegative background path distributions, and the
unit tests enumerate the paths independently. They are **observation-model
extremizers**, not physical models fitted to all spectra, coincidences and T1.

The paper uses b=HMM switching rate × window length and the upper expression.
A rate-times-window product is an expected count, not automatically the
probability of any switch. We provisionally equate them only for this sensitivity
calculation. A stationary Poisson process would instead give 1−exp(−rate × window),
but Poisson HMM switch times and negligible filtering loss have not been validated.
Likewise, mask selection could correlate signal and background, violating the
independence assumption. The displayed bounds do not cover those failures.

At the published count frequencies, moving background switches from the same
interval to a different interval changes inferred bare Q1 contrast from 0.90649
to 0.74160, a difference of 0.16488 (published contrast error 0.04). For Cu Q6
the change is only 0.00478. This comparison quantifies a calibration sensitivity;
it neither estimates the actual correction nor assigns significance to it.

## How much acceptance calibration would suffice?

Use every count frequency k/n with the same exploratory rectangle: ±4 times
its published observation error, and ±4 times the published background error.
Take the union of the physical contrast intervals above throughout that
rectangle. These are deterministic envelopes, **not** simultaneous confidence
intervals or independent binomial inference. No independence between qubits
or overlapping windows is required for the algebra.

For a specified common impact population, suppose at least a fraction ρ is
retained on each chip by **charge triggering and parity masking together**.
A bounded per-impact response h∈[0,1] with selected mean c then has population
mean in [ρc, ρc+1−ρ], allowing arbitrary selection dependence. If l_B is the
bare selected lower bound and u_C the Cu upper bound,

\[
 C_B-C_C\ge\rho(l_B-u_C)-(1-\rho).
\]

Thus a positive population difference is guaranteed when
ρ>1/(1+l_B−u_C), provided l_B>u_C. The bound is sharp over unrestricted bounded
responses: set rejected bare responses to zero and rejected Cu responses to
one. At the threshold those constructions have equal population means. These
are selection extremizers, **not full-source phonon/QP countermodels**.

| Approximate distance | Bare / Cu qubits | Selected contrast gap lower bound | Required common retention |
|---|---|---:|---:|
| 2.03 mm | Q4 / Q3 | −0.07067 | No positive gap certified |
| 4.62 mm | Q1 / Q6 | 0.36142 | >73.4528% |
| 4.06 mm | Q6 / Q1 | 0.42361 | >70.2438% |
| 5.36 mm | Q3 / Q2 | 0.31217 | >76.2095% |
| 6.66 mm | Q5 / Q4 | 0.25165 | >79.8945% |

All five non-sensor distance pairs are reported; the strongest pair is not a
newly held-out target. Distance matching is approximate and does not establish
identical geometry or event distributions. Even with equal retention, calling
C_B−C_C a *susceptibility* difference requires the same incident energy/position
population and a comparable response functional. The paper's approximately
1 mm charge-sensitive radius does not establish the required retention over
an 8 mm chip, nor a universal lower acceptance probability inside that radius.
The empirical separation therefore remains **conditional**.

## Decision and next discriminating data

The matching test failed for multiple windows, while the selected-footprint
comparison survives the limited background-timing sensitivity above. It does
not yet survive unrestricted charge/mask selection. A sufficient follow-up for
this particular contrast comparison would supply (1) an independently defined
common impact population and calibrated combined retention exceeding the
appropriate threshold, and (2) a measured signal/background-to-HMM response
kernel using synchronized windows and the same masks. Full trigger catalogs,
rejected-window flags and matched off-trigger windows would constrain the latter;
they alone cannot reveal charge-silent incident impacts. A physical mapping to
any-tunneling or relaxation additionally needs its own response calibration.

No available complete-source model has been fitted here, so there is no joint
physical-feasibility certificate, source-removal witness, or minimal-measurement
theorem for the full experiment. The increment is a corrected acquisition model
and a quantified acceptance requirement. It does not meet the review's completion
gate, and PR #19 remains draft/incomplete. Proceeding to a discovery claim would
require inventing precisely the acceptance and response sharing the review
prohibits. We preserve the conditional result without presenting it as that claim.
