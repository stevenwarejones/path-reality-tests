# Combination gate

## Target and alternatives

The preferred target is the **probability of at least one quasiparticle tunneling
event at a specified qubit after a charge-triggered radiation impact**, following
a withheld mitigation or source setting. This separates an initiating source
from chip susceptibility and avoids equating parity with coherent phase errors.
A successful prediction would calibrate the response using a different measured
intervention and propagate source, transport and readout uncertainty without
refitting the withheld outcomes.

Four primary measurement families in three inspected archives are considered.
The charge and parity families below share an archive and apparatus; they are
not separate independent experiments. Their overlap would have to be retained
in any joint likelihood. Derived spectra and figure panels are not additional
independent datasets.

## Actual records and measurement units

| Family and record | Denominators/exposure actually available | Intervention and apparatus | Clocks, geometry and calibration | Role and limitation |
|---|---|---|---|---|
| Larson charge: [16568030](https://zenodo.org/records/16568030), `Figure02/Figure2b_nonCu_Q2.csv`, `Figure2b_1umCu_Q5.csv` | 10 rate/error rows per chip; full run denominators not in these tables. Fig. 2a includes example tomography arrays, not every dose run | External Co-60 distance sweep, with no-source row; separate bare and 1 µm Cu chips | Dose coordinate 1/r² in m⁻²; rates s⁻¹; 8×8×0.525 mm Si; charge threshold and charge transport matter. No common clock with other laboratories | Initiator-sensitive response; charge efficiency is not a universal particle tag |
| Larson parity/footprint: same record, `Figure03/Figure3b_*.csv`, `Figure05/Figure5c_*.csv` | 30 bare-chip and 162 Cu-chip rate/error entries; 6 footprint summaries per chip. Paper Table 2 supplies every coincidence numerator and accepted-window denominator, with 11.78 h/5.35 h exposures; event identities, mask intersections and run grouping are not supplied | Same cooldown and cold finger, different chips; bare dose qubits Q3/Q4/Q6, Cu Q1–Q6 | Sampling varies 250 µs–10 ms in the paper; HMM and coincidence windows are protocol dependent. Source at about 59 cm for footprints. Fig. 4 traces are excerpts | Measures parity contrast and fitted spectral rates; no complete event likelihood reconstructed |
| Iaia injection and mitigation: [7249678](https://zenodo.org/records/7249678), `Figure2a.dat`, `Figure3b.dat`, `Figure4b.dat` | 30 injection-delay rows; 7 observed/extracted coincidence categories, with errors. Example traces and spectra exist, but not the complete coincidence acquisition | Junction injection on bare versus **10 µm** Cu chips; distinct from the 1 µm mitigation comparison above | Fig. 2 uses delay in µs; 1 mV injection; rate response is inferred from T1 fits. Fig. 2 errors derive from fit confidence intervals, not necessarily one standard deviation | Calibrates a particular injection response; replacing 10 µm with 1 µm copper is not a known numerical scale transformation |
| Kono mechanical controls: [11034817](https://zenodo.org/records/11034817), `Zenodo/Fig3aceFig5Fig7FigS11/data/statistics_wrapped_period` and analysis notebooks | 714 phase bins × 4 integer cells, total 12,697,600 measurements. These aggregates do not retain trial/run clustering. Larger quadrature and transition records are listed but not all acquired | Pulse-tube on/off, synchronized acceleration and qubit readout; Nb-capacitor/Al-junction transmons in a separate refrigerator | Notebook period 0.7144833333 s; binned samples 1 ms; jump protocol 3 µs. Q0/Q1 physical labels differ from some internal indices. Accelerometer is at the top plate, not a calibrated GHz phonon probe at the chip | Strong within-apparatus causal control; no measured transfer to the Larson/Iaia chips |

The handoff's gamma DOI 16568029 is the **concept record**; the API resolves it
to version 16568030. That resolution is now pinned. The gamma footprint CSV
headers say distance in mm, but contain values around 2020–6660. The 8 mm chip
geometry and the paper's 2.02–6.66 mm separations support interpreting these
numbers as µm. Results preserve the original values and label the conversion
as an inference. No probability bound uses that distance conversion.

The [matched gamma follow-up](matched-gamma.md) corrects the denominator and
observation-model omissions in the first audit, checks the proposed multiple-window
route, and quantifies conditional charge/mask-acceptance requirements.

## Feasibility calculations at the supplied precision

**1. Initiator reduction versus susceptibility reduction.** In a restricted
stationary dose model, charge and parity rates obey

\[
 C_d(z)=c_d+q_d A_d z,\qquad P_{dq}(z)=p_{dq}+h_{dq} A_d z.
\]

Both have units s⁻¹; z=1/r². A common source-only explanation can trade impact
amplitude A against charge acceptance q or poisoning/readout response h in
either family separately. Jointly, the slope ratio cancels A, but identifying
h across chips still requires a bound on q and a valid model for the rates.
The full stored curves must admit a common uncertainty construction first.

For every curve, the code searches exact three-row annihilating contrasts (and
two-row contrasts at repeated doses). If w·1=w·z=0, every affine function must
satisfy |w·y|≤k Σ|w|e for reported errors e. This is a global algebraic bound,
not a failed optimizer. The solver coefficients are converted to exact rationals and every original
decimal residual is evaluated rationally to certify a feasible upper bound.
Both exact endpoints and coefficients are recorded; no numerical cushion is
used as an exactness argument.

| Curve | Minimum multiplier of reported errors | Minimum additive rate discrepancy after a 4-error envelope (s⁻¹) |
|---|---:|---:|
| Bare charge Q2 | 1.420 | 0 |
| Cu charge Q5 | 1.916 | 0 |
| Bare parity Q3 | 100.444 | 0.76550 |
| Bare parity Q4 | 57.729 | 0.61237 |
| Bare parity Q6 | 119.534 | 0.66992 |
| Cu parity Q1 | 72.724 | 0.06708 |
| Cu parity Q2 | 30.580 | 0.03239 |
| Cu parity Q3 | 129.329 | 0.21166 |
| Cu parity Q4 | 194.873 | 0.13605 |
| Cu parity Q5 | 128.400 | 0.05805 |
| Cu parity Q6 | 97.844 | 0.05531 |

Thus the naive joint affine region is empty even before demanding a shared
source/response relation. It cannot certify useful combination gain. These
figures **do not reject radiation or the paper's qualitative dose trend**.
Fit errors, drift, nonlinear response, protocol changes or incomplete covariance
could all matter. The diagnostic says how much rate discrepancy this particular
model needs; adding it post hoc does not create a calibrated confidence region.
A Gaussian/Poisson significance or 95% simultaneous coverage is not asserted.

**2. Few intense bursts versus many weak responses.** Both can give the same
parity contrast. For contrast c=0.20, homogeneous integrated intensity
Λ=−log(0.8)/2 gives any-event probability 0.10557. A mixture of Λ=0 and Λ=5,
with active weight 0.20/(1−exp(−10)), gives exactly the same contrast and
any-event probability about 0.19866. All intensities and weights are nonnegative.
These are synthetic scalar-observation witnesses, **not** models checked against
all source-only constraints. The distinction remains larger than the deposited
0.02 error for one of the distant Cu footprint points, but neither that error nor
a four-error envelope is a complete statistical confidence statement.
The injection measurements could separate the two only if they calibrated the
same event-conditioned intensity distribution or a defensible support bound.
That bridge is not established by the inspected records.

**3. Shared periodic drive versus coincident microscopic bursts.** Mechanical
phase conditioning is directly measurable: pooled MI 0.00278066 bits versus
weighted phase-conditional MI 0.00008803 bits. The law of total covariance
assigns 84.18% of covariance to varying phase-conditioned means. The original
paper already makes the qualitative distinction. Remaining within-phase
covariance is positive; aggregated counts cannot establish its run-level
significance or assign it to QPs rather than TLS/readout effects. A gamma
experiment in another refrigerator cannot supply those missing controls.

## Candidate combinations and decision

| Combination | Potential separating information | Decision |
|---|---|---|
| Gamma charge + gamma parity, same apparatus | Source-normalized response and spatial footprint | Best empirical anchor, but no globally feasible affine model at the supplied errors and no complete footprint likelihood. Retain for a calibrated reanalysis |
| Junction injection + gamma footprint | Predict radiation response from transport calibration | Stop current fit: copper thickness differs in the downloaded injection archive, trap/boundary terms are device dependent, and footprint conversion in the closest paper is outcome-fitted |
| Mechanical phase controls + gamma + injection | Distinguish changing source strength from QP transport and detection | No justified numerical cross-chip response relation or calibrated mechanical high-frequency energy spectrum; combining likelihoods would introduce arbitrary per-chip freedom |
| Localized infrared injection as an alternative calibration source | Different positions/powers and pulsed relaxation might test transport | [11548753](https://zenodo.org/records/11548753) and its Figure 3/SM8 notebook were screened. This is another apparatus and relaxation experiment; it does not supply the missing same-chip charge/parity intensity calibration. Not used as a fourth fitted source |

The infrared archive advertises raw data, processed HDF5 and plotting notebooks.
Only its metadata and Figure 3/SM8 notebook were acquired here; no infrared
measurement reconstruction is claimed. The notebook reads already fitted
trapping, recombination and initial-density arrays. Its availability is a
promising independent validation lead, not evidence of transportability.

## Quantified requirements to reopen this gate

1. **Statistical model:** Table 2 now supplies marginal valid-window counts.
   Obtain their event identities, joint masks and run grouping, or justified simultaneous bounds on the reported means that include
   the discrepancies above. Preserve dependence between charge/parity records.
   Fit-error bars alone do not meet this requirement.
2. **Parity-to-event calibration:** to approximate any-event probability by the
   parity contrast within 5% relative error, a sufficient condition in the
   mixed-Poisson model is Λ≥log(20)≈2.996 for every active event, allowing a mass
   at zero. This has not been established for distant mitigated qubits. A known
   upper bound Λ≤L instead gives the sharp interval in [model.md](model.md).
3. **Non-ionizing fraction:** for consistently defined, detected poisoning
   events and a radiation-tag efficiency lower bound q_min>0 for *every relevant
   radiation class*, R_non≥max(0,R_total−R_tag/q_min). To establish fraction f,
   one needs q_min>R_tag/((1−f)R_total), with uncertainties included. A roughly
   1 mm charge-sensing radius on an 8 mm square gives only a geometric area
   fraction about 0.049, not a pointwise efficiency floor. Neither muon-only
   tagging nor absence of a charge jump bounds all untagged radiation.
4. **Transfer and holdout:** for the 1 µm/bare devices, establish which injection
   calibration runs belong to the radiation chips, including cooldown changes,
   and bound the transport/readout parameter drift. Freeze a complete dose,
   position or intervention as holdout before accessing its outcomes in a new
   fit. The outcomes in this exploratory gate are already inspected and cannot
   now be called a blind validation set.

No measured joint bound with full source-only countermodels is delivered.
No source removal, optimizer failure or finite model bank is presented as proof
that a full source is insufficient. The observation theorem is conditional and
standard in mathematical character; novelty remains unestablished. This compact
gate is a reproducible stopping point, not a completed discovery.
