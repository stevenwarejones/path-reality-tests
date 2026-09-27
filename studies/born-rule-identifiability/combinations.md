> Main update: the [measured scan–Viviani combination](joint-tables.md) and
> [generated joint report](results/joint-table-report.md) supersede the prospective
> quarter-turn extension as the principal candidate. The failed pairings below
> are retained to document discovery and scope.

# Expanded combination search and explicit bridges

Retrospective search 2026-09-26. Starting datasets were not a closed list.
The strongest implemented result is now a **single-qubit probability/control
nonidentifiability theorem for an entire gate-depth family**. The waveguide
analysis is the implemented runner-up: an exact mean-table ambiguity and a
prospective phase-cycling measurement. Neither gives a new empirical exclusion.
Generated quantitative results are in [combination report](results/combination-report.md).

## Candidate inventory and verified access

| Candidate combination | Actual public measurements checked | Access/sufficiency decision |
|---|---|---|
| A: processor phase sweeps + repeated-gate settings | [Zenodo 7538941](https://zenodo.org/records/7538941): `lima1-5.zip`, 290 jobs, 20,000 shots/circuit, binary counts, angle, SX count and interleaved TEST rows; [paper 2112.07567](https://arxiv.org/abs/2112.07567) and deposited circuit code read | Downloaded, hash-pinned, parsed all rows. Same hardware data contain both n=1 and n=5. No separately certified actual pulse-axis angle is supplied. |
| B: multipath shutter tables + contrast/metrology | [Zenodo 5497283](https://zenodo.org/records/5497283): two full eight-mask photodiode streams plus four-mask contrast sweep | All bytes audited. Contrast is BC, not an independently known A phase reversal. Paper discusses detector beam-combination calibration but raw curve/gain are not in these three files. |
| C: photon event statistics + detector tomography; then two detectors | [Zenodo 18483031](https://zenodo.org/records/18483031): 8-bin TMD (14 probe intensities, ~100k triggers each), PNR SNSPD (122 probability rows, event records for first 26, ~570k triggers each), supplied POVMs and nominal mean photon numbers | Downloaded both archives, inspected arrays/readmes and reconciled every available event file. Zero events retained. Different instruments/dates have separate responses. POVMs are derived from assumed coherent-probe statistics, not independent Born-free calibration. |
| D: measurement-device-agnostic tomography + processor phase sweep | [Zenodo 17302724](https://zenodo.org/records/17302724), [paper 2407.13011](https://arxiv.org/abs/2407.13011), public MIT archive `RobStarek/mda-tomo-v1.1.zip` (23,208,558 bytes) | Record/file availability verified via API; raw contents not downloaded or claimed audited. Potential methodological input, not calibration transferable to Lima. A common characterized transformation would be needed. |

Additional discovery lead: [megascale detector tomography 10810148](https://zenodo.org/records/10810148)
resolves through its API to version 10810149, with ~995 MB matrices and separate
raw time tags. The smaller 18483031 release already supplies event-level checks,
so the large release was not downloaded. No fitted quantity from it is used.

Queries: `site:zenodo.org quantum detector tomography photon number data calibration`,
`site:zenodo.org Born rule quantum tomography dataset`,
`site:zenodo.org photon counting statistics beam splitter raw data detector`,
and title-specific primary-paper searches. Archive descriptions were followed
by API metadata, readmes and byte-level inspection for A–C. Original .py files
were read, never executed. No cloud jobs or author messages were submitted.

## A — Controlled quantum gates: main implemented candidate

**Required sentence.** Dataset A (n=1 angle sweeps) cannot distinguish a
single-qubit probability deformation from an ordinary control-axis distortion
because both yield the same binary probabilities at every recorded angle.
Dataset B (n=5 and other integer gate depths) supplies repeated transformations;
under an independently certified unwarped axis law it can check depth-dependent
errors, but with a common allowed axis warp it **cannot** distinguish the two.
A second known orthogonal preparation would supply a different projection and
can distinguish this specific pair under certified preparation/measurement
consistency. That preparation is not present in the audited release.

### Explicit operational alternative and equivalence proof

For a pure qubit ray r and binary projective measurement direction m, let
p=(1+m·r)/2. Define

f_theta(p)=p+theta p(1−p)(2p−1), −1≤theta≤1.

Both outcomes use f_theta(p), f_theta(1−p); these sum to one, are in [0,1],
fix 0,1 and 1/2, and theta=0 recovers ordinary probabilities. For ensembles,
average **after** applying the pure-ray function; mixed operational states are
equivalence classes of ensembles under these effects. Do not silently identify
them with ordinary density matrices. Rotations act on every ensemble component.
There are only binary projective outcomes and their trivial coarse graining;
no arbitrary-POVM or composite-system extension is claimed. Classical assignment
errors act after this rule. θ is universal only within this proposed single-system
family, not assumed equal to the unrelated optical intensity diagnostic.

Set h_theta(x)=2f_theta((1+x)/2)−1=x+theta x(1−x²)/2.
Its derivative 1+theta(1−3x²)/2 is nonnegative on [-1,1] in the stated range,
and endpoints are ±1, proving positivity. The ideal actual circuit
SX^n Rz(phi+π) SX|0> gives
p_n(phi)=[1+s_n cos(phi)]/2, s_n=sin(nπ/2) in {0,1,−1}.
An independent complex-matrix test checks the circuit convention.

Because h is odd and fixes zero, for **every integer n**,

f_theta(p_n(phi)) = [1+s_n h_theta(cos(phi))]/2
                 = p_n(g_theta(phi)),

g_theta(phi)=signed arccos[h_theta(cos(phi))].

Thus one common one-parameter axis-distortion function g_theta mimics the
alternative for all angles and all gate depths. This is not an arbitrary
per-setting fit. The equality also survives any common classical assignment
channel q -> a+(1−a−b)q, and any number of independent replications. More shots,
more integer depths, or removing either depth changes precision but cannot
break exact equivalence. This result is conditional on ideal quarter-turn
kinematics plus this stated ordinary axis-warp family; it neither explains every
hardware residual nor proves universal untestability of Born's rule.

The deposited n=1/n=5 hardware counts are audited and summarized via a descriptive
first-harmonic least-squares fit and job-level standard errors. Those fits are
not a likelihood ratio for theta. Depth differences can reflect ordinary gate
errors; no theta estimate is inferred from them. `lima1-5s.zip` has backend
`qasm_simulator` and is kept as simulation. The `lagos_bench` CSVs omit a gate-depth
column and the description says the sequence changed after job 50; depth is not
invented from row order. It also concerns a different backend and cannot serve
as Lima calibration.

### Complementary preparation and removal analysis

A proposed calibrated ±X and ±Y preparation set, through the same unknown
quantum channel and binary measurement, gives differences Cx and Cy. Pulling the
measurement back through the channel yields a qubit effect E; positivity of
E and I−E implies Cx²+Cy²≤1. This bound allows an arbitrary *common* quantum
channel/POVM, but requires the certified antipodal orthogonal preparations.
It does not follow from four unrelated states called X and Y by software.

For the alternative at a diagonal equatorial projection, Cx=Cy=h_theta(1/sqrt2)
and radius²=(1+theta/4)². Positive theta violates the unit disk ideally. The
ordinary common-axis-warp explanation stays on it. Negative theta lies inside;
this witness cannot identify all deformations. Classical contrast loss multiplies
both C values and can remove the violation.

The executable uses four raw binomial counts and simultaneous Clopper–Pearson
intervals (Bonferroni family alpha .01); it excludes only when the minimum squared
radius over the resulting difference rectangle exceeds the certified bound.
It simulates theta=.01 at the hardware-derived 17.4 million shots per preparation,
ordinary controls, common-axis warp, 99% contrast, the hardware harmonic contrast,
and negative theta. An enlarged radius 1.01 is a conservative uncertainty
sensitivity scenario, not a measured calibration. Removing the complementary
preparation restores exact equivalence; removing the orthogonality/consistency
certificate removes the unit-disk inference. These simulations are prospective,
not public measurements of the new preparation set.

## B — Shutter data plus calibration: implemented runner-up

**Required sentence.** Dataset A (eight shutter means) cannot distinguish an
additive all-open intensity term from a conditional phase coupled to unmeasured
imaginary coherence because both affect exactly the same observable coordinate.
Dataset B (the existing BC contrast scan) supplies some coherence information but
not the required A-phase calibration; under a new independently certified A
antipodal phase exposure the combination distinguishes the constructed pair.

The observation map, nuisance domains, exact pairs, sharp regions and conservative
sample costs are fully implemented in [derivation](derivation.md). The two
existing temperatures share only a diagnostic theta, retaining distinct coherence
and phase. Raw-intensity simulations test ordinary phase drift, leakage,
nonlinear response and injections. Removing the antipodal observation restores
ambiguity; omitting calibration produces false positives. Independent detector
metrology would constrain response only for the same apparatus, range and time.
The five-path paper's calibration cannot be assigned to this detector.

## C — Photon statistics plus detector tomography: audited rejection

**Required sentence.** Dataset A (photon event histograms) cannot in general
distinguish a change in prepared photon statistics from detector response because
its probabilities factor as Y=P R. Dataset B (a response matrix reconstructed
from the same assumed coherent probes) supplies an estimate of R **conditional
on P**, not independent information to distinguish the explanations. Under
independently characterized probe statistics and independently calibrated response,
new test-state data could separate them. That independence is absent here.

A precise global surviving gauge is (P T)R_d=P(T R_d) for any stochastic latent
transformation T and independent stochastic detector responses R_d. It persists
for both detectors and after either detector is removed. `compute` verifies an
explicit finite example. A model whose “probability modification” is only such
T is an ordinary preparation/measurement change, not a genuine modified Born
rule; this route was rejected, not promoted into a physics result. More general
normalized nonlinear probability maps need a different bridge and independent
calibration. No common response or mean-photon calibration is shared between
instruments or with the waveguide experiment.

The event audit finds all TMD empirical histograms exactly reproduce P_8bin;
PNR classified events differ from the supplied outcome matrix by up to .006904
in a cell. The README points to arrival-time classification and EMG processing,
so these products must not be silently treated as identical or independent.
The TMD supplied response also has row-normalization discrepancy up to .002564.
No clipping, renormalization or fitted probability-law parameter is used to hide
these discrepancies. They require processing clarification before a calibrated
joint fit, not a claim of source error or a quantum anomaly.

## D — Trusted transformations and device-agnostic tomography

**Required sentence.** An uncalibrated processor probability table cannot
distinguish the proposed probability function from control/measurement changes
because the actual prepared rays and effects are not independently fixed.
A device-agnostic tomography dataset could supply consistency constraints under
its own trusted transformations; under demonstrably shared characterized
transformations on the target device, those constraints could break some gauge
freedom. Data from another device do not certify Lima's transformations. The
listed release is a method/data lead, not a ready physical joint inference.

## Value and novelty assessment

The main exact theorem is an observation-map specialization of known control/
measurement identifiability issues, not a claim to invent gauge freedom. The
quarter-turn circuit structure and explicit normalized deformation make the
negative result stronger than a statement that “a flexible fit could explain
anything”: a **single** ordinary control function hides the parameter across
all integer depths. Independent matrix tests and public count/metadata audits
make its scope reviewable. Literature review has not established novelty as a
physics theorem; do not market the result as a discovery or PhD-level advance.

Best next measurement: independently certified complementary ±X/±Y preparations
with a quantified common channel/readout relation. Shot count alone is not the
bottleneck: contrast/calibration can erase the ideal separation. For the optical
runner-up, the concrete missing observable is the antipodal all-open A phase.
No tested public combination currently excludes a genuine probability-law
modification or the relevant complete ordinary model class.
