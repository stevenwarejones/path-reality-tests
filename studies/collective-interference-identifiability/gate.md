# Original feasibility decision: narrow to a conditional three-atom certificate

This historical gate is retained; [revision.md](revision.md) records the
completed 2D and full-pattern checks and the still-open apparatus gate.

Question: do singleton propagation records and three-atom parity images exclude
mixtures of orthogonally labelled groups of size at most two?

Decision: **proceed with a conditional count-level certificate; do not claim the
broad apparatus-robust research goal is solved.** The useful subset is the
three-atom entry in `1D.nc` and its singleton calibration in `1D_singles.nc`,
both already bundled in `boson-1.0.4.zip` (59 MB compressed). No 3.6 GB download
is needed. The selected statistic is the original paper's full-row bunching
observable, evaluated per prepared shot instead of conditional on survival.

The preliminary calculation found 93 successes / 2,999 prepared shots. Three
singleton arrays each contain 931 shots, but are exact translated copies of
one record. Count this as one calibration family, not three independent
sources. The original notebook uses file `3/29` for all three and rolls two
copies in y. The confidence construction must respect that dependency.

| Measured family | Exact records and exposure | Intervention / date | Role / unresolved sharing |
|---|---|---|---|
| Singleton propagation | `1D_singles.nc`, keys 30, 31, 32; only 931 unique shots; initial and final binary images; 12 by 12 reference crop | Prepare one atom; notebook identifies selected run 220702/29 (July 2, 2022); per-shot timestamps absent in these NC files | Estimates one detected transition distribution. Translation to other input sites is assumed, not independently calibrated by copied arrays. |
| Three-atom interference | `1D.nc`, key 3; 2,999 initial-selected shots, 2,329 with three final occupied sites, 93 with all three in one y row | Prepare adjacent three-atom input; nominal 2.46 ms in notebook (`3NN.nc` labels duplicate subset 2.45 ms) | Tests the cluster class; final parity image must be forward-modelled. Same apparatus, different preparation. |
| Pairwise records (not used in certificate) | `2NN_2205.nc` and singleton counterparts; other settings and dates | Two-atom HOM-like evolution | Possible apparatus cross-check, not evidence for every pair in the selected three-atom preparation. No indispensable third source claimed. |

Physical bridge: fixed, noninteracting, hidden-state-independent linear
propagation and calibrated loss; exactly three prepared particles; parity
readout with no false positives. For any mixture of partitions into groups of
size at most two, permutation moments vanish outside a subgroup of size at
most two. Cauchy–Schwarz gives `P(B) <= 2 P_D(B)`. Dropping the distinct-x
requirement only enlarges the bound: `P_D(B) <= sum_y q1(y) q2(y) q3(y)`.
The translated calibration reduces this polynomial to consecutive triples of
one distribution. The cropped-out mass remains a nuisance, not a structural
zero. No interferometer phase estimate or independent-pure-state Gram model
is needed.

Individually, singletons cannot distinguish fully distinguishable from
indistinguishable internal inputs; the bunching event alone cannot distinguish
interference from stronger classical focusing. Jointly, calibration fixes a
ceiling on classical focusing. Preliminary certified-relaxation calculations
place the group-size-two ceiling just below a simultaneous lower confidence
bound on the measured event. The final report must give the resulting small
systematic-error budget. A coarse-event source-removal example is not proof of
insufficiency of every high-order count constraint.

Alternatives audited before proceeding:

* Zenodo 10829208: downloaded all four 5–8-photon CSVs, HOM histogram and code.
  These are collisionless selected outcomes, without the full rejection/trial
  record or an independently measured calibration for every photon pair. Do
  not fit a universal overlap and call it a full-class certificate.
* Zenodo 10927445: API access succeeded; downloaded the 1.33 MB ZIP. It contains
  timestamped lower-/higher-fold coincidences, normalized duplicates and
  tomography material. Entangled inputs require a different physical null.
  Not selected; no raw-record inference claimed for this archive.
* Brod et al. 2019 and Giordani et al. 2020 already witness genuine
  multiphoton indistinguishability. Menssen et al. already establishes the
  triad-phase ambiguity. Geller–Knill 2026 and Pioge–Novo–Cerf 2026 limit
  generalized-bunching interpretations. The proposed increment is a rigorous
  finite-count, crop-aware reanalysis of this atomic subset. Priority and
  consequential scientific novelty remain unestablished.

Do not expand into an expensive many-parameter fit until independently
measured propagation transfer, interactions, preparation errors and false
positives fit within the reported error budget. No absence of those controls
constitutes an impossibility theorem.
