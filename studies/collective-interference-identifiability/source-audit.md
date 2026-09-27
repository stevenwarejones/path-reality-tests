# Sources, dependencies and prior art

## Selected atomic archive

Source: Young, Geller, Eckner, Schine, Glancy, Knill and Kaufman,
[An atomic boson sampler](https://arxiv.org/html/2307.06936v2),
Nature 629, 311–316 (2024),
[DOI](https://doi.org/10.1038/s41586-024-07304-4).
Data and code: [Zenodo 10453016](https://zenodo.org/records/10453016).
Originals are hash-pinned in [sources.json](sources.json), remain external,
and were downloaded and inspected. The archive includes a NIST license notice;
our deposited JSONs are identified as derived aggregates or constructed
witnesses. The NC originals and notebook are not redistributed here.

The 59,176,561-byte `boson-1.0.4.zip` already contains the relevant NC originals,
original analysis code, processed bootstrap arrays and exports. The selected
files expand to 252,653,468 bytes (`1D.nc`) and 376,427,616 bytes
(`1D_singles.nc`). Raw original image metadata and preparation filtering are
not reconstructed from the larger separate `data.zip`. We downloaded only its
last 1 MB to inspect accessibility, not the full archive; none of that partial
file enters the inference.

The original README describes the indices as initial/final image, shot, y, x;
the selected arrays additionally have a preparation/particle-number key.
Arrays contain binary pixels and NaN padding/crops. Our parser distinguishes
padded shots from valid empty final images, checks rectangular masks, initial
preparation, binary values, no recorded excess output atoms, and exact copied
arrays. It does not assume that an empty parity pixel means a physically empty
site. Preparation selection occurred upstream; the statistical population is
that initial-selected population, not all experimental attempts.

### Crucial dependencies

* In the notebook's “N particle data” section the n=3 calibration keys 30,31,32
  all refer to the same `3/29` file. The source arrays are exactly equal after
  rolling one copy by +1 and another by -1 along y, including initial images.
  This is **one record and a transfer assumption**.
* The reference crop covers array rows 9–20; the many-body crop covers 8–19.
  Unobserved reference rows 7 and 8 can affect the latter after translation.
  Those probabilities are free variables inside the measured empty/outside
  mass. Treating all NaNs as physical zero probabilities would be unsound.
* The NC arrays retain shot ordering, but lack shot timestamps and run IDs.
  The run-map audit dates the selected experiment July 2, 2022 and specifies
  nominal n=3 evolution time 2.46 ms. `3NN.nc` lists the corresponding
  2,999-shot setting at 2.45 ms. We do not use that second packaging as a new
  independent source or infer precision from the discrepancy.
* Published Figure 2's full-row bunching is normalized by full survival.
  Our counts reproduce that selected observable as 93/2329, approximately
  0.03993, and use 93/2999 for inference. The original distinguishable
  prediction uses resampling, parity and a survival normalization; it is
  not an independent measured calibration.
* The original `inference/fullbunch.py` evaluates distinguishable 2D
  collision removal by Monte Carlo and bootstrap uncertainty. The preserved original
  row bound avoids numerical sampling. The revision certifies distinct-x 2D
  bunching directly, also without Monte Carlo.
  We do not claim byte-for-byte reproduction of the original bootstrap runs.

### Other measurement combinations actually audited

| Candidate | Access and exact files | Decision |
|---|---|---|
| Programmable multiphoton interference, [Zenodo 10829208](https://zenodo.org/records/10829208) / [paper](https://arxiv.org/abs/2305.11157) | Downloaded all four 5–8-photon CSV count distributions, HOM histogram, `Experiment.py`, `Loop.py`, notebook; CC BY 4.0 | Counts are processed collisionless events. Raw rejected events, complete trial denominators and an independently measured calibration for all pulse pairs are not established by this release. No joint certification attempted. |
| Entanglement-induced interference, [Zenodo 10927445](https://zenodo.org/records/10927445) / [paper](https://arxiv.org/html/2310.08630v2) | API accessible despite landing-page error; downloaded `ECMBIDataSet.zip` (1,332,377 bytes), 444 files: 432 coincidence files and 12 tomography files; CC BY 4.0 | Timestamped filenames, angles, singles, coincidence channels and `t_meas`/`t_actual` columns offer useful controls. Normalized `N` files are derived duplicates, not independent evidence. Internally entangled source requires a separate null and complete schema audit. Not abandoned as impossible. |
| Atomic singleton + HOM + higher-particle subsets | Bundled NC files and original code inspected | Best immediate access to trial-level parity data; selected singleton/higher-order core. HOM files do not independently certify transfer for all three selected inputs. |

The alternatives' files and checksums are pinned for reproducible discovery,
but their measured numbers are not used in the atomic certificate. No
cross-species shared parameter is invented.

## Closest primary work and scope of the increment

| Work | Already established / relationship to this study |
|---|---|
| [Shchesnovich, PRA 91, 013844 (2015)](https://arxiv.org/html/1410.1506v7) | PSD partial-distinguishability matrices and physical output-probability framework. Our moments and positivity are applications, not a new formalism. |
| [Menssen et al., PRL 118, 153603 (2017)](https://arxiv.org/abs/1609.09804) | Pairwise overlaps do not determine collective interference; triad phases are established prior art. No rediscovery claim. |
| [Shchesnovich and Bezerra (2017)](https://arxiv.org/abs/1707.03893) | Collective permutation-cycle phases and their absence from lower-order marginals. We do not equate a cluster decomposition with every definition of interference order. |
| [Brod et al., PRL 122, 063602 (2019)](https://arxiv.org/html/1804.01334v1) | Operational witnesses against mixtures of labelled groups, bounds on genuine indistinguishability, and three-photon experiments. This is especially close prior art. Our claim is a finite-count certificate for different measured atomic records with crop/loss uncertainty, not the invention of cluster witnesses. We do not assume all internal states decompose into perfectly identical/distinguishable groups. |
| [Giordani et al.](https://arxiv.org/abs/1907.01325) | Four-photon indistinguishability from measured pair overlaps. A related low-order route, not evidence of novelty for combining calibrations. |
| [van der Meer et al.](https://arxiv.org/abs/2112.00067) | Semi-device-independent photonic indistinguishability witness. Our strong stationarity/transfer assumptions must not be described as device independent. |
| Young et al. (2024), paper and bundled code above | Already reports full bunching and compares distinguishable predictions with calibration and survival corrections. Our count-level null is broader than fully distinguishable particles and has an exact confidence/optimization check, but the observable and source pairing are already used there. |
| [Geller and Knill (2026)](https://arxiv.org/html/2509.04550v4) | Visible-state and randomized generalized-bunching theory; arbitrary unrandomized bunching is not a universal monotone of indistinguishability. Our bound uses coset Cauchy–Schwarz against a specific physical class, not such a monotonicity assumption. |
| [Pioge, Novo and Cerf (2026)](https://arxiv.org/abs/2601.13792) | Time-delay-induced anomalous multimode bunching further limits naive validation criteria. Does not invalidate a direct upper bound containing all moments of the declared class. |

This is a targeted primary-literature comparison, not a completed priority
survey. A consequential new experimental conclusion would require a robust
transfer calibration, stronger count-level analysis, and comparison with
existing witnesses on the same records. The small conditional certificate
alone does not establish a PhD-level advance.

## Revision: matched-control sharing

The hash-checked notebook run maps and all 17 settings in `3NN.nc` and
`3NN_0.nc`–`3NN_2.nc` are reconstructed by `control_audit.py`. See the
[explicit sharing table and limitations](revision.md#matched-controls).
The selected 3NN and 1D images agree exactly after aligning physical x/y
coordinates; their different array dimensions do not make independent data.
The old summary quantum witness has zero probability for 68 observed patterns,
each seen once. It fails the complete parity region. The new smoothed rational
quantum witness passes every constraint of that region; the old artifact is
retained only for its original summary and singleton-source claims.
