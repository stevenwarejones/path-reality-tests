# Literature coverage and decisions

The [evidence table](evidence-table.md) records primary identifiers, versions,
methods coverage, measured quantities and access gaps. Search results were used
for discovery; scientific comparisons use the primary records listed there.
Queries were broad across years, followed by exact-title/version checks and a
365-day filter for the last current-literature pass. No exhaustive citation index
export or systematic-review completeness claim is made.

## Search record

All entries below were performed on **2026-09-26 UTC**, at successive checkpoints.
Sites covered through public search and direct retrieval included arXiv, APS,
Nature/PMC, MDPI/Preprints, J-STAGE, institutional repositories and official
SciPy documentation. General web/citation-index results were screened for primary
versions. Quoted strings below are actual submitted queries; search services
sometimes expanded them and returned irrelevant material, which was discarded.

| Checkpoint | Exact queries | Decision and affected questions |
|---|---|---|
| 1: protocol comparison | `Fermi two atom problem detector excitation signaling commutator`; `single photon precursor information velocity Zhang 2011`; `experimental no signaling random settings time tagged detection`; `retrocausal model operational no signaling future interventions measurement independence`; `"Optical Precursor of a Single Photon"`; `"Experimentally separating vacuum fluctuations from source radiation"`; `Causal Modeling Delayed Choice Experiment 1710.07323` | Separate correlations, heralded pulse fronts, Bell trials and intervention marginals. Select separated modulation for its later in-cone control; retain earlier records with explicit independence. ST01, ST03. |
| 2: channel and geometry | `"total variation" "causal influence" randomized intervention`; `"finite bandwidth" "front" causality pulse`; `spacelike separated finite duration detector switching signaling rotating wave approximation`; `optical precursor finite bandwidth pulse front causality Sommerfeld` | Keep a modest finite channel theorem plus a conditional quantum identity. Do not substitute nominal centers, group velocity or lattice/rotating-wave approximations for complete supports. ST02, ST04, ST06. |
| 3: inference/resources | `"no-signaling" martingale statistical test memory`; `"Clopper" "Pearson" exact confidence difference proportions`; `no signaling statistical test martingale device memory randomization`; `"no-signaling" "statistical" "memory" test`; `Hoeffding probability inequalities sums bounded random variables 1963 theorem 2`; `"Azuma" "1967" weighted sums martingale`; `"loophole" "Shalm" 2015 random synchronization 184` | Derive a fixed-horizon conditional score and retain IID binomial intervals only as an explicitly stronger comparison. Apparatus timing examples establish precedent, not this apparatus's capability. ST05, ST07, ST10. |
| 4: adverse cases after initial calculations | `coincidence time loophole Larsson Gill 2004 Bell timing selection`; `"no signaling" "afterpulsing" detector`; `"pre-response" filtering causality pulse`; `"Albanese" "physics8020052" correction` | Keep all clock-defined trials; demonstrate selection, predictable settings, centered-filter pre-response and sign cancellation. Synthetic samples are sensitivity checks, not standard-physics anomalies. ST03, ST10, ST12. |
| 5: closest work and final claim | `"A Minimal Operational Criterion for No-Signaling Assessment"`; `"A proof that no-signalling implies microcausality" quantum field theory`; `"spacelike" "binary" "confidence" no-signaling 2026`; `"no signaling" "randomization" "martingale" photons`; `"Bounding the Plausibility of Physical Theories" arxiv`; `"spacelike" "randomized" "source" "no-signaling" experiment 2026` | Closest binary statistical proposal is already published. Claim formal verification and a conditional resource study, not a new criterion or physical effect. Liang–Zhang strengthens the adverse randomization comparison. ST08, ST11. |
| Final focused repeat and forward/backward checks | `"no-signaling" randomized binary confidence optical modulator spacelike 2026`; `"Experimentally separating vacuum fluctuations from source radiation" correction comment`; `"A Minimal Operational Criterion for No-Signaling Assessment" correction`; `"no-signaling" "randomized" "modulator" confidence`; `"Optical Precursor of a Single Photon" comment correction`; `"A Minimal Operational Criterion" correction`; `"Experimentally separating vacuum fluctuations" correction`; `de Ramon Papageorgiou Martin Martinez 2021 causal factorization detectors relativistic causality`; `de Ramon Papageorgiou Martin Martinez 2023 noncompact detectors effective lightcone signaling` | Exact-title searches revisited precursor/vacuum/source and binary-method work; source references led backward to causal factorization and forward to noncompact-tail bounds. Recent temporal-NSIT lead screened as a distinct question. No supported new-physics or priority claim emerged. |

The final detector-paper follow-up read the causal-factorization derivation and
Appendix A of arXiv:2102.03408v2, and the perturbative estimator/exponential bounds
and Appendix B of arXiv:2305.07756v2. These sharpened the bridge: arbitrary
operations on extended detectors cannot be assumed causal merely from field
microcausality, and exponential-tail bounds are model dependent. The study thus
requires an actual localized channel/coupling calibration before a laboratory claim.

Retrieval changes matter. Liang–Zhang's HTML/publisher routes failed, but its
arXiv PDF was recovered and read. Martín-Martínez and the two de Ramón papers
were likewise read as PDFs after HTML failures. The precursor institutional
file was identified as a poster rather than mistaken for full text. Herter's
main methods were read; supplement challenge and institutional 403 responses
were respected. Albanese's full preprint was read, with publisher excerpts
confirming the published binary/interval scope. Tong's notes could not be
retrieved. These remaining gaps are ST11, not evidence of novelty or permission
to infer unseen methods. No quantitative published sensitivity is borrowed from
a source with incomplete methods/supplement coverage.

## Review-driven checkpoint, 2026-09-26 UTC

Queries: `Stenner Gauthier Neifeld 2003 speed information fast light optical medium 695`;
`Stenner Gauthier Neifeld nature02016 pdf`;
`Adenier Khrennikov 2007 Weihs fair sampling no signalling Bednorz 2017 042118`;
`Bednorz "042118" 2017`; `"Bednorz" "042118" comment response signaling`;
`Salart 2008 testing speed spooky action distance Nature 454 861`;
`"Is the Fair Sampling Assumption" "2007" "131"`; `"Stenner" "nature02587"`.

Verified all four suggested precedents against primary records. Added them to
the evidence table with source-specific coverage. Stenner publisher full text
was subscription-only and the author-group PDF failed; the public university
course PDF succeeded and its experiment/methods were read. Adenier and Bednorz
PDF methods were accessible. Salart main text and complete comment/reply PDFs
were accessible. Stenner's comment/reply was also read. Searches do not certify
absence of further responses. No source PDFs, historical counts or code were
added to this repository.

Disposition: separate the closest experimental precedent from the closest
statistical proposal. Apparent-signaling reanalyses motivate predeclared complete
records and window/multiplicity controls; they do not supply an anomaly claim
here. Control-to-primary-run transport and the assumed B values are now explicit
in the protocol, README table and ST06.
