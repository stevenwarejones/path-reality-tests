# Literature supporting the phase-intervention design

This document preserves the earlier primary-source review and adds a focused
derivation/design search on 2026-09-26 UTC. It does not claim exhaustive coverage. The comparison of Wen's reported
measurements and dataset belongs to the separate Wen study.

| Source | Material previously reviewed | Relevance and limit |
|---|---|---|
| [Lundeen et al. (2011)](https://doi.org/10.1038/nature10120) | Author manuscript: introduction, reconstruction derivation and experimental description | Optical inference uses an ensemble and a measurement model; a reconstructed quantity is not automatically an individual trajectory record. |
| [Matzkin (2020)](https://doi.org/10.1103/PhysRevResearch.2.032048) | Author manuscript v2: main weak-probe derivation and discussion | Weak responses depend on amplitudes and postselection. This informs interpretation of a probe; it does not supply the proposed phase-intervention calibration. |
| [Magaña-Loaiza et al. (2016)](https://doi.org/10.1038/ncomms13987) | Primary full-text theory, results and methods | Changing apertures changes propagation boundary conditions. Blocking a region need not be interchangeable with dephasing it under a fixed optical map. |
| [Pusey (2014)](https://doi.org/10.1103/PhysRevLett.113.200401) | Primary abstract and theorem scope, supplemented by the later paper below | Motivates a stronger contextuality objective under additional operational assumptions. The original proof has not been fully audited here. |
| [Kunjwal–Lostaglio–Pusey (2019)](https://arxiv.org/html/1812.06940v2) | Sections II.3–II.6, III.2 and IV; displayed inequalities and equivalences | The finite-pointer theorem requires joint probabilities and measurement/transformation equivalences. A full proof/instrument audit and finite-data calibration design remain necessary before applying it. |

## Search record and review limits

Prior searches included `single photon path interference local phase shift weak
measurement phase shifter contextuality experiment`, `"single photon" "phase"
Grangier Roger Aspect 1986 interference experiment`, and `"Experimental
demonstration" "contextuality" "anomalous weak values" Piacentini`, across web searches. Direct primary-source retrieval supplied the reviewed
papers above. Irrelevant or secondary results were not used as technical evidence.
The Grangier and Piacentini search hits were leads, not full protocol reviews.

The phase mechanism is established interference physics. No literature search
completed so far establishes experimental novelty for this design.

## Required later search checkpoints

1. **Before apparatus selection:** review the closest regional-phase experiment,
   its supplement, calibration and loss treatment. Compare exact preparation,
   intervention and outcome definitions; record what this design adds.
2. **Before analysis lock:** search primary work on approximate operational
   equivalences, setting-dependent loss and finite-data tests applicable to the
   selected instrument. Audit the actual chosen theorem rather than applying a
   similar-looking inequality.
3. **After pilot data:** search for mechanisms matching observed drift, phase
   cross-talk, diffraction, contamination and detector response. Carry plausible
   alternatives into the nuisance model or document why they are excluded.
4. **Before final claims:** refresh source versions/corrections and check every
   claim against the data, calibration and proof it requires.

The review record distinguishes queries, primary-source inspection, adverse
evidence and consequences for PH01–PH06. Future apparatus and pilot conclusions
depend on the corresponding empirical evidence.

## Focused derivation and design review — 2026-09-26 UTC

Queries: `phase cycling four step interferometry density matrix coherence dephasing phase shift`;
`site:arxiv.org coherence witness interferometer visibility dephasing phase shift`;
`Kunjwal Lostaglio Pusey 2019 anomalous weak values robust contextuality operational equivalences`;
`Hoeffding 1963 probability inequalities sums bounded random variables theorem 2 pdf`.
Broad phase-cycling searches included NMR and generic web material; only the
primary sources below were used for the stated technical conclusions.

| Primary source and version | Material inspected | Consequence |
|---|---|---|
| [Biswas, García Díaz, Winter, arXiv:1701.05051v3 (2017)](https://arxiv.org/html/1701.05051v3), [journal DOI](https://doi.org/10.1098/rspa.2017.0170) | Introduction, phase-unitary/POVM formulation Eqs. 1–4, two-path visibility discussion | Closest general theoretical framing. A fixed detector's fringe probes coherence; ordinary fringe detection is established physics. |
| [Hacker et al., arXiv:2305.03641v2 (2023)](https://arxiv.org/html/2305.03641v2), [journal DOI](https://doi.org/10.1088/1367-2630/ad0752) | Introduction and setup; scope of feedback/count-rate/noise analysis | Supports a practical two-mode candidate. Its feedback protocol is not a certification of independent heralded trials here. |
| [Hoeffding (1963)](https://doi.org/10.1080/01621459.1963.10500830), [primary-paper copy](https://www.cs.rpi.edu/academics/courses/spring06/random/hoefding.pdf) | Theorem 2, printed p.16, viewed in the scanned paper | Bounded independent variables give the per-setting tail bound. The four-setting union bound and power condition are derived in this protocol. |
| [Kunjwal–Lostaglio–Pusey, arXiv:1812.06940v2](https://arxiv.org/html/1812.06940v2) | Theorem 3, Eqs.29–31; Sections III.2 and IV rechecked | The stronger witness needs operational equivalences and disturbance controls absent from a phase fringe. It remains a separate research objective. |

Decision: preserve the four-phase contrast test as a test of a defined dephased
class. Add complete loss outcomes, calibration allowance and conditional power.
Keep a general contextuality conclusion out of this study. The formal contribution
is a machine-checked version of known interference algebra and its explicit
model-class boundary, not a claim of new interference physics.
