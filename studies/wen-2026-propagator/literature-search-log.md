# Staged literature and source search log

Date: 2026-09-26 UTC. This is a focused, staged
review, not a claim of exhaustive coverage of quantum foundations. Primary
papers and deposited data support technical statements; news, social discussion
and search snippets are discovery aids only. The earlier source review is
retained here with acquisition blockers updated after the verified data arrived.

## Stage 1 — paper identity, versions and acquisition

Questions: Which paper/version is being tested? Are original data and the
supplement available?

Queries included `"aeh1011" arxiv supplement`,
`"Direct experimental test of Feynman" "supplementary"`,
`"x0k6djj14" dataset version` and the exact title. Some queries returned unrelated generic code/version pages; those results
were excluded from the source review.

Primary sources read: [main paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC13510607/),
[Europe PMC XML](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13510607/fullTextXML),
[supplement API](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13510607/supplementaryFiles),
[Dryad landing page](https://datadryad.org/dataset/doi:10.5061/dryad.x0k6djj14),
[version metadata](https://datadryad.org/api/v2/versions/450489) and
[file list](https://datadryad.org/api/v2/versions/450489/files).
Full text, MathML and supplement/captions were reviewed in the earlier audit;
the current replication reread the comparison definitions and supplement
Sections 6–8. Source hashes are in provenance, originals outside git.

Historical access limitation: direct Dryad downloads returned 401/403 and some
publisher/PMC retrievals were blocked. The deposit was subsequently obtained;
all 19 originals match the pinned
manifest. Q01 is now closed. A refreshed PMC open still hit a browser challenge;
the archived primary XML remains available for equation checking. Current Dryad
landing content lists the same 19 file names. No claim that inaccessible sources
were freshly read is made.

## Stage 2 — measurement derivation and experimental lineage

Queries included `"Measuring the quantum propagator" "least action" photon`,
`"Demonstration of the quantum principle of least action with single photons"`,
`"Weak values from path integrals" Matzkin`, and the Lundeen direct-measurement
paper title. Primary papers reviewed:

- [Wen et al. 2023](https://arxiv.org/html/2305.19815v1), measurement derivation
  and paraxial mapping; refreshed during this replication.
- [Lundeen et al. 2011](https://doi.org/10.1038/nature10120), author manuscript
  reconstruction and ensemble-measurement description.
- [Matzkin 2020](https://arxiv.org/abs/2002.00832), v2 main derivation/discussion;
  refreshed version record. Its supplements have not been comprehensively audited.
- Wen 2026 main measurement equations and supplement Section 6, mapping pointer
  intensity differences and reference factors to K.

Result: preserve acquisition → calibrated K → shared-factor path products →
endpoint comparison as distinct steps. A separately acquired E is informative;
reconstruction alone does not make the comparison circular. Calibration and
conditioning remain physical assumptions. Q03/Q08/Q17 and physical-bridge.md.

## Stage 3 — deposit sufficiency, actual sheets and numerical ambiguity

Queries included `"aeh1011" "data" "code"`,
`"Wen" "propagator" "0.0518"`, the dataset identifier targeted to Dryad,
GitHub and Zenodo, and the exact paper title plus `code`.
No usable additional raw-tensor/code release was identified. The main paper's
availability statement points to the paper, supplement and Dryad; the inspected
deposit contains no executable analysis implementation. This is not proof that
the authors have no such code or that no other release exists.

Local primary-data review: every sheet, header, nonempty region and formula/hidden
state was inspected. The five image sheets have formatting dimensions larger
than their 321×101 numerical image regions. FigS4A contains a separate fixed-endpoint
x column. FigS2 histograms require care over edge counts and a nonmonotone tail.
These are documented in data-dictionary.md; originals are unchanged.

The first Q–theory calculation differed from 4.45% under our reading of the
printed comparison, triggering an additional
source check and numerical checkpoint. Queries `"aeh1011" "4.45"` and
`"aeh1011" "correction"` returned no useful clarification (one result was an
unrelated catalogue). We therefore reread the actual Results definition and
the author README rather than relying on search. Checked linear/nearest/local
cubic alignment, rounding sensitivity, unit-sum normalization, alternative
common error definitions and a conditional finite-grid model. All choices and
their limits are retained in the code and method-comparison.md.

Review feedback prompted a fresh retrieval of the primary XML; its SHA-256
matches the archived source. The exact sentence at MathML m116–m118 names Q
and theoretical Q and is quoted in method-comparison.md. Unit-sum E–Q gives
4.49%, close to the reported 4.45%; that alternative is now adjacent to
Q–theory in the comparison table. Its closeness is not evidence of the authors'
actual array choice. Both comparison identity and processing remain clarification
questions until the implementation is available.

Result: E–Q and Q–theory remain separate; E–Q versus E–C ordering survives the
tested scale convention. Exact printed metric reproduction remains open. No
subsets or fitted offsets were selected to reach a desired percentage.
Q04/Q06/Q10/Q11/Q19–Q21 updated. The per-repeat tensors and processing code are
specific remaining dependencies, not records that can be inferred from means.

## Stage 4 — alternatives, physical correspondence and ontological scope

Queries included the Feynman 1948 title, `"Anomalous Weak Values Are Proofs of
Contextuality"`, `"Anomalous weak values and contextuality" "operational
constraints"`, `"propagator" "tomography" "finite aperture" photons`, and the
looped-trajectory paper title. Primary review:

- [Feynman 1948 primary record](https://authors.library.caltech.edu/records/9h858-5hv71):
  abstract and bibliographic scope, refreshed. Finite matrix/path identity is
  independently derived, rather than claimed as a novel empirical finding.
- [Magaña-Loaiza et al. 2016](https://doi.org/10.1038/ncomms13987): main
  theory/results/methods in primary XML, previously reviewed. Changing apertures
  changes boundary conditions; spatial loops are not backward-time evidence.
- [Pusey 2014 v2](https://arxiv.org/abs/1409.1535v2): primary abstract/theorem
  scope refreshed; no new full-proof certification in this empirical PR.
- [Kunjwal–Lostaglio–Pusey 2019 v2](https://arxiv.org/html/1812.06940v2): earlier
  review of Sections II.3–II.6, III.2 and IV; current primary scope/version refresh.
  Finite-pointer Theorem 3 and its operational conditions identify a possible
  follow-up, not a witness already evaluated from this deposit.

Result: the finite formal companion and current data analysis cannot discriminate
equivalent representations. A stronger contextuality claim needs additional
operational tests. Phase intervention remains a separate study with its own
full protocol, literature and calibration obligations; no novelty claim here.

## Stage 5 — corrections and adverse evidence

Current queries: `"aeh1011" correction code`, the full title plus `erratum` or
`code`, and `"Wen" "Tian" "Feynman" correction` targeted to publisher/arXiv
sources. Search results included primary paper/deposit records and unrelated
pages; no usable correction or additional code was identified.
The [fresh Crossref record](https://api.crossref.org/works/10.1126/sciadv.aeh1011)
has an empty `relation` and no `update-to`/`updated-by` fields. Its retrieval
hash is recorded in provenance. This supports only “none identified in this
search,” not “no correction exists.”

Adverse checks retained: the broad E–Q/E–C ordering survives normalization;
the spatial-range apparent mismatch has a plausible bin-edge explanation;
rounding/interpolation do not explain the quoted error difference; the missing
repeat/covariance records prevent significance/fidelity reconstruction; a
context-dependent matching response is not itself a physical trajectory theory.
No claims about superluminal or backward-time motion follow from these tables.

## Review boundaries

This review reflects the sources and versions above. New data, code, corrections
or clarification could change the mappings, model or normalization and would
require reassessment of dependent conclusions. Empty searches do not establish
absence. Full experiment design, apparatus-specific literature review and
full-proof novelty review remain separate follow-up work.
