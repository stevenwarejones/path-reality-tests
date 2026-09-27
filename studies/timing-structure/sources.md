# Primary sources, related work and novelty boundary

Accessed 2026-09-26. This is a bounded related-work check, not an exhaustive
priority search. The mathematical ingredients below are established methods.
The contribution claimed here is a reproducible application and comparison on
this pinned archive, with explicit counterexamples and assumption accounting.
No claim of a new Bell theorem, new statistical theorem, new physics effect or
previously unstudied general use of time tags is made.

| Primary source | Inspected content and relevance | Boundary of this study |
|---|---|---|
| [NIST Bell research data](https://www.nist.gov/pml/applied-physics-division/bell-test-research-software-and-data) and the [existing provenance ledger](../nist-bell-causal-audit/sources.md) | Official archive and the audited run4 decoder, channels, overlap and pinned hashes | Same existing run; no new data collection; the catalog does not certify our new conditional assumptions |
| [Shalm et al., PRL 115, 250402 (2015)](https://doi.org/10.1103/PhysRevLett.115.250402), [author manuscript](https://pmc.ncbi.nlm.nih.gov/articles/PMC5815856/) | Apparatus, random-setting generation, memory-aware Bell-test context and excess predictability | The original Bell conclusion and its calibration are not invalidated or automatically transported into our new timing estimands |
| [Bednorz, PRA 95, 042118 (2017)](https://doi.org/10.1103/PhysRevA.95.042118), [preprint](https://arxiv.org/abs/1511.03509) | Setting dependence and assumptions in recent Bell experiments; the baseline ledger gives the detailed NIST window comparison | No claim that marginal timing-dependence searches are new; we do not reproduce his different run/window comparisons |
| [Knill et al., PRA 91, 032105 (2015)](https://doi.org/10.1103/PhysRevA.91.032105), [NIST summary](https://www.nist.gov/publications/bell-inequalities-continuously-emitting-sources) | Bell analysis on time-tag sequences with fixed observation periods and explicit random-setting assumptions | We use local marginal categories, not their distance-based joint Bell witness |
| [Christensen et al., PRA 92, 032130 (2015)](https://doi.org/10.1103/PhysRevA.92.032130), [preprint abstract](https://arxiv.org/abs/1503.07573) | Three experimental examples where timing/coincidence analysis choices change the inference | Supports treating timing analysis as substantive established territory; no implementation or replication of their distance witness is claimed |
| [Howard et al., Probability Surveys 17, 257–317 (2020)](https://doi.org/10.1214/18-PS321), [preprint](https://arxiv.org/abs/1808.03204) | Nonnegative-supermartingale framework for time-uniform concentration | We derive elementary Bernoulli finite-grid bounds and a fixed-stake mixture; these are not novel or claimed optimal |
| [Lipsitch, Tchetgen Tchetgen and Cohen, Epidemiology 21 (2010)](https://pubmed.ncbi.nlm.nih.gov/20335814/), [author repository](https://dash.harvard.edu/entities/publication/73120379-21f6-6bd4-e053-0100007fdf3b) | Negative controls and the assumptions needed for their interpretation | A conceptual analogue from another field; a future-setting diagnostic is not automatically a valid causal negative control |

The archive's original papers and the elementary probability derivations are the
sources of technical claims. Search-result summaries were used to locate papers;
no secondary review is used as a substitute for the statistical proof. Some
publisher/author PDF endpoints were unavailable; the accessible primary abstracts,
author manuscript and the existing detailed source audit delimit what was checked.

## What is established versus archive-specific

- **Established:** stratification can uncover cancellation; conditioning on a
  collider can manufacture an association; memory can create lag patterns;
  deterministic coarse-graining loses information; randomization/exogeneity
  assumptions distinguish association from intervention.
- **Implemented here:** two strictly prior-record covariates, fixed lag and
  feature families, full-source reconciliation, interval-censored covariate
  states, a pathwise missingness bound for future-setting bets, simultaneous
  weighted-TV cancellation bounds, actual occupancy stress tests, and exact
  executable witnesses tied to this archive's representation.
- **Not established:** uniqueness in the literature, new fundamental physics,
  physical sub-bin timing sensitivity, a full point-process limit, certified
  RNG freshness, or an independent-run replication.

Primary data notice and attribution remain in [the NIST notice](../nist-bell-causal-audit/NOTICE.txt).
Only derived aggregates are committed; raw archives and per-row caches stay external.
