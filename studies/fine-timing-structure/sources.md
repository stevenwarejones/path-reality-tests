# Provenance, prior work and novelty boundary

Primary-source check, 2026-09-27. This is a bounded comparison, not an exhaustive
priority search. The contribution is an executed comparison of retained versus
discarded timing information on one archive. No new theorem, optimal estimator,
new Bell witness or new physics is claimed.

| Primary source | Established idea and comparison |
|---|---|
| [NIST Bell research data](https://www.nist.gov/pml/applied-physics-division/bell-test-research-software-and-data) and [audited source ledger](../nist-bell-causal-audit/sources.md) | Same pinned HDF and raw Alice/Bob streams; exact hashes are in every aggregate and the original manifest. Attribution remains in the existing NIST notice. |
| [Knill et al. (2015), Bell inequalities for continuously emitting sources](https://www.nist.gov/publications/bell-inequalities-continuously-emitting-sources) | The primary abstract describes Bell inequalities on time-tag sequences with fixed observation periods and random independent settings. Our local marginal sub-CDFs and partitions do not reproduce that joint, distance-based Bell statistic; time-tag interrogation itself is established. |
| [Howard et al. (2020), Time-uniform Chernoff bounds via nonnegative supermartingales](https://arxiv.org/abs/1808.03204) | Exponential supermartingales and uniform concentration are established. We derive a conservative finite-grid Bernoulli inversion, inherited from #13, and do not claim its tuning is optimal. |
| [Mineiro and Howard (2023), Time-uniform confidence bands for the CDF under nonstationarity](https://arxiv.org/abs/2302.14248), [conference version](https://proceedings.neurips.cc/paper_files/paper/2023/hash/148bbc25b934211d80435b5cad5a7198-Abstract-Conference.html) | The primary abstract explicitly targets running averages of conditional distributions under arbitrary data dependence and includes an importance-weighted counterfactual extension. Thus neither a dependent-data CDF target nor counterfactual distribution bands are novel here. Our finite-grid union bound, propensity interval and missing-count enclosure are a simpler auditable construction, not a numerical superiority comparison with their algorithm. |
| [Dümbgen and Walther (2008), Multiscale inference about a density](https://arxiv.org/abs/0706.3968) | The primary abstract describes simultaneous multiscale statements using local order statistics/spacings and finite-sample significance control. Our fixed nested categorical partitions are a conservative resolution comparison, not an implementation of their density-shape procedure, a transfer of its assumptions, or a claim of multiscale optimality. |
| [Shalm et al. (2015), author manuscript](https://pmc.ncbi.nlm.nih.gov/articles/PMC5815856/) and [#13 literature ledger](../timing-structure/sources.md) | Apparatus, timing and setting assumptions remain distinct from a new archive analysis. This extension does not challenge or re-certify the original Bell conclusion. |

The accessible primary abstracts and the existing apparatus audit support this
scope comparison. The implementation's guarantee is derived in statistics.md,
not inferred from a literature citation. No stronger claim of literature novelty
or numerical superiority follows from this limited review.
