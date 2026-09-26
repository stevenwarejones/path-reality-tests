# Path reality tests

This repository tests when claims that quantum paths or trajectories are “real”
can be checked by experiment. Studies connect explicit models to observable
predictions and document the data, assumptions and reproducible calculations
needed to compare them. These are independent, reproducible checks, not claims
that any paper is wrong.

| Study | Question | Status |
|---|---|---|
| [Wen 2026 propagator](studies/wen-2026-propagator/README.md) | What can reconstructed propagators and endpoint statistics tell us about path descriptions? | In progress |

The companion [ontology-separation](https://github.com/stevenwarejones/ontology-separation)
repository contains the formal Lean framework. Its
[path-interference case study](https://github.com/stevenwarejones/ontology-separation/blob/a4d294ddb3847134ac2bdbb509c435c2eaff7a98/docs/PATH_INTERFERENCE_CASE_STUDY.md)
explains the distinction between equivalent descriptions and observable model
separation.

Original data and paper PDFs are never committed here. Each study identifies its
official source and a checksum manifest; users obtain and verify their own copies.
Repository code is distributed under the [Apache 2.0 license](LICENSE); external
data and papers retain their stated licenses. Changes are reviewed through pull
requests before merging.
