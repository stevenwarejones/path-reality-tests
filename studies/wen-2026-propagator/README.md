# Wen 2026: propagator and path reconstruction

Study of Wen et al., *Direct experimental test of Feynman's path integral
postulates with single photons*, Science Advances 12, eaeh1011 (2026),
[doi:10.1126/sciadv.aeh1011](https://doi.org/10.1126/sciadv.aeh1011).
The official dataset is [doi:10.5061/dryad.x0k6djj14](https://doi.org/10.5061/dryad.x0k6djj14),
version 3, version ID 450489, released under CC0.

This study prepares reproducible comparisons between reported endpoint
measurements, reconstructed propagators and explicitly specified models. This
initial scaffold supplies acquisition verification, read-only workbook inventory
and synthetic numerical checks. It does not report an empirical replication or
findings about the paper.

## Run the tools

Use Python 3.12 with the pinned dependencies. From this study folder:

```sh
python -m pip install -r requirements.txt
python fetch_data.py --verify-only
python extract_workbooks.py
python synthetic_checks.py
```

The pinned dependency versions were exercised on Python 3.12. The scripts resolve
paths relative to their own location, so they may also be invoked by path from
another working directory. Place user-downloaded originals in `raw/dataset/`;
that subfolder keeps the dataset's own `README.md` separate from repository
instructions. Neither originals nor complete extracted cell inventories are
committed. `--directory` selects an alternative input folder for verification and
extraction without changing the manifest.

`manifest.json` fixes all 19 filenames, byte counts and SHA-256 digests. Verification
must succeed before extraction. A missing or mismatching file yields a nonzero
exit status; do not adjust the manifest to fit a local copy. Running
`python fetch_data.py` without `--verify-only` attempts the listed public download
URLs and stops further requests after an access-control or rate-limit response.

The extractor preserves sheet names, hidden state, coordinates, formulas and
cached values in ignored `results/workbooks/` files. It inventories cells; it does
not decide scientific column mappings. The preliminary [data dictionary](data-dictionary.md)
describes source-reported roles pending sheet-by-sheet inspection.

`results/synthetic.json` and `results/synthetic-shapes.png` contain calculations
from the declared finite-grid model, with reproducible noise seeds. They are not
measured data or a fit. `provenance.json` identifies sources and hashes; original
papers and supplementary PDFs remain outside git.

## Formal companion

The Lean case study belongs in `ontology-separation` and certifies an ideal
finite example rather than this apparatus. Its [pinned source](https://github.com/stevenwarejones/ontology-separation/blob/a4d294ddb3847134ac2bdbb509c435c2eaff7a98/OntologySeparation/Experiments/PathInterference.lean)
and [case-study document](https://github.com/stevenwarejones/ontology-separation/blob/a4d294ddb3847134ac2bdbb509c435c2eaff7a98/docs/PATH_INTERFERENCE_CASE_STUDY.md)
identify the exact companion revision (PR #84, pending review).
