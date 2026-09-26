# Sources and access

Access date: 2026-09-26. This is a bounded data/processing audit, not a new literature
review or a claim of exhaustive Bell-data coverage. One run was selected from the
catalog for complete paired raw access and a post-fix experimental configuration.
The poor 03_31 training run and other independent runs are outside the analysis.

| Source | Inspected material | Use and limit |
|---|---|---|
| [NIST repository landing page](https://www.nist.gov/pml/applied-physics-division/bell-test-research-software-and-data) | Repository scope, metadata warning, redistribution notice | Official provenance; presence in the catalog alone does not certify suitability |
| [File/folder description](https://www.nist.gov/document/bell-test-data-file-folder-descriptions) | All five pages, version 2015-12-23 | Channels, formats, run quality, HDF fields, warning about spreadsheet no-click counts |
| [2017 addendum](https://s3.amazonaws.com/nist-belltestdata/belldata/File_Folder_Descriptions_Addendum_2017_02.pdf) | Three pages | Blind data later released; avoids treating the old description's absence claim as current |
| [Processing description](https://s3.amazonaws.com/nist-belltestdata/belldata/code/analysis/DataProcessingDescription.pdf) | Full one-page flow | GPS start, processing cuts and ambiguous-setting handling |
| [Analysis code ZIP](https://s3.amazonaws.com/nist-belltestdata/belldata/code/analysis/bell_analysis_code.zip) | Archive inventory; build_file_hdf5.py, settings_hdf5.py, sync_hdf5.py, hdf5_to_peter.py, hdf5_to_peter_v2.py, hfd5_to_newclicks.py, make_compressed.py, configurations and diagnostics.v2.peter.csv | Reconstruct formats and update behavior; not a claim to have run the original unavailable dependency chain |
| [Compressed raw catalog](https://www.nist.gov/pml/applied-physics-division/bell-test-research-software-and-data/repository-bell-test-research-1) | Selected Alice/Bob 03_43 archives, full streams | Independent reconstruction over the retained overlap, channel/endpoint inventory |
| [Processed catalog](https://www.nist.gov/pml/applied-physics-division/bell-test-research-software-and-data/repository-bell-test-research-2) | Full selected HDF5, every dataset and metadata group | Reconcile stored rows; note missing offsets and syncNumber rather than assuming documentation matches |
| [Shalm et al., arXiv:1511.03189v2](https://arxiv.org/abs/1511.03189v2), [PRL 115, 250402](https://doi.org/10.1103/PhysRevLett.115.250402) | Main paper apparatus, pulse aggregation and predictability discussion; archived data documentation supplies processing details | Experimental context only; no imported joint-history TV calibration or reproduction of Bell p-values |

NIST supplies these data and code for public use with notice preservation. The
small real-event fixtures and extracted diagnostic rows retain attribution and
the full notice in NOTICE.txt; their extraction date/nature is recorded. Large
archives and papers are not committed. All measured counts derive from the
manifest-pinned run; synthetic regression cases are identified as tests.
