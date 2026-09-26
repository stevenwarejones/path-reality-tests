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

## Closest marginal-dependence analyses

[Bednorz, *Physical Review A* 95, 042118 (2017)](https://doi.org/10.1103/PhysRevA.95.042118),
Section IV ([author PDF](https://www.fuw.edu.pl/~abednorz/ab17.pdf)), analyzes NIST
marginal setting dependence, principally Classical XOR3, then XOR1 and XOR2.
He considers pulses 28–800 at Alice and 37–800 at Bob, with phase regions inside
and outside 90±16 and 125±20 bins, four times the nominal radii. His Table III
caption identifies the outside-window counts (the nearby prose interchanges
Tables III/IV). He reports an adjusted p-value near 0.049 for one XOR3 comparison
and notes that the late pulses can allow ordinary subluminal communication.

This audit uses the separate 03_43 post-recovery run, primary bits [4,5,6] and
one-/five-bit sensitivities within the stored sixteen-bit windows. Its ±one-bin
radius variants do not implement Bednorz's broad phase regions or late-pulse
search. Therefore its intervals containing zero neither reproduce nor contradict
his numerical result. A direct comparison requires the same runs, windows, event
rules and inferential target. We have not performed that replication. Marginal
setting-dependence testing itself is not new here. The contributions are complete
raw/stored reconciliation, an auditable bitwise-OR comparison, conditional
predictable-mean bounds and an explicit physical-evidence ledger.

[Shalm et al. supplementary material](https://journals.aps.org/prl/supplemental/10.1103/PhysRevLett.115.250402/LHFSupplementary.pdf),
Section IV.C.2, already reports four no-signaling comparisons per pulse grouping
and run, under IID diagnostic assumptions. Table S-II gives five-pulse Classical
XOR3 counts; Table S-III summarizes 120 no-signaling tests. Those diagnostics do
not supply this audit's history-conditional assignment calibration. Section I.C
and the RNG characterization discuss predictability for the Bell test; Section II
supplies experimental timing context. Transporting either to this run's paired
endpoints and joint-history TV premise requires additional evidence.

[Howard et al., *Annals of Statistics* 49, 1055–1080 (2021)](https://doi.org/10.1214/20-AOS1991)
and [Howard et al., *Probability Surveys* 17, 257–317 (2020)](https://doi.org/10.1214/18-PS321)
develop time-uniform inference using nonnegative supermartingales. Our v2 uses
the elementary Bernoulli exponential process and a finite lambda-grid union,
derived fully in [click-statistics.md](click-statistics.md). It is not an
implementation of their empirical-Bernstein theorem or a claim of optimality.

NIST source inputs remain external and hash-pinned. No raw-event or original
count-table extracts are committed. The committed result tables are derived
analysis artifacts; offline decoder examples are explicitly synthetic. Notice
and attribution remain in NOTICE.txt. PDFs are not committed.
