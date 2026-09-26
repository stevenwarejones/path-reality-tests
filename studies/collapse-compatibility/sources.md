# Sources and access

Accessed 2026-09-26. Original files remain external. Numerical hashes and byte
sizes are in [manifest.json](manifest.json). No archive, workbook, event list,
author source file, paper PDF or FITS file is redistributed.

| Source | Evidence and audited access | Attribution / version |
|---|---|---|
| [Sodium interference](https://zenodo.org/records/17502163) | Downloaded ZIP; 95 scans of 41 measured positions/counts, acquisition workbook, mass spectra and code | Pedalino et al.; DOI 10.5281/zenodo.17502163, v1 |
| [Layered sensor](https://zenodo.org/records/3956096) | Downloaded PSD ZIP and final workbook; 14 spectra plus fit summaries/formulas | Carlesso and Vinante; DOI 10.5281/zenodo.3956096, v1 |
| [XENONnT ER subset](https://zenodo.org/records/7992017) | Downloaded ZIP; 433 selected ER energies, efficiency, background and other-model limits | XENON Collaboration; DOI 10.5281/zenodo.7992017, v2 |
| [NASA LPF DRS](https://heasarc.gsfc.nasa.gov/w3browse/all/lpffiles.html) | Directory and one 2016 Jan 10–16 FITS object downloaded; the selected small object contains HOUSEKEEPING, not angular science records | HEASARC/ST7; reformatted telemetry from ESA |
| [ESA full LPF archive](https://esdcdoi.esac.esa.int/doi/html/data/astronomy/lisa-pathfinder/Full.html) | Official metadata inspected; recommends PreProcessed because full data can contain timestamp inconsistencies | DOI 10.5270/esa-fc52vb6, v1.0; CC BY-NC 3.0 IGO |
| [LPF rotational bound](https://arxiv.org/abs/2501.08971v2) | Published numerical constraint, not a downloaded experimental likelihood: torque envelope 5.7e-34 N²m²/Hz at 3 mHz | Altamura, Vinante, Carlesso (2025), Sec. III |
| [XENONnT collapse search](https://arxiv.org/abs/2506.05507v2) | Full paper/supplement inspected; Markovian result uses 1–140 keV; its event release/likelihood is not supplied by the 2022 subset | XENON Collaboration, PRL 136, 120201 (2026) |
| [1D quantum-gas lead](https://arxiv.org/abs/2609.07195) | Direct retrieval failed; archive lead 10.5281/zenodo.22637584 not downloaded in this session | Unverified lead; no achieved limit imported |

Zenodo rendered license sections were empty in the retrieved pages; no
redistribution permission is inferred. Author code is inspected and one module is
executed only after archive verification, solely as a full-source numerical oracle.
Our implementation of the formulas is separate. The author power-scan/macroscopicity
scripts are not executed wholesale.

The NASA file is an access probe, not evidence for a CSL bound. The selected space
calculation uses the explicitly cited published torque envelope. Unrelated
observations' timestamps are never aligned.

The 2022 XENON release is for PRL 129, 161805. Its 433 events are under 30 keV.
Bosonic-dark-matter limits in that ZIP extending to 140 keV are not higher-energy
ER events. Background summaries are not multiplied into another likelihood.
