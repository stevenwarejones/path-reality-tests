# Literature and code comparison

Audit date: 2026-09-27. These are primary sources. Entries distinguish detailed
methods inspected from leads screened at publication/abstract level. A common
mechanism is prior art; the missing increment is a calibrated, unfitted prediction
or defensible full-model separation from complementary measurements.

| Work | Already established | Calibration / validation audit | Increment achieved here |
|---|---|---|---|
| [Iaia et al., Nature Communications 13, 6425 (2022)](https://doi.org/10.1038/s41467-022-33997-0); [arXiv v2](https://arxiv.org/html/2203.06586v2) | Junction injection, normal-metal mitigation and background parity coincidences | Full figure archive inspected; main and supplementary HTML methods read. Injection at selected junctions is not an arbitrary-position radiation calibration. Fit confidence intervals and extracted coincidence rates are processed quantities | Reproduced selected observed-rate ratios and injection peak delay; no novelty claim |
| [Yelton et al., Physical Review B 110, 024519 (2024)](https://doi.org/10.1103/PhysRevB.110.024519); [arXiv v2](https://arxiv.org/html/2402.15471v2) | Detailed G4CMP phonon/QP response and radiation footprint modeling, including 1 and 10 µm Cu configurations | Equations and appendices inspected. Boundary loss and trapping are fitted; recombination is fixed in the low-density analysis. The 1 µm devices are closer to the gamma experiment than the 2022 archive, but no version-pinned complete calibration run package for matching those chips/cooldowns was established here | Identified the closest transport precedent and the missing matching calibration; no G4CMP reproduction |
| [Larson et al., PRX Quantum 6, 030339 (2025)](https://doi.org/10.1103/2lyd-8swv); [arXiv v1](https://arxiv.org/html/2503.07354v1) | Controlled gamma irradiation, charge/parity correlations, mitigation footprints | Version-specific data resolved from concept DOI. Equations audited against the accessible arXiv version; published page confirms the data citation. Appendix C.5 fits the density-to-poisoning threshold using the footprint observations. Some charge-pair comparisons are unused in its parameter search, but not a complete withheld intervention | Exact envelope diagnostics and conditional parity semantics; neither is a demonstrated new physical mechanism |
| [Kono et al., Nature Communications 15, 3950 (2024)](https://doi.org/10.1038/s41467-024-48230-3) | Pulse-tube-linked correlated excitation and transition behavior; phase-conditioned correlation analysis | Main methods and selected notebooks inspected. Original code uses phase alignment and aggregate contingency counts. We independently compute plug-in MI and total-covariance decomposition, not the authors' finite-sample estimator. Processed phase counts lack run-level dependence information | Numerical replication and an explicit distinction between pooled modulation and within-phase association |
| [Benevides et al., Physical Review Letters 133, 060602 (2024)](https://doi.org/10.1103/PhysRevLett.133.060602); [data](https://zenodo.org/records/11548753) | Localized infrared irradiation with power, duration and position variation | Alternative source screened. Metadata and Figure 3/SM8 notebook inspected as text. It reads fitted decay/trapping/recombination/density arrays. Raw measurements and HDF5 were not downloaded or reconstructed here | A prospective calibration alternative, with no established cross-chip parity transfer |
| [Stress-induced phonon bursts (2024)](https://doi.org/10.1038/s41467-024-50173-8) | Mechanical/stress mechanisms and the cryogenic detector connection predate this study | Primary publication screened; source archive and full calibration not audited here | Prevents claiming the stress/radiation juxtaposition as novel |
| [Correlated phase-error bursts (2026)](https://doi.org/10.1103/1bl4-b2f7) | A gap-engineered array can retain correlated phase errors despite relaxation mitigation | Primary publication/author abstract screened; no source arrays reconstructed | Prevents treating parity/T1 mitigation as a general error-channel guarantee |

## Public code

The [G4CMP project](https://github.com/G4CMP/G4CMP) documents phonon transport,
charge-carrier physics and superconducting-film/QP response. Its current code
is not automatically the version or apparatus configuration used in the earlier
papers. No current simulator was substituted for a pinned original configuration.

The Kono archive includes executable notebooks and local helpers. We read them
to recover clocks, table meanings and processing choices but did not execute them.
The independent parser loads only the primitive phase-count object. Original
helper functions and plotted bias corrections are not silently copied into the
new analysis. Other downloaded processed files may contain numpy or sklearn
objects; they are hashed but not deserialized by the reproduction code.

The infrared notebook is a screened alternative, not a core numerical input:
`Figure_3_and_SM8.ipynb`, 949255 bytes, available from the versioned record above.
Its acquisition hash is recorded in [screened-sources.json](screened-sources.json).
The companion Figure 2 notebook request returned HTTP 504; no claim to inspect
that notebook is made. A failed request is not evidence that its data are absent.

The cosmic-ray array [Li et al. paper](https://doi.org/10.1038/s41467-025-59778-z)
and [Harrington et al. paper](https://doi.org/10.1038/s41467-025-61385-x) are distinct
leads. Neither is counted as an analyzed dataset here. Muon tagging alone would
not supply a lower tag-efficiency bound for gamma-induced events.

The sharp parity interval is a standard moment argument applied to this
observation problem. The rational affine certificates are an elementary global
consistency tool. Neither is claimed to have research priority, to establish
full physical-model identifiability, or to satisfy the mission on its own.


The mechanical follow-up uses the released acquisition-level dwell sequences
rather than only the published phase covariance. Its increment is a specified
99 µs paired-excursion forecast and a fixed-policy validation split; novelty and
full-source indispensability are not established. Kono et al. already report
phase-dependent correlations, continuous-readout backaction and mechanical
controls. Their microwave monitor is around 7 GHz and is not a rare-classification
error calibration for this new target. Trigger-synchronous mitigation is also
an established broader direction: [Tathed et al., arXiv:2606.00358](https://arxiv.org/abs/2606.00358)
use software compensation for coherent control-frame disturbances in a trapped
ion. That is a different platform and error channel; it prevents treating the
generic idea of timing-aware control as new here.
