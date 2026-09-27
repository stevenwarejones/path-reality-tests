# Bounded novelty audit

Search date: 2026-09-26. No priority, PhD-level novelty, or comprehensive literature
coverage is claimed. Searches with broad combinations of "collapse", "spectrum"
and "optimization" returned many unrelated results; those are not evidence that
a relevant result does not exist.

## Closest primary work

- [Toroš, Gasbarri, Bassi, arXiv:1601.03672](https://arxiv.org/abs/1601.03672),
  PLA 381, 3921 (2017): already compares interferometry with a declared
  macroscopic-localization requirement and colored/dissipative variants.
  Its colored approximations have explicit validity conditions; Eq. (13)
  gives a short-memory condition. Combining "coherence survives" and
  "macroscopic localization must occur" is therefore not new.
- [Carlesso, Ferialdi, Bassi, arXiv:1805.10100v3](https://arxiv.org/abs/1805.10100v3),
  EPJD 72, 159 (2018): derives spectral mechanical responses and compares
  frequency sensitivity under exponential temporal correlations.
  Low-frequency robustness relative to radiation is already established.
- [Piccione, arXiv:2511.00644v2](https://arxiv.org/abs/2511.00644v2),
  [PRA 113, 062211](https://doi.org/10.1103/h74z-zvfv), full PDF inspected:
  optimizes spatial smearing at fixed spatial variance. This differs from
  the temporal long-memory sequence here, but optimization of collapse
  phenomenology is not itself novel.
- [Altamura, Vinante, Carlesso, arXiv:2501.08971v2](https://arxiv.org/abs/2501.08971v2):
  torque/force ratios already distinguish geometry and gas-noise assumptions.
  Secs. IV–V specifically warn about subtracting phenomenological outgassing
  without independently measured pressure. Our candidate comparison retains
  that warning and does not advertise the ratio as a new discovery.
- [XENONnT, arXiv:2506.05507v2](https://arxiv.org/abs/2506.05507v2):
  the current Markovian radiation test includes electron/proton cancellation.
  Its full analysis is not reproduced by recasting an older low-energy subset.

## Baseline publications

- [Pedalino et al., Nature DOI 10.1038/s41586-025-09917-9](https://doi.org/10.1038/s41586-025-09917-9)
  and pinned Zenodo code. This study independently implements its deposited
  power-scan optics and fixed-calibration Bayesian update. The separate
  author-module oracle tests the normalization.
- [Vinante et al., arXiv:2002.09782v2](https://arxiv.org/abs/2002.09782v2),
  PRL 125, 100404 (2020): Eq. (2), Eq. (5), supplement S14 and workbook formulas
  drive the baseline. Raw fits cannot exactly match without Qprime/mask records.
- [Toroš and Bassi, arXiv:1601.02931](https://arxiv.org/abs/1601.02931):
  detailed matter-wave calculations, read as a normalization/model reference.
- [Helou et al., PRD 95, 084054](https://doi.org/10.1103/PhysRevD.95.084054):
  older LPF translation constraint is contextual; newer rotation is used for
  the implemented space-envelope candidate.

## Discovery queries and decisions

| Query / source route | Result used | Decision |
|---|---|---|
| "Narrowing the parameter space" layered force sensors | Paper, supplement and deposit | Reproduce measured baseline |
| "Colored and dissipative continuous" interferometry | 1601.03672 and primary journal | Reject unqualified long-memory reuse |
| "Diffusion minimization via optimal smearing" | 2511.00644v2 full PDF | Spatial optimization prior art |
| "CSL" "arbitrary noise spectrum"; "collapse models" "spectral optimization" | Poor relevant recall across two engines | No novelty inference from search absence |
| LISA Pathfinder rotational noise data / archive | 2501.08971v2, NASA index and ESA metadata | Implement published-envelope feasibility |
| XENONnT collapse data release likelihood | Collaboration public-data page → Zenodo 7992017 | Download and recast only declared 1–30 keV subset |
| Cold atoms collapse heating data | 1409.5388, 1605.01891; no sufficient release audited | Reject as implemented constraint |
| 2609.07195 / Zenodo 22637584 | Direct access failed | Neither measurement nor projection imported |

## Contribution actually supported

Reproducible cross-source audit, explicit source/formula discrepancies, conditional
separate-versus-joint calculations and a constructive diagnostic of the
zero-frequency gap. The finite LP dual and exact affine certificate are standard
mathematics. The OU asymptotic is not claimed as a new theorem. Its extension to
the actual moving interferometer and calibrated instrument windows is unfinished;
that extension, if successful and different from the literature, is the candidate
research contribution.

## Moving-response follow-up

The source-linked [moving derivation](moving-frame.md) explicitly compares:

- [Adler 1807.11450v3](https://arxiv.org/abs/1807.11450v3), especially Eq. (10–11):
  colored spatial noise, boosts and directional effects are prior art. A preferred
  frame or a Doppler slice is not claimed as new.
- [Carlesso et al. 1805.10100v3](https://arxiv.org/abs/1805.10100v3), Eq. (11–19):
  the stationary separable mechanical spectrum is the old baseline, whose motion
  premises are now replaced for the new response bound.
- [Toroš et al. 1601.03672](https://arxiv.org/abs/1601.03672), Appendix A2:
  the short-memory expansion's kinetic-energy condition is not extended to long
  memory. The new calculation instead declares a local *spatial-advection*
  Markov approximation, with small parameters and an independent path check;
  it does not assert equality to the full colored quantum grating map.
- [Altamura et al. 2501.08971v2](https://arxiv.org/abs/2501.08971v2):
  torque geometry and force/torque comparison are established. The present
  extension allows delayed correlations between separated masses.
- [Armano et al. 2405.05207](https://arxiv.org/abs/2405.05207): run table, calibrated
  angular acceleration, windowed/Wishart methods and nuisance decorrelation.
  This motivates the uniform finite-observation bound and the refusal to call
  housekeeping or nominal orbit vectors a calibrated joint time-series fit.

PDFs were acquired and the relevant response/acquisition derivations inspected.
No comprehensive novelty conclusion is drawn. The specific repository result is
an executed joint dephasing comparison with measured source summaries and stated
auxiliary/model limitations; priority and PhD-level significance remain unclaimed.
