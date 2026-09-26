# Theorem-to-apparatus comparison and annotated bibliography

Primary-source review dated 2026-09-26 UTC. This is a bounded literature search,
not an exhaustive priority claim. Retrieval depth and unresolved version/supplement
gaps are retained in [evidence.json](evidence.json). The calculations in
[derivation.md](derivation.md) independently multiply the specified finite instrument.

## Closest comparisons

| Work | Relevant premises and quantities | Consequence for this experiment |
|---|---|---|
| [Pusey 2014](https://arxiv.org/html/1409.1535v2), PRL 113, 200401 | Measurement noncontextuality, sharp-outcome determinism, noisy projector readout and disturbance equivalence for the final effect. Its p₋ divides joint negative-success by **bypass** success. | Do not substitute a post-probe conditional fraction; avoid relying on exact sharpness in a real detector. |
| [Kunjwal–Lostaglio–Pusey 2019](https://arxiv.org/pdf/1812.06940), PRA 100, 042116 | Theorem 3 has finite pointer outcomes and measurement/transformation equivalences. Lemma 5 allows stochastic final readout. Its p₋ is joint. | Use a≤qf+d(1−f), with q=(1+p_m)/2, and check the complete channel, not only one final marginal. |
| [Piacentini et al. 2016](https://arxiv.org/pdf/1602.02075), PRL 116, 180401 | Heralded photons, polarization system and spatial pointer; four preparation states for pointer calibration; final-effect checks and tomography. Reported anomaly 5.7 standard deviations under its analysis. | Closest experimental predecessor, but its sharpness/effect-level controls do not by themselves certify our channel equivalence. Its low pointer detection efficiency is not a demonstrated value for this proposed ancilla detector. |
| [#85 local phase](https://github.com/stevenwarejones/ontology-separation/blob/ffbbd0d3c0eb2066ea9b32b0dfd18e9c2d50bd4d/docs/PHASE_INTERVENTION.md) | Common preparation, definite region, occupation matching, invariant Q response; complete phase contrast≤occupation. | Different rejected conjunction. Occupation matching supplies neither instrument equivalence. The new bound does not require the Q-locality premise. |

Pusey's original effect-disturbance condition is weaker operationally than a full
channel equivalence but uses an ideal sharpness premise. KLP's alternative
Theorem 4 trades transformation noncontextuality for extra preparation/readout
constraints. We choose the finite transformation route because the complete
Kraus channel has a direct identity/Z implementation and no sharp-final-response
assumption. This is a choice of assumptions, not a claim of a universally stronger
test. No preparation-noncontextuality theorem is silently added to our null.

KLP Appendix F supplies the relevant finite-pointer family. The reference here
uses rational Kraus amplitudes and exact multiplication at finite strength.
Operational equivalence still must be supported independently. The derivation is
a known-theorem specialization; neither finite outcomes nor finite measurement
strength is claimed as new physics.

## Operational equivalence and finite precision

[Mazurek et al.](https://arxiv.org/pdf/1505.06244), published as Nature
Communications 7, 11780 (2016), constructs secondary procedures by convex mixing
primary procedures and supplementary controls. Its preparation/measurement
example is not an automatic construction for our channel equivalence. The present
design does not claim a validated secondary-instrument solution. It instead states
the representation premises explicitly, budgets operational audit collection, and
proves only conditional robustness. Extra preparations, transformations and a
justified tomographic space would be needed before promoting this to a calibrated
experimental generalized-noncontextuality test.

An unchanged marginal can hide disturbance, as our two-state countermodel proves.
[Ipsen 2014](https://arxiv.org/abs/1409.3538) was screened as a primary discussion
of different disturbance definitions. [Ferrie–Combes 2014](https://arxiv.org/abs/1403.2362)
and the ensuing [Brodutch comment](https://arxiv.org/abs/1410.8510) are backward
citation leads in the debate; their full comment/reply pair was not verified here.
No claim about the correctness of that exchange is needed for our explicit kernels.

[Halliwell 2017](https://arxiv.org/abs/1704.01485) distinguishes Leggett–Garg and
no-signaling-in-time tests. These use different noninvasiveness/classicality
conditions from the equivalence-based null here. Positivity of a selected
history representation, marginal NSIT and generalized noncontextuality are not
interchangeable requirements. We do not infer one from another.

## Later work and claim gate

[Zhao et al. 2026](https://www.cpsjournals.cn/en/article/doi/10.1088/1674-1056/adf4b1),
Chinese Physics B 35, 020301, studies finite disturbance and multiphoton settings.
The retrieved theorem retains sharp final measurement/outcome determinism. Its
single-photon finite-strength analysis is a relevant prior result, not evidence
of novelty here. Some displayed formula objects were missing from retrieved HTML;
its detailed corrections are not used in our calculation.

[Fan et al. 2026](https://www.nature.com/articles/s42005-026-02715-3) demonstrates
a contextuality/entanglement construction using remote preparations and structural
ensemble constraints. It does not supply the identity-plus-disturbance channel
equivalence for this sequential probe. This is a useful alternative architecture,
but importing its result would change the experiment and model class.

The 2026 search also found [Exact No Signaling in Time without Temporal
Classicality](https://arxiv.org/abs/2607.14583). It is abstract-screened only;
our conclusion about marginals already follows from the exact two-state example.
No theorem in that preprint is imported.

## Trajectory and interaction-time motivation

[Kocsis et al. 2011](https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=907599),
Science 332, 1170, reconstructs ensemble flow from weak transverse momentum and
position postselection at different planes. It does not track each photon through
all planes. The bounded path projector in our experiment is a different observable.
Its original supplementary acquisition details were not independently reanalyzed.

[Foo, Asmodelle, Lund and Ralph 2022](https://refubium.fu-berlin.de/bitstream/handle/fub188/39134/s41467-022-31608-6.pdf?sequence=1),
Nature Communications 13, 4002, proposes an operational velocity from a ratio of
real weak momentum and energy values, connected to Klein–Gordon current/density.
Its scalar single-particle description and optical regime k₀≫σ matter; density
positivity outside that regime is problematic. It is a theoretical proposal,
not a signaling demonstration. [A later experimental preprint](https://arxiv.org/html/2509.11609v1)
reconstructs relativistic Bohmian quantities. Its existence updates the original
proposal-only literature picture, but does not establish individual-worldline
ontology or the required instrument-preserving map to this projector test.

[Angulo et al. 2026](https://doi.org/10.1103/gjfq-k9dv), PRL 136, 153601, measures
a conditional excitation-time weak value; the accessible
[2024 preprint](https://arxiv.org/html/2409.03680v1) describes the cross-Kerr probe
and time-integrated conditional phase signal. This probes interaction time and
cannot fairly be dismissed as only pulse-peak reshaping. A negative value does
not establish chronological reversal. Published-version and supplement equality
were not verified, so apparatus numbers from that preprint are not adopted as
specifications for our proposal.

The obstruction is concrete: none of these momentum/energy/time instruments is
part of the allowed bounded-projector interface. Opposite unconstrained history
labels can be attached without changing that interface's probabilities. A physical
bridge must restrict those labels through dynamics and calibrated instruments.
Our experiment therefore tests the stated noncontextual representations, not
superluminal motion, backward time or all definite trajectories.
