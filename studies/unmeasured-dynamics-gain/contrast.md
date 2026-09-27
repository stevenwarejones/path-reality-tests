# Two-source gain for a discrimination resource

Both complete sources now contribute indispensable constraints to one operational resource. The result is **not** the still-open two-sided bound on a fixed-input unmeasured probability. It concerns the best binary state discrimination through a specified shared suffix and the detector. Its proof is elementary, and most of the gain already occurs for the detector alone. It is not presented as a novel dynamics theorem or a solution to the broader leakage/memory question.

Let `T = Gy Gx Gy Gy Gx Gx Gy Gx Gy Gy`. This ten-gate word is absent as a complete circuit from all 6,661 recorded words. It does occur as a suffix in both sources. Define

\[
F_T=G_T^*(E),\qquad C_T=\max_{\rho_0,\rho_1}|\operatorname{Tr}F_T(\rho_0-\rho_1)|
=\lambda_{\max}(F_T)-\lambda_{\min}(F_T).
\]

The resource is the maximum total-variation distance between the two binary output distributions, optimizing over physical input states. Equivalently, a sender encodes an equally likely bit in two freely chosen states, applies `T`, and uses the **fixed existing detector** to decode it. The optimal one-use success probability is `S_T=(1+C_T)/2`. The eigenstates of `F_T` attain it. This is a hypothetical capability with controllable inputs; the deposit contains no independently calibrated input-state test, and no such experiment was run. It is not Shannon capacity or a quantum-memory advantage.

## Certified comparison

The model class and simultaneous 95% confidence region are exactly those in [the original derivation](derivation.md): arbitrary stationary qubit CPTP native gates, shared normalized reset and binary POVM; every original cell and group retained; no statistical budget reallocation on source deletion.

| Region / physical construction | Contrast `C_T` | Optimal binary success `S_T` |
| --- | --- | --- |
| Entire joint physical region | **at least 0.7229749561** | **at least 0.86148747805** |
| Model satisfying every GST constraint | [0.6724193955, 0.6724193956] | [0.83620969775, 0.8362096978] |
| Model satisfying every RB constraint | [0.6298184453, 0.6298184454] | [0.81490922265, 0.8149092227] |
| Model satisfying every joint constraint | [0.9925826555, 0.9925826556] | [0.99629132775, 0.9962913278] |

The certified source-removal gaps in contrast are **0.0505555605** (GST-only) and **0.0931565107** (RB-only). These are existence witnesses, not source-only optima. The joint lower bound is not claimed sharp. The compatible joint model prevents an empty-region explanation of the result.

Both new witnesses satisfy all retained inequalities with minimum slack greater than **0.0000099989**, after the inherited `1e-9` numerical envelope. Each GST witness check covers 4,657 cells and six groups; each RB check covers 2,004 cells and six groups. All long GST rows are included. The interval certificate independently checks native Choi positivity, normalized Kraus physicality and the propagated effective POVM. The physical factors, complete probability vectors, and exact rational comparison margins are saved in [the certificate](results/contrast-certificate.json) and its witness files. Passing these conservative intervals still does not establish aggregate likelihood adequacy.

## Why the bound holds for the full class

The selected observations are:

| Global row | Family | Deposited expression | Plus / total | Endpoint used |
| --- | --- | --- | --- | --- |
| 1573 | GST | `GyGyGy(GxGxGyGxGyGy)^5` | 50 / 50 | lower 0.7681878105 |
| 4699 | RB | `GxGyGyGiGxGxGxGyGxGyGyGyGxGyGx(Gy)^9GxGyGyGyGxGyGyGxGxGyGxGyGy` | 0 / 285 | upper 0.0452128544 |

These are different 33- and 37-gate circuits. Both end with exactly `T`. Their preceding 23- and 27-gate words prepare potentially different unknown states `sigma_G` and `sigma_R`. Shared operations give

\[
p_G=\operatorname{Tr}(F_T\sigma_G),\quad p_R=\operatorname{Tr}(F_T\sigma_R),\quad
p_G-p_R\leq C_T.
\]

Every state expectation of an effect lies between its smallest and largest eigenvalues. Thus, on the original simultaneous coverage event,

\[
C_T\geq p_G-p_R\geq 0.7681878105-0.0452128544
=0.7229749561.
\]

This covers **every** member of the declared physical class, including boundary and singular gate sets. No optimizer coverage claim is involved. Although only two rows enter this outer bound, every constraint is imposed on all three physical constructions. A positive source-deletion gap proves that the full retained source cannot already imply this joint bound. This is the missing logical ingredient that the original one-sided target result lacked.

The spectral-range proof itself holds in any finite dimension with a common effect. It does not rely on a qubit-specific deformation, ideal gates or known prefix states. A larger Hilbert-space interpretation would bound the capability of the **entire** effective detector on that space, not necessarily access to the computational qubit alone.

## What the suffix adds—and does not add

For the bare detector, `C_empty=|u-l|`. The same two probability endpoints also imply `C_empty>=0.7229749561`, because all preceding circuits prepare valid states. The new GST-only and RB-only models have detector contrasts approximately **0.6731195677** and **0.6307008705**, already below this bound. The suffix reduces their contrasts by only approximately **0.00070017** and **0.00088243**.

Consequently **most of this combination gain is a detector-contrast constraint**. The ten-gate construction gives the task a specified common channel prefix to the detector and an unrecorded complete word, but it does not demonstrate a large new constraint arising from the repeated-gate structure. The source rows already implement two encodings through this suffix. Their probability contrast witnesses achievable discrimination; the new quantity optimizes over all inputs. This is not a new arbitrary-sequence forecast in disguise.

The resource is invariant under a physical change of Hilbert-space basis. It need not be invariant under all nonunitary gate-set gauge transformations that preserve the archived terminal probabilities. We therefore do **not** interpret a fitted contrast as a measured value: the lower bound ranges over all physical realizations in the region, and the displayed contrasts are explicitly labeled constructions. Fixing independently calibrated input preparations at the suffix entrance and performing detector tomography would resolve its effective effect (three independent Bloch axes plus the trace/bias), allowing a direct resource estimate. No matched calibration of that kind is present in these records.

## Construction and selection

The GST witness holds native gates, reset and detector axis from the saved joint model fixed, and optimizes only the two detector eigenvalues. The RB witness does the same starting from the RB-only model. With projector-response probabilities `z_i`, each prediction is `l+c*z_i`. Linear programming minimizes `c` subject to **all** retained cell and group constraints and `0<=l`, `0<=c`, `l+c<=1`, using `1e-5` inward slack. This LP is a constructive search within a subfamily, not an optimization over the complete physical class.

The retrospective screen considered all common source suffixes of lengths 1–256 that are absent as complete measured words. There were **69** such suffixes; **59** had positive endpoint bounds; **21** had numerical crossings from both new witnesses. The selected target maximized the smaller deletion gap, with declared deterministic ties. The [protocol](contrast-protocol.json), [complete screen](results/contrast-screen.json) and [discovery implementation](contrast_discover.py) preserve this selection. This is not prospective validation. The unchanged simultaneous confidence event supports the selected implication without an additional target-wise confidence allocation.

## Assumptions and adversarial checks

The essential dynamical assumption for the bound is that both occurrences of the suffix and detector implement one common effective POVM. Prefix states can differ arbitrarily; even prefix drift does not affect the spectral argument as long as those states remain valid. A context-dependent detector defeats the interpretation: effects `I` and `0` give a probability difference of one while each has zero spectral diameter. Tests explicitly exercise that failure.

If each context's effective effect differs from a common `F_T` by at most operator norm `epsilon`, then the valid conditional bound becomes `C_T>=max(0,0.7229749561-2*epsilon)`. With ten gate occurrences each within half-diamond distance `delta` of a reference operation, and fixed detector, telescoping gives `epsilon<=10*delta`. These are sensitivity hypotheses, not measured stability bounds. An enlarged memory model may invalidate identification of the same qubit-only effect after different prefixes.

The confidence label retains independent Bernoulli shots and stationary probabilities. Unrestricted within-job dependence removes the 95% guarantee. No large shot count is substituted for this assumption; the earlier dependence and composition controls remain applicable. No new injected-alternative power or held-out calibration claim is made.

## Prior art and increment

| Closest result | Relation and precise increment |
| --- | --- |
| [Oreshkov, Calsamiglia, Muñoz-Tapia and Bagan, *Optimal signal states for quantum detectors*, arXiv:1103.2365v3, equations (9)–(10)](https://arxiv.org/abs/1103.2365v3) | Already gives optimal binary success through POVM spectral spread. This formula and the eigenstate encoding are **not new**. Here it is applied to a shared suffix with simultaneous empirical constraints and both complete-source physical countermodels. |
| [Hirche et al., *Discrimination power of a quantum detector*, arXiv:1610.07644v1](https://arxiv.org/abs/1610.07644v1) | Studies optimized input-state discrimination and multi-use error rates. This study makes only a one-use binary claim and supplies no new asymptotic theorem. |
| [Original GST/RB experiment and earlier audit](../shared-quantum-dynamics/prior-art.md) | The acquisition pairing and its qualitative complementarity are established prior work. The increment here is the explicitly certified source-deletion comparison for this resource under one confidence region. |
| [PR17 measured-context contrast](../gst-rb-combination-gain/README.md) | The resource is an extremum over input states, not a fixed mixture of recorded contexts. However, its lower bound is still an elementary observed-response range, and most gain is detector-only. |

This establishes a delimited two-sided operational-resource restriction. It does not establish novelty, consequential dynamics identification, a leakage/memory tradeoff, or a two-sided prediction for one fixed reset-input unmeasured circuit. Those ambitions remain open.

## Reproduction

```sh
python studies/unmeasured-dynamics-gain/contrast_certificate.py --check
python studies/unmeasured-dynamics-gain/contrast_certificate.py --cache /tmp/unmeasured-source --check
python studies/unmeasured-dynamics-gain/contrast_discover.py --cache /tmp/unmeasured-source --check
python -m unittest discover -s tests -p 'test_discrimination_contrast.py' -v
```

Use the original hash-checked downloader to populate the external cache. Offline verification checks pinned word hashes, the exact endpoint difference, physical factors and every saved retained count constraint. Full-source verification reconstructs every original word, binds all model probabilities to both Bloch and independent density propagation, and checks target absence. The discovery route re-solves the small LPs and re-screens the full suffix family; optimizer failure never certifies exclusion.
