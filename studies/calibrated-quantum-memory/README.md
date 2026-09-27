# Calibrated quantum environmental memory

**Status: empirical feasibility gate not passed.** This study does not establish
quantum environmental memory from existing matched datasets. It supplies a
reproducible measured-data diagnostic, an exact restricted ambiguity, and a
prospective certificate tested on explicitly synthetic records.

The central question is whether dynamics require more quantum information in
the environment than a calibrated intervention can leave in the system.
Classical environmental memory may still accompany quantum system evolution.
The [feasibility gate](feasibility.md) records the candidate acquisitions,
matching evidence, alternatives checked, and missing inputs. The
[theory](theory.md) states and proves the narrower implemented result.

## Measured acquisition result

We downloaded the NMN-tomo archive at
`154235f8bbf5e70eb71c325370a67b1894490452`. All 3,240 archived rows retain
four outcome counts. Their row totals are not constant. The following inferred
inverse of outcome-dependent re-preparation regrouping restores constant
totals **for every circuit**, without losing or duplicating a count:

\[
n^{\mathrm{physical}}_{p}(0,c)=n^{\mathrm{archive}}_{p}(0,c),\qquad
n^{\mathrm{physical}}_{p}(1,c)=n^{\mathrm{archive}}_{\bar p}(1,c).
\]

Here p and bar-p are opposite eigenstates of the same axis; the other three
circuit labels are held fixed. Applying the permutation twice restores the
original archive. It matches the paper's outcome-dependent preparation
description, but no original job file or preprocessing implementation was
available to establish that this is the historical transformation. The paper
also contains different shot-number descriptions. We use the integers in the
archive and make the inference assumption explicit.

| Source | Delay settings | Archived row totals | Regrouped shots per circuit | Counts across all settings |
|---|---:|---:|---:|---:|
| IBM Perth | 9 | 7,805–8,195 overall | 8,000 | 23,328,000 |
| UQ | 1 | 8,189–10,243 | 9,216 | 2,985,984 |

Treating an archived row as an independently acquired multinomial experiment
would miss the outcome-dependent regrouping. The parser preserves all outcomes
and uses the inferred physical circuit denominators instead.

A basic shared-process requirement is that the first outcome distribution,
for fixed initial preparation and first measurement, is independent of later
re-preparation and final measurement settings. In UQ the strongest contrast is
between `xp,y,yp,x` and `xp,y,xm,x`: first-outcome-0 frequencies are 0.74641927
and 0.58018663. The difference is 0.16623264.

For M=3,240 rows, alpha=0.05 and N=9,216, the simultaneous Bernoulli Hoeffding
radius is r=sqrt(log(2M/alpha)/(2N))=0.02527218. If all actual row marginals
were within epsilon of one shared past marginal, the triangle inequality
would require

\[
\epsilon\ge\max\{0,(0.16623264-2r)/2\}=0.05784414.
\]

This selection is covered by the simultaneous region over all rows, including
all nine IBM settings; the largest contrast is not treated as prespecified.
For unequal denominators, the implementation uses the appropriate row radii.
The bound is on worst-case marginal discrepancy, also a lower bound on
worst-case joint-distribution TV discrepancy from the shared-past model.

**Interpretation:** conditional on the regrouping and sampling assumptions,
this excludes a common first-outcome distribution for the compared UQ contexts.
Possible explanations include drift, compilation/context dependence, or
incorrect record semantics. It excludes that shared description for both
classical and quantum processes; it is not backwards causation or a witness
of quantum memory. A fixed readout-error map shared across the contexts does
not by itself remove this discrepancy. Each IBM setting has a zero lower
bound with this conservative diagnostic, which does not prove a physical
model fits all its observations.

The bound concerns acquisition context dependence and is **not** the
intervention diamond-distance epsilon in the prospective theorem below.

## Prospective combination result — synthetic only

For a single-outcome qubit reset B, with trusted external probes and no
environmental coupling or extra output register, [the theory](theory.md)
proves

\[
F_6\le2/3+\tfrac12\|B-D\|_\diamond
\]

for arbitrary classical environmental memory, where D is a completely
depolarizing reset and F_6 is average six-state recovery fidelity. Eighteen
independent calibration settings bound every affine qubit-channel coefficient;
six dynamics settings estimate F_6. There is no optimizer and no finite-memory
truncation. A common simultaneous confidence construction gives a global
composite-null test for this declared protocol.

At 8,000 shots per setting (144,000 calibration shots plus 48,000 dynamics
shots), the saved synthetic example gives:

| Quantity | Value |
|---|---:|
| Calibration upper bound on half diamond distance | 0.27314757 |
| Dynamics lower bound on recovery F_6 | 0.98035664 |
| Classical-memory upper bound on F_6 | 0.93981424 |
| Separation margin | 0.04054240 |
| Minimum error required by the null | 0.31368997 |

The confidence level is at least 95% under the stated assumptions. The numbers
are conservative evaluations of analytic inequalities in ordinary double
precision, not interval-arithmetic certificates. No solver tolerances enter;
the decision uses a 1e-12 positive-margin guard. This example's margin is far
larger than floating-point residuals in the checked physical constructions
(maximum 1.12e-16).

The physical source-removal witnesses are explicit:

* Dynamics alone: identity system propagation through an uncalibrated failed
  reset, with no environmental memory, fits every dynamics constraint.
* Calibration alone: a working reset and trivial environment fit every
  calibration constraint, but predict only F_6=1/2.
* Joint quantum feasibility: two SWAPs store the state in an environment qubit
  across a working reset and fit all 24 probability constraints.

Thus the two **synthetic measurement families** are indispensable to this
compatibility decision. This does not demonstrate empirical combination gain
for NMN. The observation of the whole terminal channel alone cannot separate
the first and third constructions; one output-sensitive intervention probe
separates that pair, while all 18 settings support the full-class bound.

## Adversarial controls and sensitivity

The deterministic seeded study uses 1,000 simulated repetitions per model.
All six declared null controls give 0 rejections; each individual 95% binomial
upper bound is 0.003682. They include retained classical Z and rotated-basis
records, a four-outcome tetrahedral record, no memory, and complete/partial
reset failure. These checks do not replace the analytic coverage argument.

| Quantum SWAP return channel L_lambda | Rejections / 1,000 | Pointwise 95% interval |
|---|---:|---:|
| lambda = 0.40 | 0 | 0–0.0037 |
| lambda = 0.75 | 0 | 0–0.0037 |
| lambda = 0.90 | 18 | 0.0107–0.0283 |
| lambda = 0.97 | 981 | 0.9705–0.9885 |
| lambda = 1.00 | 1,000 | 0.9963–1 |

The conservative calibration bound has poor power for less-than-nearly-perfect
recovery at this shot budget. These are not power estimates for NMN.

Two constructed failures are kept visible: calibrating a reset but using
identity during dynamics creates a false memory claim if mismatch is ignored;
and an averaged depolarizing channel can transmit perfectly when an unmodeled
classical Pauli correction register is retained. The latter is why calibration
of the averaged map cannot justify applying this result to general flagged
instruments. Leakage and intervention–environment coupling are not covered.

## Prior art and scientific increment

See the [claim-by-claim matrix](prior-art.md). This work does not claim novelty
for quantum-memory witnesses, self-consistent tomography, the six-state EB
benchmark, or diamond-norm robustness. Its bounded contribution is the pinned
count/regrouping audit and the explicit executable prospective bridge,
source-removal constructions, and failure controls. It has not achieved the
requested new empirical certificate or established a publication-level result.

## Reproduction

From the repository root:

```bash
python -m pip install -r studies/calibrated-quantum-memory/requirements.txt
python -m unittest discover -s tests -p 'test_calibrated_quantum_memory.py' -v
python studies/calibrated-quantum-memory/certificate.py --check
python studies/calibrated-quantum-memory/audit.py --source-dir /tmp/calibrated-memory-nmn --download --check
python scripts/check_repository.py
```

The last analysis command downloads only five pinned files, about 2 MB, to an
external directory; [sources.json](sources.json) contains URLs, byte sizes,
SHA-256 hashes and licensing. `--check` compares the regenerated outputs
byte-for-byte. Reusing an existing full NMN clone is also supported. Raw data,
source notebooks, paper PDFs, and other archives are not committed. A dedicated
workflow performs source retrieval; ordinary PR checks use synthetic data.

Outputs: [measured conditional audit](results/nmn-audit.json) and
[synthetic certificate and controls](results/synthetic-certificate.json).
No author contact, new acquisition or merge was performed. The empirical
mission remains open: identify and validate a matched acquisition with
intervention calibration, then test the full declared error model.
