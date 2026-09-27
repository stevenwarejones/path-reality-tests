# Certified GST/RB combination gain

**Question:** Can the two measured families restrict predictions of shared ordinary
quantum operations more than either family alone?

**Yes, for one declared mixture of recorded contexts.** An analytic bound covers
all stationary qubit CPTP gate sets with arbitrary shared physical preparation
and binary readout. Explicit physical models satisfying every GST-only constraint
and every RB-only constraint separately exceed that joint bound. Another physical
model satisfies the complete joint region. This replaces finite-model-bank spread
with a strict, nonempty, full-class feasible-set comparison.

The observable is q=(1+μ_RB−μ_GST)/2, where each μ is the shot-weighted mean
probability over all recorded words in that family. Operationally choose a family
with probability 1/2, sample one of its recorded words with exposure weights,
and score plus for RB or minus for GST. This is **not an unseen-sequence
extrapolation**, and the bound itself is an elementary aggregate-count bound.
Novelty and a consequential dynamics-identification result remain separate goals.

| Constraint set / explicit model | Certified prediction interval for q | Smallest retained constraint slack, after numerical allowance |
|---|---:|---:|
| Joint: bound for the entire class | q ≤ 0.287814450898102 (rounded upward) | All 6,661 cells and 12 groups imposed |
| Joint-compatible point | [0.284084547181401, 0.284084549181402] | 0.0020241373 |
| GST-only witness | [0.291840997887416, 0.291840999887417] | 0.0031535448 |
| RB-only witness | [0.504854184149459, 0.504854186149460] | 0.0000999989 |

The GST-only prediction exceeds the joint upper bound by more than **0.0040265469**;
the RB-only prediction exceeds it by more than **0.2170397332**. The common
numerical allowance is 10⁻⁹ per probability or contrast. All strict claims use
exact rational endpoints in the [certificate](results/certificate.json), not
these rounded display values. No optimizer failure or source-removed relaxation
is used to establish either inequality.

| Claim | Status |
|---|---|
| 4,657 GST rows / 232,850 shots and 2,004 RB rows / 563,280 shots | Measured; immutable sources and parser inherited from PR15 |
| Shared interleaved acquisition | Supported by original paper; timestamps absent from the count files |
| Joint upper bound over all stationary qubit CPTP models | Certified conditional on the declared count region; need not be optimal |
| Both source-removal regions contain physical points above that bound | Certified by explicit normalized Kraus constructions and exact constraint checks |
| Joint region nonempty | Certified by a third physical point |
| At least 95% simultaneous confidence | Conditional on independent Bernoulli shots, fixed per-word probabilities and declared budgets |
| Complete likelihood adequacy, held-out generalization, resource identification | Not established |
| Novel scientific theorem or new GST/RB complementarity discovery | Not claimed |

## Sources and assumptions

This is a focused follow-up on current main `db1a02d`; PR11 and PR15 are preserved.
The [frozen retrospective protocol](protocol.json) fixes the observable, class,
confidence construction, length groups and deletion rule. The [derivation](theory.md)
gives the complete proof, operational semantics and numerical error budget.

Use the existing [immutable manifest](../shared-quantum-dynamics/manifest.json),
[strict parser/downloader](../shared-quantum-dynamics/sqd_sources.py) and
[source audit](../shared-quantum-dynamics/sources.md). The selected deposit is
`pyGSTio/supplemental-info-arXiv-1605.07674` at
`fc939f90f7a2b32a20aec34f36465259f38837a3`. GST SHA-256 is
`5ca509cb8e6cfe8996b8a1339457bc7f96d4efe0a3bdb8b550c3c496ce53f686`;
RB SHA-256 is
`95ec756e5d4ead57ac6ba67b9302f11badc1dc72b956f0c4cbc585c3e47026c9`.
The manifest also pins the provenance and circuit-convention files.

Both records use Gi/Gx/Gy primitive words, one preparation, and terminal binary
counts. GST denominators are 50; RB denominators are 270 or 285. The same three
maps, state and effect are shared within the reported 2015-03-30 trapped-ion
acquisition. File order is preserved, including duplicate RB words. The four-Y
word independently occurs as GST 0/50 and RB 3/285. No timestamps, per-shot order,
reset logs, independent readout calibration or leakage monitor are supplied.
This study adds no apparatus-sharing bridge to Nairobi's different epochs.

The ambiguity removed is observable, not representational: GST alone allows a
shared gate set with substantially larger randomized-family response; RB alone
allows an identity-like gate set that fails the GST responses. Adding the other
family rules out those points and every other point above the joint bound.
Stationarity and shared SPAM are required for interpreting them as one apparatus.
Independent shots are required for the confidence label. Violations of either
assumption would refute that interpretation without refuting quantum mechanics.

## Reproduction

Python 3.12, a C++ compiler, and the pinned dependencies are sufficient. No
acquisition scripts, downloaded notebooks or pickles are executed. Raw sources
remain in an external cache.

```bash
python -m pip install -r studies/gst-rb-combination-gain/requirements.txt
python studies/gst-rb-combination-gain/verify.py --check
python -m unittest discover -s tests -p 'test_gst_rb_combination_gain.py' -v
python studies/gst-rb-combination-gain/verify.py --cache /tmp/gst-rb-gain --download --check
```

The offline command rebuilds exact confidence endpoints, interval physicality and
roundoff certificates, and every retained inequality from committed counts and
propagation artifacts. The full-source command hash-checks all files, reconstructs
every word/count, and independently propagates all three models using both Bloch
prefixes and complex density superoperators. Their observed maximum differences
were 1.48×10⁻¹² (GST witness), 5.46×10⁻¹³ (RB witness), and 2.37×10⁻¹² (joint point).
The rigorous error allowance does not depend on treating those empirical
differences as error bounds.

Optional witness reconstruction, separate from certificate verification:

```bash
python studies/gst-rb-combination-gain/search.py --cache /tmp/gst-rb-gain --output /tmp/reconstructed-witnesses.json
```

The saved models are in [witnesses.json](results/witnesses.json). The joint seed
comes from PR15's all-length general CPTP fit. GST discovery scanned depolarizing
weights 0–3×10⁻⁵; RB discovery scanned identity/depolarizing weights 6×10⁻⁵–1.4×10⁻⁴
and solved two readout eigenvalues by LP. An additional GST readout LP candidate
was explored but is not needed for this certificate. The retained witness uses
the unchanged joint seed's SPAM. RB's inward LP margin was increased to 10⁻⁴
before verification. These are existence searches; no optimized supremum,
source-only training procedure, or untouched validation claim is made.

## Relation to previous results and remaining goal

The original experiment already compared GST and RB. The prior
[claim-by-claim audit](../shared-quantum-dynamics/prior-art.md) remains applicable.

| Closest result | Increment here | Limit |
|---|---|---|
| Original Blume-Kohout et al. GST/RB comparison | Explicit simultaneous region and two physical source-removal witnesses | Pairing and qualitative complementarity are not new |
| PR15 finite-bank source removal | Full-class analytic upper bound, with actual feasible points on both source-removed sides | One deliberately recorded-context contrast; no full prediction-region reconstruction |
| PR15 all-length general CPTP fit | Independently enclosed physical point certifies this region is nonempty | Region acceptance does not establish global residual adequacy |
| PR15 leakage/readout equivalence and leakage-identification literature | No new leakage claim | No population, return-rate or memory lower bound |

The immediate requested combination-gain criterion is met within the stated
region and class. The more ambitious unresolved problem is an informative bound
for an **unmeasured primitive sequence**, with physical source-removal attainers,
or an operational resource boundary robust to justified dependence/drift.
This certificate does not establish that such a result is impossible. The
required next ingredient is a nontrivial dynamics constraint for that target,
rather than its own measured aggregate. No companion formal PR is proposed for
the elementary linear-bound algebra.
