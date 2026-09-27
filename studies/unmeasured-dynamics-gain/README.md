# Can combined records predict a new quantum experiment?

**A certified one-sided result, with the requested two-sided result still open.** For the unmeasured word `Gi Gi Gy Gy`, every stationary qubit CPTP model in the unchanged simultaneous GST/RB confidence region has plus-outcome probability at least **0.1113640705**. An explicit physical model satisfying **all** RB constraints predicts at most **0.002342601**. Adding GST therefore removes physically attainable RB-only predictions. We have not found a model satisfying every GST constraint that crosses the same bound. This is not a certificate that RB contributes indispensable information to this prediction.

The scientific question remains whether shared ordinary dynamics across complementary measured families restrict an informative unmeasured experiment more than either complete source does alone. The result below advances that question in one direction; it does not establish the requested combined-source payoff or novelty.

| Claim | Status |
| --- | --- |
| GST and RB counts and primitive words from the same reported interleaved March 2015 ion acquisition | Measured; source matching inherited and independently reconstructed |
| `Gi Gi Gy Gy` absent from both complete expanded circuit sets | Verified from all 6,661 word hashes and original files |
| Lower bound 0.1113640705 over **all** stationary qubit states, effects and CPTP gates satisfying the joint region | Certified analytic implication; statistical coverage conditional on independent shots |
| Physical RB-only point below the bound, and nonempty joint region | Certified, including every retained cell/group and independent full-source propagation |
| Physical GST-only point below the same bound | Unavailable; two local searches failed to find one |
| Detailed counts add information beyond two overall means | Certified physical aggregate-only counterexample |
| All-word-equivalent physical models disagree on entanglement-breaking status of `Gy^10749` | Certified restricted operational identifiability obstruction; a small boundary crossing |
| Useful two-source restriction on this resource, novelty, or major scientific advance | Not established |

This builds on [PR17's study](../gst-rb-combination-gain/README.md), preserving its measured-context mixture result. The new target is a literal unrecorded sequence, not a relabeling of recorded frequencies. In the nominal apparatus interpretation it is two idles followed by a two-Y readout rotation; no ideal rotation assumption is used in the proof.

The [two-sided follow-up](followup.md) now ranks targets by both deletion margins and proves a limitation of the entire transfer-inequality method using all GST constraints. It does not supply the missing physical GST-only witness or establish that GST alone implies the present joint bound.

## Decisive comparison

All predictions use the same model class and the same original 95% confidence construction. Removing a source deletes all its constraints without reallocating statistical budgets. Intervals below include a common rigorous `1e-9` propagation envelope and outward quantization to `1e-8`.

| Region / physical point | Certified target interval or bound | Minimum retained inequality slack | Interpretation |
| --- | --- | --- | --- |
| Full joint class | **[0.1113640705, 1]** (outer interval) | — | Global lower bound; not an optimal interval |
| Joint feasible point | [0.989130599, 0.989130611] | GST 0.0031338120; RB 0.0020241373 | Prevents vacuity; not a validation experiment |
| Complete RB-only witness | [0.002342589, 0.002342601] | 0.0000999990 | Violates joint lower bound by **at least 0.1090214695** |
| Existing complete GST-only point | [0.989071589, 0.989071601] | 0.0031535448 | Does not separate |
| GST search, seed 0 | [0.771414139, 0.771414151] | 0.00000111086199 | Feasible construction, no certified optimum |
| GST search, seed 260927 | [0.771559299, 0.771559311] | 0.00000107145458 | Feasible construction, no certified optimum |
| Only the two overall means, physical model | [0.059440759, 0.059440771] | GST 0.0038399730; RB 0.0024689032 | Violates joint lower bound by **at least 0.0519232995** |

The certificate checks 4,657 GST cells plus six GST groups, and 2,004 RB cells plus six RB groups. For the aggregate control only the two explicitly declared overall means are retained. It is not a source-only witness. All other relevant retained constraints are checked, including the original long GST sequences. Joint compatibility with these conservative intervals is not a likelihood goodness-of-fit claim and does not undo the narrower-model validation failures in the earlier study.

The analytic bound uses just three of the joint region's rows: empty word (0/50), `GiGi` (RB 0/285), and `GyGy` (GST 50/50). It says that, in dimension two with shared preparation/readout, two low-probability preparations cannot become arbitrarily distinguishable under the same suffix. The complete RB witness shows that RB alone leaves the readout-flip response unconstrained enough to evade this prediction. **Which ambiguity RB removes beyond complete GST is unresolved.** Merely observing an RB row in the proof is not evidence of indispensable information.

## What the controls say

* **Remove composition:** keeping the measured probabilities at the joint point while assigning the new word any probability in [0,1] satisfies all count constraints in the independent-word relaxation. Shared qubit dynamics, not a direct frequency constraint on the new word, supplies the bound.
* **Keep only aggregate means:** an explicit stationary qubit gate set matches both overall centers but predicts approximately 0.05944. Detailed circuit information matters. The construction uses opposite small offsets about a common rotation axis, not the intended ideal gate geometry.
* **Allow a third level:** a stationary qutrit cyclic-shift example violates the composition inequality. This is an assumption counterexample, not a qutrit fit to the full acquisition or a measured leakage bound.
* **Allow bounded context variation:** hypothetical half-diamond deviations of 0, 0.0001, 0.001 and 0.01 per gate give lower bounds approximately 0.111364, 0.110402, 0.101780 and 0.018862. No such drift cap is measured here. See [derivation](derivation.md).
* **Allow unrestricted shot dependence:** 95% coverage is lost. Algebraic implications of confidence-region membership remain, but the observations no longer justify that membership with the stated coverage. Recorded row order is preserved; it is not a timestamp or evidence of shot independence.
* **Prospective decision:** conditional on the training region containing the true stationary qubit model, 41 independent new shots and rejection on at most one plus outcome have size below 0.048497 and power above 0.995765 against the displayed RB-only point. This is a designed discrimination experiment, not measured validation. Combining 95% training coverage with 5% conditional test size gives a 10% unconditional union bound, not a 5% one.

## A precisely delimited resource obstruction

Let a gate set have Bloch matrices `G`, state `rho` and terminal effect `E`. The transformation `G -> S G S^-1`, `rho -> S rho`, `E -> E S^-1`, with `S=diag(1,s,s,s)`, preserves every terminal word probability exactly. For this acquisition's joint-compatible point, both `s=1` and `s=1009/1000` yield physical gate sets. Their native unnormalized Choi matrices have certified minimum eigenvalue at least `1e-6`.

For the unmeasured block `Gy^10749`, however, a trusted Bell-pair input would yield a normalized Choi state whose partial transpose is at least `3e-6 I` at the first point, and has a Rayleigh quotient below `-3e-6` at the second. In two-by-two dimensions this certifies respectively entanglement-breaking and non-entanglement-breaking behavior. This is a difference in an **additional operational task with a fixed external reference**, not an observable distinction between two representations of the archived binary experiment.

No amount of additional terminal-word counting of the archived type can identify that binary property at these two realizations. This does **not** prove that the joint data cannot tighten a resource interval, or that all useful resource boundaries are impossible. The crossing is small and deliberately selected near the entanglement-breaking boundary; we have not demonstrated a robust practical memory advantage. The block is 31% longer than the longest deposited word, with no available conversion to calibrated storage time.

An independently trusted Z measurement of the reset preparation distinguishes this particular scale orbit by a plus-probability gap between 0.0044564 and 0.0044565. Such a calibration would enter as its own measured likelihood with a declared trusted measurement, not as a fixed nuisance value inferred from a Born-based GST fit. It removes this scale ambiguity, not every gate-set gauge. Alternatively, a trusted entangled input and reference-assisted output measurement directly addresses the resource task. Neither calibrated data product has been verified in the matched archive.

## Limits and next concrete step

The missing mathematical object for the primary target is a fully GST-compatible physical model with `p(GiGiGyGy) < 0.1113640705`, or a different unmeasured target with a certified joint interval and **both** complete-source crossing models. The local GST searches reaching 0.771 supply neither a lower bound on that optimum nor evidence that certification is impossible. The simple geometry relaxation ignores almost all GST composition constraints and is likely too loose; that is a limitation of this method.

For the resource task, terminal observations cannot break the displayed exact equivalence. A matched trusted-state/readout or reference-assisted calibration is a concrete missing observation. The [source and prior-art audit](audit.md) distinguishes the public deposits actually reconstructed from alternative experimental leads for which calibration likelihoods were not verified. No unrelated hardware calibration was imported.

## Reproduce

Python 3.12, a C++11 compiler, and [the inherited pinned dependencies](../gst-rb-combination-gain/requirements.txt) suffice. No acquisition script is executed; raw archives stay outside the repository.

```sh
python -m pip install -r studies/gst-rb-combination-gain/requirements.txt
OPENBLAS_NUM_THREADS=1 python studies/unmeasured-dynamics-gain/verify.py --check
python studies/unmeasured-dynamics-gain/controls.py --check
python -m unittest discover -s tests -p test_unmeasured_dynamics_gain.py -v
OPENBLAS_NUM_THREADS=1 python studies/unmeasured-dynamics-gain/verify.py --cache /tmp/unmeasured-source --download --check
python studies/unmeasured-dynamics-gain/aggregate_discovery.py --cache /tmp/unmeasured-source
```

The offline route verifies saved count probabilities, physicality, analytic bounds and the resource certificate. The **separate full-source route is essential**: it reconstructs every primitive word and denominator, binds all saved probabilities to their physical parameters, cross-checks density-matrix propagation, verifies target absence and the entire screen hash, and expands every alternate Clifford-format RB row. Both routes run in CI. To reconstruct the optional literature supplement audit, run `python studies/unmeasured-dynamics-gain/source_audit.py --cache /tmp/unmeasured-source --download --literature`; this PDF is not part of the likelihood. Discovery optimizers are not rerun by certification:

```sh
python studies/unmeasured-dynamics-gain/discover.py --cache /tmp/unmeasured-source
python studies/unmeasured-dynamics-gain/discover.py --cache /tmp/unmeasured-source --source GST --seed 0 --maxiter 80
python studies/unmeasured-dynamics-gain/discover.py --cache /tmp/unmeasured-source --source GST --seed 260927 --maxiter 80
```

[Protocol](protocol.json), [discovery history](discovery-history.md), [derivation](derivation.md), [machine certificate](results/certificate.json), [resource certificate](results/resource-certificate.json), and [controls](results/controls.json) separate retrospective discovery from deterministic verification. No source or target was reserved for confirmatory validation in this follow-up.
