# Dynamics after removing detector contrast

PR21 completed a two-source operational discrimination-resource result. Its scope and certificates are preserved. This follow-up asks whether the combination identifies **dynamics**, rather than detector spectral range, and whether the exhibited physical models are adequate aggregate explanations.

**The new dynamics endpoint is not established.** The implemented results identify three specific obstacles: detector retuning cancels exactly from the proposed normalized resource; a stronger finite-sample aggregate test excludes the saved PR21 constructions; and the matched repeated-word records alone cannot certify an operator-level stability allowance. None is a theorem excluding the full stationary CPTP model class or proving that dynamics combination gain is impossible.

| Claim | Status |
| --- | --- |
| PR21 two-source discrimination-resource separation | Completed baseline, unchanged |
| Detector-normalized ratio removes affine detector bias/scale for every sequence | Exact identity, including nonunital channels |
| Previous GST contrast witness can separate this ratio from its joint seed | Impossible: ratios are exactly equal for every word |
| New normalized two-sided gain | Not established; no crossings in the inherited 59 positive-bound candidates |
| Saved joint and two contrast witnesses pass stronger aggregate adequacy check | No; finite-binomial interval certificates reject these particular tables on their retained records |
| General stationary CPTP class is excluded | **Not established**; local refitting is not a global certificate |
| Matched timestamps or feedback/calibration records establish a context allowance | Unavailable in the audited deposit |
| Eight common words certify detector stability | No; an explicit two-context detector control fits all 23 corresponding count cells with operator-norm difference one |
| New publication-level dynamics inference or novelty | Not claimed |

## A target that removes detector scale

For a nonconstant binary effect `E`, define

\[
R_T=\frac{\operatorname{diam}(G_T^*(E))}{\operatorname{diam}(E)}\in[0,1].
\]

This is the fraction of the detector's maximum contrast retained after a sequence. In Bloch coordinates, if `E=aI+e.sigma` and the channel's linear part is `A_T`, then `R_T=||A_T^T e||/||e||`. Channel translations and detector bias cancel. Detector orientation remains a nuisance. The ratio is undefined for a constant detector; the positive observed response range rules that case out in the joint region.

For any physically admissible affine detector change `E'=aI+bE`, with `b` nonzero, trace preservation gives

\[
G_T^*(E')=aI+bG_T^*(E),\qquad
\operatorname{diam}(E')=|b|\operatorname{diam}(E).
\]

Thus **`R_T(E')=R_T(E)` for every finite word**, without assuming unital channels, isotropic noise or ideal gates. The old GST contrast witness differs from the joint seed only by this affine detector transformation. Exact rational factors verify the identity, so this witness construction cannot supply a normalized-dynamics crossing for any target, regardless of further search or shots.

For PR21's ten-gate suffix, the normalized contrasts are approximately:

| Physical construction | Normalized contrast |
| --- | --- |
| Joint seed | 0.9989598102 |
| GST-only detector-retuned witness | Exactly the same as the joint seed |
| RB-only detector-retuned witness | 0.9986008817 |

The directed intervals and exact affine coefficients are in [normalized.json](results/normalized.json). All 59 earlier positive response-range bounds were checked against these physical points; none crosses on both sides. This is not a search over the full physical region.

The ratio removes one nuisance but is not a general gate-error metric. Every unitary sequence has `R_T=1`, even when coherent errors radically change a fixed-reset prediction. It also need not be invariant under arbitrary nonunitary gate-set gauges. Any inference must range over all physical realizations or use independently fixed input calibration; a fitted value alone is not an observed resource.

## Why response range alone is insufficient for the next goal

The earlier two rows imply `R_T>=0.7229749561`, using `diam(E)<=1`. This is still a response-range argument and does not meet the new standard.

It is tight in the relaxation that keeps an arbitrary CPTP block, effect and two unknown input states but drops the repeated native-gate composition. Set `E=diag(0,1)` and use a measure-and-prepare channel that sends input `|0>` to a diagonal state with plus probability 0.0452128544, and input `|1>` to one with plus probability 0.7681878105. Its dual effect is `diag(0.0452128544,0.7681878105)`, so the normalized contrast is exactly 0.7229749561. Explicit Kraus tests verify complete positivity and trace preservation.

This is a block-level tightness witness, **not** a native gate-set model satisfying every source constraint. It shows precisely what additional information the next certificate must use: shared composition across more recorded circuits. The old all-word transfer-relaxation obstruction remains applicable to the earlier fixed-input target. More instances of that inequality cannot restore the information it omits.

## Aggregate adequacy: exact finite-sample test

Cell-region compatibility is unchanged. A separate diagnostic now uses the complete product-binomial sampling model. For counts `K_i~Bin(n_i,p_i)`, define the saturated binomial deviance

\[
D_p(K)=\sum_i2\left[K_i\log\frac{K_i}{n_ip_i}+(n_i-K_i)\log\frac{n_i-K_i}{n_i(1-p_i)}\right],
\]

with zero-count terms defined by continuity. For a fixed candidate `p` and any saved rational `0<t<1/2`, independence and Markov's inequality give

\[
\Pr_p\{D_p(K)\ge d\}\le\min\left(1,e^{-td}\prod_i M_i(t,p_i)\right),
\]

where the **finite** sum is

\[
M_i(t,p)=\sum_{k=0}^{n_i}{n_i\choose k}
(k/n_i)^{2tk}(1-k/n_i)^{2t(n_i-k)}
p^{(1-2t)k}(1-p)^{(1-2t)(n_i-k)}.
\]

No chi-square approximation, estimated degrees of freedom or known-nuisance recovery argument is used. Positive-coefficient Horner evaluation with 50-digit directed interval arithmetic encloses the sums. The actual physical probabilities lie within the inherited independently checked `1e-9` envelope of each saved value; the interval score and moment bounds cover that whole envelope. Native physical factors are checked and full-source reproduction independently binds every probability to density-matrix propagation. Tail bounds are rounded **up** on a `1e-12` grid, including extreme underflow cases.

Candidate-specific tilts were selected retrospectively. This remains a valid candidatewise tail bound: each Chernoff expression bounds the same exact survival probability at the observed score. Inverting these tests defines a pointwise confidence region; searching a nonconvex physical class is not a certificate that this region is empty. No statement below excludes all stationary qubit dynamics.

| Candidate | Retained records | Deviance (approx.) | Certified tail upper bound |
| --- | --- | ---: | ---: |
| PR21 joint-compatible model | All 6,661 rows | 7288.335 | 0.000057704437 |
| PR21 GST contrast witness | All 4,657 GST rows | 22368.973 | 0.000000000001 |
| PR21 RB contrast witness | All 2,004 RB rows | 2534.809 | 0.000000000003 |
| New local general-CPTP refit | All 6,661 rows | 7271.541 | 0.000114550079 |

All four are rejected at the declared candidatewise cutoff **0.025**. The [generated adequacy artifact](results/adequacy.json) also reports a new general-CPTP refit, its exact diagnostic and original-region membership. The stronger test does **not** retroactively change PR21's simultaneous 95% cell/group region or its valid separation theorem. Intersecting that old region with this separate 97.5% region would guarantee only at least 92.5% coverage by the union bound, not 95%; no such stronger joint-region certificate is claimed here.

The new numerical refit holds only three redundant factor scales fixed and varies the remaining general CPTP/state/effect coordinates. It uses all recorded sequences. The final saved attempt reaches deviance approximately 7271.541 after 150 evaluations, without convergence; the earlier 120-evaluation attempt using coarser finite differences reached approximately 7273.327. These improve the point objective but provide no global optimization result. The final point also fails the aggregate diagnostic. It is retained as a reproducible local attempt, not a stationary-class exclusion or a successful explanation.

This test still assumes independent Bernoulli shots and independent rows. Under perfect within-row correlation, a count can be either zero or the entire denominator while individual shots retain the same marginal probability. A test explicitly demonstrates that the binomial guarantee then fails. No recorded chronology or block structure is invented to calibrate that dependence away.

## Shared-operation stability: what the matched records say

The pinned deposit inventory contains the March 30 GST and RB counts, their alternate Clifford encoding, analysis notebooks and earlier experimental epochs. Its Figure 5 notebook creates the alternating-error data with `generate_fake_data`; these are **simulations**, not additional measured calibration. The paper describes interleaved GST/RB and active pulse-time/frequency feedback, but the audited count files contain no timestamps, feedback outcomes or calibrated probe-state descriptions. Earlier “drift” files are different acquisition dates and cannot supply a numerical allowance for March 30.

Exactly **eight expanded words** occur in both selected families, spanning **23 rows**. The [stability audit](results/stability.json) preserves every index, count and hash. Their agreements test those responses, not an operator-norm bound on the detector or diamond-distance bound on a suffix.

An explicit restricted control makes the distinction sharp. Take shared identity native gates, reset `|0><0|`, and two family-specific effects

\[
E_{GST}=\operatorname{diag}(0.01,0),\qquad
E_{RB}=\operatorname{diag}(0.01,1).
\]

Every word then has probability 0.01 in either family. This satisfies **all 23 matched-word cell constraints** under the unchanged endpoints, with minimum slack 0.0093200413, while `||E_GST-E_RB||=1`. These are two physical detectors and common ordinary gates, not one operation per circuit or arbitrary memory. The construction does **not** fit the complete GST acquisition or source-wide groups. Its exact scope is the inability of these duplicate probes alone to certify detector stability.

A bound on response differences is not an upper bound on operation distance without a trusted, sufficiently well-conditioned informationally complete probe frame. The concrete missing observation is a time-linked set of calibrated preparations spanning the qubit operator space, measured through each relevant suffix/detector context, with counts and a declared sampling model. Matching feedback/control logs would also establish what stability law is justified. We found no such measured records in this pinned deposit; we do not claim they do not exist anywhere or that all combined-data stability inference is impossible.

## Prior art and next decision

The [original experiment](https://doi.org/10.1038/ncomms14485), especially “Quantifying non-Markovianity,” “Comparison to randomized benchmarking,” and “Experimental details,” already reported residual non-Markovianity, a GST/RB discrepancy and an alternating-error explanation. The new point rejections are not a discovery of those phenomena. The increment is finite-sample, interval-checked adequacy bounds for the specific constructions used in this repository, alongside a detector-normalization obstruction and a matched-probe stability audit. Chernoff bounds and affine-channel algebra are established methods; no novelty claim is attached.

The next physical inference needs all three ingredients together: a stronger composition-preserving certificate for the normalized ratio or a fixed-reset probability; credible candidate constructions checked against aggregate adequacy as well as cell constraints; and empirically justified context allowances. The current calculations do not supply that combination. A better fit alone will not be enough, and neither an unsuccessful local fit nor missing calibration metadata proves that the requested dynamics result is impossible.

## Reproduction

Use Python 3.12, a C++ compiler, and the unchanged pinned dependencies in `../gst-rb-combination-gain/requirements.txt`. All downloaded files stay in external caches; notebooks and serialized computed objects are never executed.

```sh
python studies/detector-normalized-dynamics/adequacy.py --check
python studies/detector-normalized-dynamics/normalized.py --check
python studies/detector-normalized-dynamics/stability.py --check
python -m unittest discover -s tests -p 'test_detector_normalized_dynamics.py' -v

python studies/shared-quantum-dynamics/sqd_sources.py --cache /tmp/dynamics-counts --download
python studies/detector-normalized-dynamics/adequacy.py --cache /tmp/dynamics-counts --check
python studies/detector-normalized-dynamics/stability.py --source-cache /tmp/dynamics-counts --audit-cache /tmp/dynamics-stability --download --check
```

The separate local-search route is `python studies/detector-normalized-dynamics/refit.py --cache /tmp/dynamics-counts`. Its saved physical point is certified by the full-source route; CI does not treat a repeated optimizer run as a proof. [Protocol](protocol.json) and [stability manifest](stability-manifest.json) pin the scope and audited bytes. The original count manifest, parsers, physical certificates and PR21 artifacts are reused without modification.
