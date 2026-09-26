# Records needed for an end-to-end reconstruction

The current deposit supports a descriptive Figure 3/4 replication. The following
records would allow the missing steps to be checked. This is a technical input
specification, not a statement that the records do not exist elsewhere.

| Record | Minimum information | Questions affected |
|---|---|---|
| Complete complex propagator estimates | Explicit repeat, slice, input-mode and output-mode indices; real/imaginary components; units; slice coordinates; missing-value mask; uncertainty definition | Q04, Q05, Q14 |
| Pointer and reference acquisitions | Original frames or integer counts; polarization setting; exposure/gain; dark/background frames; reference amplitudes/phases; detector calibration; matching trial/repeat identifiers | Q03, Q04, Q08 |
| Independent endpoint acquisition | Original counts/frames per repeat; preparation identifier; settings; complete successes/failures or a documented calibration of their probabilities | Q04, Q08 |
| Acquisition dependence | Timestamp/order/block identifiers, shared references and normalizations, paired repeats or a justified full covariance model | Q04, Q05, Q14 |
| Figure 3 evaluation | Exact E, Q and theoretical arrays used for 4.45%; comparison identity; coordinate matching; normalizations; averaging order; executed code/version | Q06 |
| Fidelity evaluation | Actual complex path vectors or reproducible generation inputs; normalization convention; exact implemented expression; code/version | Q07, Q09 |
| Spatial and aperture model | Camera-to-object mapping, centers and edges, mode profiles/overlap, apertures at each plane, integration weights, retained/omitted amplitudes and discretization controls | Q10, Q17 |
| Action/length grouping | Individual labels or a reproducible path-index rule; bin boundaries and endpoint treatment; counts per bin; phase reference/wrapping; magnitude weighting; arithmetic/circular averaging | Q11, Q20 |
| Noise and summary generation | Noise distribution and cross-entry/repeat correlations, seed if reproducibility is intended, amplitude/phase normalization, histogram and simulation plotting code | Q14, Q19, Q21 |

Self-describing array files plus a machine-readable metadata table would suffice;
a particular proprietary format is not required. Every array needs axis labels,
units and a source identifier. Original frames should be linked to derived K
entries so shared references remain visible. Independent repeat counts are distinct from path counts or pixels.

An ingestion check would first validate file hashes and dimensions, compare the
mapping against the published method, and reproduce the already deposited
summaries. Only then could covariance-preserving resampling or a justified
likelihood be used. Any new uncertainty analysis would state its trial and
calibration assumptions. The existing manifest continues to identify the same
version-3 files; new material would have separate provenance.

[Exact synthetic witnesses](results/identifiability.json) show why marginal
means/variances, separate magnitude/phase summaries and success-conditioned
shapes do not generally identify these records. More decimal places in the
current outputs would not repair that loss of information.
