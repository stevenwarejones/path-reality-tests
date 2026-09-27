# Collapse compatibility: source audit and spectral feasibility

This study reproduces sodium-cluster interference and layered-force-sensor
baselines, audits an additional XENONnT event release and a LISA Pathfinder archive
probe, and implements two cross-experiment feasibility calculations. **It does
not establish a new exclusion of objective-collapse models.**

The [moving-field report](results/moving-report.md) extends the stationary
comparison with one shared noise frame, the same-frame moving benchmark,
archived mass/velocity/optics, JPL trajectory metadata, published acquisition
windows, and conservative assembly/orientation bounds. Motion produces a finite
response; a conditional ensemble-dephasing survivor remains. This is not a
physical outcome-selection model or a calibrated joint experimental acceptance.

Read [the moving derivation and acquisition audit](moving-frame.md) for the
precise local-Markov approximation and unresolved empirical limitations. The
original stationary calculations remain as explicit comparison baselines.

## Reproduce

Use Python 3.12 and install `requirements.txt`. Keep originals outside this repo.

```sh
python studies/collapse-compatibility/fetch_data.py --data-dir /tmp/collapse-sources
python studies/collapse-compatibility/moving_analysis.py --data-dir /tmp/collapse-sources
python studies/collapse-compatibility/analyze.py --data-dir /tmp/collapse-sources
python studies/collapse-compatibility/analyze.py --data-dir /tmp/collapse-sources --check
python studies/collapse-compatibility/analyze.py --check
python studies/collapse-compatibility/moving_analysis.py --check
python -m unittest discover -s tests -p 'test_collapse_*.py' -v
```

The full-source command verifies eight external objects (including normalized Horizons responses), reproduces fits and
Bayesian updating, compares 45 evaluations against the inspected author module,
and regenerates derived artifacts. Offline checks use the committed small
response summaries, re-evaluate theoretical calculations and certificates, and
check the source/configuration/code hashes. They do not pretend to reread originals.

Read [the generated report](results/report.md), [candidate inventory](combinations.md),
[source sufficiency](data-dictionary.md), [literature audit](literature.md), and
[derivation](derivation.md). All results associate with input/code hashes rather
than a circular commit or report hash. The ordinary CI is offline; a separate
manual workflow tests complete public-source reproduction.

No Lean companion is opened: the current general result is elementary conditional
spectral mathematics and does not justify a separate formal PR.
