# Collapse compatibility: source audit and spectral feasibility

This study reproduces sodium-cluster interference and layered-force-sensor
baselines, audits an additional XENONnT event release and a LISA Pathfinder archive
probe, and implements two cross-experiment feasibility calculations. **It does
not establish a new exclusion of objective-collapse models.**

The strongest completed result is a conditional limitation: a finite collection
of positive-frequency spectral envelopes cannot bound finite-time stationary
dephasing over an unrestricted OU correlation-time family. A constructive sequence
keeps the declared macroscopic suppression fixed while the measured-frequency
responses tend to zero. Applying that sequence to an actual universal collapse
field additionally requires the motion, noise-frame and instrument-response
premises listed in [the derivation](derivation.md). Those premises are not proved
by the datasets. Novelty is not established.

The original interference/sensor pairing remains valuable as a baseline. It is
not converted into an arbitrary-spectrum LP: mass/velocity mixing is nonlinear,
the author posterior fixes empirical calibration, and a colored Talbot-Lau map
has not been established here.

## Reproduce

Use Python 3.12 and install `requirements.txt`. Keep originals outside this repo.

```sh
python studies/collapse-compatibility/fetch_data.py --data-dir /tmp/collapse-sources
python studies/collapse-compatibility/analyze.py --data-dir /tmp/collapse-sources
python studies/collapse-compatibility/analyze.py --data-dir /tmp/collapse-sources --check
python studies/collapse-compatibility/analyze.py --check
python -m unittest discover -s tests -p 'test_collapse_*.py' -v
```

The full-source command verifies five downloaded objects, reproduces fits and
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
