# Born-rule identifiability across measured datasets

**Main measured result:** the December 2022 Nairobi preparation scan and October
2023 Viviani acquisition supply different constraints on a common probability
deformation. Both full mean tables also admit explicit ordinary qutrit models
with small preparation-dependent leakage. No Born-rule violation is identified.
See the [measured joint report](results/joint-table-report.md) and
[observation models, proofs and prior-art comparison](joint-tables.md).

The earlier [rotation/portfolio report](results/combination-report.md) remains a
secondary feasibility result: its exact common angle warp explains why combining
quarter-turn depths cannot identify the deformation. The new measured combination
replaces that uninformative pairing as the main analysis.

The optical runner-up is also complete: The waveguide archive's nonzero Sorkin means
admit a specified ordinary optical model and an indistinguishable additive
all-open intensity diagnostic. This is mean-table compatibility, not a full
acquisition likelihood fit or Born-rule violation. Combining temperatures does
not exclude ordinary theory with separate 0.03-radian phase regions.

A proposed extra all-open exposure, advancing one path by π, separates the
constructed pair under explicit calibration. Phase cycling and the underlying
algebra are known; no new fundamental theorem or priority is claimed.

- [Generated optical report](results/report.md)
- [Exact derivation](derivation.md) and [theory dossier](theory.md)
- [Expanded dataset combinations](combinations.md)
- [Literature](literature.md), [sources](sources.md), [data dictionary](data-dictionary.md)
- [Manifest](manifest.json), [protocol](protocol.json), [sufficiency](sufficiency.json)

## Reproduce

From the repository root, Python 3.12:

```sh
python -m pip install -r studies/born-rule-identifiability/requirements.txt
python studies/born-rule-identifiability/fetch_data.py --cache /tmp/born-rule-sources
python studies/born-rule-identifiability/analyze.py --cache /tmp/born-rule-sources --check
python studies/born-rule-identifiability/fetch_data.py --cache /tmp/born-rule-sources --portfolio
python studies/born-rule-identifiability/combinations.py --cache /tmp/born-rule-sources --check
python studies/born-rule-identifiability/fetch_data.py --cache /tmp/born-rule-sources --tables
python studies/born-rule-identifiability/joint_tables.py --cache /tmp/born-rule-sources --check
# Refit the nonconvex models and rerun count simulations (several minutes):
python studies/born-rule-identifiability/joint_tables.py --cache /tmp/born-rule-sources --refit --output-dir /tmp/born-table-refit
python -m unittest discover -s tests -p 'test_born_rule_*.py' -v
```

Sources must be outside the checkout; all originals are SHA-256 pinned and
excluded from git. `analyze.py --check` without `--cache` regenerates the
constructions, designs, simulations and report from the committed aggregate
summary but does not reread originals. Use `--cache` for source-to-report
verification. Without `--check`, outputs go to `results/` or `--output-dir`.
Ordinary CI is offline for this study; manual source reproduction is separate.

Unresolved: Peres mean differences around .0001, original autocorrelation
processing, absent gain/response calibration records, unmeasured leakage phases,
and no independently calibrated antipodal exposure. Published calibration
statements are not numerical calibration records. No source author was contacted.
The conditional optical result does not warrant an elementary-only Lean PR.

The table `--check` path verifies saved physical constructions, analytic bounds,
source reduction and report. It does not claim optimizer identity. `--refit`
repeats the full exploratory optimization; different equally feasible nuisance
parameters need not match bit-for-bit. All probability predictions are reviewable.
