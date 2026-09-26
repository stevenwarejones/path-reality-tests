# Born-rule identifiability against optical imperfections

**Limited result; draft research.** The expanded search now has a stronger
single-qubit result: one ordinary control-axis warp exactly mimics a normalized
probability deformation at every integer gate depth. See the
[main combination report](results/combination-report.md).

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
