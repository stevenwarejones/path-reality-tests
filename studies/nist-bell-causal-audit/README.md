# NIST Bell-record causal audit

This study reconstructs **107,109,596 paired records** from one public NIST
Bell-test run and asks what they can establish about operational influence.
It produces measured marginal contrasts and assumption-indexed bounds, with
explicit barriers to a spacelike causal interpretation.

The independent raw decoder reproduces every stored setting and the archived
click-update behavior. Preserving every eligible pulse bit corrects 105 Alice
and 94 Bob words, adding 31 and 27 primary-window clicks respectively. All
no-click rows remain. Ambiguous settings and timestamp excursions receive
worst-case score treatment. This is not a reanalysis of the published Bell
p-value or evidence of a signaling anomaly.

The [numerical report](results/report.md) gives both directions, each receiver
setting, chronological drift, pulse and phase-window sensitivity, and assumed
randomization/leakage budgets. With **assumed** ideal joint assignment and zero
ordinary leakage, the largest primary conditional upper limit is 0.2652%.
That covers an absolute average signed effect under the statistical and causal
premises; it is not an apparatus-certified bound or a bound on every trial's
absolute effect. Every displayed interval includes zero.

## What is established and what is missing

- [Reconstruction](reconstruction.md): full raw overlap, source reconciliation,
  multiclick correction, explicit excluded prefixes/suffixes and censored tails.
- [Statistical derivation](statistics.md): history-conditional assignment,
  receiver-specific score, simultaneous retrospective-interval coverage,
  arbitrary completion of unknown scores and nuisance sensitivity.
- [Machine-readable sufficiency ledger](sufficiency.json): required evidence,
  external premises and disabled causal/spacelike claim flags.
- [Source inventory](sources.md): exact files and documents inspected.

Remaining inputs for an apparatus-certified causal/spacelike interpretation are
an independent cross-site timing/endpoint calibration, evidence of complete
physical acquisition, justified history-conditional assignment calibration and
a transported ordinary-leakage budget. The raw files improve the reconstruction
substantially; they do not themselves supply these premises. A model preserving
operational no-signaling is not excluded by these tests.

## Reproduce

Python 3.12, with the pinned [requirements](requirements.txt). Full inputs are
about 2.3 GB without the optional 0.8 GB source-code archive. Store them outside
this repository. The decoder streams ZIP members without extracting multi-GB
raw files. Its working memory is bounded by a chunk, not the entire trial table.

```sh
python -m pip install -r studies/nist-bell-causal-audit/requirements.txt
python studies/nist-bell-causal-audit/fetch_data.py --directory /tmp/nist-inputs
python studies/nist-bell-causal-audit/reconstruct.py \
  --hdf5 /tmp/nist-inputs/hdf5.hdf5 \
  --alice /tmp/nist-inputs/alice.zip --bob /tmp/nist-inputs/bob.zip \
  --output /tmp/nist-output/reconstruction.json
python studies/nist-bell-causal-audit/audit.py \
  --hdf5 /tmp/nist-inputs/hdf5.hdf5 \
  --reconstruction /tmp/nist-output/reconstruction.json \
  --output-dir /tmp/nist-output
python studies/nist-bell-causal-audit/verify_reproduction.py /tmp/nist-output
```

The last command compares full regenerated count tables, analysis and raw
reconstruction summary to committed artifacts. It fails on discrepancies.
The optional `--keys code` download provides the archived processing sources.
The original input hashes are checked before reconstruction.

Offline CI runs real-byte decoder fixtures, source-diagnostic reconciliation,
mathematical/adversarial tests and regeneration from the committed count tables:

```sh
python -m unittest discover -s tests -v
python studies/nist-bell-causal-audit/audit.py --check
python studies/nist-bell-causal-audit/report.py --output-dir /tmp/nist-report
```

Offline CI does not re-download or re-decode the 2.3 GB dataset. Full source
reproduction was run locally; its provenance is in the result summary. The
small fixtures are narrow regression evidence, not a substitute for that run.
No synthetic measurements are mixed with the measured result tables.
