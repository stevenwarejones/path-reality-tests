# NIST Bell-record causal audit

This study reconstructs **107,109,596 paired records** from one public NIST
Bell-test run and asks what they can establish about operational influence.
It produces measured marginal contrasts and assumption-indexed bounds, with
explicit barriers to a spacelike causal interpretation.

The independent raw decoder reproduces every stored setting and the archived
click-update behavior. Preserving every eligible pulse bit differs from 105 Alice
and 94 Bob words, adding 31 and 27 primary-window clicks respectively. All
no-click rows remain. Ambiguous settings and timestamp excursions receive
worst-case score treatment. This is not a reanalysis of the published Bell
p-value or evidence of a signaling anomaly.

The [numerical report](results/report.md) gives both directions, each receiver
setting, chronological drift, pulse and phase-window sensitivity, and assumed
randomization/leakage budgets. With **assumed** ideal joint assignment and zero
ordinary leakage, and the additional **assumed detector-event completeness and
record-order premise**, the revised largest primary conditional upper limit is
**0.005482% (54.82 per million)**. Limits across the four comparisons are
30.95–54.82 per million, about 13–29% of each comparison's smaller observed
click rate (130–429 per million across arms). These are descriptive scale
comparisons, not confidence bounds on a relative effect.

The original centered-score result, 0.2652%, is retained as a conservative
reference. It is 6–20 times the click rate and cannot exclude complete suppression
of a detector at these rates. The click-only, count-sensitive revision also
reports an unrestricted unknown-row envelope; its 217–237 per million limits
remain above the lower-rate baselines. The tighter result therefore depends
materially on the detector-record premise. Neither version is an unconditional
experimental limit on faster-than-light influence.

The revision was specified after seeing the original results. It is retrospective,
not preregistered or independently confirmatory. The target is the **absolute
average signed effect**, so sign-changing effects can cancel. All causal and
spacelike claim gates remain disabled. See the [derivation](click-statistics.md)
and [closest prior analyses](sources.md#closest-marginal-dependence-analyses).

## What is established and what is missing

- [Reconstruction](reconstruction.md): full raw overlap, source reconciliation,
  multiclick comparison, explicit excluded prefixes/suffixes and censored tails.
- [Click-only derivation](click-statistics.md) and [v1 reference](statistics.md): history-conditional assignment,
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
python studies/nist-bell-causal-audit/click_audit.py --directory /tmp/nist-output
python studies/nist-bell-causal-audit/verify_reproduction.py /tmp/nist-output
```

The last command compares full regenerated count tables, analysis and raw
reconstruction summary to committed artifacts. It fails on discrepancies.
The optional `--keys code` download provides the archived processing sources.
The original input hashes are checked before reconstruction.

Offline CI runs synthetic decoder fixtures, count-conservation checks,
mathematical/adversarial tests and regeneration from the committed count tables:

```sh
python -m unittest discover -s tests -v
python studies/nist-bell-causal-audit/audit.py --check
python studies/nist-bell-causal-audit/click_audit.py --check
python studies/nist-bell-causal-audit/report.py --output-dir /tmp/nist-report
```

The separate `nist-reproduce.yml` workflow downloads the three pinned inputs,
re-decodes both full raw streams, regenerates both analyses and requires exact
artifact equality. It provides a separate execution environment for the same
implementation, not an independent decoder or an independent scientific review.
The ordinary offline job remains network-free. No source event/count extracts
are committed, and synthetic examples never enter measured result tables.
