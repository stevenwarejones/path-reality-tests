# Fine timing versus four coarse outcomes

This follow-up to [timing structure](../timing-structure/README.md) asks what the
same first-event records buy when pulse identity and one-tag phase cells are
retained. It adds no experiment, lag search or history search. Start with the
[measured report](report.md), [paired recovery](recovery-report.md) and
[conclusions](conclusions.md). The [frozen protocol](protocol.md),
[derivation](statistics.md), [provenance and prior work](sources.md), and
[validation](validation.md) distinguish the assumptions and resolution limits.

## Reproduce offline

```sh
python -m pip install -r studies/synthetic-timing-shifts/requirements.txt -r studies/wen-2026-propagator/requirements.txt
python studies/fine-timing-structure/analyze.py --check
python studies/fine-timing-structure/recovery.py --check
python studies/fine-timing-structure/plot.py --check
python -m unittest discover -s tests -p test_fine_timing.py -v
```

The 216 paired empirical-background cases use 400 repetitions each, plus six
explicit row-model controls. Recovery takes a few minutes; no raw archive is
needed for offline checks. Generated JSON and figures are committed. All bands
and all sensitivity cases are available, not only selected plots.

## Re-extract pinned originals

```sh
python studies/nist-bell-causal-audit/fetch_data.py --directory /tmp/nist-inputs
for side in alice bob; do
  python studies/fine-timing-structure/extract.py \
    --hdf5 /tmp/nist-inputs/hdf5.hdf5 \
    --archive "/tmp/nist-inputs/$side.zip" --side "$side" \
    --cache-dir /tmp/fine-events \
    --output "studies/fine-timing-structure/$side-counts.json" --check
done
```

This validates all source hashes, all 107,109,596 overlapping settings/legacy
click words and timestamp excursions per receiver. The final open raw interval
is censored. A sparse external cache preserves each selected first event's row,
pulse bit and original decoder float64 phase. It includes guarded rows solely
for provenance; those phases are never trusted by the inference. Raw streams
and event caches are not committed. The aggregates include no-event exposure,
uncertainty, half counts for exact #13 comparison, and 16 setting-blind block
histograms for simulations. Halves are not a new inferential branch.

`--chunk-records 200003` independently reproduces the default 1,000,000-record
extraction, including the cache checksum. The same script runs in full-source CI.
No fine-phase clipping or empirical support selection is permitted.
