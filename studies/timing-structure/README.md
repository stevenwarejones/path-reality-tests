# Timing structure in an existing Bell archive

Three complementary interrogations of the same pinned NIST run: cancellation
behind pooled averages, temporal setting fingerprints, and constructive limits
on what the recorded data identify. This extends the [timing feasibility baseline](../synthetic-timing-shifts/README.md).

**Nothing was detected; this is not strong evidence against structured effects.**
The four occupancy-matched toy cases each resolve positive hidden TV in 0/400
repetitions. The displayed future-setting examples require 873–1,247 events of
imbalance to reject, versus 11–133 observed (lag +1, any event, eta=0,
event-supported completion; known-event total and uncertain rows held fixed).
The weighted state-contribution contrast includes state prevalence: a nonzero
contrast alone would establish neither different conditional effects nor
cancellation. Cancellation is assessed by the hidden-TV result.

Start with the [results report](report.md) and [completed conclusions](conclusions.md).
The [protocol](protocol.md), [statistical derivation](statistics.md),
[related-work ledger](sources.md) and [recovery report](recovery-report.md)
explain what each result means. All stated comparisons, including halves,
uncertainty envelopes and assumption sensitivity, are in `results.json`.

## Reproduce offline

```sh
python -m pip install -r studies/synthetic-timing-shifts/requirements.txt
python studies/timing-structure/analyze.py --check
python studies/timing-structure/recovery.py --check
python -m unittest discover -s tests -p test_timing_structure.py -v
```

Python 3.12 and the pinned baseline NumPy/SciPy/h5py versions are used in CI.
`--check` compares generated outputs without overwriting them. Omitting it
regenerates derived artifacts. The recovery simulations take tens of seconds.

## Reproduce from pinned originals

```sh
python studies/nist-bell-causal-audit/fetch_data.py --directory /tmp/nist-inputs
for side in alice bob; do
  python studies/timing-structure/extract.py \
    --hdf5 /tmp/nist-inputs/hdf5.hdf5 \
    --archive "/tmp/nist-inputs/$side.zip" --side "$side" \
    --cache-dir /tmp/nist-structure \
    --output "studies/timing-structure/$side-counts.json" --check
done
python studies/timing-structure/analyze.py --check
```

Raw archives and per-row caches remain outside the repository. Budget roughly
3 GB for inputs plus 0.65 GB for the two caches. The extractor reconciles all raw
settings, legacy click words and timestamp excursions, then retains a hashed
external cache with one coarse outcome, detector-presence flag and prior sync
interval flag per row. It records multiplicity and censored raw endpoints.
The aggregate JSONs contain no per-event timestamp stream.

`--chunk-records` controls raw decoding and `--chunk-rows` controls aggregation.
Local validation uses 1,000,000/1,000,000 and 200,003/173,011 and requires identical
aggregate files, including identical cache hashes. `--from-cache` skips decoding
but verifies the HDF and cache hashes before aggregation; it is not a substitute
for full-source CI. Both receivers are regenerated from source in CI.

## Interpretation

Part A bounds current-setting effects under an explicit per-history assignment
model. Part B reports descriptive lag scores and a separate conditional test of
future-setting freshness. Part C proves exact equivalences under specified
observation maps. None uses the sender's outcome, coincidence selection or
postselected detected-photon denominators. These are different questions and
must not be combined into a claim of superluminal or backwards influence.

The fixed local-history states are uneven and the measured records are sparse.
The occupancy-matched toy study exposes the resulting power loss. No new
experimental collection is necessary to reproduce or finish this study.

## Figures

With the repository's pinned plotting dependencies installed, run:

```sh
python studies/timing-structure/plot.py --output-dir /tmp/timing-figures
```

The regular CI review artifact includes the cancellation and lag figures.
Lag bars are missing-record completion ranges, explicitly not confidence intervals.
