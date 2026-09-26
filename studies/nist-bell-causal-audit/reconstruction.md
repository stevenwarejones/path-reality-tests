# Reconstruction and provenance

Inputs are the NIST 2015-09-18 `03_43` post-modelock-recovery run, its two
compressed raw streams, the corresponding HDF5, and the archived processing
code. URLs, lengths and SHA-256 hashes are in [manifest.json](manifest.json).
The source catalog describes refrigerator warming near the end of this run.
The 03_31 training run is documented as lacking mode lock, so it is not silently
used as a control. See [sources.md](sources.md) for the exact documents inspected.

## Binary records and recovered overlap

Compressed events are little-endian packed (uint8 channel, uint64 time tag,
uint16 transfer index), 11 bytes. Time tags in this run are below 2^63, so the
signed arithmetic used for differences preserves their values. Channel 0 is
detection, 2 setting 0, 4 setting 1, 5 GPS PPS, and 6 sync. Channel 64 is counted
as an additional recorded marker; no physical interpretation is supplied here.
Do not treat it as a detector or a new trial. Streaming keeps a sync interval
across buffer boundaries and refuses malformed records or unsupported sync gaps.
The final open interval on each side is recorded as censored, never filled with 0.

The HDF5 omits the advertised offsets and per-click syncNumber. We locate its
first 64 settings uniquely in the first 500,000 raw records, then verify EVERY
setting over the full retained interval. The recovered raw sync offsets are
122,626 (Alice) and 143,611 (Bob). Alice has 107,232,222 closed raw intervals;
Bob has 107,264,281. The retained overlap has 107,109,596 rows. Alice has no
closed suffix, Bob has 11,074. These prefix/suffix counts are explicit and are
not represented as no-detection rows. No interior sync deletion is required to
reconcile this run. There are two code-3 settings at Alice; Bob's are unambiguous.

This is reconstruction of the archive's pairing, not an independent physical
synchronization certificate. Matching local setting sequences cannot prove that
the original cross-site alignment was correct. A GPS/propagation audit and
calibration of the hardware setting-generation/detection endpoints remain needed
for a spacelike interpretation. The parser supplies complete records in the
retained overlap, not a proof that every physical clock trial was recorded.

## Pulse assignment and multiclick correction

The archived builder estimates a local laser period from the previous corrected
sync interval: divide its duration by round(duration/(129102/800)). For enormous
timestamp jumps it reuses the preceding normal interval duration. We reproduce
that convention, then divide each detector's delay from sync into a pulse number
and phase. A click is eligible when abs(phase−peak)<radius. Sixteen pulse bits
begin at the configured bit offset. The nominal configurations are:

| Side | Peak (time-tag bins) | Radius (bins) | Pulse bit offset |
|---|---:|---:|---:|
| Alice | 90 | 4 | 28 |
| Bob | 125 | 5 | 37 |

A time-tag bin is 78.125 ps. Our narrow/wide variants decrease/increase the radius
by one bin, with fixed centers. The pulse groups use zero-based bits [5], [4,5,6]
and [3,4,5,6,7]. They are never coincidence-selected. These are distinct from
changing the phase radius and do not certify spacelike separation.

The archived `build_file_hdf5.py` updates `results[syncidx] += 2**bitidx` with
repeated advanced indices. NumPy buffered advanced indexing can overwrite earlier
updates for the same sync. The diagnostic implementation reproduces all archived
click words exactly. The independent implementation uses `np.bitwise_or.at` so
all distinct eligible pulse bits survive. It changes 105 Alice and 94 Bob words.
The primary three-bit outcome gains 31 Alice clicks and 27 Bob clicks.

Two short real excerpts reproduce this mechanism in offline tests. They contain
raw bytes, raw and HDF row indices, nominal configuration and expected corrected
and archived words. They are attributed NIST extracts under [NOTICE.txt](NOTICE.txt),
not synthetic measurements. The original full inputs are never edited.

The archived `diagnostics.v2.peter.csv` primary 16 counts reconcile exactly if
the two ambiguous code-3 Alice no-click rows are assigned to code 2. This is a
numerical reconciliation convention, not evidence for which physical setting
occurred. The audit retains ambiguity. The main documentation's suggested
mapping differs, so we do not infer a unique physical mapping from the table.
We do not use old spreadsheet no-click formulas, which NIST flags as potentially
incorrect. This audit is not a reproduction of a published Bell p-value.

## Timestamp uncertainty and provenance boundary

Ten Alice and eight Bob huge timestamp jumps occur in opposite-signed pairs.
`results/reconstruction.json` records their locations and the reconstruction
checks. The score interval treats all rows from each first jump through the
interval following its return as arbitrary. It also accounts for ambiguous
settings. This is a conservative robustness exercise, not an assertion that
these are all possible hardware failures or an estimate of their probability.

The available archive lacks `build_file_txt.py`, which the archived Python 2
scripts import. Thus rerunning the original chain unchanged is not possible from
this ZIP. Our independent decoder is instead validated by full setting equality,
full equality to the archived update behavior, explicit correction fixtures and
published-repository count reconciliation. The standalone HDF5's missing fields
are repaired for this run by reading raw events, not silently invented.
