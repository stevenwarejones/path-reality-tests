# Source audit and selection

Four combinations were considered. “Available” below means downloaded and parsed;
statements only supported by a paper are labeled separately. All inspected bytes
are in the selected [manifest](manifest.json) or [candidate manifest](candidate-manifest.json).
Downloads remain in an external cache. No acquisition code, notebook, or pickle
from a deposit is executed. [Candidate audit output](results/candidates.json)
records exact inspected fields and count reconstruction.

| Combination | Ambiguity the extra family could remove | Audit and decision |
|---|---|---|
| Nairobi preparation tables + single-repeated-operation dimension records | Independent preparation/effect fits need not arise from shared gates; a length series could constrain repetition | Prior #11's 2022/2023 epochs remain separate. Downloaded Belem May 2023 repetition records are another device; deposit also lists Lima/Quito, not a matched Nairobi acquisition. No bridge verified. Stop this cross-acquisition candidate; no shared hardware parameters invented. |
| Quantinuum H1-1 2025-05-02 SPAM + single-qubit RB with leakage gadget | SPAM could constrain binary assignment; c/l outputs could separate survival from flagged population | Actual JSON/QASM and all counts verified. Same dated machine group, but no persistent physical-ion IDs or interleaved chronology/calibration map in these files. RB includes a two-qubit ancilla gadget absent from SPAM; its response cannot be identified from the two supplied binary SPAM circuits. A richer conditional analysis is possible, but a calibrated population bridge is unverified. Do not transfer its gadget results to Nairobi/GST. |
| Superconducting leakage-control Fig. 3d RB + Fig. 3e leakage RB | Ground-state survival and third-level population could separate gate loss from computational error | Downloaded v2 archive has experimental probability arrays, error bars, and fitted curves, with simulations explicitly named separately. Sequence indices/lengths are present, but the inspected figure records do not supply per-shot counts and primitive words. Useful method-control source, not a circuit-level matched count analysis here. This is not a claim that no useful reanalysis is possible. |
| Trapped-ion 2015-03-30 GST + RB | Periodic germ sequences constrain coherent gate error; randomized orderings constrain effects suppressed by repetition | Selected. Both files contain measured primitive words and counts. The primary paper explicitly states the two families were interleaved during one experiment. Use one apparatus-sharing scope, while retaining the absence of timestamps and shot chronology. |

## Selected immutable records

Repository: [supplemental-info-arXiv-1605.07674](https://github.com/pyGSTio/supplemental-info-arXiv-1605.07674/tree/fc939f90f7a2b32a20aec34f36465259f38837a3).
Commit `fc939f90f7a2b32a20aec34f36465259f38837a3`, associated deposit
[Zenodo 231329](https://doi.org/10.5281/zenodo.231329).
Every URL uses the commit; SHA-256 and size independently pin bytes.

| Field | GST | RB |
|---|---|---|
| File | `2015_03_30-GST_BB1_XYXY_8192_condensed.txt` | `2015_03_30-RB_0320_condensed.txt` |
| Measured versus simulated | Repository README identifies ExperimentalData as experimental runs; these are measured counts | Same; no simulated notebook results used as observations |
| Rows / distinct words | 4,657 / 4,657 | 2,004 / 1,996; repeats retained, never double-counted as separate evidence sources |
| Counts / denominators | Integral “plus count”, total exactly 50; 232,850 total shots | Integral “plus count”, total 270 or 285; 563,280 total shots |
| Sequence | Literal Gi/Gx/Gy with parenthesized repeated germs; lengths 0–8,198 | Expanded primitive Gi/Gx/Gy words; lengths 2–1,970 |
| Gate alphabet | Declared I, X(π/2), Y(π/2), including dynamically corrected implementation in paper | Same primitive alphabet compiled from Cliffords; no per-Clifford rescaling used |
| Reset/readout | One preparation and terminal binary measurement per word in the documented analysis convention | Same declared scope; no intermediate instrument in these words |
| Qubit / epoch | Single trapped-Yb-ion qubit, 2015-03-30 per source/paper | Same interleaved acquisition per paper, not inferred solely from filename |
| Timing actually recorded | File row order only; no timestamps, slot durations, shot order or per-shot reset log | Same; adjacent file rows are not claimed to be consecutive hardware events |
| Calibration | No independent calibrated leakage population/response in selected files | No leakage monitor or clock-phase intervention |
| Shared parameters proposed | Reset contrast, binary effect, three primitive CPTP operations; optional uniform transfer rates or two-state error process | Exactly the same parameters under stationarity within the interleaved acquisition |
| Missing fields | Chronological jobs/blocks, independent reset calibration, computational-subspace monitor, leakage-level response | Same; compiled scheduling details beyond primitive word are unavailable |

The source README and analysis notebook are also pinned. The notebook is read as
JSON/text to establish target signs and germ conventions, never executed.
The nominal plus effect is |1><1|. The empty word has 0/50 plus counts; the
four-Y word independently expands to four quarter turns and has 0/50 in GST,
3/285 in the first RB row. This is the Gate A prototype, not a fitted calibration.
The parser checks the exact column header, integer counts, denominators, trailing
zero fields, bounded expansion, and the complete gate grammar. Counts can be
reconstructed without pyGSTi or unpickling any computed objects.

Paper-supported apparatus match: [Blume-Kohout et al., 2017, “Comparison to
randomized benchmarking”](https://www.nature.com/articles/ncomms14485).
The deposited files alone do not independently establish physical timestamps.
Stationarity across that interleaving remains a tested model assumption, not a fact
proved by the acquisition description. A mislabeled gate, sequence-dependent
compilation, or reset variation would refute the common-operation interpretation.

## Runner-up schema details

The repeated-operation record [8208198 v1](https://zenodo.org/records/8208198)
contains Belem `data.json` specifying 200 jobs, 12 repetitions, 20,000 shots,
backend `ibmq_belem`, start `2023-05-29 13:08:34`. The audit reconstructs every
CSV denominator. Source text uses logical qubit 0, triple resets, SX/RZ and terminal
measurement, and a submitted repetition delay. These are source-code intentions;
no execution logs verify pulse timing. Its preparation/settings differ from Nairobi.
The presence of a Nairobi backend name in a generic script is not a Nairobi dataset.

The Quantinuum [hardware download documentation](https://docs.quantinuum.com/systems/user_guide/hardware_user_guide/benchmarks/hardware_data.html)
links the dated H1-1 SPAM and SQ_RB files. The documented schema is not blindly
copied: **`leakage_postselect` actually contains unflagged counts**, not leakage
probabilities, in the inspected file. Independent reconstruction from c/l strings
(with reversed classical-bit order) verifies all survival and unflagged counts.
SPAM has 2 circuits × 10,000 shots; RB has 16 circuits × 100 shots, ten logical
qubits per shot, lengths 2/128/512/2048 and four repetitions. QASM stores slot-ordering
operations and logical q/a registers, but neither stable physical-ion identities nor
per-circuit acquisition timestamps. The gadget includes ancilla preparation, ZZ
operations and measurement. These are measured shot records; expected_output is a
target label, not another measured dataset. Simultaneous ten-qubit observations are
not assumed independent. Dated download URLs have no immutable-version contract;
SHA-256 pins inspected bytes and fails closed on upstream change.

The leakage-control [13808824 v2](https://zenodo.org/records/13808824) archive
is paired with [Hyyppä et al., PRX Quantum 5, 030353](https://doi.org/10.1103/PRXQuantum.5.030353).
Fig. 3d/e READMEs describe 25 random sequences per length, experimental ground/third
level probabilities and fitted curves. The complete inspected file inventory is
reproducible. Gate words, raw shot denominators, time ordering and calibration-to-job
links needed here are not in those inspected figure arrays. Simulated curves are
explicitly tagged `sim`; no experimental population limit from this device is used
for the target ion.

## Baseline preservation

Current main at selection was `95d209e`, after PR #11 merge `4f3b825`.
The reviewed PR head `a2a2da6e3dd036f50061aa87af69670d8a07f815` is retained unchanged
through the existing [Born-rule study](../born-rule-identifiability/README.md).
Its table parsers, qutrit certificates and epoch separation remain the baseline.
This study does not reinterpret its tested .1%/.2% preparation-population caps as
measured or minimal leakage, and supplies no temporal explanation of Nairobi.
