# Sharp full-table compatibility

This study extends the path-contextuality design with necessary bounds and exact
attaining models for the complete ideal reference table. See [derivation.md](derivation.md)
and [evidence.md](evidence.md). It is not a measured violation.

At q=16/25, the positive-success bound forces d≥297/1225, versus 13/320
from the negative-success witness alone. At d=1/50 its absolute gap is smaller:
109/6250 versus 297/15625. This is improved disturbance tolerance, not a demonstrated
power advantage. The joint decision and power comparison are implemented in
[the contextuality design](../path-contextuality/design.py).

Run `python studies/path-compatibility/analyze.py --check`.
The exact-fraction attainers cover the boundary; floating-point LP tests are
independent numerical corroboration. Formal companion: ontology-separation #87.

The earlier Wen reconstruction and illustrative slope-family design are parked at
[parked/wen-reconstruction](https://github.com/stevenwarejones/path-reality-tests/tree/parked/wen-reconstruction).
They are outside this contextuality study. Their practical recommendation is to
measure the reference wavefront phase against an independently characterized
external reference to test the near-collimation approximation. Calibration and
total trial cost remain unidentified; the illustrative 22-trial number is not
an experimental budget.
