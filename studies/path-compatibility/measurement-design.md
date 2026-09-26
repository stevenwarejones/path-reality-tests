# Which extra measurement distinguishes which ambiguity?

The target is the spatial reference-phase ambiguity in the deposited example,
not a general distinction between quantum interpretations. A new acquisition's
response is a prospective model prediction, not a property inferred from an
unmeasured intervention.

## A finite, conditional decision problem

Specify two alternative phase families before data collection:

- H+: θ(x)=κx+φ with κ/π in [1/20,1/16].
- H−: κ/π in [−1/16,−1/20].
- In both, φ/π in [−1/40,1/40] relative to a separately characterized reference.

The nonzero slope restriction is a substantive hypothesis, not a conclusion
from the images or a confidence interval. Assume reference-fringe visibility
v≥9/10 and all-trial detection efficiency η≥4/5. Both require independent
calibration; no supplied image certifies them.

An interferometric quadrature readout has two detected outcomes and no detection.
Assign Z=+1,−1,0 respectively. For Y quadrature its mean is ηv sin θ. Thus losses
are retained, not conditioned away. X quadrature has mean ηv cos θ.

| Candidate | Guaranteed distinction | Main-trial consequence |
|---|---|---|
| Repeat local polarization images | Phase-gauge invariant | No finite uniform test |
| Reference intensity | Phase-gauge invariant | No finite uniform test |
| X reference quadrature at x=8 | The allowed φ=0 and opposite slopes give identical distributions | Worst-case separation zero |
| Y reference quadrature at x=4 | Signed mean magnitude ≥(4/5)(9/10)(1/2)=9/25 | 72 eligible trials suffice under the stated bounds |
| Y reference quadrature at x=8 | Signed mean magnitude ≥(4/5)(9/10)(9/10)=81/125 | 22 eligible trials suffice under the stated bounds |

At x=4, H+ phases lie in [7π/40,11π/40], so sin θ≥1/2. At x=8 they lie in
[3π/8,21π/40], where sin θ≥cos(π/8)>9/10; the negative family is reflected.
These are analytic conservative envelopes, not minima inferred from a sampled
grid or two best-fit models. The inequality cos(π/8)>9/10 follows from
cos²(π/8)=(2+sqrt(2))/4>81/100. The x=8 upper endpoint has sine cos(π/40),
which is larger than cos(π/8).

For independent stationary eligible trials, use the sign of the mean Z. With a
signed mean margin m, each error is at most exp(−n m²/2) by Hoeffding's inequality
for a variable in [−1,1]. We require type-I ≤0.01 and type-II ≤0.1, taking
n=ceil(2 log(100)/m²). The rounded n also guarantees the stricter 0.01 bound for
both errors. These guarantees are uniform over the declared restricted families.
The cost is **eligible main trials**, including erasures. No photon-count
likelihood is applied to the deposited grayscale images.

Recommendation: the Y reference quadrature at x=8 has the smallest certified
main-trial budget in this declared menu. This is not an apparatus-wide optimum.
The dataset supplies neither a calibrated implementation of this extra readout
nor the trials needed to certify η,v,φ. Therefore calibration-trial cost and total
end-to-end cost are explicitly `null`, not zero. A device proposal must supply
that calibration protocol and its cost before this becomes an executable total
budget.

## What remains impossible without the extra restrictions

If the two allowed slope families include κ=0, they contain an identical
completion. Any test then has type-I plus type-II error equal to one at that
shared completion. No number of trials, and no new measurement, gives uniform
separation of those overlapping classes. The Lean companion proves this
error-sum identity for any finite distribution and randomized decision rule.

For the two nonzero completions specifically, all existing local pointer and
reference-intensity settings remain invariant. Their inability to distinguish
is exact and survives arbitrarily precise repetition. The proposed reference
quadrature breaks that measurement symmetry only after its reference phase is
independently constrained. This is an assumption-versus-evidence result, not an
assumption-free exclusion of contextual or trajectory models.
