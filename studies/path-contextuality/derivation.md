# Instrument, null model and adversarial checks

## The finite instrument

Let 0≤a≤b and a²+b²=1. On the two-dimensional path space define
K₋=diag(b,a), K₊=diag(a,b), E=diag(0,1), Z=diag(1,−1).
Direct multiplication gives ΣK†K=I and

- K₋†K₋ = (b²−a²)(I−E)+a²I;
- M(X)=ΣKXK†=(1−d)X+dZXZ, with d=(a−b)²/2;
- q=b²=(1+p_m)/2 and p_m=b²−a².

The first identity is the noisy path-measurement equivalence; the second is a
channel identity for every complex input X. Each branch is completely positive
because it is a Kraus map. Completeness supplies trace preservation. An explicit
isometric implementation is V|Q⟩=|Q⟩(b|−⟩+a|+⟩),
V|P⟩=|P⟩(a|−⟩+b|+⟩). The orthogonal path labels and a²+b²=1 prove V†V=I.
It can be extended to a unitary on path and a qubit probe; this specifies the
instrument without attributing capabilities to an existing device.

For source (x,y), success (z,−w), failure (w,z), each vector normalized and real,
the complete probabilities are

| Probe | Success | Failure |
|---|---|---|
| − | (bxz−ayw)² | (bxw+ayz)² |
| + | (axz−byw)² | (axw+byz)² |

Bypass success is f=(xz−yw)². The reference parameters are a=3/5,b=4/5,
x=z=4/5,y=w=3/5. The path-projector weak value is −9/7, but the test uses
the exact table, not that ratio or its weak-coupling approximation. Measured
post-probe success is 1513/15625, different from bypass f=49/625.
Pusey's p₋ uses division by the bypass probability; the later KLP p₋ is our joint
a. Neither is automatically the observed pointer frequency conditional on
post-probe success. Here that conditional frequency is 1369/1513.

## An arbitrary finite ontic model

Choose a normalized μ, a normalized K(m,λ′|λ), and 0≤r(λ)≤1. Preparation is the
same in probe and bypass trials; the final response is the same after either
procedure. The constraints are Σλ′K(−,λ′|λ)≤q and
ΣmK(m,λ′|λ)=(1−d)δλ′λ+dD(λ′|λ), with normalized D and 0≤d≤1.
They do not assume deterministic r, deterministic path or preparation
noncontextuality over an additional family of preparations.

For each input λ, separate λ′=λ from transitions to other states. Since
K(−,λ′|λ)≤dD(λ′|λ) off the diagonal, and r(λ′)≤1,

K(−,λ′|λ)r(λ′) ≤ K(−,λ′|λ)r(λ)+dD(λ′|λ)(1−r(λ)).

The same inequality holds on the diagonal by nonnegativity. Sum over λ′ and
then μ. This proves a≤qf+d(1−f). Also a≤q. This direct proof works for every
finite cardinality and never asserts that a deterministic enumeration is complete.
In an appropriate measurable kernel model the same pointwise bound can be
integrated, provided the diagonal/identity representation is measurable; that
extension is not a Lean theorem here.

The full parameterized ontic class is not presented as an LP: μ, K and r occur
multiplicatively. Fixing a guessed ontic cardinality and solving an LP would not
prove completeness. No such reduction or failed optimization is used for exclusion.

## Nonempty null and exact attacks

A one-state model, fair negative pointer and success response 1/4 is normalized,
undisturbed and saturates a=qf=1/8. At the reference parameters q=16/25,
d=1/50 and f=49/625, a fair one-state probe and stochastic final response
r=49/625 are also admissible: choose D=identity, so (1−d)I+dD=I. This gives
a=49/1250, below the bound. The reduced null at those parameters is therefore
nonempty; this is not a model of the full quantum calibration family.
The following complete two-state models are
independently evaluated with Fractions in `exact_countermodels()`; the first two
are also formalized in the companion.

1. **Hidden disturbance:** prepare a fair bit λ. Toss a fair pointer m and replace
   λ by m; final success means λ′=0. Bypass and final success marginals are both
   1/2, but a=1/2. Thus marginal change zero does not imply d=0: that mistaken
   substitution would give a gap 1/4. The pointer cap is correctly 1/2.
2. **Contextual pointer:** keep λ unchanged, output m=λ, final success iff λ=0.
   The entire joint table is the same as case 1 and the pointer marginal is fair,
   but the negative response at λ=0 is one. Operational fairness on one
   preparation does not imply the pointwise cap. The channel is exactly identity.
3. **Selection amplification within the null:** prepare failure state 0. Negative
   branches transition 0→1 with weight d=1/100 and stay 0 with weight q−d;
   positive branches stay 0 with weight 1−q, q=1/2. At state 1 both branches stay
   there, with weights q and 1−q. Final success iff state 1. Then f=0,a=d,
   and **every postselected pointer is negative** although the null is satisfied.
4. **Setting-dependent preparation:** even an identity transition with a fair
   pointer gives a=1/2 if probe trials prepare an always-success state, whereas
   f=0 if bypass trials prepare an always-failure state. Both preparations are
   normalized; comparing them as one μ invents a violation at d=0.

The first two examples establish non-identifiability from the *displayed* marginal
checks; they are not models for the entire tomographically complete quantum
calibration family. A claimed equivalence across all those procedures is a
stronger hypothesis. Generic approximate noncontextuality is an additional
assumption: small operational distance does not constrain representation distance.

## Derived robustness and losses

If actual a and f are within ε_a and ε_f of a model satisfying the exact
representation premises, the null obeys
W≤ε_a+|q−d|ε_f. This follows by subtracting the ideal score and bounding its two
linear differences. The companion proves it directly. Also a≤q+ε_a under
the same near-model premise. Combining the two bounds gives
 a≤min(q,qf+d(1−f))+ε_a+|q−d|ε_f,
which justifies using this allowance with the capped confidence rule. If an independent premise
instead allows response cap q+u and disturbance weight d+v, the bound becomes
qf+d(1−f)+uf+v(1−f); it follows by the same row proof. Operational tomography
alone supplies neither u nor v. We report a conditional premise, not a validated
secondary-procedure construction or a calibration-to-ontic-distance theorem.

To expose all failures, append independent probe erasure with probability 1−β,
final erasure with probability 1−η, and symmetric final-bit flip e. For ideal joint
matrix T, let T′ be T followed by that bit-flip channel. The nine outcomes are

- P(m,f)=βη T′mf for m∈{−,+}, f∈{success,failure};
- P(m,no detection)=β(1−η)Σf Tmf;
- P(lost probe,f)=(1−β)ηΣm T′mf;
- P(lost probe,no detection)=(1−β)(1−η).

Their sum is one. Bypass f′=η[(1−2e)f+e], a′=βηT′−,success,
q′=βq, while d is unchanged by erasing a *record after* the complete instrument.
Thus W′=a′−q′f′−d(1−f′) is directly evaluable without fair sampling. Unknown
physical absorption before the instrument is not this erasure model and requires
recalibration of the channel. The main theorem tolerates arbitrary stochastic
final response; the prospective η/e model is used only to predict power.

## Why the motivating trajectory bridge is absent

The experiment measures a bounded projector on a two-mode system. A velocity
ratio needs momentum/energy instruments, a position effect and dynamics relating
conditional moments to a worldline law. An integrated excitation time needs a
system–medium coupling and integration over time. None of those maps is determined
by this two-mode table. Mathematically, appending an unconstrained velocity or
history label to any λ without coupling it to K or r leaves every probability
unchanged; opposite labels are observationally indistinguishable on this interface.
Consequently there is no identifiable trajectory label from this instrument
family alone. A concrete dynamical bridge could add constraints, but it is not
an unstated consequence of projector contextuality.

## Finite calibration cannot establish equality by tolerance alone

For 0<p,p′<1, any finite Bernoulli count k out of n has positive probability
under both parameters. Hence the same finite record is possible under an exact
equivalence p=p′ and under an arbitrarily small violation of that equivalence.
A residual acceptance region alone cannot promote the equality to a logical
fact. Noncontextuality constrains equal procedures, and does not give a modulus
of continuity for unequal ones. The explicit hidden-variable countermodels above
show the separate problem with inferring representation constraints from a small
set of operational marginals. Secondary procedures must address the actual
procedure family and its physical availability; citing the existence of a
preparation-mixing method does not discharge that instrument obligation.
