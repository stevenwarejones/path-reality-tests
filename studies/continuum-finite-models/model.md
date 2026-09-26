# Model and mathematical derivations

These are mathematical derivations, with numerical checks in `model.py`. They
are **not yet certified by the companion Lean development**.

## Common experiment

A free nonrelativistic particle has mass m>0 and moves on a stationary circle of
circumference L>0. The Fourier Hilbert space is l²(Z;C), with normalized
coefficients c and E_j=hbar²(2πj/L)²/(2m). Evolution multiplies c_j by
exp(-itE_j/hbar). Moduli, normalization and inner products are preserved;
U(t+s)=U(t)U(s), U(-t)=U(t) inverse. Its position representation is the usual
unitary Fourier-series identification with L² of the circle. Continuous position
does not require a continuous energy spectrum.

Preparations and effects are calibrated in this common Fourier convention.
A complete finite POVM, including failure, is transported through the same
isometric embedding under both hypotheses. Controls fix time, mode labels and
readout phase; these are not separately fitted under each hypothesis. The core
witness uses |0>+|j> normalized by sqrt(2), exact spectral evolution, and a
balanced coherent recombiner followed by mode detection.

| Description | Finite objects | What it implies |
|---|---|---|
| Finite amplitude sum / matrix product | route index, intermediate bases | Identical amplitude by distributivity |
| Finite time slices, continuous positions | intermediate times only | Each position integration remains infinite |
| Numerical spatial grid | sites | Approximation, unless asserted fundamental |
| Spectral effective model | accessible Fourier modes | Exact on its invariant accessible subspace |
| Fundamental nearest-neighbor ring | sites and specified couplings | Distinct dispersion; continuous time still permits infinitely many histories |
| Finite ontic model | hidden states | Neither a spatial lattice nor a count of occupied trajectories |

## Two alternatives and aliasing

Set a=L/N. With cyclic shift S on C^N, H_a=hbar²(2I-S-S*)/(2ma²).
For F_nj=N^(-1/2)exp(2πinj/N), the geometric sum proves F*F=I.
Acting with S and S* multiplies a column by exp(∓2πij/N), so

    E_a(j)=hbar²(1-cos(k_j a))/(ma²).

Both shift terms must be retained for N=2; for N=1 H_a=0. Comparisons require
N>2J for accessible labels -J,…,J: labels then remain distinct modulo N.
The spectral alternative H_spec=F diag(E_j) F* (using one representative of each
residue) assigns exactly the continuum energies to the accessible modes.
Embedding them into l² is isometric, intertwines evolution and preserves every
transported POVM probability, including failure. Its position couplings are
generally long range. Finite dimensionality alone is therefore not excluded.

## Approximation and the full state

The Taylor integral remainder, together with 1-cos x=2sin²(x/2) and
|sin u|≤|u|, gives globally

    0 ≤ x²/2-(1-cos x) ≤ x⁴/24.

For |k|≤K, the energy gap lies in [0,hbar²a²K⁴/(24m)]. The scalar inequality
|exp(ix)-exp(iy)|≤|x-y| and the l² norm give uniformly for |t|≤T

    ||U(t)c - U_a(t)c|| ≤ epsilon = T hbar a² K⁴/(24m).

For normalized pure states, trace distance D=sqrt(1-|<u,v>|²)≤||u-v||.
A complete POVM contracts it to TV=half the outcome L1 distance. Thus TV≤epsilon.
No claim that the bound is optimal is made; a common energy shift can improve it.

Let P_J retain |j|≤J, with tail tau=||(I-P_J)c||². If tau<1, normalized projection
c_J=P_Jc/sqrt(1-tau) satisfies

    ||c-c_J|| = sqrt(2(1-sqrt(1-tau))),  D(c,c_J)=sqrt(tau).

Consequently TV≤min(1,sqrt(tau)+epsilon). If tau=1 normalized projection is
undefined: use any normalized accessible state and the trivial TV≤1 instead.
For every fixed normalized c, T<infinity, independent n<infinity and margin d>0,
choose J so sqrt(tau_J)<d/(2n), then an integer N>2J so epsilon<d/(2n).
This proves the fixed-state bounded-time approximation quantifiers. Uniformity
over a class requires a common tail envelope; an energy/bandwidth condition
cannot be silently omitted. Finite-rank projections cannot converge in operator
norm to the identity on all of l²: a unit vector outside their range has error 1.

For independent trials, telescoping the product laws and using the L1 norm gives
TV(P^n,Q^n)≤min(1,n TV(P,Q)); the same sum bound covers fixed nonidentical settings.
For any randomized rejection test f with 0≤f≤1, alpha+beta=1-(Qf-Pf)≥1-TV.
Thus no fixed finite resource experiment separates the continuum uniformly from
all sufficiently fine approximants in this access class. A fixed N may be testable.
This does not cover adaptive settings, correlated trials, ancillary channel
access or optional stopping.

## Implementable witness and its failures

For phase theta=(E_j-E_0)t/hbar-reference, the complete readout is

    p_plus=eta(1+v cos theta)/2,
    p_minus=eta(1-v cos theta)/2, p_failure=1-eta.

The ideal eta=v=1 form follows directly by projecting the evolved balanced
superposition onto (|0>±exp(-i reference)|j>)/sqrt(2). Dephasing multiplies the
cross term by v; equal efficiency multiplies both detector effects by eta.
Relative lattice phase is Delta=t[(E_a(j)-E_a(0))-(E_j-E_0)]/hbar.
At a continuum-locked zero reference the cosine channel difference is
eta*v*(1-cos Delta)/2; at the quadrature reference it is eta*v*|sin Delta|/2.
The two-phase menu has a nonzero contrast unless Delta is an integer multiple
of 2π. For 0<|Delta|<π it necessarily separates if eta*v>0.
For example L=2π, m=hbar=1, N=8, j=1 has E_a=16(1-sqrt(2)/2)/π²<1/2;
t=1 gives a nonzero contrast. The gap is a consequence of the specified kinetic
operator, not a measurement of how many paths exist.

Symmetric j,-j have equal energies in both models: zero signal at every time.
Also zero time, zero visibility and full loss give exact equality. A nonzero
dispersion correction wraps to equality at t=2πhbar/|delta E|.
An unconstrained per-setting reference offset absorbs every Delta exactly.
At one nonzero momentum, a free kinetic scale can match the energy ratio for
all times; multiple times alone do not fix that degeneracy. Two positive
momenta distinguish curvature if a common kinetic scale and phase origin are
calibrated: E_a(2)/E_a(1)=4cos²(π/N)<4=E(2)/E(1), for N>4.

`outcome_box` analytically encloses every point in a continuous nuisance box,
including phase extrema between endpoints. It relaxes cross-setting correlations.
A strictly positive disjoint coordinate gap is a sufficient robust separation
in full-outcome TV. Zero certified gap is inconclusive, unless an explicit
shared prediction (such as zero visibility) is constructed.
