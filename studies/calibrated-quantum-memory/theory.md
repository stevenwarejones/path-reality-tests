# Null, global bound, and exact restricted ambiguity

## Scope and physical assumptions

This prospective protocol is narrower than general multi-outcome NMN
tomography. It uses a **single-outcome**, memoryless qubit intervention B,
intended to implement D(X) = Tr(X) I/2. There is no accessible classical output
register from B, no leakage, and no action of B on the environment. The same B
is used in calibration and dynamics. The six Pauli eigenstates and the output
Pauli measurements are trusted and fixed in one physical representation.
Preparation labels are not leaked to the process. A fresh input is initially
independent of the environment. Final measurement choices do not change the
earlier output state. Counts have independent trials within a fixed
setting and stable probabilities; cross-setting independence is not needed for
the union bound. These are assumptions, not properties established for NMN.

Let the classical-memory null contain every pre-intervention quantum
instrument {A_l}, with an arbitrarily large classical register l, and every
post-intervention conditional CPTP channel C_l. Thus

\[
\Lambda_B(X)=\sum_l C_l(B(A_l(X))),\qquad
A_l\text{ CP},\quad \sum_l A_l\text{ TP}.
\]

Integrals may replace the sum. Initial classical randomness and classical
history are absorbed into l. This allows input-dependent classical records,
not merely random unitary noise or a bounded two-state memory. It also permits
quantum system propagation through B. This is the classical-memory class
at the specified cut; it is not a fully classical theory of preparations and
measurements. For the broader comb problem, conditional instruments imply
separable temporal Choi factors, and separability implies PPT. A PPT relaxation
therefore outer-contains the null, but neither converse is assumed. No PPT
optimization is used by this implementation.

## Calibration-to-recovery bound

Define epsilon = (1/2)||B-D||_diamond. For an ideal D,

\[
\Lambda_D(X)=\sum_l \operatorname{Tr}[A_l(X)] C_l(I/2)
             =\sum_l \operatorname{Tr}(M_lX)\sigma_l.
\]

This is an entanglement-breaking (EB) channel, whatever the classical memory
dimension. For the six-state recovery score

\[
F_6(\Lambda)=\tfrac16\sum_{a\in\{x,y,z\},s=\pm1}
\operatorname{Tr}[\rho_{a,s}\Lambda(\rho_{a,s})],\quad
\rho_{a,s}=(I+s\sigma_a)/2,
\]

write M_l = m_l(I+v_l·sigma) and sigma_l = (I+w_l·sigma)/2.
Positivity gives |v_l|, |w_l| <= 1 and normalization gives sum m_l = 1.
Summing over the six states yields

\[
F_6(\Lambda_D)=\sum_l m_l(1/2+v_l\cdot w_l/6)\le 2/3.
\]

The pre-instrument with its retained classical flag and the post-processing
are CPTP channels. Contractivity of the diamond norm under these compositions
implies (1/2)||Lambda_B-Lambda_D||_diamond <= epsilon. Each output success
probability, and hence their average, changes by at most epsilon. Therefore

\[
\boxed{F_6(\Lambda_B)\le\min\{1,2/3+\epsilon\}.}
\]

This covers the entire declared null without fitting its memory dimension or
optimizing over a finite bank. It is a valid upper bound, not asserted to be the
sharp tradeoff at every epsilon. The ingredient EB benchmark and norm
composition are standard; no novelty claim is made.

## Operationally anchored calibration from 18 tables

Write the Bloch action of B as r -> t+Tr. Measure each of three output Paulis on
each of six input eigenstates. Let m_(a,b,s) be that +/-1 expectation. Then

\[
m_{a,b,\pm}=t_a\pm T_{ab}.
\]

These are physical probabilities using trusted probes, not distances between
gauge-dependent GST estimates. For each interval [l,u] on m, intersect the
three constraints on t_a from the +/- pairs, and bound T_ab by their difference.
If the t intervals do not intersect, return an acquisition/model warning,
not a memory rejection. Other physicality inconsistencies can remain: this
outer region is not advertised as an independent feasibility solver.
Likewise, a positive `infer` margin alone excludes the declared null but does
not certify joint quantum feasibility. The saved demonstration separately
checks an explicit quantum model against every observed table before giving
the memory interpretation. A new empirical application must do the same.

The linear difference map has the exact expansion

\[
(B-D)(X)=\tfrac12\sum_a t_a\sigma_a\operatorname{Tr}(X)
 +\tfrac12\sum_{a,b}T_{ab}\sigma_a\operatorname{Tr}(\sigma_bX).
\]

Each rank-one superoperator in this sum has diamond norm at most one:
the output operator sigma_a/2 has trace norm one, and the input trace
functional has completely bounded norm ||I||_infinity or ||sigma_b||_infinity,
both one. The triangle inequality gives the conservative certificate

\[
\epsilon\le\epsilon_U=
\min\{1,\tfrac12(\sum_a |t_a|_U+\sum_{ab}|T_{ab}|_U)\}.
\]

This leaves every unconstrained direction free. It does not replace partial
calibration by an assumed isotropic norm ball. All 18 settings are required;
the parser rejects incomplete tables.

## Finite-data coverage

Use alpha/2 for calibration and alpha/2 for dynamics, retaining these budgets
in source-removed regions. At n shots per setting the 18 expectation intervals
have radius sqrt(2 log(36/(alpha/2))/n); the six dynamics probability intervals
have radius sqrt(log(12/(alpha/2))/(2n)). Clip to their valid ranges.
Hoeffding and the union bound give simultaneous coverage at least 1-alpha.
Let F_L be the average of the six lower endpoints. Reject the null only if
F_L > min(1,2/3+epsilon_U). This follows on the same simultaneous event for
every null member, so it is a composite-null guarantee under the assumptions.
No post-selection of a witness or bootstrap coverage approximation is used.

For a separately justified calibration-to-dynamics channel drift bound d in
half diamond norm, replace epsilon_U by min(1,epsilon_U+d). Such a bound is
not inferred from the available NMN records. Independent calibration alone
does not establish transferability.

## Exact source-removed physical constructions

Define L_lambda(X) = lambda X+(1-lambda)Tr(X)I/2, 0 <= lambda <= 1.

* **Dynamics-only null:** the environment is trivial, both evolution segments
  are identity, and the supposed reset is B=L_lambda. This fits every terminal
  preparation/effect probability of L_lambda; information can remain in the
  system. A single dummy outcome is reported by B.
* **Quantum alternative:** swap the input system into an environment qubit,
  apply B=D on the system, then swap back and apply L_lambda. Its complete
  terminal channel is also L_lambda, and B has the same single dummy outcome.
  Equality holds for every input and final effect, not just the six sampled
  probabilities. For lambda>1/3, F_6=(1+lambda)/2>2/3; with a trusted D this
  cannot have classical-only environmental memory in the declared class.
* **Calibration-only null:** B=D with identity evolutions and no environment.
  It fits every calibration probability of the quantum construction but gives
  F_6=1/2. It is excluded by the complete dynamics source in the demonstration.

For lambda=1 the dynamics-only null has F_6=1 and is excluded by calibration.
The calibration-only null is excluded by dynamics. The quantum construction
fits all tables jointly, establishing nonempty joint quantum feasibility.
These are complete source-only witnesses for this **synthetic** experiment,
not numerical ablations or a witness for NMN. In the generated example all
24 probability constraints are checked, with the same confidence budgets.

An input |+x> and output X calibration measurement already separates the
two exact mechanisms: the success probabilities are (1+lambda)/2 and 1/2.
This is a minimal distinguishing test for this pair, not enough on its own to
upper-bound arbitrary instrument error. Full 18-setting calibration is used
for the composite-null certificate.

## Adversarial limits

The simulations include three-/four-direction classical measure-and-prepare
memory, a retained Z record, no memory, and imperfect resets. They test code,
not uniform coverage (the proof supplies coverage). Power intervals are
pointwise exact binomial intervals over simulated repetitions.

Three critical exclusions must remain visible:

1. Calibrate D, then silently use identity in dynamics: this memoryless
   adversary passes the simulated calibration source and yields perfect
   recovery. The procedure rejects because acquisition matching is false.
2. A random Pauli operation with its classical Pauli label retained has
   averaged channel D, yet a decoder with the label perfectly recovers the
   system. Calibrating only the averaged channel is insufficient if that
   register is available later. The single-outcome/no-side-channel assumption
   is essential; a multi-outcome extension must calibrate the flagged channel.
3. Hidden quantum leakage or direct intervention on the environment changes
   the system boundary or the intervention semantics. No simulated pass can
   certify these away. Likewise, imperfect unanchored preparations/readout
   are outside the present theorem and cannot inherit its coverage guarantee.

This is an exact obstruction for a specified restricted observation set, not a
general impossibility theorem about existing process-tomography datasets.
