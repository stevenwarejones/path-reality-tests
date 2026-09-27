# Frozen extension: fine digitized timing

Frozen on 2026-09-27 UTC, before this extension's extraction or new
setting-dependent results. PR #13 and its measured coarse results were already
known. This is retrospective, not preregistered. No simulation-informed boundary
selection is used. The first commit records this protocol and constants; subsequent
commits implement and execute it. This error budget applies to this fixed extension,
not to the sequence of exploratory studies in the repository.

## Population and retained variable

Same pinned NIST run, both receivers, original HDF correspondence, rows
[64,n-64), full population only. No new history/lag/subset searches. Select the
first channel-0 event in raw record order whose decoded pulse minus bitoffset
is 4, 5 or 6, without a phase-radius cut. Count multiplicities before selection.
No selected event is an outcome. Keep all uncertain interior rows.

Reuse reconstruct.py: raw timestamps are signed 64-bit integer tags; sync channel
6, detector 0. Pulse=floor(delay/period); phase=delay modulo period. The period
uses the preceding corrected sync interval divided by 800. Rounding a normal
interval to 800 pulses implies integer duration 129022..129182, hence phase
in [0,161.4775), with a safe fixed enclosing domain [0,162). Huge paired jumps
use the audited preceding-duration substitution; guards remain uncertain. Refuse
unhandled intervals, out-of-domain phases, or changed config. Peaks are 90/125
and pulse offsets 28/37 for Alice/Bob. No wrap across pulses is added.

Preserve the lossless ordered pair (pulse bit, decoder float64 phase) for each
first event in an EXTERNAL sparse cache. Phase is the decoder's floating modulo,
not a newly rounded integer or an inferred analog time. Report boundary-rounding
limitations. No-event rows are represented by absence from this selected-event
cache, while raw detector presence is evaluated separately for completion rules.
The committed aggregates use cell=162*(pulse_bit-4)+floor(phase), 486 cells.
Coarse category is exactly phase<pk, pk<=phase<pk+2, or phase>=pk+2; it must
reproduce #13 counts from these fine aggregates, including halves and uncertainty.

## Frozen families and error allocation

A: all 486 lexicographic grid thresholds, at the upper edge of each integer-phase
cell in pulse bits 4,5,6. Event sub-CDF F_xy(t) is on the ALL-interior-row scale;
no detected-photon normalization. The final threshold is event incidence.
Bound the grid supremum K_grid. Monotone interpolation gives a conservative upper
bound on K_full between grid points; retaining the cache does not turn the grid
into a full-resolution test.

B: nested partitions: four original outcomes; then separate pulse identity;
then split each pulse/category by floor(phase/8), floor(phase/2), floor(phase).
Always retain the no-event atom, deriving its contrast from total event mass.
Bins intersect the legal enclosing domain and the original categories, so every
level refines its predecessor even at peak boundaries. No empty cell is removed.
The finest partition still discards sub-tag decoder phase. TV_full >= TV_width1;
a conservative full-record upper bound uses event incidence. J_level is TV_level
minus TV_coarse, not a conditional-on-event effect or a causal attribution.

Each family receives alpha=.005. For each binary feature, receiver and setting
context pay both signs and all 10 lambda values (.005,.01,.02,.04,.08,.16,.32,
.64,1,2) in an exponential adapted-Bernoulli count inversion. One fixed start and
endpoint; no half/history inferential branch. A has 2*4*486=3888 count labels.
B pays 4*(1 + number of event bins across all five levels), summed over receivers;
common.py computes the exact count from the frozen maps, including duplicate
features. The extra 1 is event incidence. Total fixed-family error <=.01 under
the assignment and correspondence assumptions. No IID KS test or permutation.

q_xy,i in [.25-epsilon,.25+epsilon], epsilon=0,.001,.01, even conditional on
relevant latent/pretrial history. Infer time-averaged CURRENT intervention laws
over those histories; not sustained interventions on past settings. Propensities
are not certified by marginal frequencies. Both missingness envelopes use the
same covered latent counts: unrestricted allows an event in every compatible
uncertain row; event-supported requires a raw detector record and assumes its
completeness/order. Known settings restrict compatible contexts. No coordinate
restriction is trusted on uncertain detector events: allow any legal cell or
no event. Independent binary enclosures are conservative, not a sharp joint
missing-data optimization. Impose probability normalization and nested sum
consistency where feasible; flag empty confidence/model intersections.

## Recovery fixed before measured associations

Use 400 repetitions per case, seed 20260927: worst-case Monte Carlo standard
error <=.025, roughly five percentage points 95% half-width. Report exact
pointwise binomial intervals; zero/400 does not prove a 1% tail guarantee.
Compare every resolution on each identical realization. Use setting-blind
16 chronological block histograms (pooled over remote settings, retaining local
setting) as an empirical-background law, explicitly independent categorical
sampling within blocks, not an arbitrary-memory replica. Sample sizes .25,1,4
times archive size; sparse original event rates. Scale missingness with sample
size; include ideal no-missingness controls to separate statistical cost.

Within-category redistribution: mix a fraction a=0,.1,.25,.5,1 of event mass
into opposite fixed integer-phase cells pk and pk+1 in pulse bit 5 (core),
transferring only existing core mass, preserving every coarse probability.
Digital translation: shift the recorded phase-cell distribution by +1 or +2
integer cells in one remote arm within each pulse, saturating at phase cell161;
report boundary mass. Event-rate-only changes multiply incidence by 1+a while
keeping timing shape, a=.1,.5. Ordinary drift is retained in block rates. A
separate predictable-history stress control uses fresh settings with outcome
memory; assignment-bias and intentional one-row misalignment controls distinguish
violated assumptions from bugs. No exact analog-shift power claim. Preserve
#13's exact ideal-quantizer invisibility construction.

Freeze all choices above. Complete the family, quantify gains and limitations,
and stop even if no difference or positive J is resolved. No boundary tuning,
extra lags, favorable subset selection, or repeated error-budget optimization.
