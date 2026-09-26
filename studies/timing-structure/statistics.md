# Estimands, filtrations and error accounting

## Current-setting strata: an extension of the baseline count bound

Fix a receiver, local setting y, feature f and predictable local state G_i=g.
Let C_i=1{X_i=x,Y_i=y,G_i=g,f(O_i)=1}. Under the stated exogeneity model,
its predictable mean is q_xy,i * 1{G_i=g} * p_xy,i(f), where p is the response
law conditional on the pretrial history and intervention on the CURRENT settings.
The history may contain arbitrary past device states and earlier settings.
Averaging over realized histories is not a sustained intervention on past settings.
For the pooled group omit G. Define

    a_gxyf = (1/N) sum_i 1{G_i=g} p_xy,i(f),
    d_gyf = a_g1yf - a_g0yf.

These are weighted stratum contributions; do NOT divide them by the empirical
number of clicks or by the number of observed assignments in that stratum.
The no-event contrast is -d_gy,any because the state mass is common to both
current-setting interventions. Thus T_g=.5 sum_{four categories}|d_gyf|.
The sum over both states is the mean conditional TV within those states;
T_pool is TV after forgetting the state. Convexity gives
0 <= H=T_0+T_1-T_pool <=1. Our lower bound max(0,L_0+L_1-U_pool) is conservative
and simultaneous. H>0 means cancellation at this chosen state resolution,
not unbounded instantaneous signalling. A zero lower bound does not show H=0.
The interaction interval covers d_1-d_0 and can reveal opposite contributions.

For arbitrary adapted Bernoulli C_i with mean mu_i,
E[exp(lambda*C_i)|F_(i-1)] = 1+mu_i*(exp(lambda)-1)
<= exp(mu_i*(exp(lambda)-1)). Product exponentials are nonnegative
supermartingales. Ville plus the fixed finite family gives b=
log(2*10*160*2/.005). On that event,

    L = max(0, max_lambda (lambda*k-b)/(exp(lambda)-1)),
    U = min(N, min_lambda (-lambda*u-b)/(exp(-lambda)-1)).

Here k<=sum C_i<=u are pathwise count completions. With q_xy,i in
[.25-epsilon,.25+epsilon], bound a by [L/(N*(.25+epsilon)),
U/(N*(.25-epsilon))] intersect [0,1], then subtract the arms. Both deterministic
starts (64 and floor(n/2)) are paid for. The first covers both the first-half
and full endpoints; the second covers the second half. Five overlapping groups,
four contexts, four features and both receivers are all in the family.
An unknown state can contribute to either state upper bound. We deliberately
ignore some cross-feature consistency constraints, making intervals wider.
Empty model/confidence intersections are flagged, never counted as a discovery.

## Why lagged comparisons need a different argument

For a past remote setting, X_(i-d) is already in the pre-outcome history.
One cannot replace its conditional probability by 1/2 at time i. Ordinary
apparatus memory can create a dependence without any spacelike influence.
The fixed-lag score 4*(k1-k0)/N is therefore a DESCRIPTIVE count contrast.
Its completion range covers the finite recorded-data score under missingness,
not a population parameter. Unequal setting frequencies can affect the score.
No error bar on that score is labelled as statistical confidence.

For future lag d>0, reindex at the assignment t=i+d. Let
A_t=1{Y_(t-d)=y,f(O_(t-d))=1}. Assume all of this earlier local record is in
F_(t-1), including correct row ordering across stations. Let Z_t=2X_t-1 and
assume |E[Z_t|F_(t-1)]|<=2*eta. Then for fixed -1<lambda<1,

    M_t/M_(t-1) = (1+lambda*Z_t*A_t)/(1+2*eta*abs(lambda)*A_t)

has conditional expectation <=1. No IID assumption, no independence among lags,
and no hypothesis about how earlier outcomes were generated is required.
For each signed stake the terminal log product, on complete data, is
k1*log(1+lambda)+k0*log(1-lambda)-(k0+k1)*log(1+2*eta*abs(lambda)).
An uncertain row has factor at least
(1-abs(lambda))/(1+2*eta*abs(lambda)); a certified no-event row has factor one.
Multiply that lower bound once per compatible uncertain source row, not once
per possible remote label. This remains a PATHWISE lower bound even if the
observed missingness mechanism is not predictable. The completed underlying
process must still satisfy the null and ordering assumptions.

Average the 16 terminal stake e-values equally. A lower bound on an e-value
still has expectation <=1. Reject at E_lower >=64/.005 by Markov/union bound.
The JSON includes min(1,64/E_lower), a conservative family-adjusted p upper
bound. Only the full interior endpoint is tested. Halves are descriptive at
nonzero lags. Stake directions and magnitudes are mixed, not selected for free.
The two completion envelopes bound the same underlying process, so do not
require two independent alpha allocations when their assumptions both hold.
As eta grows the e-value decreases: report nested sensitivity nulls, not the
smallest p across an unpenalized eta search. No inference certifies eta from
observed frequencies. A failure may indicate setting predictability, timing
misalignment, record errors or a failed ordering assumption; it is not evidence
of retrocausality on its own. A non-rejection is not an RNG certification.

## Joint scope

The count family spends .005 and the future family .005: total <=.01 under
BOTH sets of assumptions. This budget applies to this fixed study, not to every
past or future exploratory analysis in the repository. Confidence for multiple
assumed epsilon values is read as a sensitivity curve: if a bound epsilon_0
holds, the same underlying count event covers all epsilon>=epsilon_0.
If an external calibration can fail with probability gamma, add gamma.
No apparatus calibration, physical spacelike interpretation, independence of
Alice/Bob or independent-run replication is asserted here.
