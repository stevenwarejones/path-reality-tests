"""Two distinct inferential families; see statistics.md for their assumptions."""
from math import expm1, isfinite, log, log1p
import numpy as np
from common import LAGS, GROUPS, FEATURES

COUNT_ALPHA = .005
FUTURE_ALPHA = .005
COUNT_LABELS = 2*len(GROUPS)*4*len(FEATURES)
COUNT_STARTS = 2
FUTURE_TESTS = 2*2*len(FEATURES)*sum(l > 0 for l in LAGS)
LAMBDAS = (.005, .01, .02, .04, .08, .16, .32, .64, 1., 2.)
STAKES = tuple(s*v for v in LAMBDAS[:8] for s in (-1, 1))
EPSILONS = (0., .001, .01)
ETAS = (0., .0001, .001, .01)


def count_interval(k, u, n):
    if any(type(v) is not int for v in (k, u, n)) or not 0 <= k <= u <= n or not n:
        raise ValueError('invalid binary count completion')
    b = log(2*len(LAMBDAS)*COUNT_LABELS*COUNT_STARTS/COUNT_ALPHA)
    return [max([0.] + [(v*k-b)/expm1(v) for v in LAMBDAS]),
            min([float(n)] + [(-v*u-b)/expm1(-v) for v in LAMBDAS])]


def contrast(arms, n, epsilon):
    if not isfinite(epsilon) or not 0 <= epsilon < .25 or len(arms) != 2:
        raise ValueError('invalid assignment model')
    bounds = [count_interval(int(k), int(u), n) for k, u in arms]
    p = [(l/(n*(.25+epsilon)), min(1., u/(n*(.25-epsilon)))) for l, u in bounds]
    if any(l > u for l, u in p):
        return None
    return [max(-1., p[1][0]-p[0][1]), min(1., p[1][1]-p[0][0])]


def difference(a, b):
    return None if a is None or b is None else [a[0]-b[1], a[1]-b[0]]


def tv(intervals):
    if len(intervals) != 4 or any(c is None or c[0] > c[1] for c in intervals):
        return None
    near = [max(0., l, -u) for l, u in intervals]
    far = [max(abs(l), abs(u)) for l, u in intervals]
    lo, hi = max(max(near), .5*sum(near)), min(1., .5*sum(far))
    if lo > hi+1e-14:
        return None
    return [min(lo, hi), hi]


def cancellation(pooled, state0, state1):
    """TV weighted by stratum mass, and its excess over pooled TV."""
    if any(v is None for v in (pooled, state0, state1)):
        return dict(stratified_tv=None, hidden_tv=None)
    stratified = [state0[0]+state1[0], min(1., state0[1]+state1[1])]
    if stratified[0] > stratified[1] or stratified[1] < pooled[0]:
        return dict(stratified_tv=None, hidden_tv=None)
    hidden = [max(0., stratified[0]-pooled[1]), min(1., stratified[1]-pooled[0])]
    return dict(stratified_tv=stratified, hidden_tv=hidden)


def descriptive_score(k0, k1, u0, u1, n):
    if n <= 0 or min(k0, k1, u0, u1) < 0:
        raise ValueError('invalid descriptive counts')
    return dict(trusted_score=4*(k1-k0)/n,
                completion_range=[4*(k1-k0-u0)/n, 4*(k1+u1-k0)/n])


def future_log_e(k0, k1, uncertain, eta):
    """Pathwise lower bound on a mixture of fixed-stake freshness e-values."""
    if any(type(v) is not int or v < 0 for v in (k0, k1, uncertain)):
        raise ValueError('invalid event counts')
    if not isfinite(eta) or not 0 <= eta < .5:
        raise ValueError('invalid freshness sensitivity')
    values = []
    for stake in STAKES:
        penalty = log1p(2*eta*abs(stake))
        values.append(k1*log1p(stake)+k0*log1p(-stake)-(k0+k1)*penalty+
                      uncertain*(log1p(-abs(stake))-penalty))
    maximum = max(values)
    return maximum+log(float(np.exp(np.array(values)-maximum).mean()))


def future_result(k0, k1, uncertain, eta):
    value = future_log_e(k0, k1, uncertain, eta)
    threshold = log(FUTURE_TESTS/FUTURE_ALPHA)
    return dict(log_e_lower=value, reject_freshness=value >= threshold,
                family_adjusted_p_upper=float(np.exp(min(0., log(FUTURE_TESTS)-value))),
                assumed_eta=eta, log_rejection_threshold=threshold)


def minimum_imbalance(known, uncertain, eta):
    """Smallest majority-minus-minority count crossing the fixed family threshold."""
    if type(known) is not int or known < 0:
        raise ValueError('invalid known event total')
    threshold = log(FUTURE_TESTS/FUTURE_ALPHA)
    if future_log_e(0, known, uncertain, eta) < threshold:
        return None
    lo, hi = (known+1)//2, known
    while lo < hi:
        middle = (lo+hi)//2
        if future_log_e(known-middle, middle, uncertain, eta) >= threshold:
            hi = middle
        else:
            lo = middle+1
    return 2*lo-known
