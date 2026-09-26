"""Bounded-score sensitivity for a fixed family; physical assumptions are external."""
import math


def score_interval(total, unknown, n, *, alpha, comparisons, setting_tv, leakage, start=0):
    """Cover a mean signed effect, conditional on the premises in statistics.md.

    Known scores are in {-2,0,2}. Unknown rows may have ANY score in [-2,2].
    setting_tv bounds the history-conditional joint setting law's TV from uniform.
    leakage bounds the absolute shift of the average signed target effect.
    Neither parameter can be estimated just by inspecting setting balance.
    """
    if not isinstance(n, int) or not isinstance(unknown, int) or n <= 0 or not 0 <= unknown <= n:
        raise ValueError('invalid trial or unknown count')
    if not isinstance(comparisons, int) or comparisons < 1:
        raise ValueError('comparisons must be a positive integer')
    if not all(math.isfinite(v) for v in (total, alpha, setting_tv, leakage)):
        raise ValueError('nonfinite argument')
    if not 0 < alpha < 1 or not 0 <= setting_tv <= 1 or not 0 <= leakage <= 2:
        raise ValueError('invalid probability or gap allowance')
    if abs(total) > 2*(n-unknown):
        raise ValueError('impossible score total')
    if not isinstance(start, int) or start < 0:
        raise ValueError('start must be a nonnegative row index')
    # Sum_{s>=0} 1/((s+1)(s+2)) = Sum_{n>=1} 1/(n(n+1)) = 1.
    # Union over every contiguous interval permits retrospective endpoints.
    log_penalty = (math.log(2*comparisons/alpha) + math.log(start+1)
                   + math.log(start+2) + math.log(n) + math.log(n+1))
    sampling = math.sqrt(8*log_penalty/n)
    missing = 2*unknown/n
    radius = sampling + missing + 4*setting_tv + leakage
    center = total/n
    # Intersect with [-1,1]; an empty interval is an incompatibility with premises,
    # never silently converted into a point bound.
    lower, upper = max(-1., center-radius), min(1., center+radius)
    empty = lower > upper
    return dict(score_mean=center, start=start, sampling_radius=sampling,
                unknown_radius=missing, setting_bias_allowance=4*setting_tv,
                leakage_allowance=leakage, lower=lower, upper=upper,
                empty=empty, absolute_upper=None if empty else max(abs(lower), abs(upper)),
                excludes_zero=empty or lower > 0 or upper < 0)


def missing_outcome_difference(k0, n0, missing0, k1, n1, missing1):
    """Deterministic identification interval with known arm denominators; no CI."""
    for k,n,m in ((k0,n0,missing0),(k1,n1,missing1)):
        if any(not isinstance(v,int) for v in (k,n,m)) or n <= 0 or min(k,m)<0 or k+m>n:
            raise ValueError('invalid arm counts')
    return (k1/n1-(k0+missing0)/n0, (k1+missing1)/n1-k0/n0)
