"""Count-sensitive, time-uniform predictable-mean bounds; see click-statistics.md."""
from math import expm1, isfinite, log

# Fixed family; optimizing over this grid is paid for in the error allocation.
LAMBDAS = (.005, .01, .02, .04, .08, .16, .32, .64, 1., 2.)
EVENT_LABELS = 2*4*3*3  # receiver side, setting context, phase, pulse group


def count_interval(lower_count, upper_count, n, *, start=0, alpha=.01):
    """Cover sum E[C_i|history] for binary C, even with interval-censored counts."""
    if any(not isinstance(v,int) or isinstance(v,bool) for v in (lower_count,upper_count,n,start)):
        raise ValueError('counts, n and start must be integers')
    if not (0 <= lower_count <= upper_count <= n and n > 0 and start >= 0 and isfinite(alpha) and 0 < alpha < 1):
        raise ValueError('invalid count interval arguments')
    budget = log(2*len(LAMBDAS)*EVENT_LABELS/alpha)+log(start+1)+log(start+2)
    lower = max([0.] + [(v*lower_count-budget)/expm1(v) for v in LAMBDAS])
    upper = min([float(n)] + [(-v*upper_count-budget)/expm1(-v) for v in LAMBDAS])
    return dict(lower=lower,upper=upper,empty=lower>upper,log_threshold=budget)


def effect_interval(arms, n, *, start=0, alpha=.01, per_history_tv=0., leakage=0.):
    """Invert q_xy in [1/4-epsilon,1/4+epsilon], then subtract arm means."""
    if len(arms)!=2 or not all(isfinite(x) for x in (per_history_tv,leakage)) or not (0<=per_history_tv<.25 and 0<=leakage<=1):
        raise ValueError('invalid nuisance budget or arm count')
    bounds=[count_interval(k,u,n,start=start,alpha=alpha) for k,u in arms]
    lo=[b['lower']/(n*(.25+per_history_tv)) for b in bounds]
    hi=[min(1.,b['upper']/(n*(.25-per_history_tv))) for b in bounds]
    lower=max(-1.,lo[1]-hi[0]-leakage)
    upper=min(1.,hi[1]-lo[0]+leakage)
    empty=any(b['empty'] for b in bounds) or any(a>b for a,b in zip(lo,hi)) or lower>upper
    return dict(lower=lower,upper=upper,absolute_upper=None if empty else max(abs(lower),abs(upper)),
                empty=empty,arm_probability_bounds=[list(pair) for pair in zip(lo,hi)],
                assumed_per_history_joint_tv=per_history_tv,assumed_leakage_gap=leakage)
