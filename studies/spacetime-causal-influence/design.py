#!/usr/bin/env python3
"""Prospective, synthetic causal-influence study. No experimental inputs.

Primary inference: a fixed-horizon randomized score allowing device memory.
Optional comparison: simultaneous exact binomial intervals under an IID model.
Probability coverage is derived in statistics.md; floating point is not a proof.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.stats import beta, binom

HERE = Path(__file__).resolve().parent
C_M_PER_NS = 0.299792458


def integer(x, name, minimum=0):
    if isinstance(x, (bool, np.bool_)) or not isinstance(x, (int, np.integer)) or x < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return int(x)


def probability(x, name, open_interval=False):
    if isinstance(x, (bool, np.bool_)) or not isinstance(x, (int, float, np.number)):
        raise ValueError(f"{name} must be a real probability")
    x = float(x)
    if not math.isfinite(x) or not (0 < x < 1 if open_interval else 0 <= x <= 1):
        raise ValueError(f"invalid {name}")
    return x


def absolute_interval(lo, hi):
    if not -1 <= lo <= hi <= 1:
        raise ValueError("invalid difference interval")
    return max(0.0, lo, -hi), min(1.0, max(-lo, hi))


def score_radius(n, alpha):
    n = integer(n, "n", 1)
    alpha = probability(alpha, "alpha", True)
    return math.sqrt(2 * math.log(2 / alpha) / n)


def score_interval(matches, n, alpha=0.01, setting_bias=0.0):
    """Bounds mean conditional signed effect; constant effect gives binary TV.

    setting_bias bounds |P(X_i=1|history)-1/2|, including relevant hidden history.
    No optional stopping, outcome deletion, or average absolute-effect guarantee.
    """
    n = integer(n, "n", 1)
    matches = integer(matches, "matches")
    if matches > n:
        raise ValueError("matches exceeds n")
    setting_bias = probability(setting_bias, "setting_bias")
    if setting_bias > 0.5:
        raise ValueError("setting bias exceeds 1/2")
    estimate = 2 * matches / n - 1
    radius = score_radius(n, alpha) + 2 * setting_bias
    lo, hi = max(-1.0, estimate-radius), min(1.0, estimate+radius)
    lower, upper = absolute_interval(lo, hi)
    return dict(estimate=estimate, signed_lower=lo, signed_upper=hi,
                absolute_lower=lower, absolute_upper=upper, radius=radius)


def cp_interval(k, n, error):
    n = integer(n, "n")
    k = integer(k, "k")
    error = probability(error, "error", True)
    if k > n:
        raise ValueError("k exceeds n")
    if n == 0:
        return 0.0, 1.0
    lo = 0.0 if k == 0 else float(beta.ppf(error/2, k, n-k+1))
    hi = 1.0 if k == n else float(beta.isf(error/2, k+1, n-k))
    return lo, hi


def iid_difference(k0, n0, k1, n1, alpha=0.01):
    """Bonferroni CP difference interval; all counts include failure outcomes."""
    alpha = probability(alpha, "alpha", True)
    l0, u0 = cp_interval(k0, n0, alpha/2)
    l1, u1 = cp_interval(k1, n1, alpha/2)
    lo, hi = l1-u0, u1-l0
    lower, upper = absolute_interval(lo, hi)
    return dict(signed_lower=lo, signed_upper=hi,
                absolute_lower=lower, absolute_upper=upper)


def classify(interval, nuisance):
    nuisance = probability(nuisance, "nuisance")
    return dict(reject=interval['absolute_lower'] > nuisance,
                clean_absolute_upper=min(1.0, interval['absolute_upper']+nuisance))


def coupling_budget(setting0, setting1):
    """Union within each coupling, triangle across settings; no independence."""
    def total(xs):
        return min(1.0, sum(probability(x, "failure bound") for x in xs))
    e0, e1 = total(setting0), total(setting1)
    return min(1.0, e0+e1)


def geometry(distance_m, radius_a_m, radius_b_m, distance_error_m,
             a_earliest_ns, a_latest_ns, b_earliest_ns, b_latest_ns):
    vals = [distance_m, radius_a_m, radius_b_m, distance_error_m,
            a_earliest_ns, a_latest_ns, b_earliest_ns, b_latest_ns]
    if not all(isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)
               for x in vals):
        raise ValueError("nonfinite geometry")
    if min(distance_m, radius_a_m, radius_b_m, distance_error_m) < 0:
        raise ValueError("negative distance budget")
    if a_latest_ns < a_earliest_ns or b_latest_ns < b_earliest_ns:
        raise ValueError("reversed time interval")
    dmin = distance_m-radius_a_m-radius_b_m-distance_error_m
    dt = max(b_latest_ns-a_earliest_ns, a_latest_ns-b_earliest_ns)
    margin = dmin-C_M_PER_NS*dt
    return dict(minimum_distance_m=dmin, maximum_time_difference_ns=dt,
                margin_m=margin, margin_ns=margin/C_M_PER_NS, spacelike=margin > 0)


def sufficient_trials(observed_gap, nuisance, alpha=0.01, beta_miss=0.1):
    """Uniform per-history observed signed gap; strict, fixed-horizon power bound."""
    d = probability(observed_gap, "observed_gap")
    nuisance = probability(nuisance, "nuisance")
    alpha = probability(alpha, "alpha", True)
    beta_miss = probability(beta_miss, "beta_miss", True)
    if d <= nuisance:
        return None
    threshold = (math.sqrt(2*math.log(2/alpha)) +
                 math.sqrt(2*math.log(1/beta_miss)))**2 / (d-nuisance)**2
    n = math.floor(threshold)+1
    # Guard float rounding at the strict inequality boundary.
    while not d > nuisance+score_radius(n, alpha)+math.sqrt(2*math.log(1/beta_miss)/n):
        n += 1
    return n


def exact_score_power(n, gap, nuisance, alpha=0.01):
    """Numerical binomial sum for stationary independent score trials only."""
    n = integer(n, "n", 1)
    gap = probability(gap, "gap")
    nuisance = probability(nuisance, "nuisance")
    t = nuisance+score_radius(n, alpha)
    if t >= 1:
        return 0.0
    # Strict inequalities; recheck endpoints using the actual decision statistic.
    high = math.floor(n*(1+t)/2)+1
    while high <= n and not 2*high/n-1 > t:
        high += 1
    while high > 0 and 2*(high-1)/n-1 > t:
        high -= 1
    low = math.ceil(n*(1-t)/2)-1
    while low >= 0 and not 2*low/n-1 < -t:
        low -= 1
    while low < n and 2*(low+1)/n-1 < -t:
        low += 1
    p = (1+gap)/2
    return float(binom.cdf(low, n, p)+binom.sf(high-1, n, p))


def results(config):
    alpha, miss = config['alpha'], config['beta_miss']
    rows = []
    for d in config['observed_gaps']:
        for b in config['nuisance_budgets']:
            n = sufficient_trials(d, b, alpha, miss)
            rows.append(dict(observed_gap=d, nuisance=b, total_trials=n,
                certified_power_lower=None if n is None else 1-miss,
                stationary_binomial_power=None if n is None else exact_score_power(n,d,b,alpha),
                acquisition_seconds_at_hypothetical_rate=None if n is None else n/config['trial_rate_hz']))
    rng = np.random.Generator(np.random.PCG64(integer(config['seed'], 'seed')))
    n, reps = config['simulation_trials'], config['simulation_repeats']
    synthetic = []
    for d in (0.0, 0.03, 0.1):
        matches = rng.binomial(n, (1+d)/2, reps)
        t = config['simulation_budget']+score_radius(n,alpha)
        rejection_count = int(np.count_nonzero(np.abs(2*matches/n-1) > t))
        synthetic.append(dict(observed_gap=d, rejection_count=rejection_count,
            repeats=reps, empirical_rejection_rate=rejection_count/reps,
            stationary_binomial_power=exact_score_power(n,d,config['simulation_budget'],alpha)))
    null_example = score_interval(50000,100000,alpha)
    return dict(label='SYNTHETIC / PROSPECTIVE; no apparatus observations',
        estimand='mean conditional signed effect; TV only for a constant binary effect',
        config=config, geometry=geometry(**config['geometry']), power=rows,
        synthetic=synthetic, null_example=dict(**null_example, **classify(null_example,0.002)))


def compare(actual, expected, path='$'):
    """Float roundoff only; discrete values, structure and labels stay exact."""
    if type(actual) is not type(expected):
        raise AssertionError(f'{path}: type mismatch')
    if isinstance(actual, dict):
        if actual.keys() != expected.keys():
            raise AssertionError(f'{path}: key mismatch')
        for k in actual:
            compare(actual[k], expected[k], path+'.'+k)
    elif isinstance(actual, list):
        if len(actual) != len(expected):
            raise AssertionError(f'{path}: length mismatch')
        for i, (a,e) in enumerate(zip(actual,expected)):
            compare(a,e,f'{path}[{i}]')
    elif isinstance(actual, float):
        if not math.isfinite(actual) or not math.isclose(actual,expected,rel_tol=1e-10,abs_tol=1e-12):
            raise AssertionError(f'{path}: {actual} != {expected}')
    elif actual != expected:
        raise AssertionError(f'{path}: {actual!r} != {expected!r}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check',action='store_true')
    args = parser.parse_args()
    output = results(json.loads((HERE/'config.json').read_text()))
    target = HERE/'results.json'
    if args.check:
        compare(output,json.loads(target.read_text()))
        print('Spacetime prospective snapshot agrees.')
    else:
        target.write_text(json.dumps(output,indent=2,allow_nan=False)+'\n')


if __name__ == '__main__':
    main()
