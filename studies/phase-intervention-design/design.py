#!/usr/bin/env python3
"""Conservative fixed-sample precision and power for a prospective four-phase test.

No measured data enter these calculations. Trials at each of the four settings
are independent Bernoulli with fixed probabilities; calibration is an external
bound on the true null contrast. Selection and stopping are fixed in advance.
"""
import argparse
import json
import hashlib
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def probability(value, name):
    if not math.isfinite(value) or not 0 < value < 1:
        raise ValueError(f'{name} must be finite and strictly between 0 and 1')
    return value


def contrast_radius(n, failure):
    """Twice the simultaneous per-setting Hoeffding radius (four union terms)."""
    probability(failure, 'failure probability')
    if isinstance(n, bool) or not isinstance(n, int) or n < 1:
        raise ValueError('n must be a positive integer')
    return math.sqrt(2 * math.log(8 / failure) / n)


def sample_plan(contrast, nuisance, alpha=.01, beta=.1):
    """Sufficient n for power >= 1-beta at size <= alpha, conditional on premises.

On the simultaneous beta event the empirical maximum contrast is at least
contrast-radius_beta. The rejection rule requires it to exceed
nuisance+radius_alpha. Strict inequality is retained at integer rounding.
"""
    probability(alpha, 'alpha')
    probability(beta, 'beta')
    if not (math.isfinite(contrast) and math.isfinite(nuisance)
            and 0 <= nuisance < contrast <= 1):
        raise ValueError('need 0 <= nuisance < anticipated contrast <= 1')
    numerator = math.sqrt(2 * math.log(8 / alpha)) + math.sqrt(2 * math.log(8 / beta))
    n = math.floor((numerator / (contrast - nuisance)) ** 2) + 1
    # Explicitly preserve the strict test in case floating-point rounding moves n.
    while contrast <= nuisance + contrast_radius(n, alpha) + contrast_radius(n, beta):
        n += 1
    return dict(anticipated_contrast=contrast, null_contrast_budget=nuisance,
                alpha=alpha, beta=beta, trials_per_setting=n, total_trials=4*n,
                rejection_contrast_threshold=nuisance + contrast_radius(n, alpha),
                contrast_radius_alpha=contrast_radius(n, alpha),
                contrast_radius_beta=contrast_radius(n, beta),
                power_lower_bound=1-beta)


def classify(counts, totals, nuisance, alpha=.01):
    """Fixed-sample rule on complete herald-denominated selected-bin counts.

Returns only reject-calibrated-null or inconclusive. Unequal totals are accepted
using separate simultaneous radii; they are never pooled or detection-conditioned.
"""
    probability(alpha, 'alpha')
    if not math.isfinite(nuisance) or not 0 <= nuisance <= 1:
        raise ValueError('null contrast budget must be in [0,1]')
    if len(counts) != 4 or len(totals) != 4:
        raise ValueError('exactly four phase settings are required')
    for k, n in zip(counts, totals):
        if any(isinstance(x, bool) or not isinstance(x, int) for x in (k, n)):
            raise ValueError('counts and totals must be integers')
        if n < 1 or not 0 <= k <= n:
            raise ValueError('require 0 <= count <= positive total')
    estimates = [k/n for k, n in zip(counts, totals)]
    radii = [contrast_radius(n, alpha)/2 for n in totals]
    lower = [abs(estimates[i] - estimates[i+2]) - radii[i] - radii[i+2] for i in (0, 1)]
    return dict(estimates=estimates, per_setting_radius=radii,
                lower_absolute_contrasts=lower, null_contrast_budget=nuisance,
                decision='reject-calibrated-null' if max(lower) > nuisance else 'inconclusive')



def binomial_interval(k, n, alpha, family=4, method='cp'):
    """Two-sided interval, Bonferroni allocation alpha/family to each setting."""
    from scipy.stats import beta as beta_dist
    probability(alpha, 'alpha')
    if isinstance(family, bool) or not isinstance(family, int) or family < 1:
        raise ValueError('family must be a positive integer')
    if any(isinstance(x, bool) or not isinstance(x, int) for x in (k, n)) or n < 1 or not 0 <= k <= n:
        raise ValueError('require integer 0 <= k <= positive n')
    if method == 'hoeffding':
        r = math.sqrt(math.log(2*family/alpha)/(2*n))
        return max(0., k/n-r), min(1., k/n+r)
    if method != 'cp':
        raise ValueError('method must be cp or hoeffding')
    tail = alpha/(2*family)
    return (0. if k == 0 else float(beta_dist.ppf(tail, k, n-k+1)),
            1. if k == n else float(beta_dist.isf(tail, k+1, n-k)))


def interval_classify(counts, totals, nuisance, alpha=.01, method='cp', occupation=None):
    """Simultaneous interval test; optional occupation is (P-count, total).

An occupation trial has complete P/Q outcomes or an independently justified
calibration bound; missing region detections cannot simply be called Q.
"""
    if len(counts) != 4 or len(totals) != 4:
        raise ValueError('exactly four phase settings are required')
    if not math.isfinite(nuisance) or not 0 <= nuisance <= 1:
        raise ValueError('nuisance must be in [0,1]')
    family = 5 if occupation is not None else 4
    intervals = [binomial_interval(k, n, alpha, family, method) for k, n in zip(counts, totals)]
    lower = [max(intervals[i][0]-intervals[i+2][1],
                 intervals[i+2][0]-intervals[i][1], 0.) for i in (0, 1)]
    occ = None if occupation is None else binomial_interval(*occupation, alpha, family, method)
    bound = nuisance + (0. if occ is None else occ[1])
    return dict(intervals=intervals, occupation_interval=occ,
                lower_absolute_contrasts=lower, null_upper_bound=bound,
                decision='reject-calibrated-null' if max(lower) > bound else 'inconclusive')


def certified_power(probabilities, totals, nuisance, alpha=.01, beta=.1,
                    method='cp', occupation=None):
    """A sufficient power certificate, not a normal approximation or exact power.

Under fixed independent binomial sampling, beta-tail count envelopes hold
simultaneously with probability >=1-beta. Monotonic confidence endpoints give
an explicit worst case over that envelope. Occupation is (true w, sample size).
"""
    from scipy.stats import binom
    probability(beta, 'beta')
    if len(probabilities) != 4 or len(totals) != 4:
        raise ValueError('exactly four probabilities and totals are required')
    family = 5 if occupation is not None else 4
    ps = list(probabilities) + ([] if occupation is None else [occupation[0]])
    ns = list(totals) + ([] if occupation is None else [occupation[1]])
    if any(not math.isfinite(p) or not 0 <= p <= 1 for p in ps):
        raise ValueError('probabilities must be in [0,1]')
    for n in ns:
        binomial_interval(0, n, alpha, family, method)
    if not math.isfinite(nuisance) or not 0 <= nuisance <= 1:
        raise ValueError('nuisance must be in [0,1]')
    tail = beta/(2*family)
    low = [int(binom.ppf(tail, n, p)) for n, p in zip(ns, ps)]
    high = [int(binom.isf(tail, n, p)) for n, p in zip(ns, ps)]
    lows = [binomial_interval(k, n, alpha, family, method)[0] for k, n in zip(low, ns)]
    highs = [binomial_interval(k, n, alpha, family, method)[1] for k, n in zip(high, ns)]
    contrast_lower = max(lows[0]-highs[2], lows[2]-highs[0],
                         lows[1]-highs[3], lows[3]-highs[1])
    null_upper = nuisance + (highs[4] if occupation is not None else 0.)
    certified = contrast_lower > null_upper
    return dict(certified=certified, power_lower_bound=1-beta if certified else 0.,
                count_lower=low, count_upper=high, contrast_lower=contrast_lower,
                null_upper=null_upper, method=method, alpha=alpha, beta=beta)


def interval_plan(probabilities, nuisance, alpha=.01, beta=.1, method='cp', w=None):
    """Return a verified sufficient equal-arm sample size; no minimality claim.

Integer binomial quantiles are discrete. Binary refinement only reduces the
search bracket; the final returned integer is checked against the certificate.
"""
    # Validate the whole design before indexing or searching.
    certified_power(probabilities, [1]*4, nuisance, alpha, beta, method,
                    None if w is None else (w, 1))
    d = max(abs(probabilities[0]-probabilities[2]), abs(probabilities[1]-probabilities[3]))
    if d <= nuisance + (0 if w is None else w):
        raise ValueError('anticipated contrast does not exceed null budget')
    def test(n):
        return certified_power(probabilities, [n]*4, nuisance, alpha, beta, method,
                               None if w is None else (w, n))
    upper = 1
    while not test(upper)['certified']:
        upper *= 2
        if upper > 2**31:
            raise ValueError('no design found within computational search limit')
    lower = upper//2
    while lower+1 < upper:
        middle = (lower+upper)//2
        if test(middle)['certified']:
            upper = middle
        else:
            lower = middle
    result = test(upper)
    if not result['certified']:
        raise ArithmeticError('final integer sample size failed its strict certificate')
    return dict(trials_per_setting=upper, total_trials=upper*(5 if w is not None else 4),
                probabilities=list(probabilities), occupation=w, nuisance=nuisance, **result)


def lossy_probabilities(eta, visibility):
    if any(not math.isfinite(x) or not 0 <= x <= 1 for x in (eta, visibility)):
        raise ValueError('efficiency and visibility must be in [0,1]')
    return [eta*(1+visibility)/2, eta/2, eta*(1-visibility)/2, eta/2]


def examples():
    from scipy import __version__ as scipy_version
    return {
        'kind': 'prospective conditional designs, not observed data or achieved power',
        'source': {'filename': 'design.py',
                   'sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                   'scipy_version': scipy_version},
        'premises': ['independent fixed-probability Bernoulli trials within each setting',
                     'complete herald-denominated outcomes',
                     'predeclared settings, selected bin, sample sizes and stopping',
                     'externally justified null-contrast and anticipated-signal bounds'],
        'plans': [sample_plan(d, b) for d, b in ((.8, .02), (.1, .02), (.05, .02), (.03, .02))],
        'interval_comparison': [dict(contrast=d, hoeffding_baseline=sample_plan(d, .02),
            cp=interval_plan([(.5+d/2), .5, (.5-d/2), .5], .02),
            hoeffding_binomial_power=interval_plan([(.5+d/2), .5, (.5-d/2), .5], .02, method='hoeffding'))
            for d in (.8, .1, .05, .03)],
        'low_probability_example': dict(center=.1, contrast=.03,
            cp=interval_plan([.115, .1, .085, .1], .02),
            hoeffding=interval_plan([.115, .1, .085, .1], .02, method='hoeffding')),
        'occupation_examples': [dict(eta=e, visibility=v,
            plan=interval_plan(lossy_probabilities(e, v), .02, w=.5))
            for e, v in ((1., .6), (.9, .9), (.8, .8), (.8, .7))],
        'occupation_grid': [dict(eta=e, visibility=v, nuisance=b, trials_per_setting=10000,
            **certified_power(lossy_probabilities(e, v), [10000]*4, b, occupation=(.5, 10000)))
            for b in (0., .01, .02, .05) for e in (.6, .7, .8, .9, 1.)
            for v in (.6, .7, .8, .9, 1.)],
        'calibration': 'If calibration can fail with probability alpha_cal, total type-I error is at most alpha+alpha_cal.',
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = examples()
    target = ROOT / 'results/design.json'
    if args.check:
        if json.loads(target.read_text()) != result:
            raise SystemExit('Prospective design snapshot differs; inspect before updating')
        print('Prospective precision/power snapshot passed')
    else:
        target.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')


if __name__ == '__main__':
    main()
