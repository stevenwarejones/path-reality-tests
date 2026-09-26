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


def examples():
    return {
        'kind': 'prospective conditional designs, not observed data or achieved power',
        'source': {'filename': 'design.py',
                   'sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
        'premises': ['independent fixed-probability Bernoulli trials within each setting',
                     'complete herald-denominated outcomes',
                     'predeclared settings, selected bin, sample sizes and stopping',
                     'externally justified null-contrast and anticipated-signal bounds'],
        'plans': [sample_plan(d, b) for d, b in ((.8, .02), (.1, .02), (.05, .02), (.03, .02))],
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
