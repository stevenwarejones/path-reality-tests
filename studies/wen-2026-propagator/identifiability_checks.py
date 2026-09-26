#!/usr/bin/env python3
"""Exact synthetic counterexamples to recovering missing joint data from summaries.

These examples are not fits, reconstructed repeats, or estimates of experimental error.
No Dryad data enter this script. Rational arithmetic checks each equality.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def moments(values):
    n = len(values)
    mean = sum(values) / n
    return mean, sum((x-mean)**2 for x in values) / n


def analyze():
    # Marginal means AND variances coincide; only cross-repeat pairing changes.
    x = [F(9, 10), F(11, 10)]
    y_aligned, y_reversed = x[:], x[::-1]
    assert moments(y_aligned) == moments(y_reversed)
    aligned_products = [a*b for a, b in zip(x, y_aligned)]
    reversed_products = [a*b for a, b in zip(x, y_reversed)]
    assert moments(aligned_products)[0] == F(101, 100)
    assert moments(reversed_products)[0] == F(99, 100)
    # Magnitude and phase marginals coincide; their pairing changes |sum amplitudes|^2.
    # Complex values represented as exact rational (real, imaginary) pairs.
    a = [(F(1), F(0)), (F(0), F(2)), (F(-3), F(0))]
    b = [(F(2), F(0)), (F(0), F(3)), (F(-1), F(0))]
    intensity = lambda values: sum(z[0] for z in values)**2 + sum(z[1] for z in values)**2
    assert sorted(z[0]**2 + z[1]**2 for z in a) == sorted(z[0]**2 + z[1]**2 for z in b)
    assert intensity(a) == 8 and intensity(b) == 10
    # Same recorded success-conditioned shape admits different absolute losses.
    shape = [F(1, 4), F(3, 4)]
    efficiency = [F(1, 5), F(4, 5)]
    full = [[eta*x for x in shape] + [1-eta] for eta in efficiency]
    assert all(sum(row) == 1 for row in full)
    assert all([x/sum(row[:2]) for x in row[:2]] == shape for row in full)
    return {
        'kind': 'exact synthetic non-identifiability witnesses; no experimental input',
        'source': {'filename': 'identifiability_checks.py',
                   'sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
        'shared_factor_pairing': {'factor_mean': str(moments(x)[0]),
                                 'factor_population_variance': str(moments(x)[1]),
                                 'aligned_product_mean': str(moments(aligned_products)[0]),
                                 'reversed_product_mean': str(moments(reversed_products)[0])},
        'phase_magnitude_pairing': {'same_magnitudes': [1, 2, 3],
                                    'same_phases_in_pi_units': ['0', '1/2', '1'],
                                    'coherent_intensities': [str(intensity(a)), str(intensity(b))]},
        'loss_identifiability': {'same_conditional_shape': [str(x) for x in shape],
                                'full_distributions_including_failure': [[str(x) for x in row] for row in full]},
        'scope': 'These witnesses disprove general recovery from such marginal summaries. They do not identify the missing experimental records or bound their actual effect.'
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = analyze()
    target = ROOT / 'results/identifiability.json'
    if args.check:
        if json.loads(target.read_text()) != result:
            raise SystemExit('Identifiability snapshot differs; inspect before updating')
        print('Exact identifiability witnesses and source digest passed')
    else:
        target.write_text(json.dumps(result, indent=2) + '\n')


if __name__ == '__main__':
    main()
