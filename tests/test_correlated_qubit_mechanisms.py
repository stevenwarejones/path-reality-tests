"""Independent checks of observation identities, sharpness and certificates."""
from fractions import Fraction
import hashlib
import io
import math
from pathlib import Path
import pickle
import sys
import unittest
import zipfile

import numpy as np
from scipy.stats import poisson

STUDY = Path(__file__).resolve().parents[1] / "studies/correlated-qubit-mechanisms"
sys.path.insert(0, str(STUDY))
from cqm_response import (any_bounds, any_probability, corrected_contrast,
                          odd_probability, two_point_witness)
from cqm_analysis import affine_certificate, mutual_information, PrimitiveUnpickler
from cqm_sources import decode_member


class CorrelatedQubitMechanisms(unittest.TestCase):
    def test_poisson_counting_and_independent_background(self):
        # Direct summation of count probabilities is independent of the PGF formula.
        for intensity in (0, .001, .2, 3.2, 15):
            p = poisson.pmf(np.arange(160), intensity)
            self.assertAlmostEqual(p[1::2].sum(), odd_probability(intensity), places=13)
            self.assertAlmostEqual(p[1:].sum(), any_probability(intensity), places=13)
            for background in (0, .1, .49):
                odd = odd_probability(intensity)
                observed = odd*(1-background)+(1-odd)*background
                self.assertAlmostEqual(corrected_contrast(observed, background), 2*odd, places=12)
        for bad in (-1, float("nan"), float("inf")):
            with self.assertRaises(ValueError):
                odd_probability(bad)

    def test_sharp_mixed_poisson_bounds(self):
        with self.assertRaises(ValueError):
            any_bounds(1, 100)  # A finite cap cannot attain contrast exactly one.
        for contrast in (.001, .2, .8, .999):
            fixed = -math.log1p(-contrast)/2
            low, high = any_bounds(contrast)
            self.assertAlmostEqual(any_probability(fixed), low)
            for cap in (fixed, fixed+1, fixed+5):
                witness = two_point_witness(contrast, cap)
                reproduced = sum(w*2*odd_probability(l) for w, l in zip(witness['weights'], witness['intensities']))
                p_any = sum(w*any_probability(l) for w, l in zip(witness['weights'], witness['intensities']))
                self.assertAlmostEqual(reproduced, contrast)
                self.assertAlmostEqual(p_any, any_bounds(contrast, cap)[1])
                self.assertTrue(low-1e-14 <= p_any <= high+1e-14)
            with self.assertRaises(ValueError):
                any_bounds(contrast, fixed/2)
        # Random heterogeneous distributions test an entire region, not just endpoints.
        rng = np.random.default_rng(20260927)
        for _ in range(100):
            lam = rng.uniform(0, 5, 7)
            weight = rng.dirichlet(np.ones(7))
            c = sum(w*2*odd_probability(l) for w, l in zip(weight, lam))
            a = sum(w*any_probability(l) for w, l in zip(weight, lam))
            lo, hi = any_bounds(c, 5)
            self.assertLessEqual(lo, a+1e-14)
            self.assertLessEqual(a, hi+1e-14)

    def test_global_affine_infeasibility_certificate(self):
        data = [['0', '0', '.1'], ['1', '1', '.1'], ['2', '0', '.1']]
        r = affine_certificate(data)
        self.assertEqual(Fraction(r['required_error_multiplier']['exact_lower_bound']), 5)
        self.assertEqual(Fraction(r['additive_discrepancy_s_minus1']['exact_lower_bound']), Fraction(1, 10))
        self.assertLess(r['feasible_multiplier_upper']-5, 2e-8)
        dup = affine_certificate([['1', '0', '.01'], ['1', '.2', '.01']])
        self.assertEqual(Fraction(dup['required_error_multiplier']['exact_lower_bound']), 10)
        # Fits each of three observations exactly: no falsely positive obstruction.
        exact = affine_certificate([['0','1','.1'], ['1','3','.1'], ['2','5','.1']])
        self.assertTrue(exact['fits_four_reported_errors'])

    def test_common_phase_can_create_pooled_dependence(self):
        low = np.array([81, 9, 9, 1])
        high = np.array([1, 9, 9, 81])
        self.assertAlmostEqual(mutual_information(low), 0)
        self.assertAlmostEqual(mutual_information(high), 0)
        self.assertGreater(mutual_information(low+high), .3)

    def test_source_integrity_and_nonprimitive_pickle_rejection(self):
        stream = io.BytesIO()
        with zipfile.ZipFile(stream, 'w', compression=zipfile.ZIP_DEFLATED) as z:
            z.writestr('test.txt', b'pinned data')
        block = stream.getvalue()
        with zipfile.ZipFile(io.BytesIO(block)) as z:
            info = z.getinfo('test.txt')
        m = {'name':'test.txt', 'size':11, 'compressed_size':info.compress_size,
             'sha256':hashlib.sha256(b'pinned data').hexdigest()}
        self.assertEqual(decode_member(block, m), b'pinned data')
        with self.assertRaises(ValueError):
            decode_member(block, dict(m, sha256='0'*64))
        self.assertEqual(PrimitiveUnpickler(io.BytesIO(pickle.dumps([[1,2],[3,4]]))).load(), [[1,2],[3,4]])
        with self.assertRaises(ValueError):
            PrimitiveUnpickler(io.BytesIO(pickle.dumps(Path('/tmp')))).load()


if __name__ == '__main__':
    unittest.main()
