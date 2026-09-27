"""Independent checks of the run-preserving readout forecast machinery."""
import io
from pathlib import Path
import pickle
import sys
import unittest

import numpy as np

STUDY = Path(__file__).resolve().parents[1]/'studies/correlated-qubit-mechanisms'
sys.path.insert(0, str(STUDY))
from cqm_mechanical import (NumericUnpickler, count_run, dwell_states,
                            fit_marginals, window_events, PERIOD, PHASE_ORIGIN)


class MechanicalForecast(unittest.TestCase):
    def test_dwell_partition_and_type_validation(self):
        self.assertEqual(dwell_states([(0,0,2),(1,2,5),(2,5,6)], samples=6).tolist(), [0,0,1,1,1,2])
        bad = [[(0.,0,2)], [(True,0,2)], [(0,0,2**64)], [(0,0,0)],
               [(3,0,2)], [(0,1,2)], [(0,0,1),(1,2,3)],
               [(0,0,2),(1,1,3)], [(0,0,1),(0,1,3)]]
        for rows in bad:
            with self.subTest(rows=rows), self.assertRaises(ValueError):
                dwell_states(rows, samples=3)
        with self.assertRaises(ValueError):
            window_events(np.array([0., .5, 1., 0.]), window=3)

    def test_persistent_window_counts_against_direct_oracle(self):
        a = np.zeros(101, dtype=np.int8)
        b = np.zeros(101, dtype=np.int8)
        a[2:5] = 1
        b[4:7] = 2
        a[34:36] = 1  # A two-readout excursion must NOT qualify.
        b[40:43] = 1
        a[68:71] = [1,2,1]  # Changing E/F does not interrupt non-G persistence.
        b[69:72] = 1
        lag = PERIOD-.00005-PHASE_ORIGIN
        bins = 715
        expected = np.zeros((bins, 9), dtype=np.int64)
        for start in range(0, len(a)-33, 33):
            ga, gb = a[start] == 0, b[start] == 0
            def event(s):
                return any(all(s[start+j+k] != 0 for k in range(3)) for j in range(1,32))
            ea, eb = event(a), event(b)
            index = int(((start*3e-6+lag+PHASE_ORIGIN) % PERIOD)/PERIOD*bins)
            expected[index] += [1,ga,gb,ga and gb,ga and ea,gb and eb,
                                ga and gb and ea and eb,ga and gb and ea,ga and gb and eb]
        np.testing.assert_array_equal(count_run([a,b], lag, bins), expected)
        self.assertEqual(expected[:, 6].sum(), 2)
        with self.assertRaises(ValueError):
            count_run([a, b[:-1]], lag, bins)

    def test_policy_does_not_fit_joint_target(self):
        train = np.zeros((5,9), dtype=np.int64)
        train[:, 0:4] = 100
        train[:, 4] = [1,2,5,9,3]
        train[:, 5] = [3,2,6,8,4]
        first, start = fit_marginals(train, policy_bins=2)
        train[:, 6:9] = 99  # All joint future outcomes change; inputs do not.
        second, other = fit_marginals(train, policy_bins=2)
        np.testing.assert_array_equal(first, second)
        self.assertEqual(start, other)
        self.assertEqual(start, 0)
        train[0, 1] = 0
        with self.assertRaises(ValueError):
            fit_marginals(train, policy_bins=2)

    def test_numeric_pickle_scope_and_native_scalars(self):
        value = [np.arange(4, dtype=np.float64), np.int64(7), np.float64(.2)]
        restored = NumericUnpickler(io.BytesIO(pickle.dumps(value, protocol=4))).load()
        np.testing.assert_array_equal(restored[0], value[0])
        self.assertIs(type(restored[1]), int)
        self.assertIs(type(restored[2]), float)
        with self.assertRaises(ValueError):
            NumericUnpickler(io.BytesIO(pickle.dumps(Path('/tmp')))).load()
        with self.assertRaises(ValueError):
            NumericUnpickler(io.BytesIO(pickle.dumps(np.complex128(1j), protocol=4))).load()


if __name__ == '__main__':
    unittest.main()
