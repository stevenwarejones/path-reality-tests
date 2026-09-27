"""Calibration controls: conventions, channel joins and acquisition blind spots."""
import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

STUDY = Path(__file__).resolve().parents[1] / 'studies/moving-noise-boundary'
sys.path.insert(0, str(STUDY))
from calibration_gate import (affine_wls, analyze, complex_periodogram,
                              inverse_susceptibility, load_lockin, lockin_power_gain)


class CalibrationGateTests(unittest.TestCase):
    def test_weighted_baseline_and_factor_two_trap(self):
        inputs = json.loads((STUDY/'results/calibration-inputs.json').read_text())
        result = analyze(inputs)
        for axis, energy, error in [('x', 231, 9), ('z', 720, 79)]:
            fit = result['cal_width_reconstruction'][axis]
            self.assertAlmostEqual(fit['reconstructed_caption_energy_pK'], energy, delta=1)
            self.assertAlmostEqual(fit['formal_energy_error_pK'], error, delta=1)
            self.assertGreater(abs(fit['literal_header_energy_pK']-energy), 5*error)
        self.assertFalse(result['cal_2020_to_2019_calibration_transfer_established'])
        self.assertFalse(result['new_empirical_spectrum_rejection'])

    def test_wls_recovers_affine_signal_and_weight_scaling(self):
        t = np.array([0., .1, .3, .7, 1.])
        y = 2.3 + 7.2*t
        sd = np.array([.3, .1, .4, .2, .7])
        p, c, chi = affine_wls(t, y, sd)
        q, d, _ = affine_wls(t, y, 3*sd)
        np.testing.assert_allclose(p, [2.3, 7.2], rtol=1e-12)
        np.testing.assert_allclose(q, p, rtol=1e-12)
        np.testing.assert_allclose(d, 9*c, rtol=1e-12)
        self.assertLess(chi, 1e-23)

    def test_export_blocks_and_clock_join(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'export.txt'
            t = np.arange(8)*.25-2
            data = np.concatenate([np.c_[t, np.full(len(t), j)] for j in range(9)])
            np.savetxt(path, data, delimiter=';')
            _, times, channels, dt = load_lockin(path)
            self.assertEqual(channels.shape, (9, 8))
            np.testing.assert_allclose(times, t-t[0])
            self.assertEqual(dt, .25)
            data[10, 0] += .01
            np.savetxt(path, data, delimiter=';')
            with self.assertRaisesRegex(ValueError, 'clock grids'):
                load_lockin(path)

    def test_rms_demodulation_psd_and_off_grid_line(self):
        # Independent real-signal Fourier normalization: RMS demodulation
        # requires no extra factor two when mapping back to positive RF.
        n, dt, amplitude = 4096, .01, 2.1
        t = np.arange(n)*dt
        fc, offset = 500/(n*dt), 17/(n*dt)
        z = amplitude*np.exp(2j*np.pi*offset*t)
        x = np.sqrt(2)*np.real(z*np.exp(2j*np.pi*fc*t))
        rf = 2*dt/n*abs(np.fft.rfft(x))**2
        baseband = complex_periodogram(z, dt, np.ones(n))
        self.assertAlmostEqual(rf[517], baseband[17], delta=1e-8)
        for df in [offset, offset+.173/(n*dt)]:
            z = amplitude*np.exp(2j*np.pi*df*t)
            for w in [np.ones(n), np.hanning(n)]:
                psd = complex_periodogram(z, dt, w)
                self.assertAlmostEqual(psd.sum()/(n*dt), amplitude**2, places=12)

    def test_force_injection_and_missing_filter_cannot_be_ignored(self):
        n, dt, f0, mass, gamma = 8192, .25, 26.7, .43e-6, 1.84e-5
        offset = 7/(n*dt)
        force_rms, gain = 1e-16, .16e6
        for tau in [0., 10., 100.]:
            # Analytic steady linear-oscillator injection into an explicitly
            # declared four-pole acquisition filter. Not an archive fit.
            h = (1+2j*np.pi*offset*tau)**-4
            inv = inverse_susceptibility(f0+offset, f0, mass, gamma)
            z = force_rms/inv*gain*h*np.exp(2j*np.pi*offset*np.arange(n)*dt)
            p = complex_periodogram(z/gain, dt, np.ones(n))
            recovered = np.sqrt(p[7]/(n*dt))*abs(inv)
            self.assertAlmostEqual(recovered/force_rms, abs(h), places=12)
            self.assertAlmostEqual((recovered/force_rms)**2,
                                   lockin_power_gain(offset, tau, 4), places=12)
        self.assertLess(lockin_power_gain(.004, 1000, 8), 4e-23)
        # A finite scan over tau is not a lower response certificate if the
        # archive does not bound tau: its infimum at nonzero offset is zero.
        self.assertLess(lockin_power_gain(.004, 1e8, 1), 1e-12)


if __name__ == '__main__':
    unittest.main()
