"""Reproduce source-level quantities before admitting CAL or levitation data.

The output is deliberately a failed calibration gate, not a joint noise limit.
Raw archives live outside this repository. No fit here has a joint confidence
interpretation: the input workbook lacks covariance and lock-in samples overlap.
"""
import argparse
import hashlib
import json
from pathlib import Path
import urllib.request

import numpy as np
import openpyxl
from scipy.optimize import curve_fit

HERE = Path(__file__).resolve().parent
MANIFEST = HERE / 'gate-manifest.json'
INPUTS = HERE / 'results/calibration-inputs.json'
RESULT = HERE / 'results/calibration-gate.json'
RB87_MASS_KG = 1.443160e-25
KB = 1.380649e-23


def encode(obj):
    return json.dumps(obj, indent=2, sort_keys=True, allow_nan=False) + '\n'


def verify(path, entry):
    raw = path.read_bytes()
    if len(raw) != entry['bytes'] or hashlib.sha256(raw).hexdigest() != entry['sha256']:
        raise ValueError('Source changed: ' + path.name)


def fetch(directory):
    directory.mkdir(parents=True, exist_ok=True)
    for entry in json.loads(MANIFEST.read_text())['files']:
        path = directory / entry['name']
        if not path.exists():
            temporary = path.with_suffix(path.suffix + '.part')
            temporary.write_bytes(urllib.request.urlopen(entry['url'], timeout=180).read())
            verify(temporary, entry)
            temporary.replace(path)
        verify(path, entry)


def affine_wls(t, y, sd):
    """Reported diagonal weights only; covariance is NOT established empirically."""
    t, y, sd = map(np.asarray, (t, y, sd))
    if np.any(sd <= 0):
        raise ValueError('Positive uncertainties required')
    design = np.column_stack((np.ones_like(t), t))
    weighted = design / sd[:, None]
    covariance = np.linalg.inv(weighted.T @ weighted)
    coefficients = covariance @ weighted.T @ (y / sd)
    residual = (y - design @ coefficients) / sd
    return coefficients, covariance, float(residual @ residual)


def read_workbook(path):
    book = openpyxl.load_workbook(path, data_only=True, read_only=True)
    widths = [[float(v) for v in row[1:6]]
              for row in list(book['Figure 4g'].values)[6:12]]
    positions = []
    for row in list(book['Figure 4f'].values)[5:]:
        if isinstance(row[1], (float, int)):
            positions.append(list(row[1:6]))
    book.close()
    if len(widths) != 6 or len(positions) != 23:
        raise ValueError('Unexpected published source table')
    return {'widths': widths, 'positions': positions,
            'columns': ['TOF_ms', 'x_um', 'x_error_um', 'z_um', 'z_error_um'],
            'width_header': 'x-rms-width / z-rms-width (µm)',
            'reported_error_definition': 'figure caption: SD, n >= 3; full covariance unavailable'}


def load_lockin(path):
    """LabOne plotter exports nine sequential (time,value) blocks, not 18 columns."""
    header = []
    with path.open() as stream:
        for line in stream:
            if not line.startswith('%'):
                break
            header.append(line.rstrip())
    data = np.loadtxt(path, delimiter=';', comments='%')
    if data.ndim != 2 or data.shape[1] != 2 or not np.isfinite(data).all():
        raise ValueError('Malformed lock-in export')
    starts = np.r_[0, np.flatnonzero(np.diff(data[:, 0]) < 0) + 1]
    if len(starts) != 9:
        raise ValueError('Expected 3 amplitude, 3 phase, 3 frequency blocks')
    ends = np.r_[starts[1:], len(data)]
    blocks = [data[lo:hi] for lo, hi in zip(starts, ends)]
    t = blocks[0][:, 0]
    if any(b.shape != blocks[0].shape or not np.array_equal(b[:, 0], t) for b in blocks):
        raise ValueError('Channel clock grids do not match')
    dt = float(np.median(np.diff(t)))
    # Exported decimal times have finite precision; do not silently resample.
    if dt <= 0 or np.max(np.abs(np.diff(t) - dt)) > 2e-6:
        raise ValueError('Nonuniform sample grid')
    return header, t - t[0], np.array([b[:, 1] for b in blocks]), dt


def complex_periodogram(z, dt, window):
    """Two-sided complex PSD. For RMS demodulation this is real one-sided RF PSD.

    x(t)=sqrt(2) Re[z(t) exp(i omega_c t)] gives S_x^(1)(fc+df)=S_z^(2)(df)
    for non-overlapping sidebands. No extra factor two is applied.
    """
    z, window = np.asarray(z), np.asarray(window)
    return dt / np.sum(window**2) * np.abs(np.fft.fft(window * z))**2


def inverse_susceptibility(f, resonance_hz, mass_kg, gamma):
    return mass_kg * ((2*np.pi*resonance_hz)**2 - (2*np.pi*f)**2
                      + 1j*gamma*2*np.pi*f)


def lockin_power_gain(offset_hz, time_constant_s, order):
    """Sensitivity family, NOT inferred settings of the archived lock-in."""
    if time_constant_s < 0 or order < 1:
        raise ValueError('Invalid low-pass response')
    return (1 + (2*np.pi*np.asarray(offset_hz)*time_constant_s)**2)**(-order)


def lockin_summary(path):
    header, t, channels, dt = load_lockin(path)
    # Demod 2 is the near-26.7-Hz oscillator with decaying amplitude; demod 3
    # at the same reference has the detuned wheel phase. This role assignment
    # is supported by Fig S6, but the export has no wiring/configuration bundle.
    amp = channels[1]
    phase = np.unwrap(np.deg2rad(channels[4]))
    oscillator = channels[7]
    if np.ptp(oscillator) > 1e-8:
        raise ValueError('Reference frequency changed within acquisition')
    coefficients = np.polyfit(t, phase, 1)
    fc = float(oscillator[0] + coefficients[0] / (2*np.pi))
    parameters, _ = curve_fit(lambda t, A, tau: A*np.exp(-t/tau), t, amp,
                              p0=(amp[0], 1e5), bounds=([0, 1], [1, 1e8]),
                              x_scale=[.01, 1e5], ftol=1e-12, xtol=1e-12, gtol=1e-12)
    envelope = parameters[0] * np.exp(-t/parameters[1])
    residual_phase = phase - np.polyval(coefficients, t)
    residual = amp*np.exp(1j*residual_phase) - envelope
    offsets = np.fft.fftfreq(len(t), dt)
    band = (np.abs(offsets) >= .002) & (np.abs(offsets) <= .004)
    # Nominal published calibration. Filter transfer, run gain, mode coupling
    # and nonlinear response are NOT inferred by matching the best noise floor.
    gain = .16e6  # V/m, Supplement C, stated relative uncertainty 7%
    mass = .43e-6  # kg
    gamma = 1.84e-5  # s^-1; do not use the dimensionally wrong Q/f as bandwidth
    transfer = np.abs(inverse_susceptibility(fc+offsets, fc, mass, gamma))**2
    spectra = {}
    for name, window in [('rectangular', np.ones(len(t))), ('Hann', np.hanning(len(t)))]:
        psd = complex_periodogram(residual/gain, dt, window) * transfer
        spectra[name] = {
            'mean_band_ASD_N_per_sqrtHz': float(np.sqrt(np.mean(psd[band]))),
            'median_log2_diagnostic_ASD_N_per_sqrtHz': float(np.sqrt(np.median(psd[band])/np.log(2))),
            'interpretation': 'nominal reconstruction; median correction is not a confidence interval'}
    return {
        'header': header, 'samples_per_channel': len(t), 'sample_interval_s': dt,
        'duration_s': float(t[-1]), 'time_scale': 'export time label; UTC mapping not established',
        'reference_frequency_Hz': float(oscillator[0]), 'phase_slope_frequency_Hz': fc,
        'amplitude_fit_initial_V': float(parameters[0]), 'amplitude_fit_tau_s': float(parameters[1]),
        'published_ringdown_lower_bound_s': 109000,
        'fit_to_published_lower_bound_ratio': float(parameters[1]/109000),
        'fitted_flux_at_0_and_5_hours_mPhi0': [float(parameters[0]*np.exp(-u/parameters[1])/.43*1000)
                                              for u in [0, 18000]],
        'published_displacement_gain_V_per_m': gain, 'published_relative_gain_uncertainty': .07,
        'nominal_mass_kg': mass, 'nominal_damping_s_inv': gamma,
        'band_offset_Hz': [.002, .004], 'band_bins_both_sides': int(band.sum()),
        'window_reconstructions': spectra,
        'filter_settings_in_export': False,
        'pipeline_status': 'exploratory reconstruction, not a reproduced calibrated noise likelihood'}


def extract(directory):
    manifest = json.loads(MANIFEST.read_text())
    for entry in manifest['files']:
        verify(directory/entry['name'], entry)
    return {'sources': manifest['files'], 'cal': read_workbook(directory/'cal2020-source.xlsx'),
            'levitation': {name: lockin_summary(directory/name) for name in
                           ['fuchs-minus3cm.txt', 'fuchs-plus3p5cm.txt']}}


def analyze(inputs):
    a = np.array(inputs['cal']['widths'])
    fits = {}
    for name, j, reference, error in [('x', 1, 231., 9.), ('z', 3, 720., 79.)]:
        p, cov, chi2 = affine_wls(a[:, 0]*1e-3, a[:, j], a[:, j+1])
        slope, slope_error = float(p[1]), float(np.sqrt(cov[1, 1]))
        energy = RB87_MASS_KG / KB * slope**2 / 2  # (µm/s)^2 to pK factors cancel
        energy_error = RB87_MASS_KG / KB * slope*slope_error
        fits[name] = {
            'intercept_um': float(p[0]), 'source_width_slope_um_per_s': slope,
            'formal_slope_error_um_per_s': slope_error, 'diagonal_chi2': chi2, 'residual_dof': 4,
            'inferred_rms_conversion': float(1/np.sqrt(2)),
            'slope_div_sqrt2_um_per_s': slope/np.sqrt(2),
            'reconstructed_caption_energy_pK': energy, 'formal_energy_error_pK': energy_error,
            'caption_energy_pK': reference, 'caption_error_pK': error,
            'energy_difference_pK': energy-reference, 'error_difference_pK': energy_error-error,
            'literal_header_energy_pK': 2*energy,
            'baseline_numerically_reproduced_with_inferred_conversion':
                bool(abs(energy-reference) < 1 and abs(energy_error-error) < 1)}
    gain_scan = {str(tau): {str(order): float(lockin_power_gain(.004, tau, order))
                           for order in [1, 4, 8]} for tau in [0, 1, 10, 100, 1000]}
    return {
        'status': 'Gate A remains blocked; no new joint constraint or source-removal gain claimed',
        'cal_width_reconstruction': fits,
        'cal_raw_2019_52pK_baseline_reproduced': False,
        'cal_2020_to_2019_calibration_transfer_established': False,
        'cal_error_correlation_matrix_obtained': False,
        'levitation_nominal_to_best_published_ASD_ratios': {
            name: data['window_reconstructions']['Hann']['mean_band_ASD_N_per_sqrtHz']/5e-16
            for name, data in inputs['levitation'].items()},
        'filter_sensitivity_power_gain_at_4mHz': gain_scan,
        'maximum_tau_for_at_least_90percent_power_gain_at_4mHz_s': {
            str(order): float(np.sqrt(.9**(-1/order)-1)/(2*np.pi*.004)) for order in [1, 4, 8]},
        'empirical_filter_time_constant_bound_s': None,
        'original_SSB_witness_resolved': False,
        'new_empirical_spectrum_rejection': False}


def compare(expected, actual, location='root'):
    if isinstance(expected, dict):
        if expected.keys() != actual.keys():
            raise AssertionError('Keys differ at ' + location)
        for key in expected:
            compare(expected[key], actual[key], location+'.'+key)
    elif isinstance(expected, list):
        if len(expected) != len(actual):
            raise AssertionError('Lengths differ at ' + location)
        for i, (x, y) in enumerate(zip(expected, actual)):
            compare(x, y, location+f'[{i}]')
    elif isinstance(expected, float):
        if not np.isclose(expected, actual, rtol=3e-6, atol=1e-25):
            raise AssertionError(f'Numeric mismatch at {location}: {expected} != {actual}')
    elif expected != actual:
        raise AssertionError('Mismatch at ' + location)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-dir', type=Path)
    parser.add_argument('--fetch', action='store_true')
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    if args.fetch:
        if args.data_dir is None:
            parser.error('--fetch needs --data-dir')
        fetch(args.data_dir)
    if args.data_dir:
        inputs = extract(args.data_dir)
        if args.check:
            compare(json.loads(INPUTS.read_text()), inputs)
        else:
            INPUTS.write_text(encode(inputs))
    else:
        inputs = json.loads(INPUTS.read_text())
    result = analyze(inputs)
    if args.check:
        compare(json.loads(RESULT.read_text()), result)
    else:
        RESULT.write_text(encode(result))
    print(result['status'])


if __name__ == '__main__':
    main()
