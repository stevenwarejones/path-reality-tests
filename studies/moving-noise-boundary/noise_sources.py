"""Reproduce compact measured-source summaries; external originals stay external."""
import hashlib
import json
from pathlib import Path
import sys
import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parent/'collapse-compatibility'
sys.path.insert(0,str(BASE))
from io_data import sodium, force, verify_sources
from baselines import harmonic_fit, weighted_line
from response import D


def summarize(directory):
    manifest = verify_sources(directory)
    scans,masses,weights = sodium(directory)
    rows,spectra,cells,formulas = force(directory)
    old = json.loads((BASE/'results/baselines.json').read_text())['data']
    old_scans = {x['scan']:x for x in old['sodium']['scans']}
    result=[]
    for scan in scans:
        fit=harmonic_fit(scan.position_m,scan.counts)
        design=np.column_stack([np.ones(len(scan.counts)),np.cos(2*np.pi*scan.position_m/D),np.sin(2*np.pi*scan.position_m/D)])
        # Independent normal-equation reproduction of the deposited baseline.
        beta=np.linalg.solve(design.T@design,design.T@scan.counts)
        visibility=np.hypot(*beta[1:])/beta[0]
        if abs(visibility-old_scans[scan.name]['visibility']) > 2e-10:
            raise ArithmeticError('Sodium baseline differs')
        result.append({'scan':scan.name, 'mean_counts':fit['mean_counts'],
                       'visibility':fit['visibility'], 'visibility_se':fit['visibility_se_gaussian'],
                       'optical_mean_probability':old_scans[scan.name]['transmission'],
                       'optical_visibility':old_scans[scan.name]['quantum_visibility'],
                       'ols_row_l1':np.abs(np.linalg.pinv(design)).sum(axis=1).tolist(),
                       'incident_rate_hz':float(scan.incident_rate_hz), 'integration_s':float(scan.integration_s),
                       'count_total':int(scan.counts.sum())})
    high=rows[:,0]>=100
    beta,cov,chi=weighted_line(rows[high,2]*1e7,rows[high,3]*1e19,rows[high,4]*1e19)
    slope=float(beta[1]*1e-12)
    if not np.isclose(slope,old['force']['weighted_y_only_fit']['B1'],rtol=2e-10):
        raise ArithmeticError('Thermal slope baseline differs')
    sensor=[]
    for temperature,table in spectra.items():
        f,y,error=table.T
        nav=int(round(np.median((y/error)**2)))
        if not np.allclose((y/error)**2,nav,rtol=2e-6):
            raise ValueError('PSD error convention differs')
        sensor.append({'temperature_mk':temperature, 'rows':len(f), 'minimum_hz':float(f.min()),
                       'maximum_hz':float(f.max()), 'minimum_psd_phi0_squared_per_hz':float(y.min()),
                       'averages':nav, 'minimum_step_hz':float(np.diff(f).min()),
                       'maximum_step_hz':float(np.diff(f).max()),
                       'peak_hz':float(f[np.argmax(y)]),
                       'first_row':table[0].tolist(), 'last_row':table[-1].tolist()})
    # B1 = C_phiF (4 kB k / omega0); includes source's Phi0-normalized units.
    calibration=cells['D23']*2*np.pi*cells['N20']/(4*cells['K20']*cells['L20'])
    return {'source_files':[i for i in manifest['files'] if i['name'] in
                            ['Na-Cluster-Interference.zip','PSD.zip','Final_data.xlsx']],
            'sodium':{'scans':result,'count_bins':sum(len(s.counts) for s in scans),
                      'mass_model_minmax_u':[float(masses.min()),float(masses.max())],
                      'mass_model_mean_u':float(masses@weights),
                      'support_is_measured_tail_bound':False},
            'sensor':{'spectra':sensor,'workbook_cells':cells,
                      'thermal_regression_slope':slope,'thermal_regression_chi_squared':chi,
                      'dc_flux_force_calibration_phi0_squared_per_N2':calibration,
                      'baseline_tolerance_relative':2e-10,
                      'measured_timestamps_available':False},
            'limitations':['optical baseline has unreconciled residuals',
                           'incident rate is not independently calibrated Bernoulli trial count',
                           'exposure and mass tails are not bounded by the archive',
                           'mechanical mode/control/DC calibration is not independently measured']}


def canonical(value):
    """Stable precision for compact results, not an uncertainty assertion."""
    if isinstance(value,dict): return {k:canonical(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)): return [canonical(v) for v in value]
    if isinstance(value,(float,np.floating)): return float(f'{value:.11g}')
    if isinstance(value,np.integer): return int(value)
    return value


def encode(value):
    return json.dumps(canonical(value),indent=2,sort_keys=True,allow_nan=False)+'\n'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()
