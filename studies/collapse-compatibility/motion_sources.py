"""Acquisition metadata and semantic normalization for external Horizons inputs."""
import json
from pathlib import Path
import numpy as np
from io_data import sodium, force


def normalize_horizons(raw):
    """Remove ONLY the query-generation timestamp, retaining model/version headers."""
    return ('\n'.join(x for x in raw.decode('utf-8').splitlines() if not x.startswith('Ephemeris / API_USER '))+'\n').encode()


def horizons(path,target):
    text=Path(path).read_text()
    if f'({target})' not in text or 'Reference frame : ICRF' not in text or 'Output units    : KM-S' not in text or 'JDUT' not in text:
        raise ValueError('Unexpected trajectory target/frame/units/time')
    rows=[]
    for line in text.split('$$SOE')[1].split('$$EOE')[0].strip().splitlines():
        c=line.split(','); rows.append([float(c[0])]+[1000*float(x) for x in c[2:8]])
    a=np.array(rows)
    if not np.isfinite(a).all() or np.any(np.diff(a[:,0])<=0): raise ValueError('Invalid trajectory')
    return {'jd_utc':a[:,0].tolist(),'position_m':a[:,1:4].tolist(),'velocity_m_s':a[:,4:].tolist(),
            'source':[x.strip() for x in text.splitlines() if x.startswith(('Target body name:','Center body name:'))],
            'prediction_after_2016_04_14':'nominal' in text and target=='-141043'}


def source_summary(directory):
    directory=Path(directory)
    scans,mass,mw=sodium(directory); rows,spectra,_,_=force(directory)
    return {'evidence':'observed_metadata_and_derived_summaries',
            'sodium':{'mass_quadrature_u':mass.tolist(),'normalized_mass_weights':mw.tolist(),
                      'powers_by_scan':{s.name:s.powers_w.tolist() for s in scans},
                      'integration_times_s':sorted(set(s.integration_s for s in scans)),
                      'velocity_mean_m_s':158.,'velocity_sigma_m_s':9.,
                      'velocity_evidence':'author analysis parameters; paper describes TOF measurement; event-resolved velocities not deposited'},
            'sensor':{'frequency_grid_hz':spectra[100][:,0].tolist(),
                      'reported_sample_rate_hz':100000.,'reported_fft_length':2**22,
                      'window':'Blackman; symmetric/periodic convention not specified',
                      'averages_by_temperature':{str(int(t)):int(round(np.median((s[:,1]/s[:,2])**2))) for t,s in spectra.items()},
                      'epoch_and_attitude':'not in released PSD/workbook; no time-modulation claim'},
            'trajectory':{'lpf':horizons(directory/'lpf_201702.txt','-141043'),
                          'earth_lpf':horizons(directory/'earth_201702.txt','399'),
                          'earth_sodium':horizons(directory/'earth_sodium.txt','399')}}
