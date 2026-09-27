#!/usr/bin/env python3
"""A common observed quadrature flag in deposited 1 ms intervention records."""
import argparse
import json
from pathlib import Path

import numpy as np
from sklearn.decomposition import PCA

from cqm_analysis import close
from cqm_mechanical import HERE, intervention_clocks
from cqm_sources import acquire, digest

NORMS = np.array([2.3888843269295563, 2.7012968453752717])


def calibrate(ground, excited):
    if ground.shape != excited.shape or ground.shape != (2,32768) or not np.iscomplexobj(ground):
        raise ValueError('unexpected calibration schema')
    params=[]
    for g,e in zip(ground,excited):
        xy=lambda c: np.column_stack([c.real,c.imag])
        pca=PCA(n_components=2,svd_solver='full').fit(xy(np.r_[g,e]))
        if pca.transform(xy(e))[:,0].mean() < pca.transform(xy(g))[:,0].mean():
            pca.components_[0] *= -1
        g_projected=pca.transform(xy(g))[:,0]
        e_projected=pca.transform(xy(e))[:,0]
        mu_g,mu_e=float(g_projected.mean()),float(e_projected.mean())
        params.append({'ground_mean':mu_g, 'excited_mean':mu_e,
                       'threshold':(mu_g+mu_e)/2,
                       'prepared_G_flag_fraction':float(np.mean(g_projected>(mu_g+mu_e)/2)),
                       'prepared_E_flag_fraction':float(np.mean(e_projected>(mu_g+mu_e)/2))})
    return params


def normalized_flags(projected, params, norms):
    if projected.ndim != 3 or projected.shape[0] != 2 or projected.shape[-1] != 32768 or not np.isfinite(projected).all():
        raise ValueError('invalid projected readout record')
    output=[]
    for q in range(2):
        g,e=params[q]['ground_mean'],params[q]['excited_mean']
        if e<=g or norms[q]<=0:
            raise ValueError('invalid calibration orientation/normalization')
        u=(projected[q]*norms[q]-g)/(e-g)
        output.append(u>.5)
    return np.array(output)


def counts(flags):
    if flags.ndim!=3 or flags.shape[0]!=2 or flags.dtype!=bool:
        raise ValueError('expected paired boolean acquisitions')
    rows=[]
    for run in range(flags.shape[1]):
        a,b=flags[:,run]
        initial=~a[:-1]&~b[:-1]
        observed=int((initial&a[1:]&b[1:]).sum())
        rows.append({'acquisition':run,'samples':len(a), 'A_flags':int(a.sum()), 'B_flags':int(b.sum()),
                     'paired_flags':int((a&b).sum()), 'initial_both_unflagged':int(initial.sum()),
                     'next_1ms_paired_flags':observed})
    keys=[k for k in rows[0] if k!='acquisition']
    total={k:sum(row[k] for row in rows) for k in keys}
    total['paired_flag_fraction']=total['paired_flags']/total['samples']
    total['conditional_1ms_pair_fraction']=total['next_1ms_paired_flags']/total['initial_both_unflagged']
    return {'total':total,'acquisitions':rows}


def analyze(cache):
    base=cache/'mechanical/Zenodo'
    verified=[]
    for number in [617,618,619]:
        a=(base/f'Fig2/{number}.npy').read_bytes()
        b=(base/f'Fig6/{number}.npy').read_bytes()
        if a!=b:
            raise ValueError('control notebooks do not share identical calibration bytes')
        verified.append({'calibration':number,'sha256':digest(a)})
    params=calibrate(np.load(base/'Fig2/617.npy',allow_pickle=False), np.load(base/'Fig2/618.npy',allow_pickle=False))
    result={}
    for label,path,norms in [('PT_off_774','Fig2/774_r.npy',np.ones(2)),
                             ('PT_on_782','Fig2/782_r.npy',np.ones(2)),
                             ('controlled_shock_833','Fig6/833_r.npy',NORMS)]:
        result[label]=counts(normalized_flags(np.load(base/path,allow_pickle=False),params,norms))
    accel=np.load(base/'Fig6/832.npy',allow_pickle=False)
    if accel.shape!=(3,700000) or not np.isfinite(accel).all():
        raise ValueError('unexpected shock accelerometer record')
    time,voltage,trigger=accel
    edges=np.flatnonzero(trigger>1)
    if not len(edges) or np.any(np.diff(time)<=0):
        raise ValueError('missing shock record trigger')
    return {'status':'Common calibrated observed flag exists; physical intervention transfer has not been fitted',
            'definition':'u=(saved projected quadrature times notebook normalization - prepared G mean)/(prepared E mean - prepared G mean); flag u>0.5',
            'calibrations_identical':verified, 'calibration':params,
            'calibration_scope':'Prepared-reference flag fractions mix preparation, dynamics and separation errors; not pure readout-error calibration.',
            'clocks':intervention_clocks(cache), 'records':result,
            'shock_accelerometer':{'samples':len(time),'trigger_time_s':float(time[edges[0]]),
                                  'sample_interval_s':float(np.median(np.diff(time)))},
            'limitations':['No raw IQ for the Fig2 projected records; notebook PCA orientation is taken as deposited.',
                           'Matching projection and sampling is not a bound on state-dependent backaction or run drift.',
                           'One off record, one on record and fifty shock records are not independent randomized interventions.',
                           'The 3 us continuous record is not downsampled or declared physically equivalent.',
                           'The shock accelerometer record is not an online phase log for the separate 3100 dwell acquisitions.']}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache',type=Path,required=True)
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    acquire(args.cache,verify_only=True,mechanical_prediction=True,mechanical_interventions=True)
    result=analyze(args.cache)
    output=HERE/'results/intervention-observable.json'
    if args.check:
        close(result,json.loads(output.read_text()))
        print('Common intervention-observable audit reproduced')
    else:
        output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
        print(output)
