#!/usr/bin/env python3
"""Conditional forcing-to-readout transfer and exact calibration-support ambiguity."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.signal import lfilter

from cqm_analysis import close
from cqm_mechanical import HERE
from cqm_sources import acquire

TAUS=(.001,.003,.01,.03,.1,.3,1.)
PHASE_EDGES=np.array([0,.05,.1,.2,.4,.8,1.5])
DT=.001


def causal_filter(values, dt, tau):
    if dt<=0 or tau<=0 or not np.isfinite(values).all():
        raise ValueError('invalid filter input')
    r=np.exp(-dt/tau)
    return lfilter([1-r],[1,-r],values)


def read_forcing(path, center=False):
    a=np.load(path,allow_pickle=False)
    if a.ndim!=2 or a.shape[0]!=3 or not np.isfinite(a).all():
        raise ValueError('invalid acceleration/trigger schema')
    t,v,trigger=a
    dt=float(t[1]-t[0])
    if dt<=0 or not np.allclose(np.diff(t),dt,rtol=0,atol=1e-10):
        raise ValueError('nonuniform forcing clock')
    edges=np.flatnonzero(trigger>1)
    if not len(edges):
        raise ValueError('no acquisition trigger')
    if center:
        if edges[0] == 0:
            raise ValueError('no pre-trigger baseline')
        v=v-v[:edges[0]].mean()
    return t,5*v,float(t[edges[0]]),dt


def reservoir(forcing, qubit_times, tau):
    t,a,trigger,dt=forcing
    filtered=causal_filter(a,dt,.001)
    energy=causal_filter(filtered**2,dt,tau)
    target=qubit_times+trigger
    # Last available accelerometer sample; never interpolate with future input.
    index=np.searchsorted(t,target,side='right')-1
    if np.any(index<0) or np.any(target>t[-1]):
        raise ValueError('forcing record does not cover prediction')
    return energy[index]


def transition(energy, parameters):
    b,k,d=parameters
    if min(b,k,d)<0 or d<=0 or np.any(energy<0):
        raise ValueError('nonphysical rate')
    up=b+k*energy
    rate=up+d
    decay=np.exp(-rate*DT)
    steady=up/rate
    p01=steady*(-np.expm1(-rate*DT))
    p11=steady+(1-steady)*decay
    return p01,p11


def transition_loss(flags, energy, parameters):
    p01,p11=transition(energy,parameters)
    p=np.where(flags[:-1],p11,p01)
    p=np.clip(p,1e-15,1-1e-15)
    return float(-np.sum(flags[1:]*np.log(p)+(1-flags[1:])*np.log1p(-p)))


def fit(controls):
    grid=[]
    for tau in TAUS:
        energies=[reservoir(r['forcing'],r['time'][:-1],tau) for r in controls]
        # Unit rescaling from training covariates only, shared across both channels.
        scale=float(max(np.max(e[:4095]) for e in energies))
        energies=[e/scale for e in energies]
        params=[];diagnostics=[]
        for q in range(2):
            def objective(theta):
                return sum(transition_loss(r['flags'][q,0,:4096],e[:4095],np.exp(theta)) for r,e in zip(controls,energies))
            bounds=[(-12,9),(-12,11),(0,11)]
            candidates=[]
            for b,k,d in [(1,100,1000),(.1,10,5000),(10,1000,10000),(.1,1000,500)]:
                solved=minimize(objective,np.log([b,k,d]),method='L-BFGS-B',bounds=bounds,
                                options={'maxiter':500,'ftol':1e-12,'gtol':1e-6})
                candidates.append(solved)
            best=min(candidates,key=lambda c:c.fun)
            params.append(np.exp(best.x))
            diagnostics.append({'success':bool(best.success),'message':str(best.message),
                'training_nll':float(best.fun),'at_search_boundary':bool(any(abs(v-lo)<1e-4 or abs(v-hi)<1e-4 for v,(lo,hi) in zip(best.x,bounds)))})
        grid.append({'tau_s':tau,'scale':scale,'params':np.array(params),'diagnostics':diagnostics,
                     'nll':sum(d['training_nll'] for d in diagnostics)})
    best=min(grid,key=lambda r:r['nll'])
    return best, [{'tau_s':r['tau_s'],'training_nll':r['nll']} for r in grid]


def forward(initial, energy, params):
    if np.shape(initial)!=(2,):
        raise ValueError('need one paired initial state')
    transitions=[transition(energy,p) for p in params]
    probabilities=np.empty((2,len(energy)+1));probabilities[:,0]=initial
    for q,(p01,p11) in enumerate(transitions):
        for j in range(len(energy)):
            probabilities[q,j+1]=(1-probabilities[q,j])*p01[j]+probabilities[q,j]*p11[j]
    onset=np.array([(1-probabilities[q,:-1])*transitions[q][0] for q in range(2)])
    return probabilities,onset,transitions


def run_duration_counts(flags):
    result=np.zeros(4,dtype=int);censored=0
    changes=np.diff(np.r_[False,flags,False].astype(int))
    for start,end in zip(np.flatnonzero(changes==1),np.flatnonzero(changes==-1)):
        if start==0 or end==len(flags):
            censored+=1
            continue
        length=end-start
        result[0 if length==1 else 1 if length==2 else 2 if length<=5 else 3]+=1
    return result,censored


def duration_predictions(onset, transitions):
    p01,p11=transitions
    n=len(onset)
    result=np.zeros(4)
    # An onset at sample j+1 uses transition j; subsequent survival uses j+1,... .
    for length in range(1,6):
        count=n-length
        if count<=0:
            continue
        prob=onset[:count].copy()
        for offset in range(1,length):
            prob*=p11[offset:offset+count]
        prob*=1-p11[length:length+count]
        result[0 if length==1 else 1 if length==2 else 2]+=prob.sum()
    # 6+ complete runs: at least five survivals, then at least one exit by the end.
    if n>6:
        count=n-6
        prob=onset[:count].copy()
        for offset in range(1,6):
            prob*=p11[offset:offset+count]
        log_tail=np.r_[np.cumsum(np.log(np.clip(p11,1e-300,1))[::-1])[::-1],0.]
        result[3]=np.sum(prob*(-np.expm1(log_tail[6:6+count])))
    return result


def summarize(flags, energy, params, time, origin=None):
    if flags.ndim!=3 or flags.shape[0]!=2 or flags.shape[-1]!=len(energy)+1:
        raise ValueError('acquisition shape mismatch')
    observed_onsets=(~flags[:,:,:-1])&flags[:,:,1:]
    lag_values=np.arange(-10,11)
    observed_lag=np.zeros(21,dtype=int);predicted_lag=np.zeros(21)
    observed_duration=np.zeros((2,4),dtype=int);predicted_duration=np.zeros((2,4));censored=np.zeros(2,dtype=int)
    pred=np.zeros((2,flags.shape[-1]));joint=np.zeros(flags.shape[-1]);rows=[]
    for run in range(flags.shape[1]):
        p,onset,tr=forward(flags[:,run,0],energy,params)
        pred+=p;joint+=p[0]*p[1]
        rows.append({'acquisition':run,'samples_after_initial':flags.shape[-1]-1,
                     'observed_flags':flags[:,run,1:].sum(1).tolist(),
                     'predicted_flags':p[:,1:].sum(1).tolist(),
                     'observed_paired':int(np.all(flags[:,run,1:],axis=0).sum()),
                     'predicted_paired':float((p[0,1:]*p[1,1:]).sum())})
        for q in range(2):
            c,edge=run_duration_counts(flags[q,run]);observed_duration[q]+=c;censored[q]+=edge
            predicted_duration[q]+=duration_predictions(onset[q],tr[q])
        for i,lag in enumerate(lag_values):
            if lag>=0:
                a=slice(0,len(energy)-lag);b=slice(lag,len(energy))
            else:
                a=slice(-lag,len(energy));b=slice(0,len(energy)+lag)
            observed_lag[i]+=int(np.sum(observed_onsets[0,run,a]&observed_onsets[1,run,b]))
            predicted_lag[i]+=float(onset[0,a]@onset[1,b])
    n=flags.shape[1]*(flags.shape[-1]-1)
    total={'exposure_samples':n,'observed_flags':flags[:,:,1:].sum((1,2)).tolist(),
           'predicted_flags':pred[:,1:].sum(1).tolist(),
           'observed_paired':int(np.all(flags[:,:,1:],axis=0).sum()),'predicted_paired':float(joint[1:].sum())}
    total['observed_paired_per_million']=total['observed_paired']/n*1e6
    total['predicted_paired_per_million']=total['predicted_paired']/n*1e6
    phase=[]
    if origin is not None:
        start=int(np.argmin(abs(time-origin)))
        cycles=(len(time)-start)//1500
        index=np.arange(start,start+cycles*1500)
        coordinates=(index-start)%1500/1000
        for low,high in zip(PHASE_EDGES[:-1],PHASE_EDGES[1:]):
            selected=index[(coordinates>=low)&(coordinates<high)]
            observed=int(np.all(flags[:,:,selected],axis=0).sum())
            expected=float(joint[selected].sum());exposure=len(selected)*flags.shape[1]
            phase.append({'phase_s':[float(low),float(high)],'exposure_samples':exposure,
                          'observed_flags':flags[:,:,selected].sum((1,2)).tolist(),
                          'predicted_flags':pred[:,selected].sum(1).tolist(),
                          'observed_paired':observed,'predicted_paired':expected,
                          'observed_paired_per_million':observed/exposure*1e6,
                          'predicted_paired_per_million':expected/exposure*1e6})
    return {'total':total,'acquisitions':rows,'protocol_phase':phase,
            'onset_lags_ms':lag_values.tolist(),'observed_onset_pairs':observed_lag.tolist(),
            'predicted_onset_pairs':predicted_lag.tolist(),
            'onset_pair_exposures':(flags.shape[1]*(len(energy)-abs(lag_values))).tolist(),
            'observed_onset_pairs_per_million':(observed_lag/(flags.shape[1]*(len(energy)-abs(lag_values)))*1e6).tolist(),
            'predicted_onset_pairs_per_million':(predicted_lag/(flags.shape[1]*(len(energy)-abs(lag_values)))*1e6).tolist(),
            'duration_labels':['1 ms','2 ms','3-5 ms','6+ ms'],
            'observed_complete_flag_runs':observed_duration.tolist(),
            'predicted_complete_flag_runs':predicted_duration.tolist(),
            'observed_duration_fractions':(observed_duration/np.maximum(observed_duration.sum(1,keepdims=True),1)).tolist(),
            'predicted_duration_fractions':(predicted_duration/np.maximum(predicted_duration.sum(1,keepdims=True),np.finfo(float).tiny)).tolist(),
            'boundary_censored_observed_runs':censored.tolist()}


def analyze(cache, center=False):
    base=cache/'mechanical/Zenodo'
    controls=[]
    for num,label in [(562,'PT_on'),(567,'PT_off')]:
        time=np.load(base/f'Fig2/{num}_t.npy',allow_pickle=False)
        r=np.load(base/f'Fig2/{num}_r.npy',allow_pickle=False)
        if r.shape!=(2,1,8192) or not np.isfinite(r).all() or not np.array_equal(time,np.arange(8192)*DT):
            raise ValueError('unexpected matched control')
        controls.append(dict(label=label,flags=r>0,time=time,forcing=read_forcing(base/f'Fig2/{num}_v.npy',center=center)))
    best,grid=fit(controls);tau,scale,params=best['tau_s'],best['scale'],best['params']
    control_energy=[reservoir(r['forcing'],r['time'][:-1],tau)/scale for r in controls]
    cap=max(float(e.max()) for e in control_energy)
    validations={}
    for r,e in zip(controls,control_energy):
        validation=r['flags'][:,:,4096:]
        validations[r['label']]=summarize(validation,e[4096:],params,r['time'][4096:])
        validations[r['label']]['conditional_transition_nll']=[transition_loss(validation[q,0],e[4096:],params[q]) for q in range(2)]
    time=np.load(base/'Fig6/833_t.npy',allow_pickle=False)
    flags=np.load(base/'Fig6/833_r.npy',allow_pickle=False)>0
    forcing=read_forcing(base/'Fig6/832.npy',center=center)
    energy=reservoir(forcing,time[:-1],tau)/scale
    linear=summarize(flags,energy,params,time,origin=1.7-forcing[2])
    saturated=summarize(flags,np.minimum(energy,cap),params,time,origin=1.7-forcing[2])
    return {'pretrigger_baseline_removed':center,'status':'Conditional transfer test, not an adequate physical countermodel certificate',
            'model':'Nonnegative birth/death generator conditional on reference forcing; saved quadrature flags treated as states without independently calibrated emission/backaction law.',
            'selected_reservoir_tau_s':tau,'training_forcing_scale':scale,'parameters_b_k_d_per_s':params.tolist(),
            'fit_diagnostics':best['diagnostics'],'training_grid':grid,
            'support':{'control_energy_cap':cap,'shock_energy_max':float(energy.max()),
                       'shock_fraction_above_control_cap':float(np.mean(energy>cap)),
                       'control_transition_laws_identical':all(np.array_equal(e,np.minimum(e,cap)) for e in control_energy)},
            'calibration_validation':validations,'shock_predictions':{'linear':linear,'saturated_extension':saturated},
            'limitations':['Shock forcing 832 is a separate reference run, not a measured input for every 833 acquisition.',
                           'Top-plate acceleration is not a calibrated absorbed-power measurement at the chip.',
                           'No intrinsic emission, preparation or backaction calibration establishes a physical-state interpretation.',
                           'No global optimizer certificate or dependence-calibrated adequacy region.',
                           'Both extensions have the exact same PT calibration laws, but adequacy must be assessed separately.']}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache',type=Path,required=True)
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    acquire(args.cache,verify_only=True,mechanical_transfer=True)
    result={'status':'Initial uncentered test and subsequent pre-trigger-centered diagnostic; neither is a full physical adequacy certificate',
            'initial_uncentered':analyze(args.cache),
            'pretrigger_centered_diagnostic':analyze(args.cache,center=True)}
    output=HERE/'results/transfer.json'
    if args.check:
        close(result,json.loads(output.read_text()))
        print('Conditional intervention transfer reproduced')
    else:
        output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
        print(output)
