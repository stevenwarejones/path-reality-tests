#!/usr/bin/env python3
"""Joint July-7 Hamiltonian diagnostic, with explicit limits on calibration transfer.

Numerical fitted-point diagnostics are not full-class exclusion certificates.
The exact conditional Hamiltonian envelope is in control_design.py.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,io,json,re,zipfile
import numpy as np
from numpy.polynomial.hermite import hermgauss
from scipy.special import jv,xlogy
from scipy.optimize import minimize
from control_audit import audit
HERE=Path(__file__).resolve().parent
PARAMETERS=[121.94618326893601,106.58118068399486,-0.374281912946175,0.08116996089253452,0.9627199765171192]


def prediction(parameters,times,order=61):
    """Infinite separable nearest-neighbor walk, Gaussian shared scale and fixed loss.

Parameters: Jx/h,Jy/h in Hz; additive time offset in ms; relative scale SD;
per-atom survival. The input sites within the observed 12x12 crop are (4,5),
(5,5),(6,5), independently reconstructed by source_geometry().
"""
    jx,jy,off,sig,eta=parameters;z,w=hermgauss(order);w/=np.sqrt(np.pi)
    t=(np.array(times)[:,None]+off)*(1+sig*np.sqrt(2)*z[None,:])/1000
    px=jv(np.arange(12)[None,None,:]-5,4*np.pi*jx*t[:,:,None])**2
    py=jv(np.arange(12)[None,None,None,:]-np.array([4,5,6])[None,None,:,None],4*np.pi*jy*t[:,:,None,None])**2
    q=eta*np.einsum('tzjy,tzx,z->tjyx',py,px,w).reshape(len(times),3,144)
    q=np.concatenate([q,1-q.sum(axis=2,keepdims=True)],axis=2)
    dx=px.sum(axis=2)**3-3*px.sum(axis=2)*(px**2).sum(axis=2)+2*(px**3).sum(axis=2)
    D=eta**3*np.einsum('tzy,tzy,tzy,tz,z->t',py[:,:,0],py[:,:,1],py[:,:,2],dx,w)
    assert np.all(q>=0) and np.allclose(q.sum(axis=2),1)
    return q,D


def cases():
    cs=json.loads((HERE/'results/matched-regions.json').read_text())['cases']
    return [c for c in cs if c['time_ms']<4]


def fit(cs):
    cts=np.array([[s['counts'] for s in c['singletons']] for c in cs]);ts=[c['time_ms'] for c in cs]
    def loss(p):return -np.sum(xlogy(cts,np.maximum(prediction(p,ts,15)[0],1e-250)))
    fits=[minimize(loss,[128,119,-.4,s,.95],bounds=[(50,220),(50,220),(-.49,.5),(0,.3),(.8,.999)],method='Nelder-Mead',options={'maxiter':3000,'xatol':1e-8}) for s in [0,.02,.1]]
    best=min(fits,key=lambda r:r.fun)
    return {'parameters':best.x.tolist(),'solver_success':bool(best.success),'meaning':'Local point proposal only; fit uses singleton records, not many-body outcomes.'}


def evaluate(parameters=PARAMETERS):
    cs=cases();assert len(cs)==15
    cts=np.array([[s['counts'] for s in c['singletons']] for c in cs]);ts=[c['time_ms'] for c in cs]
    q,D=prediction(parameters,ts);q2,D2=prediction(parameters,ts,81)
    convergence=float(max(np.max(abs(q-q2)),np.max(abs(D-D2))))
    assert convergence<1e-10 # numerical convergence check, not an interval certificate
    N=cts.sum(axis=2,keepdims=True);dev=2*np.sum(xlogy(cts,cts/np.maximum(N*q,1e-250)))
    records=[]
    for i,c in enumerate(cs):
        failures=[]
        for j,s in enumerate(c['singletons']):
            for k,(count,pr) in enumerate(zip(s['counts'],q[i,j])):
                lo,hi=map(lambda x:float(F(x)),s['intervals'][str(count)])
                if not lo<=pr<=hi:failures.append({'input':j,'output':k,'count':count,'probability':float(pr),'lower':lo,'upper':hi})
        records.append({'time_ms':c['time_ms'],'singleton_cell_interval_failures':failures,
            'C2_bunch_probability':float(2*D[i]),'boson_bunch_probability':float(6*D[i]),'observed_bunch_fraction':c['bunch_count']/c['many_shots']})
    return {'model':'One shared infinite separable NN Hamiltonian and Gaussian scale law for July 7 only.',
        'parameters_Jx_Hz_Jy_Hz_offset_ms_scale_sd_survival':parameters,'singleton_shots':int(N.sum()),
        'singleton_deviance_descriptive':float(dev),'quadrature_61_vs_81_max_difference':convergence,
        'singleton_cell_interval_failure_count':sum(len(c['singleton_cell_interval_failures']) for c in records),
        'all_records_compatible':False,'cases':records,
        'interpretation':'This fitted point fails calibration constraints. Neither its C2 prediction deficit nor local optimizer termination excludes the Hamiltonian class. No aggregate-adequate joint model is claimed.'}


def source_geometry(cache):
    import xarray as xr
    audit(cache) # pins original archive, notebook and all input members
    result=[]
    for j in range(3):
        d=xr.open_dataarray(cache/'controls'/f'3NN_{j}.nc')
        for t in [c['time_ms'] for c in cases()]:
            a=d.sel(key=t).values;ys,xs=np.where(~np.isnan(a[1,0]));site=np.argwhere(np.nan_to_num(a[0,0]))-[ys.min(),xs.min()]
            assert site.tolist()==[[4+j,5]]
        result.append([4+j,5])
    return result


def apparatus_audit(cache):
    geometry=source_geometry(cache);inventory={}
    with zipfile.ZipFile(cache/'boson-1.0.4.zip') as z:
        for name in ['H_resampling.py','H_spectroscopic.npy','2212_H.npy','U_2ms_calibrated.npy','U_late_time_calibrated.npy','UtoJ_conversion_Hz.npy']:
            b=z.read('boson-1.0.4/inference/'+name)
            inventory[name]={'sha256':hashlib.sha256(b).hexdigest(),'bytes':len(b)}
            if name.endswith('.npy') and name!='UtoJ_conversion_Hz.npy':
                a=np.load(io.BytesIO(b),allow_pickle=False);inventory[name]['shape']=list(a.shape)
        script=z.read('boson-1.0.4/inference/H_resampling.py').decode()
        assert '220822 to 220823' in script and 'coordinates on 220805' in script
        for token in ['popt = np.array','perr = np.array','t_offset_s_fixed = -0.0005366067918397108']:assert token in script
    notebook=(cache/'atomic-notebook.ipynb').read_text()
    assert '220812' in notebook and 'had to adjust mirrors on this day' in notebook
    return {'independent_input_sites_in_crop':geometry,'archive_members':inventory,
        'spectroscopic_model_available':True,'published_resampling_reference':'August 5 coordinate frame and August 22-23 alignment comments; not a July transfer certificate.',
        'can_share_across_July2_July7_August_without_extra_assumptions':False,
        'missing_bridge':'A confidence set for the July lattice potential, alignment, effective time and shot-scale noise, or a validated transfer of the later calibration. Published fit errors are not a measured shot-noise law.'}


def close(a,b):
    if isinstance(a,float):assert np.isclose(a,b,rtol=2e-9,atol=2e-11),(a,b)
    elif isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:close(a[k],b[k])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for x,y in zip(a,b):close(x,y)
    else:assert a==b,(a,b)


def main():
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');p.add_argument('--fit',action='store_true');p.add_argument('--cache',type=Path);a=p.parse_args()
    if a.fit:print(json.dumps(fit(cases()),indent=2));return
    r=evaluate();out=HERE/'results/dynamics-diagnostic.json'
    if a.check:close(r,json.loads(out.read_text()))
    else:out.write_text(json.dumps(r,sort_keys=True,indent=2)+'\n')
    if a.cache:
        ar=apparatus_audit(a.cache);out=HERE/'results/apparatus-audit.json'
        if a.check:assert ar==json.loads(out.read_text())
        else:out.write_text(json.dumps(ar,sort_keys=True,indent=2)+'\n')
    print(json.dumps({k:v for k,v in r.items() if k!='cases'},indent=2))
if __name__=='__main__':main()
