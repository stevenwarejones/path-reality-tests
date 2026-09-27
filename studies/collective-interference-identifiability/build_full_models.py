#!/usr/bin/env python3
"""Reconstruct full histogram and find rational physical feasibility witnesses."""
from collections import Counter
from pathlib import Path
from math import isqrt,comb
import argparse,json
import numpy as np
import xarray as xr
from scipy.optimize import least_squares
from scipy.stats import beta
from fractions import Fraction as F
from sources import reconstruct
from certificate import check_endpoint
HERE=Path(__file__).resolve().parent

def full_histogram(cache):
    reconstruct(cache)
    a=xr.open_dataarray(cache/'selected/1D.nc').sel(key=3).isel(image_id=1).values
    a=np.nan_to_num(a[~np.all(np.isnan(a),axis=(1,2))]);imgs=a[:,8:20,8:20].reshape(-1,144)
    assert imgs.sum()==a.sum() and np.all(imgs.sum(axis=1)<=3)
    c=Counter(tuple(map(int,np.flatnonzero(v))) for v in imgs)
    return {'schema':1,'mode_count':144,'patterns':[[list(p),k] for p,k in sorted(c.items())]},imgs

def make_model(base,kind,seed):
    Q=10**7;R=1105;unit=[]
    for a in range(-R,R+1):
        b=isqrt(R*R-a*a)
        if a*a+b*b==R*R:
            unit.append(complex(a,b))
            if b:unit.append(complex(a,-b))
    unit=np.array(unit);nums=np.rint(np.sqrt(base)*Q).astype(np.int64)
    assert not nums[0].any() and not nums[-1].any()
    amp=np.stack([nums,np.roll(nums,1,axis=0),np.roll(nums,-1,axis=0)],axis=-1)/Q
    def res(v):
        T=(amp*np.exp(1j*np.c_[np.zeros(28),v.reshape(28,2)][:,None,:])).reshape(-1,3)
        g=T.conj().T@T;z=np.array([g[0,1],g[0,2],g[1,2]])
        return np.r_[z.real,z.imag]
    r=least_squares(res,np.random.default_rng(seed).uniform(-np.pi,np.pi,56),xtol=1e-13,ftol=1e-13,gtol=1e-13,max_nfev=2000)
    assert r.success
    phase=np.c_[np.ones(28)*R,unit[np.argmin(abs(np.exp(1j*r.x)[:,None]-unit[None,:]/R),axis=1)].reshape(28,2)]
    return {'kind':kind,'amplitude_denominator':Q,'phase_radius':R,'base_amplitude_numerators':nums.tolist(),
        'row_phase_numerators':[[[int(z.real),int(z.imag)] for z in row] for row in phase]}

def endpoint(n,k,tail,side):
    scale=10**18
    if side=='lower':p=F(0) if k==0 else F(int(np.floor(beta.ppf(float(tail),k,n-k+1)*scale)),scale)
    else:p=F(1) if k==n else F(int(np.ceil(beta.isf(float(tail),k+1,n-k)*scale)),scale)
    # Inverse tails propose only; move outward until integer arithmetic certifies.
    direction=-1 if side=='lower' else 1
    for _ in range(1000):
        try:check_endpoint(n,k,p,tail,side);return p
        except AssertionError:p+=F(direction,scale)
    raise AssertionError('Endpoint proposal failed')

def main():
    p=argparse.ArgumentParser();p.add_argument('--cache',type=Path,required=True);p.add_argument('--check-source',action='store_true');a=p.parse_args()
    h,imgs=full_histogram(a.cache);out=HERE/'results';text=json.dumps(h,separators=(',',':'))+'\n'
    if a.check_source:
        assert (out/'parity-histogram.json').read_text()==text
        print('Complete 2,999-shot parity histogram matches the pinned source.');return
    (out/'parity-histogram.json').write_text(text)
    K=sum(comb(144,k) for k in range(4));N=len(imgs);tail=F(1,640*K)
    region={'schema':1,'alphabet_size':K,'pattern_intervals':{str(k):[str(endpoint(N,k,tail,s)) for s in ['lower','upper']] for k in sorted({v for _,v in h['patterns']}|{0})},
        'bunch_lower':str(endpoint(N,93,F(1,320),'lower'))}
    (out/'full-region.json').write_text(json.dumps(region,indent=2,sort_keys=True)+'\n')
    for kind in ['boson','cluster']:
        base=np.zeros((28,12))
        if kind=='boson':
            w=json.loads((out/'physical-witness.json').read_text())
            for y,x,n in w['nonzero_cell_counts']:base[y,x-8]=n
            base[9:21]+=.1;base/=931+144*.1
        else:
            prob=(imgs.mean(axis=0)/3+1e-6).reshape(12,12);x=prob.sum(axis=0);x/=x.sum()
            g=np.exp(-.5*((np.arange(12)-6)/1.2)**2);g/=g.sum()
            base[8:20]=.96*(.8*g[:,None]*x[None,:]+.2*prob/prob.sum())
        m=make_model(base,kind,728 if kind=='boson' else 730)
        (out/('full-'+kind+'.json')).write_text(json.dumps(m,indent=2,sort_keys=True)+'\n')
    print('Proposals saved; run full_models.py for the exact acceptance decision.')
if __name__=='__main__':main()
