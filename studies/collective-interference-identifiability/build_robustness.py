#!/usr/bin/env python3
"""Propose robust-region endpoints and physical models; exact checks are separate."""
from pathlib import Path
from collections import Counter
from fractions import Fraction as F
from math import comb,isqrt,floor,ceil,log
import argparse,json
import numpy as np
import xarray as xr
from scipy.optimize import brentq,minimize,least_squares
from scipy.stats import beta
from build_full_models import make_model
from certificate import check_endpoint
from control_audit import audit
from robustness import check_chernoff,generic_transport
HERE=Path(__file__).resolve().parent


def endpoint(n,k,tail,side):
    if side=='lower' and k==0:return F(0)
    if side=='upper' and k==n:return F(1)
    value=beta.ppf(float(tail),k,n-k+1) if side=='lower' else beta.isf(float(tail),k+1,n-k)
    grid=10**10 if value>1e-6 else 10**18
    p=F(floor(value*grid) if side=='lower' else ceil(value*grid),grid)
    for _ in range(1000):
        try:check_endpoint(n,k,p,tail,side);return p
        except AssertionError:p+=F(-1 if side=='lower' else 1,grid)
    raise AssertionError('Could not certify an outward endpoint')

def chernoff_endpoints(n,k,tail):
    def bound(p):
        value=n*log(n)
        if k:value+=k*(log(p)-log(k))
        if n-k:value+=(n-k)*(np.log1p(-p)-log(n-k))
        return value-log(float(tail))
    grid=10**6;q=k/n
    lo=F(0) if k==0 else F(floor(brentq(bound,1e-100,q,xtol=1e-14)*grid),grid)
    hi=F(1) if k==n else F(ceil(brentq(bound,q,1-1e-15,xtol=1e-14)*grid),grid)
    check_chernoff(n,k,lo,tail,'lower');check_chernoff(n,k,hi,tail,'upper')
    return [str(lo),str(hi)]

def independent_model(prob,seed=728):
    Y,W,_=prob.shape;Q=10**8;R=1105;unit=[]
    for a in range(-R,R+1):
        b=isqrt(R*R-a*a)
        if a*a+b*b==R*R:
            unit.append(complex(a,b))
            if b:unit.append(complex(a,-b))
    unit=np.array(unit);nums=np.rint(np.sqrt(prob)*Q).astype(np.int64);amp=nums/Q
    def res(v):
        T=(amp*np.exp(1j*np.c_[np.zeros(Y),v.reshape(Y,2)][:,None,:])).reshape(-1,3);g=T.conj().T@T;z=np.array([g[0,1],g[0,2],g[1,2]])
        return np.r_[z.real,z.imag]
    r=least_squares(res,np.random.default_rng(seed).uniform(-np.pi,np.pi,2*Y),xtol=1e-12,ftol=1e-12,gtol=1e-12,max_nfev=1000)
    assert np.isfinite(r.x).all() # Only the exact contraction check decides admissibility.
    ph=np.c_[np.ones(Y)*R,unit[np.argmin(abs(np.exp(1j*r.x)[:,None]-unit[None,:]/R),axis=1)].reshape(Y,2)]
    model={'kind':'cluster','width':W,'amplitude_denominator':Q,'phase_radius':R,'amplitudes':nums.tolist(),
        'phases':[[[int(z.real),int(z.imag)] for z in row] for row in ph]}
    generic_transport(model);return model

def source_cases(cache):
    a=audit(cache) # validates hashes, schema and the non-independent duplicate record
    many=xr.open_dataarray(cache/'controls/3NN.nc');singles=[xr.open_dataarray(cache/f'controls/3NN_{j}.nc') for j in range(3)]
    def read(d,k):
        v=d.sel(key=k).values;v=v[:,~np.all(np.isnan(v),axis=(0,2,3))];mask=~np.isnan(v[1,0]);return np.nan_to_num(v),mask
    cases=[]
    for target in [float(k) for k in many.key.values if k!=2.45]:
        key=float(many.key.values[np.argmin(abs(many.key.values-target))]);v,mask=read(many,key)
        ys,xs=np.where(mask);ya,yb=ys.min(),ys.max()+1;xa,xb=xs.min(),xs.max()+1;v=v[1,:,ya:yb,xa:xb];W=v.shape[-1];assert v.shape[1]==W
        N=len(v);M=W*W;pats=Counter(tuple(map(int,np.flatnonzero(im))) for im in v.reshape(N,M));ss=[]
        for ds in singles:
            ar,ma=read(ds,key);assert np.array_equal(ma,mask);ar=ar[1,:,ya:yb,xa:xb];n=len(ar);ct=ar.sum(axis=0).astype(int)
            assert np.all(ar.sum(axis=(1,2))<=1)
            ss.append({'shots':n,'counts':list(map(int,np.r_[ct.reshape(M),n-ct.sum()]))})
        cases.append({'time_ms':key,'width':W,'many_shots':N,'bunch_count':sum(k for p,k in pats.items() if len(p)==3 and len({o//W for o in p})==1),
            'patterns':[[list(p),k] for p,k in sorted(pats.items())],'singletons':ss})
    return cases

def drift_models(read):
    c=read('cell-certificate.json');m=read('full-boson.json');base=np.array(m['base_amplitude_numerators'],float)**2/m['amplitude_denominator']**2
    q0=np.r_[base[9:21].sum(axis=1),1-base.sum()];r=base[9:21]/q0[:12,None]
    bounds=np.array([[float(F(v)) for v in p] for p in c['row_intervals']]);lower=float(F(read('full-region.json')['bunch_lower']))
    h=read('parity-histogram.json')['patterns'];hr=np.array([sum(k for p,k in h if any(o//12==y for o in p)) for y in range(12)])
    # Search in stricter CP row intervals; final exact region uses Chernoff envelopes.
    hlo=beta.ppf(1/15360,hr,3000-hr);hhi=beta.isf(1/15360,hr+1,2999-hr)
    pc=np.stack([base,np.roll(base,1,axis=0),np.roll(base,-1,axis=0)],axis=-1);qc=pc.sum(axis=1)
    rc=np.divide(pc,qc[:,None,:],out=np.ones_like(pc)/12,where=qc[:,None,:]>0)
    a,b,z=[rc[8:20,:,j] for j in range(3)];s=np.sqrt(a*b);hh=s.sum(axis=1);ab=(a*b).sum(axis=1);ac=(a*z).sum(axis=1);bc=(b*z).sum(axis=1);sz=(s*z).sum(axis=1)
    w=(z*((1-a)*(1-b)+(hh[:,None]-s)**2-2*ab[:,None]+2*a*b)).sum(axis=1)
    def prob(q):
        v=np.zeros(28);v[9:21]=q[:12];a,b,c=v[8:20],v[7:19],v[9:21]
        pb=(w*a*b*c).sum()
        empty=(1-c)*((1-a)*(1-b)+a*b*(hh**2+2*ab))+a*c*ac*(1-b)+b*c*bc*(1-a)-2*a*b*c*hh*sz
        return pb,1-empty
    def constraints(v):
        q=v[:13];d=v[13:26];pb1,h1=prob(q+d);pb2,h2=prob(q-d);hit=(h1+h2)/2
        return np.r_[q+d-1e-6,q-d-1e-6,v[26:]-d,v[26:]+d,(pb1+pb2)/2-lower-.00005,hit-hlo-1e-6,hhi-hit-1e-6]
    def eq(v):return np.array([v[:13].sum()-1,v[13:26].sum()])
    bnd=[(l+1e-6,u-1e-6) for l,u in bounds]+[(-.4,.4)]*12+[(0,0)]+[(0,1)]*13;bnd[12]=(.045,float(bounds[-1,1]))
    best=None
    for seed in range(3):
        v=np.random.default_rng(seed).normal(0,.03,13);v[-1]=0;v[:-1]-=v[:-1].mean();x=np.r_[q0,v,np.abs(v)]
        res=minimize(lambda x:1.5*x[26:].sum(),x,bounds=bnd,constraints=[{'type':'ineq','fun':constraints},{'type':'eq','fun':eq}],method='SLSQP',options={'ftol':1e-11,'maxiter':1500})
        if min(constraints(res.x))>-1e-8 and max(abs(eq(res.x)))<1e-8 and (best is None or res.fun<best.fun):best=res
    assert best is not None
    result=[]
    for sign in [-1,1]:
        q=best.x[:13]+sign*best.x[13:26];p=np.zeros((28,12));p[9:21]=q[:12,None]*r
        result.append(make_model(p,'cluster',728))
    return result

def main():
    p=argparse.ArgumentParser();p.add_argument('--cache',type=Path,required=True);p.add_argument('--check-source',action='store_true');p.add_argument('--models',action='store_true');a=p.parse_args()
    out=HERE/'results';read=lambda n:json.loads((out/n).read_text());cases=source_cases(a.cache)
    if a.check_source:
        old=read('matched-regions.json')['cases']
        for raw,saved in zip(cases,old):
            assert len(cases)==len(old)
            for key in ['time_ms','width','many_shots','bunch_count','patterns']:assert raw[key]==saved[key]
            for r,s in zip(raw['singletons'],saved['singletons']):assert r['shots']==s['shots'] and r['counts']==s['counts']
        print('All 16 complete independent-calibration settings reproduce from pinned originals.');return
    def save(name,value):(out/name).write_text(json.dumps(value,sort_keys=True,indent=None if name.startswith('matched-') else 2,separators=(',',':') if name.startswith('matched-') else None)+'\n')
    h=read('parity-histogram.json')['patterns'];N=read('counts.json')['many_shots'];K=sum(comb(144,k) for k in range(4))
    reg={'schema':1,'alphabet_size':K,'bunch_lower':read('full-region.json')['bunch_lower'],
        'pattern_intervals':{str(k):[str(endpoint(N,k,F(1,1280*K),side)) for side in ['lower','upper']] for k in sorted({k for p,k in h}|{0})},
        'row_hit_intervals':[chernoff_endpoints(N,sum(k for p,k in h if any(o//12==y for o in p)),F(1,15360)) for y in range(12)]}
    save('robust-region.json',reg)
    for case in cases:
        W=case['width'];M=W*W;N=case['many_shots'];K=sum(comb(M,k) for k in range(4))
        case['model_file']='matched-'+str(round(case['time_ms'],6))+'.json'
        case['bunch_lower']=str(endpoint(N,case['bunch_count'],F(1,1280),'lower'))
        case['pattern_intervals']={str(k):[str(endpoint(N,k,F(1,2560*K),side)) for side in ['lower','upper']] for k in sorted({k for p,k in case['patterns']}|{0})}
        ps=[]
        for record in case['singletons']:
            n=record['shots'];record['intervals']={str(k):[str(endpoint(n,k,F(1,7680*(M+1)),side)) for side in ['lower','upper']] for k in sorted(set(record['counts']))}
            ct=np.array(record['counts'][:M]).reshape(W,W)
            rowcts=list(map(int,ct.sum(axis=1)))+[record['counts'][-1]]
            record['row_intervals']=[[str(endpoint(n,k,F(1,7680*(W+1)),side)) for side in ['lower','upper']] for k in rowcts]
            ps.append((ct+.05)/(ct.sum()+.05*M)*.95)
        if a.models:
            pr=np.stack(ps,axis=-1);par=read('matched-search-parameters.json')[str(round(case['time_ms'],6))];mix,power=par['mix'],par['row_power']
            pr=(1-mix)*pr+mix*pr.mean(axis=2)[:,:,None];pr*=pr.sum(axis=1)[:,None,:]**(power-1);pr*=par['survival']/pr.sum(axis=(0,1))
            save(case['model_file'],independent_model(pr,par['seed']))
    save('matched-regions.json',{'schema':1,'search_family_size':16,'cases':cases})
    if a.models:
        for name,model in zip(['drift-minus.json','drift-plus.json'],drift_models(read)):save(name,model)
    print('Proposals saved; run robustness.py for independent exact acceptance.')
if __name__=='__main__':main()
