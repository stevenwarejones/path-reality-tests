"""Exploratory training-only continuation; fixes three redundant factor scales."""
import argparse,json
import numpy as np
from scipy.optimize import least_squares
from sqd_sources import HERE,sources,split
from sqd_general import GeneralEvaluator,bounds,certificate
from sqd_inference import residual,deviance,simultaneous_intervals
from format_artifacts import formatted


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--cache',required=True);ap.add_argument('--check',action='store_true');a=ap.parse_args();rows=sources(a.cache);train=split(rows);fam=np.array([r['family'] for r in rows]);rr=[r for r,t in zip(rows,train) if t];E=GeneralEvaluator(rr);full=GeneralEvaluator(rows);k=np.array([r['k'] for r in rr]);n=np.array([r['n'] for r in rr]);path=HERE/'results/general-polish.json'
    if a.check:

        for saved in [path,path.with_name('general-polish-initial.json')]:
            r=json.loads(saved.read_text());np.testing.assert_allclose(full(r['x']),r['probabilities'],atol=2e-8,rtol=0)
        print('both training-only continuations reproduce');return
    base=np.array(json.loads((HERE/'results/general-fits.json').read_text())['joint']['x']);free=np.array([i for i in range(52) if i not in (0,16,32)]);lo,hi=bounds();low=lo[free]+1e-11;high=hi[free]-1e-11
    def expand(v):
        x=base.copy();x[free]=v;return x
    def fun(v):return residual(k,n,E(expand(v)))
    def jac(v):
        cols=[]
        for j in range(len(v)):
            u=v.copy();w=v.copy();u[j]=max(low[j],v[j]-1e-4);w[j]=min(high[j],v[j]+1e-4)
            cols.append((fun(w)-fun(u))/(w[j]-u[j]))
        return np.array(cols).T
    res=least_squares(fun,np.clip(base[free],low,high),bounds=(low,high),jac=jac,x_scale=1.,max_nfev=1200,ftol=1e-9,xtol=1e-9,gtol=1e-5)
    x=expand(res.x);p=full(x);out={'status':'exploratory training-only numerical continuation, not a new split or exclusion certificate','x':x.tolist(),'deviance':float(2*res.cost),'nfev':int(res.nfev),'success':bool(res.success),'optimality':float(res.optimality),'physical_certificate':certificate(x),'probabilities':p.tolist(),'scores':{}}
    for family in ['GST','RB']:
        for part,tt in [('train',train),('validation',~train)]:
            mask=tt&(fam==family);ks=np.array([r['k'] for r,m in zip(rows,mask) if m]);ns=np.array([r['n'] for r,m in zip(rows,mask) if m]);ll,hh=simultaneous_intervals(ks,ns)
            out['scores'][family+'_'+part]={'deviance':float(deviance(ks,ns,p[mask]).sum()),'cell_violations':int(np.sum((p[mask]<ll)|(p[mask]>hh)))}
    path.write_text(formatted(out)+'\n');print({k:v for k,v in out.items() if k not in ['x','probabilities']},flush=True)

if __name__=='__main__':main()
