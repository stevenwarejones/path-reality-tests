"""Local full-CPTP likelihood refinement; no model-class exclusion from status."""
import argparse,json,sys,time
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'unmeasured-dynamics-gain'))
from dynamics_common import sources,GeneralEvaluator,bounds,save,baseline_models
from sqd_general import certificate
from sqd_inference import residual

def run(cache,max_nfev=150):
    rows=sources(cache);engine=GeneralEvaluator(rows);k=np.array([r['k'] for r in rows]);n=np.array([r['n'] for r in rows])
    start=baseline_models()['joint']['x'];x=np.array(start);lo,hi=bounds();free=np.array([i for i in range(52) if i not in [0,16,32]])
    def full(z):
        out=x.copy();out[free]=z;return out
    t=time.time();fit=least_squares(lambda z:residual(k,n,engine(full(z))),x[free],bounds=(lo[free],hi[free]),jac='3-point',diff_step=1e-5,x_scale='jac',max_nfev=max_nfev,ftol=2e-8,xtol=2e-8,gtol=2e-5)
    params=full(fit.x);out={'x':params.tolist(),'probabilities':engine(params).tolist(),'deviance':float(2*fit.cost),'success':bool(fit.success),'nfev':int(fit.nfev),'optimality':float(fit.optimality),'elapsed_seconds':time.time()-t,'physical_certificate':certificate(params),'scope':'Local likelihood refinement, no global optimum or model-class exclusion','recipe':{'seed':'PR21 joint','jacobian':'three point','relative_step':1e-5,'fixed_redundant_scales':[0,16,32],'max_nfev':max_nfev}}
    save(HERE/'results/stationary-refit.json',out);print({k:v for k,v in out.items() if k not in ['x','probabilities']},flush=True)
if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--cache',required=True);a.add_argument('--max-nfev',type=int,default=150);v=a.parse_args();run(v.cache,v.max_nfev)
