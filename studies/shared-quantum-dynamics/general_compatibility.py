"""Constructive in-sample compatibility, not held-out prediction or global MLE."""
import argparse,json
from functools import lru_cache
import numpy as np
from sqd_sources import HERE
from sqd_general import certificate
from sqd_repeat import exact_upper,GRID
from format_artifacts import formatted

BINS=[(0,16),(17,64),(65,256),(257,1024),(1025,8198)]


def calculate():
    rows=json.loads((HERE/'results/observations.json').read_text());fit=json.loads((HERE/'results/general-fits.json').read_text())['all_joint'];p=np.array(fit['probabilities']);k=np.array([r['k'] for r in rows]);n=np.array([r['n'] for r in rows]);N=len(rows)
    # Total alpha .05: .025 over two-sided cells and .025 over ten means.
    # Tail denominator is exactly 2*N/.025 = 80*N.
    upper=lru_cache(None)(lambda k,n:exact_upper(int(k),int(n),80*N))
    lo=np.array([GRID-upper(nn-kk,nn) for kk,nn in zip(k,n)])/GRID
    hi=np.array([upper(kk,nn) for kk,nn in zip(k,n)])/GRID
    cell_slack=float(min(np.min(p-lo),np.min(hi-p)));groups=[]
    for family in ['GST','RB']:
        for low,high in BINS:
            mask=np.array([r['family']==family and low<=r['length']<=high for r in rows]);exposure=int(n[mask].sum());observed=float(k[mask].sum()/exposure);prediction=float(n[mask]@p[mask]/exposure)
            radius=float(np.sqrt(np.log(2*10/.025)/(2*exposure)))
            groups.append({'family':family,'length_min':low,'length_max':high,'shots':exposure,'observed_mean':observed,'predicted_mean':prediction,'hoeffding_radius':radius,'slack':radius-abs(observed-prediction)})
    phys=certificate(fit['x']);margin=min(cell_slack,min(g['slack'] for g in groups))
    return {'status':'retrospective in-sample constructive ordinary compatibility in a declared conservative simultaneous count region; no held-out claim',
        'fit':'all_joint','confidence':.95,'cell_alpha':.025,'group_alpha':.025,'cell_count':N,'group_count':10,
        'exact_cell_tail_denominator':80*N,'cell_violations':int(np.sum((p<lo)|(p>hi))),
        'minimum_cell_slack':cell_slack,'groups':groups,'minimum_all_constraint_slack':margin,
        'numerical_verification_tolerance':2e-8,'physical_certificate':phys,
        'compatible':bool(margin>2e-8 and phys['trace_preservation_max_error']<1e-10 and phys['choi_min_eigenvalue']>=-1e-10),
        'resource_point':{'leakage_population':0,'classical_memory_states':1,'gate_time_dependence':False},
        'scope':'No strictly positive resource lower bound follows from this region, since it contains an explicit stationary qubit model. This does not certify the likelihood optimum, adequacy under stronger tests, or prediction of untouched data.',
        'sampling':'Independent Bernoulli shots for Hoeffding means; binomial within each row for cells. Unrestricted dependence voids these confidence claims.'}


if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--check',action='store_true');args=a.parse_args();r=calculate();path=HERE/'results/general-compatibility.json'
    if args.check:
        from certificates import close
        close(r,json.loads(path.read_text()))
    else:path.write_text(formatted(r)+'\n')
    print({k:v for k,v in r.items() if k!='groups'})
