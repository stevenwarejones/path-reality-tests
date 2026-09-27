"""Discovery only: screening and complete-source local witness searches."""
import argparse,json,time,hashlib
import numpy as np
from scipy.optimize import minimize
from dynamics_common import *

def screen(rows,R,write=True):
    targets=candidates(rows,R);E=GeneralEvaluator(targets);models=baseline_models();ps={name:E(m['x']) for name,m in models.items()}
    summaries=[]
    for kind in ['repeat','append_xx','prepend_xx','append_yy','prepend_yy']:
        ids=[j for j,t in enumerate(targets) if t['kind']==kind]
        summaries.append({'kind':kind,'count':len(ids),'prediction_ranges':{n:[float(p[ids].min()),float(p[ids].max())] for n,p in ps.items()},
            'nonvacuous_bounds':int(sum(targets[j].get('lower',0)>0 or targets[j].get('upper',1)<1 for j in ids))})
    eligible=[j for j,t in enumerate(targets) if t['kind'].startswith('append') and t['lower']>ps['RB_only'][j]+1e-6]
    # Primary: largest guaranteed gap against the complete-RB witness, then shortest
    # word/hash. This is outcome-dependent retrospective selection, reported in full.
    chosen=max(eligible,key=lambda j:(targets[j]['lower']-ps['RB_only'][j],-targets[j]['length'],targets[j]['word_sha256']))
    t=targets[chosen];t={**t,'word':list(t['word']),'baseline_predictions':{n:float(p[chosen]) for n,p in ps.items()}}
    artifact={'summaries':summaries,'primary':t,'search_count':len(targets),'complete_family_sha256':hashlib.sha256('\n'.join(t['kind']+':'+str(t['source_index'])+':'+t['word_sha256'] for t in targets).encode()).hexdigest(), 'eligible_one_sided_targets':len(eligible), 'family_definition':'Every RB row, five transformations in protocol order, recorded expanded words excluded; duplicates retained; all source rows and target hashes reconstructed by --check full-source verification'}
    if write:save(HERE/'results/screen.json',artifact);print('PRIMARY',t,flush=True)
    return artifact

def local_search(rows,R,target,source,seed,maxiter):
    baseline=baseline_models()['joint']['x'];x=np.array(baseline);rng=np.random.default_rng(seed)
    if seed:
        x[:48]+=rng.normal(0,.015,48)
    lo,hi=bounds();x=np.clip(x,lo+1e-10,hi-1e-10)
    # Remove only the three redundant process-factor overall scales.
    free=np.array([i for i in range(52) if i not in [0,16,32]])
    allrows=rows+[{'word':tuple(target['word'])}];E=GeneralEvaluator(allrows);cache={};best=None;history=[];start=time.time();iterations=0
    def expand(v):
        z=x.copy();z[free]=v;return z
    def probabilities(v):
        nonlocal cache,best
        if 'v' not in cache or not np.array_equal(cache['v'],v):
            z=expand(v);p=E(z);cache={'v':v.copy(),'p':p}
            slack=float(residual_constraints(p[:-1],R,source).min())
            if slack>2e-7 and (best is None or p[-1]<best['target_probability']):
                best={'x':z.tolist(),'probabilities':p[:-1].tolist(),'target_probability':float(p[-1]),'minimum_slack':slack,'retained':source}
        return cache['p']
    def callback(v):
        nonlocal iterations
        iterations+=1
        if iterations%10==0:
            entry={'iteration':iterations,'seconds':round(time.time()-start,3),'best_feasible':None if best is None else best['target_probability']};history.append(entry);print(entry,flush=True)
    p=probabilities(x[free])
    result=minimize(lambda v:probabilities(v)[-1],x[free],method='SLSQP',bounds=list(zip(lo[free]+1e-10,hi[free]-1e-10)),constraints=[{'type':'ineq','fun':lambda v:residual_constraints(probabilities(v)[:-1],R,source)-1e-6}],callback=callback,options={'maxiter':maxiter,'ftol':1e-9,'eps':1e-4})
    artifact={'source':source,'seed':seed,'maxiter':maxiter,'status':str(result.message),'success':bool(result.success),'iterations':int(result.nit),'elapsed_seconds':time.time()-start,'history':history,'best':best,'scope':'local search only; no exclusion or optimum certificate'}
    save(HERE/('results/search-'+source+'-'+str(seed)+'.json'),artifact);print(artifact['status'],None if best is None else (best['target_probability'],best['minimum_slack']),flush=True)

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--cache',required=True);a.add_argument('--source',choices=['GST','RB']);a.add_argument('--seed',type=int,default=0);a.add_argument('--maxiter',type=int,default=80);args=a.parse_args()
    rows=rows_from_cache(args.cache);R=region(rows)
    target=json.loads((HERE/'results/screen.json').read_text())['primary'] if args.source else screen(rows,R)['primary']
    if args.source:local_search(rows,R,target,args.source,args.seed,args.maxiter)
