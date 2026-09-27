"""Exploratory physical witness bank, never an outer confidence region.

Simultaneous count intervals exclude individual observable vectors globally;
local fits supply attainable candidates only. Selection may be adaptive because
coverage is simultaneous over the prediction vector.
"""
import argparse,json
import numpy as np
from sqd_sources import HERE,sources,split
from sqd_models import Evaluator
from sqd_inference import fit,simultaneous_intervals


def acceptance(rows,probabilities,cell_budget=3785,group_budget=8):
    k=np.array([r['k'] for r in rows]);n=np.array([r['n'] for r in rows])
    # Cell constraints plus predeclared coarse length-bin averages.
    lo,hi=simultaneous_intervals(k,n,alpha=.025*len(rows)/cell_budget)
    cell=np.maximum(np.maximum(lo-probabilities,probabilities-hi),0)
    groups=[]
    for fam in ('GST','RB'):
        for low,high in [(0,16),(17,64),(65,256),(257,1024)]:
            mask=np.array([r['family']==fam and low<=r['length']<=high for r in rows])
            if mask.any():groups.append(mask)
    gaps=[]
    for mask in groups:
        N=n[mask].sum();f=k[mask].sum()/N;p=np.dot(n[mask],probabilities[mask])/N
        radius=np.sqrt(np.log(2*group_budget/.025)/(2*N))
        gaps.append(max(0,abs(f-p)-radius))
    return {'accepted':bool(max(float(cell.max()),max(gaps,default=0))<=1e-12),
            'cell_violations':int(np.sum(cell>1e-12)),'max_cell_gap':float(cell.max()),'max_group_gap':float(max(gaps,default=0))}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--cache',required=True);a=ap.parse_args();rows=sources(a.cache);train=split(rows);families=np.array([r['family'] for r in rows]);E=Evaluator(rows)
    fits=json.loads((HERE/'results/fits.json').read_text());bank={}
    for label,f in fits.items():bank[label]={'x':f['x'],'model':f['model'],'origin':f['source'],'probabilities':f['probabilities']}
    # Explicit exploratory extension after the frozen primary fits. No validation
    # outcome used in grid choice or nuisance fitting.
    for source in ('GST','RB','joint'):
        mask=train if source=='joint' else train&(families==source);rr=[r for r,m in zip(rows,mask) if m]
        for theta in (.5,1.,2.,3.,4.,5.):
            f=fit(rr,'alternating',starts=1,fixed={15:theta})
            p=E(f['x'],'alternating')
            bank[f'alternating_{source}_theta{theta}']={'x':f['x'],'model':'alternating','origin':source,'probabilities':p.tolist(),'local_success':f['success']}
        print('profile bank',source,flush=True)
    for label,b in bank.items():
        p=np.array(b['probabilities']);b['acceptance']={}
        for source in ('GST','RB','joint'):
            mask=train if source=='joint' else train&(families==source)
            b['acceptance'][source]=acceptance([r for r,m in zip(rows,mask) if m],p[mask])
    ranges={}
    for source in ('GST','RB','joint'):
        accepted=[(name,b) for name,b in bank.items() if b['acceptance'][source]['accepted']]
        if not accepted:ranges[source]={'accepted':0};continue
        ps=np.array([b['probabilities'] for _,b in accepted]);spread=np.ptp(ps,axis=0)
        ranges[source]={'accepted':len(accepted),'labels':[name for name,_ in accepted],
            'validation_max_attainable_spread':float(spread[~train].max()),
            'GST_validation_max_attainable_spread':float(spread[(~train)&(families=='GST')].max()),
            'RB_validation_max_attainable_spread':float(spread[(~train)&(families=='RB')].max()),
            'status':'finite-bank attainable spread; lower bound on diameter of interval-compatible physical predictions, NOT an outer bound'}
    result={'status':'exploratory witness bank; not a global model confidence set','bank':bank,'ranges':ranges,
            'coverage':'One simultaneous 95% joint prediction region under independent Bernoulli shots; fixed budget of 3785 CP cells alpha=.025 plus 8 Hoeffding groups alpha=.025. Source removal drops constraints from that same region.'}
    (HERE/'results/removal.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps(ranges,indent=2),flush=True)


if __name__=='__main__':main()
