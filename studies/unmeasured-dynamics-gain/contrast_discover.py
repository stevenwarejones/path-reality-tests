"""Retrospective detector LP and suffix selection; no optimization exclusion claim."""
import argparse
from scipy.optimize import linprog
from dynamics_common import *


def witness(rows,R,source):
    model=baseline_models()['joint' if source=='GST' else 'RB_only']
    x=np.array(model['x']);x[49:51]=[0,1];engine=GeneralEvaluator(rows);p=engine(x)
    # p(l,c)=l+c*p(0,1); l and c are lower eigenvalue and spectral diameter.
    A=np.c_[np.ones(len(p)),p];W=R['W']@A
    cells=R['family']==source;groups=np.array([g['family']==source for g in R['groups']])
    C=np.r_[A[cells],-A[cells],W[groups],-W[groups],[[1,1]]]
    b=np.r_[R['hi'][cells],-R['lo'][cells],R['group_hi'][groups],-R['group_lo'][groups],1]-1e-5
    result=linprog([0,1],A_ub=C,b_ub=b,bounds=[(0,1),(0,1)],method='highs')
    if not result.success:raise ValueError('No LP witness found; no physical exclusion follows')
    low,contrast=result.x;x[49:51]=[low,low+contrast];p=engine(x)
    if residual_constraints(p,R,source).min()<9e-6:raise ValueError('LP inward margin lost')
    return {'x':x.tolist(),'probabilities':p.tolist(),'retained':[source],
            'recipe':'Detector eigenvalue LP, fixed gates/reset/detector axis, all retained cells and groups, 1e-5 inward slack'}


def screen(rows,R,models):
    known={r['word'] for r in rows};gst={};rb={}
    for i,r in enumerate(rows):
        for length in range(1,min(256,r['length'])+1):
            w=r['word'][-length:]
            if w in known:continue
            dest=gst if r['family']=='GST' else rb
            value=R['lo_num'][i] if r['family']=='GST' else -R['hi_num'][i]
            if w not in dest or value>dest[w][0]:dest[w]=(value,i)
    ops={k:operations(v['x']) for k,v in models.items()};ranked=[];common=0
    for w,(low,i) in gst.items():
        if w not in rb:continue
        common+=1;bound=(low+rb[w][0])/GRID
        if bound<=0:continue
        contrasts={}
        for name,(g,initial,effect,channels) in ops.items():
            e=effect.copy()
            for gate in reversed(w):e=e@g[0,gate]
            contrasts[name]=float(2*np.linalg.norm(e[1:4]))
        ranked.append({'word':list(w),'word_sha256':target_hash(w),'length':len(w),'GST_row':i,'RB_row':rb[w][1],
                       'bound_numerator':low+rb[w][0],'bound_denominator':GRID,'contrasts':contrasts,
                       'two_sided_margin':min(bound-contrasts['GST'],bound-contrasts['RB'])})
    ranked.sort(key=lambda r:(-r['two_sided_margin'],-r['length'],r['word']))
    return {'status':'Retrospective discovery, not a physical region approximation','common_unmeasured_suffixes':common,
            'positive_bounds':len(ranked),'positive_two_sided_margins':sum(r['two_sided_margin']>1e-6 for r in ranked),
            'primary':ranked[0],'candidates':ranked}


def run(cache,check=False):
    rows=sources(cache);R=region(rows);models={s:witness(rows,R,s) for s in ['GST','RB']}
    selected=screen(rows,R,{**models,'joint':baseline_models()['joint']})
    if check:
        old=json.loads((HERE/'results/contrast-screen.json').read_text())
        if old['primary']['word_sha256']!=selected['primary']['word_sha256'] or old['positive_bounds']!=selected['positive_bounds']:raise ValueError('selection changed')
        for s in models:
            saved=json.loads((HERE/f'results/contrast-{s}.json').read_text())
            if np.max(abs(np.array(saved['x'])-np.array(models[s]['x'])))>1e-7:raise ValueError('witness recipe changed')
    else:
        for s,m in models.items():save(HERE/f'results/contrast-{s}.json',m)
        save(HERE/'results/contrast-screen.json',selected)
    print(selected['primary'])

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--cache',required=True);p.add_argument('--check',action='store_true');a=p.parse_args();run(a.cache,a.check)
