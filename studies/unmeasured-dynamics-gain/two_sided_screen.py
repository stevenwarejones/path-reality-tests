"""Retrospective two-sided selection; physical point pools are discovery aids only."""
import argparse
from dynamics_common import *

def run(cache):
    rows=sources(cache);R=region(rows);targets=candidates(rows,R);E=GeneralEvaluator(targets);models=baseline_models()
    for seed in [0,260927]:models['GST_search_'+str(seed)]=json.loads((HERE/f'results/search-GST-{seed}.json').read_text())['best']
    ps={k:E(v['x']) for k,v in models.items()};gst=np.stack([v for k,v in ps.items() if k.startswith('GST')]);rb=ps['RB_only'];ranked=[]
    for j,t in enumerate(targets):
        for side,bound in [('lower',t.get('lower',0)),('upper',t.get('upper',1))]:
            if (side=='lower' and bound<=0) or (side=='upper' and bound>=1):continue
            if side=='lower':gp=float(gst[:,j].min());gm=bound-gp;rm=bound-rb[j]
            else:gp=float(gst[:,j].max());gm=gp-bound;rm=rb[j]-bound
            ranked.append({'word_sha256':t['word_sha256'],'kind':t['kind'],'source_index':t['source_index'],'length':t['length'],'side':side,'bound':float(bound),'GST_point':gp,'RB_point':float(rb[j]),'joint_point':float(ps['joint'][j]),'two_sided_margin':float(min(gm,rm))})
    ranked.sort(key=lambda r:(-r['two_sided_margin'],r['length'],r['word_sha256'],r['side']))
    return {'status':'Discovery screen; no physical region inferred from finite point spread','candidate_count':len(targets),'nonvacuous_directed_bounds':len(ranked),'positive_two_sided_margins':sum(r['two_sided_margin']>1e-6 for r in ranked),
        'objective':'max min(GST crossing margin, RB crossing margin), using the bound direction; fixed point pool is not an optimum certificate',
        'point_pool':list(models),'best':ranked[0],'top_ten':ranked[:10]}

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--cache',required=True);a.add_argument('--check',action='store_true');args=a.parse_args();out=run(args.cache);p=HERE/'results/two-sided-screen.json'
    if args.check:
        old=json.loads(p.read_text())
        for k in ['candidate_count','nonvacuous_directed_bounds','positive_two_sided_margins','point_pool']:
            if old[k]!=out[k]:raise ValueError('screen changed')
        if old['best']['word_sha256']!=out['best']['word_sha256'] or abs(old['best']['two_sided_margin']-out['best']['two_sided_margin'])>1e-8:raise ValueError('best target changed')
    else:save(p,out)
    print(out['positive_two_sided_margins'],out['best'])
