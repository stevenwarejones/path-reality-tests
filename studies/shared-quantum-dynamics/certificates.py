"""Regenerate source geometry and all-word equivalence certificates."""
import argparse,json
import numpy as np
from sqd_sources import HERE,sources,split
from sqd_models import parameters,density_probability,Evaluator,certificate
from sqd_inference import diagnostic
from sqd_equivalence import synthetic_certificate,invariants,fiber_interval,fiber_point,leakage_population


def generate(rows,fits):
    train=split(rows);fams=np.array([r['family'] for r in rows]);x,_,_=parameters('leakage');x[:9]=np.arange(9)*.05;x[15:]=[1,2,.8]
    d={'status':'Source match and restricted-model obstruction; no structural leakage-vs-flag gain possible','generic_parameters':x.tolist(),'prototype':{},'diagnostics':{}}
    for f in ['GST','RB']:
        first=next(r for r in rows if r['family']==f and r['word']==(2,2,2,2))
        d['prototype'][f]={k:first[k] for k in ['expression','k','n','length','word_sha256']}
        d['prototype'][f]['independent_expansion']=[2,2,2,2]
        d['prototype'][f]['density_probability']=density_probability((2,2,2,2),x,'leakage')
    for f in ['GST','RB','joint']:
        mask=train if f=='joint' else train&(fams==f);rr=[r for r,m in zip(rows,mask) if m]
        d['diagnostics'][f]=diagnostic(Evaluator(rr),x,'leakage',[r['n'] for r in rr])
    E=Evaluator(rows);out={'synthetic':synthetic_certificate(),'measured_fit_fibers':{}}
    for source in ['GST','RB','joint']:
        x=fits['leakage_'+source]['x'];lo,hi=fiber_interval(x);a=fiber_point(x,lo+.05*(hi-lo));b=fiber_point(x,lo+.95*(hi-lo));p=E(a,'leakage');q=E(b,'leakage')
        out['measured_fit_fibers'][source]={'invariants':invariants(x),'loss_interval':[lo,hi],'parameters_a':a.tolist(),'parameters_b':b.tolist(),
            'max_probability_gap_all_measured_sequences':float(np.max(abs(p-q))), 'physical_a':certificate(a,'leakage'),'physical_b':certificate(b,'leakage'),
            'population_at_1000_a':leakage_population(1000,a),'population_at_1000_b':leakage_population(1000,b),
            'status':'exact fiber conditional on one feasible fitted observable vector; not a data confidence interval'}
    return {'gate-a.json':d,'equivalence.json':out}


def close(actual,expected,path='root'):
    if isinstance(expected,dict):
        if set(actual)!=set(expected):raise ValueError('keys differ '+path)
        for k in expected:close(actual[k],expected[k],path+'.'+k)
    elif isinstance(expected,list):
        if len(actual)!=len(expected):raise ValueError('length differs '+path)
        for i,(a,b) in enumerate(zip(actual,expected)):close(a,b,path+f'[{i}]')
    elif isinstance(expected,float):
        # Derivative SVD null values are platform dependent; absolute 2e-7 is
        # tiny relative to the explicitly recorded rank threshold (~.002).
        if not np.isclose(actual,expected,atol=2e-7,rtol=2e-5):raise ValueError('numeric mismatch '+path)
    elif actual!=expected:raise ValueError('value mismatch '+path)


if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--cache',required=True);a.add_argument('--check',action='store_true');args=a.parse_args()
    fits=json.loads((HERE/'results/fits.json').read_text())
    for name,obj in generate(sources(args.cache),fits).items():
        p=HERE/'results'/name
        if args.check:close(obj,json.loads(p.read_text()))
        else:p.write_text(json.dumps(obj,indent=2)+'\n')
    print('Source geometry and equivalence certificates reproduce')
