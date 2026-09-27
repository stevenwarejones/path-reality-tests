"""Full-class analytic bound, exact retained constraints, independent propagation."""
import argparse, importlib.util, math
from fractions import Fraction as F
from dynamics_common import *
from interval_physics import certify
from resource_certificate import verify_resource
spec=importlib.util.spec_from_file_location('prior_gain_verify',GAIN/'verify.py')
prior=importlib.util.module_from_spec(spec);spec.loader.exec_module(prior)
ERROR=F(1,10**9)
TARGET=(0,0,2,2)

def record(q):return {'numerator':q.numerator,'denominator':q.denominator,'decimal':float(q)}
def ceilroot(n):
    r=math.isqrt(n);return r+(r*r<n)
def lower_bound(a,b,f):
    if min(a,b,f)<0 or max(a,b,f)>GRID or a+b>GRID:raise ValueError('bound outside monotonic domain')
    return F(max(0,f-ceilroot(a*(GRID-b))-ceilroot(b*(GRID-a))),GRID)

def audit(rows,m,R,retained,aggregate=False):
    p=[F.from_float(float(v)) for v in m['probabilities']]
    if len(p)!=len(rows) or any(v<0 or v>1 for v in p):raise ValueError('invalid probabilities')
    cert=certify(m['x'],max(r['length'] for r in rows))
    if F(cert['probability_error_numerator'],cert['probability_error_denominator'])>ERROR:raise ValueError('numerical envelope exceeded')
    slacks={f:[] for f in ['GST','RB']}
    if not aggregate:
        for i,(r,v) in enumerate(zip(rows,p)):
            slacks[r['family']]+=[v-F(R['lo_num'][i],GRID)-ERROR,F(R['hi_num'][i],GRID)-v-ERROR]
    for g in R['groups']:
        if aggregate and (g['min_length'],g['max_length'])!=(0,8198):continue
        ids=[i for i,r in enumerate(rows) if r['family']==g['family'] and g['min_length']<=r['length']<=g['max_length']]
        mean=sum(rows[i]['n']*p[i] for i in ids)/g['n']
        slacks[g['family']]+=[mean-F(g['lower_numerator'],g['lower_denominator'])-ERROR,F(g['upper_numerator'],g['upper_denominator'])-mean-ERROR]
    for f in retained:
        if min(slacks[f])<=0:raise ValueError('retained constraint violated: '+f)
    return {'retained':retained,'aggregate_only':aggregate,'minimum_slack':{f:record(min(v)) for f,v in slacks.items()},'checked_two_sided_constraints':{f:len(v)//2 for f,v in slacks.items()},'violations':{f:sum(s<0 for s in v) for f,v in slacks.items()}}

def target_probability(x):
    row={'expression':'GiGiGyGy','word':TARGET}
    p=float(prior.independent_probabilities([row],x)[0])
    # The saved decimal center is recomputed; outward quantization stabilizes artifacts.
    # Independent evaluation is cross-checked with the certified Bloch propagation.
    z=float(GeneralEvaluator([row])(x)[0])
    if abs(p-z)>float(ERROR)/10:raise ValueError('target density/Bloch disagreement')
    return (F(math.floor(z*10**8),10**8)-ERROR,F(math.ceil(z*10**8),10**8)+ERROR)

def binding(rows,m):
    p=GeneralEvaluator(rows)(m['x']);saved=np.array(m['probabilities'])
    cert=certify(m['x'],max(r['length'] for r in rows))
    err=F(cert['probability_error_numerator'],cert['probability_error_denominator'])
    diff=max(abs(F.from_float(float(a))-F.from_float(float(b))) for a,b in zip(p,saved))
    if err+diff>ERROR:raise ValueError('artifact/model binding failed')
    independent=prior.independent_probabilities(rows,m['x'])
    if np.max(abs(independent-p))>float(ERROR):raise ValueError('independent propagation failed')

def run(cache=None,download=False):
    rows=json.loads((BASE/'results/observations.json').read_text());rows=rows['rows'] if isinstance(rows,dict) else rows
    if cache:
        raw=sources(cache,download)
        if len(raw)!=len(rows) or any(any(a[k]!=b[k] for k in ['family','row','k','n','length','word_sha256']) for a,b in zip(rows,raw)):raise ValueError('source binding mismatch')
        rows=raw
    if target_hash(TARGET) in {r['word_sha256'] for r in rows}:raise ValueError('target is measured')
    for i,w,k,n in [(0,(),0,50),(10,(2,2),50,50),(4662,(0,0),0,285)]:
        if (rows[i]['word_sha256'],rows[i]['k'],rows[i]['n'])!=(target_hash(w),k,n):raise ValueError('bound input mismatch')
    R=region(rows);L=lower_bound(R['hi_num'][0],R['hi_num'][4662],R['lo_num'][10])
    models=baseline_models();models['GST_search_0']=json.loads((HERE/'results/search-GST-0.json').read_text())['best'];models['GST_search_260927']=json.loads((HERE/'results/search-GST-260927.json').read_text())['best'];models['aggregate']=json.loads((HERE/'results/aggregate-model.json').read_text())
    result={'status':'one-sided unmeasured prediction; two-sided goal unresolved','target_word':list(TARGET),'target_sha256':target_hash(TARGET),'joint_lower':record(L),'inputs':{str(i):{'k':rows[i]['k'],'n':rows[i]['n'],'lower_numerator':R['lo_num'][i],'upper_numerator':R['hi_num'][i],'denominator':GRID} for i in [0,10,4662]},'models':{}}
    for name,m in models.items():
        retained=['GST','RB'] if name in ['joint','aggregate'] else ['RB'] if name=='RB_only' else ['GST']
        a=audit(rows,m,R,retained,name=='aggregate');low,high=target_probability(m['x']);a.update({'target_lower':record(low),'target_upper':record(high),'separation_below_joint_lower':record(L-high)})
        if name in ['RB_only','aggregate'] and L-high<=F(1,100):raise ValueError('separation lost')
        if cache:binding(rows,m)
        result['models'][name]=a
    resource=json.loads((HERE/'results/resource-certificate.json').read_text())
    checked=verify_resource(resource['models']['1009/1000']['rayleigh_vector'])
    if checked!=resource:raise ValueError('resource certificate mismatch')
    if resource['word_sha256'] in {r['word_sha256'] for r in rows}:raise ValueError('resource word is recorded')
    if cache:
        # Numerical cross-check of the exact similarity identity on every source word.
        engine=GeneralEvaluator(rows);g,initial,effect,_=operations(models['joint']['x']);scale=np.array([1.,1.009,1.009,1.009,1.])
        transformed=np.ascontiguousarray(g*scale[None,None,:,None]/scale[None,None,None,:]);out=np.empty(engine.count)
        engine.fast.qubit_prefix(engine.offsets,engine.words,engine.reuse,engine.count,transformed,np.ascontiguousarray(initial*scale),np.ascontiguousarray(effect/scale),out)
        if np.max(abs(out[engine.inverse]-np.array(models['joint']['probabilities'])))>float(ERROR):raise ValueError('similarity/source binding failed')
        from source_audit import audit_extra
        audit_extra(cache,rows,download)
        from discover import screen
        old=json.loads((HERE/'results/screen.json').read_text());new=screen(rows,R,write=False)
        for k in ['complete_family_sha256','search_count','eligible_one_sided_targets']:
            if old[k]!=new[k]:raise ValueError('screen mismatch')
        if new['primary']['word_sha256']!=old['primary']['word_sha256']:raise ValueError('selection changed')
    return result

def main():
    a=argparse.ArgumentParser();a.add_argument('--cache');a.add_argument('--download',action='store_true');a.add_argument('--check',action='store_true');args=a.parse_args()
    if args.download and not args.cache:a.error('--download requires --cache')
    out=run(args.cache,args.download);path=HERE/'results/certificate.json'
    if args.check:
        if json.loads(path.read_text())!=out:raise ValueError('artifact mismatch')
    else:save(path,out)
    print(out['status'],out['joint_lower']['decimal'])
    for name,v in out['models'].items():print(name,v['target_lower']['decimal'],v['target_upper']['decimal'],v['minimum_slack'][v['retained'][0]]['decimal'])
if __name__=='__main__':main()
