"""Finite-binomial deviance Chernoff bounds with directed interval enclosures.

These are candidatewise tests, never an optimizer-based model-class exclusion.
"""
import argparse,json,math,sys,time
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
import numpy as np
from mpmath import iv
HERE=Path(__file__).resolve().parent
PREVIOUS=HERE.parent/'unmeasured-dynamics-gain'
sys.path.insert(0,str(PREVIOUS))
from dynamics_common import BASE,baseline_models,save,sources,region
from dynamics_verify import binding,ERROR,record,audit
from interval_physics import certify

def rational(q):
    q=F(q);return iv.mpf(q.numerator)/q.denominator

def enclosure(center,error=ERROR):
    q=F.from_float(float(center));a=q-error;b=q+error
    if not 0<a<b<1:raise ValueError('probability envelope reaches boundary; special treatment required')
    return iv.mpf([rational(a).a,rational(b).b])

@lru_cache(None)
def coefficients(n,tilt):
    a=2*rational(tilt);out=[]
    for k in range(n+1):
        c=iv.mpf(math.comb(n,k))
        if k:c*= (iv.mpf(k)/n)**(a*k)
        if n-k:c*= (iv.mpf(n-k)/n)**(a*(n-k))
        out.append(c)
    return out

def log_mgf(n,p,tilt):
    if n<1 or not F(0)<tilt<F(1,2):raise ValueError('invalid n or tilt')
    # Positive-coefficient Horner evaluation, exact finite binomial sum.
    s=1-2*rational(tilt);a=p**s;b=(1-p)**s;y=a/b;c=coefficients(n,tilt);v=c[-1]
    for z in reversed(c[:-1]):v=v*y+z
    return n*iv.log(b)+iv.log(v)

def deviance(k,n,p):
    if not 0<=k<=n:raise ValueError('invalid count')
    out=iv.mpf(0)
    if k:out+=k*(iv.log(iv.mpf(k)/n)-iv.log(p))
    if n-k:out+=(n-k)*(iv.log(iv.mpf(n-k)/n)-iv.log(1-p))
    return 2*out

def outward(v,grid=10**6):
    low=F.from_float(float(np.nextafter(float(v.a),-np.inf)));high=F.from_float(float(np.nextafter(float(v.b),np.inf)))
    return [record(F(math.floor(low*grid),grid)),record(F(math.ceil(high*grid),grid))]

def evaluate(rows,probabilities,tilt):
    previous=iv.dps
    try:
        iv.dps=50
        return evaluate_at_precision(rows,probabilities,tilt)
    finally:
        iv.dps=previous

def evaluate_at_precision(rows,probabilities,tilt):
    if len(rows)!=len(probabilities):raise ValueError('row/probability mismatch')
    D=iv.mpf(0);M=iv.mpf(0)
    for row,p in zip(rows,probabilities):
        P=enclosure(p);D+=deviance(row['k'],row['n'],P);M+=log_mgf(row['n'],P,tilt)
    B=M-rational(tilt)*D
    # Always positive outward rounding, including tiny underflowed binary64 values.
    if float(B.b)>=0:upper=F(1)
    else:
        z=outward(iv.exp(B),10**12)[1];upper=min(F(1),F(z['numerator'],z['denominator']))
    return {'rows':len(rows),'tilt':record(tilt),'deviance_enclosure':outward(D),'log_mgf_enclosure':outward(M),
            'log_chernoff_factor_enclosure':outward(B),'tail_upper':record(upper),'excluded_at_0.025':upper<F(1,40)}

def models():
    out={'PR21_joint':(baseline_models()['joint'],'joint',F(38,1024))}
    for family,t in [('GST',399),('RB',101)]:
        out['PR21_contrast_'+family]=(json.loads((PREVIOUS/f'results/contrast-{family}.json').read_text()),family,F(t,1024))
    path=HERE/'results/stationary-refit.json'
    out['stationary_refit']=(json.loads(path.read_text()),'joint',F(38,1024))
    return out

def run(cache=None):
    rows=json.loads((BASE/'results/observations.json').read_text())
    if cache:
        raw=sources(cache)
        if len(rows)!=len(raw) or any(any(a[k]!=b[k] for k in ['family','row','k','n','length','word_sha256']) for a,b in zip(rows,raw)):raise ValueError('source binding mismatch')
        rows=raw
    R=region(rows)
    out={'sampling':'Independent binomial rows under fixed candidate; no fitted degrees of freedom',
         'scope':'Candidatewise rejection, not exclusion of stationary CPTP class; original PR21 confidence region unchanged',
         'interval_digits':50,'probability_envelope':record(ERROR),'models':{}}
    for name,(m,source,t) in models().items():
        c=certify(m['x'],max(r['length'] for r in rows))
        if F(c['probability_error_numerator'],c['probability_error_denominator'])>ERROR:raise ValueError('kernel envelope exceeded')
        if cache:binding(rows,m)
        ids=[i for i,r in enumerate(rows) if source=='joint' or r['family']==source]
        start=time.time();result=evaluate([rows[i] for i in ids],[m['probabilities'][i] for i in ids],t)
        original=audit(rows,m,R,[])
        families=['GST','RB'] if source=='joint' else [source]
        result['original_region_membership']={f:{'minimum_slack':original['minimum_slack'][f],'violations':original['violations'][f]} for f in families}
        result['retained']=source;out['models'][name]=result
        print(name,result['tail_upper']['decimal'],result['excluded_at_0.025'],'seconds',round(time.time()-start),flush=True)
    return out

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--cache');p.add_argument('--check',action='store_true');a=p.parse_args();out=run(a.cache);path=HERE/'results/adequacy.json'
    if a.check:
        if out!=json.loads(path.read_text()):raise ValueError('adequacy certificate mismatch')
    else:save(path,out)
