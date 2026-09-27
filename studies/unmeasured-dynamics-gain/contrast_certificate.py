"""Full-class two-sided operational contrast certificate, not a fixed-input bound."""
import argparse,math
from fractions import Fraction as F
from mpmath import iv
from dynamics_common import *
from dynamics_verify import audit,binding,record
from resource_certificate import exact_gates,matmul,positive_definite,choi

TARGET=(2,1,2,2,1,1,2,1,2,2)
INPUTS={1573:('GST','GyGyGy(GxGxGyGxGyGy)^5',50,50),
        4699:('RB','GxGyGyGiGxGxGxGyGxGyGyGyGxGyGx(Gy)^9GxGyGyGyGxGyGyGxGxGyGxGyGy',0,285)}

def enclosing_grid(value,denominator=10**10):
    # Outward binary64 conversion followed by exact rational directed rounding.
    low=F.from_float(float(np.nextafter(float(value.a),-np.inf)))
    high=F.from_float(float(np.nextafter(float(value.b),np.inf)))
    return F(math.floor(low*denominator),denominator),F(math.ceil(high*denominator),denominator)

def effective_contrast(x,word=TARGET):
    # Independently certify native Choi positivity (TP also follows structurally).
    G,P=exact_gates(x)
    for g in G:positive_definite(choi(g,P))
    c,l,u,theta=[iv.mpf(float(z)) for z in x[48:]]
    if not 0<=x[48]<=1 or not 0<=x[49]<=1 or not 0<=x[50]<=1:raise ValueError('nonphysical SPAM')
    e=[(l+u)/2,(l-u)*iv.sin(theta)/2,iv.mpf(0),(l-u)*iv.cos(theta)/2]
    for g in reversed(word):e=matmul([e],G[g])[0]
    radius=iv.sqrt(sum(z.real**2 for z in e[1:]))
    if float((e[0].real-radius).a)<0 or float((1-e[0].real-radius).a)<0:raise ValueError('effective POVM not certified')
    lo,hi=enclosing_grid(2*radius)
    return {'contrast_lower':record(lo),'contrast_upper':record(hi),
            'binary_success_lower':record((1+lo)/2),'binary_success_upper':record((1+hi)/2),
            'effective_eigenvalue_intervals':[[float(a),float(b)] for a,b in [enclosing_grid(e[0].real-radius),enclosing_grid(e[0].real+radius)]],
            'method':'60-digit directed normalized-Kraus reconstruction; native Choi LDL; dual-effect propagation and outward spectral-radius enclosure'}

def validate_inputs(rows):
    for i,(family,expr,k,n) in INPUTS.items():
        w=parse_word(expr);r=rows[i]
        if (r['family'],r['word_sha256'],r['k'],r['n'])!=(family,target_hash(w),k,n):raise ValueError('resource input mismatch')
        if w[-len(TARGET):]!=TARGET:raise ValueError('shared suffix changed')
    if target_hash(TARGET) in {r['word_sha256'] for r in rows}:raise ValueError('target is recorded')

def run(cache=None):
    rows=json.loads((BASE/'results/observations.json').read_text())
    if cache:
        raw=sources(cache)
        if len(rows)!=len(raw) or any(any(a[k]!=b[k] for k in ['family','row','word_sha256','k','n','length']) for a,b in zip(rows,raw)):raise ValueError('source binding mismatch')
        rows=raw
    validate_inputs(rows);R=region(rows);lower=F(R['lo_num'][1573]-R['hi_num'][4699],GRID)
    out={'status':'Certified two-sided gain for an operational discrimination resource; fixed-input sequence gain remains open',
         'target_word':list(TARGET),'target_sha256':target_hash(TARGET),'target_absent_from_complete_source_words':True,
         'resource':'C=diameter(G_TARGET^*(E)); optimal equal-prior binary success=(1+C)/2',
         'joint_contrast_lower':record(lower),'joint_binary_success_lower':record((1+lower)/2),
         'inputs':{str(i):{'family':f,'expression':ex,'k':k,'n':n,'lower_numerator':R['lo_num'][i],'upper_numerator':R['hi_num'][i],'denominator':GRID} for i,(f,ex,k,n) in INPUTS.items()},'models':{}}
    models={s:json.loads((HERE/f'results/contrast-{s}.json').read_text()) for s in ['GST','RB']}
    models['joint']=baseline_models()['joint']
    for name,m in models.items():
        retained=['GST','RB'] if name=='joint' else [name]
        membership=audit(rows,m,R,retained)
        if cache:binding(rows,m)
        contrast=effective_contrast(m['x']);q=contrast['contrast_upper'];upper=F(q['numerator'],q['denominator'])
        gap=lower-upper
        if name!='joint' and gap<=F(1,20):raise ValueError('two-sided resource margin lost')
        detector=abs(F.from_float(m['x'][50])-F.from_float(m['x'][49]))
        out['models'][name]={'membership':membership,**contrast,'separation_below_joint_lower':record(gap),
                             'detector_only_contrast':record(detector),'detector_only_separation':record(lower-detector)}
    out['scope']='All stationary qubit CPTP models with shared SPAM in unchanged simultaneous confidence region. Resource assumes freely controllable external input states; no matching calibration was observed. Spectral-range bound itself extends to any finite dimension.'
    return out

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--cache');p.add_argument('--check',action='store_true');a=p.parse_args();out=run(a.cache);path=HERE/'results/contrast-certificate.json'
    if a.check:
        if out!=json.loads(path.read_text()):raise ValueError('contrast certificate mismatch')
    else:save(path,out)
    print(out['status'],out['joint_contrast_lower']['decimal'])
    for name,m in out['models'].items():print(name,m['contrast_upper']['decimal'],m['separation_below_joint_lower']['decimal'])
