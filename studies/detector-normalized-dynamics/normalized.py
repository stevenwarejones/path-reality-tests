"""Remove detector affine scale and certify why the old LP witnesses cannot help."""
import json,sys,argparse
from pathlib import Path
from fractions import Fraction as F
from mpmath import iv
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'unmeasured-dynamics-gain'))
from dynamics_common import baseline_models,save
from contrast_certificate import TARGET,effective_contrast
from dynamics_verify import record

def ratio_interval(x,word):
    c=abs(F.from_float(x[50])-F.from_float(x[49]))
    if c==0:raise ValueError('normalized contrast undefined for constant detector')
    result=effective_contrast(x,word)
    return [F(result[k]['numerator'],result[k]['denominator'])/c for k in ['contrast_lower','contrast_upper']]

def affine_relation(base,changed):
    # Necessary for using the exact identity for every word, rather than a finite comparison.
    if any(base[i]!=changed[i] for i in list(range(49))+[51]):raise ValueError('gates, reset or detector axis changed')
    l,u=map(F.from_float,base[49:51]);ll,uu=map(F.from_float,changed[49:51])
    if l==u:raise ValueError('constant base detector')
    b=(uu-ll)/(u-l);a=ll-b*l
    return a,b

def run():
    joint=baseline_models()['joint']['x'];gst=json.loads((HERE.parent/'unmeasured-dynamics-gain/results/contrast-GST.json').read_text())['x']
    rb=json.loads((HERE.parent/'unmeasured-dynamics-gain/results/contrast-RB.json').read_text())['x'];a,b=affine_relation(joint,gst)
    values={name:[record(z) for z in ratio_interval(x,TARGET)] for name,x in [('joint',joint),('GST_contrast',gst),('RB_contrast',rb)]}
    # All suffix candidates from the earlier fully specified screen; bounds remain response-range bounds.
    screen=json.loads((HERE.parent/'unmeasured-dynamics-gain/results/contrast-screen.json').read_text())
    margins=[]
    for t in screen['candidates']:
        L=F(t['bound_numerator'],t['bound_denominator'])
        highs=[ratio_interval(x,tuple(t['word']))[1] for x in [gst,rb]]
        margins.append(L-max(highs))
    return {'scope':'Normalized dynamics target defined; no two-source dynamics gain established',
            'primary_word':list(TARGET),'normalized_contrasts':values,
            'exact_all_word_identity':{'a':record(a),'b':record(b),'equation':'E_new=a I+b E; TP gives G_w^*(E_new)=a I+b G_w^*(E), so every nonzero normalized contrast is unchanged'},
            'response_range_screen':{'positive_bounds':len(margins),'two_sided_crossings':sum(m>0 for m in margins),'largest_two_sided_margin':record(max(margins))},
            'response_only_tightness':{'E':'diag(0,1)','F':'diag(0.0452128544,0.7681878105)',
              'channel':'Measure input in Z; prepare diagonal output with plus probability 0.0452128544 on input 0 and 0.7681878105 on input 1',
              'normalized_contrast':record(F(7229749561,10**10)),
              'not_claimed':'This block-level response witness does not realize the entire recorded native gate set or satisfy every count constraint'},
            'constant_detector_case':'Undefined; the positive observed response gap excludes it from the joint region'}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();out=run();path=HERE/'results/normalized.json'
    if a.check:
        if out!=json.loads(path.read_text()):raise ValueError('normalized artifact mismatch')
    else:save(path,out)
    print(out['response_range_screen'])
