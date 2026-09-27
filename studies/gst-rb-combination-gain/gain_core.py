"""One confidence region and source deletion rule for full-class gain witnesses."""
from pathlib import Path
import sys,json,math
from fractions import Fraction
from functools import lru_cache
import numpy as np

HERE=Path(__file__).resolve().parent
BASE=HERE.parent/'shared-quantum-dynamics'
sys.path.insert(0,str(BASE))
from sqd_repeat import exact_upper,GRID
from sqd_general import GeneralEvaluator,operations,certificate,bounds
from sqd_sources import sources
from format_artifacts import formatted

BINS=[(0,16),(17,64),(65,256),(257,1024),(1025,8198)]


def exp_lower(x,terms=50):
    total=term=Fraction(1)
    for j in range(1,terms+1):term=term*x/j;total+=term
    return total


def radius_numerator(n):
    a=math.ceil(math.sqrt(math.log(960)/(2*n))*GRID)
    while exp_lower(Fraction(2*n*a*a,GRID*GRID))<960:a+=1
    return a


def region(rows):
    if len(rows)!=6661:raise ValueError('unexpected observation count')
    k=np.array([r['k'] for r in rows]);n=np.array([r['n'] for r in rows]);family=np.array([r['family'] for r in rows])
    up=lru_cache(None)(lambda k,n:exact_upper(int(k),int(n),80*6661))
    lo_num=[GRID-up(nn-kk,nn) for kk,nn in zip(k,n)]
    hi_num=[up(kk,nn) for kk,nn in zip(k,n)]
    lo=np.array(lo_num)/GRID;hi=np.array(hi_num)/GRID
    groups=[];means={};W=[]
    for f in ['GST','RB']:
        for low,high in [(0,8198)]+BINS:
            mask=np.array([r['family']==f and low<=r['length']<=high for r in rows]);N=int(n[mask].sum());K=int(k[mask].sum());a=radius_numerator(N)
            center=Fraction(K,N);radius=Fraction(a,GRID);lower=max(Fraction(0),center-radius);upper=min(Fraction(1),center+radius)
            w=np.zeros(len(rows));w[mask]=n[mask]/N;W.append(w)
            g={'family':f,'min_length':low,'max_length':high,'k':K,'n':N,'radius_numerator':a,'radius_denominator':GRID,
                'lower_numerator':lower.numerator,'lower_denominator':lower.denominator,'upper_numerator':upper.numerator,'upper_denominator':upper.denominator}
            groups.append(g)
            if (low,high)==(0,8198):means[f]=(lower,upper,w)
    upper_q=(1+means['RB'][1]-means['GST'][0])/2
    return {'lo_num':lo_num,'hi_num':hi_num,'lo':lo,'hi':hi,'family':family,'groups':groups,'W':np.array(W),
        'group_lo':np.array([g['lower_numerator']/g['lower_denominator'] for g in groups]),
        'group_hi':np.array([g['upper_numerator']/g['upper_denominator'] for g in groups]),
        'q_weights':(means['RB'][2]-means['GST'][2])/2,'RB_weights':means['RB'][2],
        'joint_upper_q':float(upper_q),'joint_upper_q_exact':upper_q,'joint_upper_RB':float(means['RB'][1])}


def residual_constraints(p,R,source='joint'):
    cell=np.ones(len(p),bool) if source=='joint' else R['family']==source
    grp=np.ones(len(R['groups']),bool) if source=='joint' else np.array([g['family']==source for g in R['groups']])
    means=R['W']@p
    return np.r_[p[cell]-R['lo'][cell],R['hi'][cell]-p[cell],means[grp]-R['group_lo'][grp],R['group_hi'][grp]-means[grp]]


def evaluate(p,R):
    p=np.array(p)
    return {'q':float(.5+R['q_weights']@p),'mu_RB':float(R['RB_weights']@p),
        'joint_upper_q':R['joint_upper_q'],'q_separation':float(.5+R['q_weights']@p-R['joint_upper_q']),
        'sources':{s:{'minimum_constraint_slack':float(residual_constraints(p,R,s).min()),'violations':int(np.sum(residual_constraints(p,R,s)<-2e-8))} for s in ['GST','RB','joint']}}


def save(path,obj):path.write_text(formatted(obj)+'\n')
