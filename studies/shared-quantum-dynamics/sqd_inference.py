"""Descriptive likelihood fits and deliberately conditional uncertainty calculations."""
import numpy as np
from scipy.optimize import least_squares
from scipy.special import xlogy
from scipy.stats import beta
from sqd_models import Evaluator, parameters


def deviance(k,n,p):
    k,n,p=map(np.asarray,(k,n,p))
    p=np.clip(p,1e-14,1-1e-14)
    f=k/n
    return 2*(xlogy(k,f/p)+xlogy(n-k,(1-f)/(1-p)))


def residual(k,n,p):
    return np.sign(k/n-p)*np.sqrt(np.maximum(0,deviance(k,n,p)))


def fit(rows,model,starts=2,max_nfev=100,seed=260927,counts=None,initial=None,fixed=None):
    evaluator=Evaluator(rows)
    k=np.array([r['k'] for r in rows]) if counts is None else np.asarray(counts)
    n=np.array([r['n'] for r in rows])
    x,lo,hi=parameters(model)
    if initial is not None:x=np.asarray(initial)
    rng=np.random.default_rng(seed)
    fixed={} if fixed is None else fixed
    free=np.array([i for i in range(len(x)) if i not in fixed],int)
    for i,v in fixed.items():x[i]=v
    def expand(v):
        full=x.copy();full[free]=v;return full
    best=None
    for start in range(starts):
        xx=x.copy()
        if start:
            xx[:9]+=rng.normal(0,.1,9)
            xx[9:12]*=rng.uniform(.5,1.5,3)
            if model!='qubit':xx[15]*=2
        xx=np.clip(xx,lo+1e-10,hi-1e-10)
        result=least_squares(lambda v:residual(k,n,evaluator(expand(v),model)),xx[free],bounds=(lo[free],hi[free]),
                             max_nfev=max_nfev,ftol=2e-7,xtol=2e-7,gtol=1e-5,x_scale='jac',diff_step=1e-4)
        score=float(2*result.cost)
        if best is None or score<best['deviance']:
            best={'x':expand(result.x).tolist(),'deviance':score,'success':bool(result.success),
                  'nfev':int(result.nfev),'optimality':float(result.optimality)}
    return best


def diagnostic(evaluator,x,model,n):
    """Finite-difference observable Jacobian. Zero singular values retain gauges."""
    x=np.asarray(x);p=evaluator(x,model)
    J=[]
    _,lo,hi=parameters(model)
    for i in range(len(x)):
        h=1e-4*max(1,abs(x[i]));a=x.copy();b=x.copy()
        a[i]=max(lo[i],x[i]-h);b[i]=min(hi[i],x[i]+h)
        J.append((evaluator(b,model)-evaluator(a,model))/(b[i]-a[i]))
    J=np.array(J).T
    weighted=J*np.sqrt(np.asarray(n)/np.maximum(p*(1-p),1e-8))[:,None]
    U,s,V=np.linalg.svd(weighted,full_matrices=False)
    return {'singular_values':s.tolist(),'rank_relative_1e-6':int(np.sum(s>s[0]*1e-6)),
            'parameter_count':len(x),'threshold':float(s[0]*1e-6),'note':'local diagnostic; parameter units matter; not a global identification proof'}


def simultaneous_intervals(k,n,alpha=.05):
    k,n=np.asarray(k),np.asarray(n)
    tail=alpha/(2*len(k))
    low=np.where(k==0,0,beta.ppf(tail,k,n-k+1))
    high=np.where(k==n,1,beta.ppf(1-tail,k+1,n-k))
    return low,high


def context_certificate(rows,alpha=.05):
    """All shared-word family contrasts; any dimension/shared-reset model obeys equality.

    Bonferroni CP bounds require iid shots for each word/family pool. Return a
    lower bound on half the cross-family probability discrepancy, NOT leakage.
    """
    pools={}
    for r in rows:
        key=(r['family'],r['word_sha256'])
        k,n=pools.get(key,(0,0));pools[key]=(k+r['k'],n+r['n'])
    common=sorted({w for f,w in pools if f=='GST'} & {w for f,w in pools if f=='RB'})
    counts=[pools[(f,w)] for w in common for f in ('GST','RB')]
    k,n=np.array(counts).T
    lo,hi=simultaneous_intervals(k,n,alpha)
    witnesses=[]
    for j,w in enumerate(common):
        a,b=2*j,2*j+1
        witnesses.append({'word_sha256':w,'GST_counts':counts[a],'RB_counts':counts[b],
                          'GST_interval':[float(lo[a]),float(hi[a])], 'RB_interval':[float(lo[b]),float(hi[b])],
                          'half_context_gap_lower':float(max(0,lo[a]-hi[b],lo[b]-hi[a])/2)})
    return {'alpha':alpha,'words':witnesses,'max_half_context_gap_lower':max(w['half_context_gap_lower'] for w in witnesses),
            'assumptions':'iid Bernoulli shots within each word/family; no cross-pool independence needed for union bound',
            'unrestricted_within_acquisition_dependence_lower':0.}
