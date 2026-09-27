#!/usr/bin/env python3
"""Find a certificate with SciPy; certificate.py verifies it independently."""
from fractions import Fraction as F
from pathlib import Path
import json
import math
import numpy as np
from scipy.optimize import linprog
from scipy.stats import beta
from certificate import constraints, dual_bound, verify

HERE=Path(__file__).resolve().parent
GRID=10**8

def endpoints(n,k,tail):
    l=F(0) if k==0 else F(math.floor(beta.ppf(float(tail),k,n-k+1)*GRID),GRID)
    u=F(1) if k==n else F(math.ceil(beta.ppf(1-float(tail),k+1,n-k)*GRID),GRID)
    return l,u

def main():
    counts=json.loads((HERE/'results/counts.json').read_text())
    cats=counts['reference_row_counts']+[counts['single_shots']-sum(counts['reference_row_counts'])]
    intervals=[endpoints(counts['single_shots'],k,F(1,4160)) for k in cats]
    lower=endpoints(counts['many_shots'],counts['bunch_shots'],F(1,160))[0]
    loss=intervals[-1]
    lo=[F(0),F(0)]+[x[0] for x in intervals[:-1]]+[F(0)]
    hi=[loss[1],loss[1]]+[x[1] for x in intervals[:-1]]+[loss[1]]
    target=F(1166,100000)
    def build(l,u,depth=0):
        assert depth < 40, 'No certificate within depth limit'
        rows,rhs,bounds=constraints(l,u,loss)
        A=np.zeros((len(rows),39))
        for i,row in enumerate(rows):
            for j,v in row.items(): A[i,j]=float(v)
        eq=np.zeros((1,39));eq[0,:15]=1
        c=np.r_[np.zeros(27),-np.ones(12)]
        res=linprog(c,A_ub=A,b_ub=list(map(float,rhs)),A_eq=eq,b_eq=[1],
                    bounds=[(float(a),float(b)) for a,b in bounds],method='highs')
        assert res.success, res.message
        lam=[F(round(max(0,-v)*10**10),10**10) for v in res.ineqlin.marginals]
        mu=F(round(-res.eqlin.marginals[0]*10**10),10**10)
        ub=dual_bound(l,u,loss,lam,mu)
        if ub <= target:
            return {'lambda':{str(i):str(v) for i,v in enumerate(lam) if v},'mu':str(mu)}
        q=res.x[:15];t=res.x[27:]
        gaps=[t[i]-q[i]*q[i+1]*q[i+2] for i in range(12)]
        term=int(np.argmax(gaps));j=max(range(term,term+3),key=lambda j:u[j]-l[j])
        assert u[j]>l[j]
        mid=(l[j]+u[j])/2
        lu=u.copy();lu[j]=mid;rl=l.copy();rl[j]=mid
        return {'split':j,'at':str(mid),'left':build(l,lu,depth+1),'right':build(rl,u,depth+1)}
    cert={'schema':1,'category_intervals':[[str(v) for v in p] for p in intervals],
          'bunch_lower':str(lower),'polynomial_upper':str(target),'tree':build(lo,hi)}
    result=verify(cert,counts)
    (HERE/'results/certificate.json').write_text(json.dumps(cert,indent=2,sort_keys=True)+'\n')
    (HERE/'results/bounds.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(result)

if __name__=='__main__':main()
