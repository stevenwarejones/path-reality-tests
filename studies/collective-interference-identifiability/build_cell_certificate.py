#!/usr/bin/env python3
"""Find rational 2D bounds; the standard-library checker supplies the proof."""
from fractions import Fraction as F
from pathlib import Path
import json,time
import numpy as np
from scipy.optimize import linprog
from build_certificate import endpoints
from cell_certificate import conditional_problem,row_problem,repaired_bound,verify
HERE=Path(__file__).resolve().parent

def solve(problem):
    rows,rhs,bounds,eq,c,constant=problem
    A=np.zeros((len(rows),len(c)));E=np.zeros((len(eq),len(c)))
    for i,row in enumerate(rows):
        for j,v in row.items():A[i,j]=float(v)
    for i,row in enumerate(eq):
        for j,v in row.items():E[i,j]=float(v)
    r=linprog(-np.array(list(map(float,c))),A_ub=A,b_ub=list(map(float,rhs)),
      A_eq=E,b_eq=np.ones(len(eq)),bounds=[(float(l),float(u)) for l,u in bounds],method='highs')
    assert r.success,r.message
    lam=[F(round(max(0,-v)*10**10),10**10) for v in r.ineqlin.marginals]
    mu=[F(round(-v*10**10),10**10) for v in r.eqlin.marginals]
    d={'lambda':{str(i):str(v) for i,v in enumerate(lam) if v},'mu':list(map(str,mu))}
    return d,repaired_bound(problem,d),r.x

def upper_grid(x):return F(-(-x.numerator*10**8//x.denominator),10**8)

def main():
    start=time.perf_counter()
    read=lambda n:json.loads((HERE/'results'/n).read_text())
    c=read('counts.json');w=read('physical-witness.json');C=np.zeros((12,12),int)
    for y,x,n in w['nonzero_cell_counts']:C[y-9,x-8]=n
    cats=c['reference_row_counts']+[c['single_shots']-int(C.sum())]
    intervals=[endpoints(c['single_shots'],k,F(1,8320)) for k in cats]
    cond=[[endpoints(int(C[y].sum()),int(C[y,x]),F(1,92160)) for x in range(12)] for y in range(12)]
    ds=[];weights=[]
    for y in range(8,20):
        ps=[cond[z-9] if 9<=z<=20 else [(F(0),F(1))]*12 for z in [y-1,y,y+1]]
        d,ub,_=solve(conditional_problem([[p[0] for p in r] for r in ps],[[p[1] for p in r] for r in ps]))
        weight=min(F(1),upper_grid(ub));d['upper']=str(weight);weights.append(weight);ds.append(d)
    loss=intervals[-1]
    lo=[F(0),F(0)]+[p[0] for p in intervals[:-1]]+[F(0)]
    hi=[loss[1],loss[1]]+[p[1] for p in intervals[:-1]]+[loss[1]]
    def tree(l,u,weights,target,depth=0):
        assert depth<20,'Requested certificate unresolved'
        d,ub,x=solve(row_problem(l,u,loss,weights))
        if ub<=target:return d
        gaps=[float(weights[i])*(x[27+i]-x[i]*x[i+1]*x[i+2]) for i in range(12)]
        i=int(np.argmax(gaps));assert gaps[i]>1e-10,'Target below feasible relaxation point'
        j=max(range(i,i+3),key=lambda j:u[j]-l[j]);mid=(l[j]+u[j])/2
        lu=u.copy();lu[j]=mid;rl=l.copy();rl[j]=mid
        return {'split':j,'at':str(mid),'left':tree(l,lu,weights,target,depth+1),'right':tree(rl,u,weights,target,depth+1)}
    row_upper=F(119,10000);cell_upper=F(112,10000)
    cert={'schema':1,'row_intervals':[[str(v) for v in p] for p in intervals],
      'conditional_intervals':[[[str(v) for v in p] for p in r] for r in cond],
      'conditional_duals':ds,'row_upper':str(row_upper),'cell_upper':str(cell_upper),
      'bunch_lower':str(endpoints(c['many_shots'],c['bunch_shots'],F(1,160))[0]),
      'row_tree':tree(lo,hi,[F(1)]*12,row_upper),'cell_tree':tree(lo,hi,weights,cell_upper)}
    result=verify(cert,c,w)
    (HERE/'results/cell-certificate.json').write_text(json.dumps(cert,indent=2,sort_keys=True)+'\n')
    (HERE/'results/cell-bounds.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps(result,indent=2));print('generation+verification seconds',round(time.perf_counter()-start,3))

if __name__=='__main__':main()
