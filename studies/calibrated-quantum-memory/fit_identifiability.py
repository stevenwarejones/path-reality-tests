#!/usr/bin/env python3
"""Optional convex discovery of IBM witnesses. Verification is independent.

Writes a candidate instrument catalog; a solver result alone is not evidence.
Run identifiability.py on the output before interpreting it.
"""
import argparse
import json
from pathlib import Path
import cvxpy as cp
import numpy as np
from scipy.stats import beta
import identifiability as verify


def discover(source_dir, mapping):
    weight=.003
    counts,labels,_=verify.load_counts(source_dir,mapping)
    n=counts.sum(axis=-1)[...,None];tail=1/(40*counts.size)
    lo=np.where(counts==0,0,beta.ppf(tail,np.maximum(counts,1),n-counts+1))
    hi=np.where(counts==n,1,beta.ppf(1-tail,counts+1,np.maximum(n-counts,1)))
    settings=sorted({','.join(k.split(',')[1:3]) for k in labels})
    js=[[cp.Variable((4,4),hermitian=True) for _ in range(2)] for _ in settings]
    slack=cp.Variable();constraints=[]
    def pt(j):
        return cp.bmat([[cp.trace(j[2*i:2*i+2,2*k:2*k+2]) for k in range(2)] for i in range(2)])
    for yi,pair in enumerate(js):
        constraints += [pt(sum(pair))==np.eye(2)]
        ref=next(i for i,k in enumerate(settings) if k[0]==settings[yi][0])
        for b,j in enumerate(pair):
            q=cp.real(cp.trace(j))/2
            constraints += [j-weight*q*verify.IDENTITY_CHOI >> 1e-7*np.eye(4),pt(j-js[ref][b])==0]
    for r,label in enumerate(labels):
        a,m,p,z=label.split(',');yi=settings.index(m+','+p)
        for b in range(2):
            for c in range(2):
                k=np.kron(verify.STATES[a].T,(verify.I+(1-2*c)*verify.PAULIS['xyz'.index(z)])/2)
                probability=cp.real(cp.trace(js[yi][b]@k))
                constraints += [probability>=float(lo[:,r,2*b+c].max())+slack,
                                probability<=float(hi[:,r,2*b+c].min())-slack]
    problem=cp.Problem(cp.Maximize(slack),constraints)
    problem.solve(solver='CLARABEL',tol_gap_abs=1e-10,tol_feas=1e-10,tol_gap_rel=1e-10,max_iter=200)
    if slack.value is None or slack.value<=0:
        raise RuntimeError('No candidate found; this is not a classical-null exclusion')
    denominator=10**9
    values=np.array([[j.value for j in pair] for pair in js])
    re=np.rint(values.real*denominator).astype(np.int64)
    im=np.rint(values.imag*denominator).astype(np.int64)
    for y in range(18):
        for b in range(2):
            re[y,b]=(re[y,b]+re[y,b].T)//2
            im[y,b]=(im[y,b]-im[y,b].T)//2
    # Repair exact common first effects and TP after rounding, before verification.
    for axis in 'xyz':
        ref=next(i for i,k in enumerate(settings) if k[0]==axis)
        er=verify.partial_trace_output(re[ref,0]);ei=verify.partial_trace_output(im[ref,0])
        for y,setting in enumerate(settings):
            if setting[0]!=axis:continue
            for b in range(2):
                tr=er if b==0 else denominator*np.eye(2,dtype=int)-er
                ti=ei if b==0 else -ei
                dr=tr-verify.partial_trace_output(re[y,b]);di=ti-verify.partial_trace_output(im[y,b])
                for i in range(2):
                    for j in range(2):
                        re[y,b,2*i+1,2*j+1]+=int(dr[i,j]);im[y,b,2*i+1,2*j+1]+=int(di[i,j])
    witness=dict(mapping=mapping,denominator=denominator,memory_weight=[3,1000],settings=settings,
                 choi=[[dict(real=re[y,b].reshape(-1).tolist(),imag=im[y,b].reshape(-1).tolist()) for b in range(2)] for y in range(18)])
    verify.check_instruments(witness)
    print(mapping,'candidate slack',slack.value,'— full count verification still required')
    return witness


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source-dir',required=True,type=Path)
    p.add_argument('--output',required=True,type=Path)
    args=p.parse_args()
    result=dict(kind='rational_physical_witnesses_not_raw_data',
                discovery='Optional SDP; optimizer status is not used by verifier',
                witnesses=[discover(args.source_dir,m) for m in ['regroup_first_1','literal_rows']])
    args.output.write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
