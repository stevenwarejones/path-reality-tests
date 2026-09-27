"""Reconstruct the discovered witnesses; no optimizer output is an outer bound."""
import argparse
import json
import numpy as np
from scipy.optimize import linprog
from gain_core import HERE,BASE,region,sources,operations,GeneralEvaluator,save
from sqd_general import pack,unitary
from sqd_models import PAULI,IDEAL
from sqd_transfer import reset_kraus,TAU


def pack_channels(channels,spam):
    x=[]
    for g,kk in enumerate(channels):
        U=unitary(IDEAL[g]);coeff=np.array([[np.trace(P@K@U.conj().T)/2 for P in PAULI] for K in kk]).T
        chi=coeff@coeff.conj().T;L=np.linalg.cholesky(chi+1e-16*np.eye(4));x.extend(pack(L))
    return np.r_[x,spam]


def optimize_readout(x,E,R,source,margin=1e-4):
    z=np.array(x);z[49]=0;z[50]=1;B=E(z);A=np.c_[1-B,B]
    cell=R['family']==source;groups=np.array([g['family']==source for g in R['groups']]);GA=R['W']@A
    C=np.r_[A[cell],-A[cell],GA[groups],-GA[groups]]
    b=np.r_[R['hi'][cell],-R['lo'][cell],R['group_hi'][groups],-R['group_lo'][groups]]-margin
    sol=linprog(-(R['q_weights']@A),A_ub=C,b_ub=b,bounds=[(0,1),(0,1)],method='highs')
    if not sol.success:raise RuntimeError('readout discovery failed; not an exclusion certificate')
    z[49:51]=sol.x
    return z


def main():
    a=argparse.ArgumentParser();a.add_argument('--cache',required=True);a.add_argument('--output',required=True);args=a.parse_args()
    rows=sources(args.cache);R=region(rows);E=GeneralEvaluator(rows)
    base=np.array(json.loads((BASE/'results/general-fits.json').read_text())['all_joint']['x'])
    _,_,_,channels=operations(base);delta=3e-5
    kk=[[np.sqrt(1-delta)*K for K in ks]+[np.sqrt(delta)*K for K in reset_kraus(TAU)] for ks in channels]
    gst=pack_channels(kk,base[48:]);rate=.00014
    kk=[[np.sqrt(1-rate)*np.eye(2)]+[np.sqrt(rate)*K for K in reset_kraus(TAU)] for _ in range(3)]
    rb=optimize_readout(pack_channels(kk,[.995,0,1,0]),E,R,'RB')
    from pathlib import Path
    save(Path(args.output),{name:{'retained':retained,'x':x.tolist(),'probabilities':E(x).tolist()}
        for name,retained,x in [('joint','joint',base),('GST_only','GST',gst),('RB_only','RB',rb)]})

if __name__=='__main__':main()
