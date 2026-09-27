"""Reconstruct aggregate-only physical countermodel; no full-source feasibility claim."""
import argparse
from dynamics_common import *
from search import pack_channels
from sqd_general import unitary

def construct(rows,delta=.005189999999999999):
    R=region(rows);nx=np.array([r['word'].count(1) for r in rows]);ny=np.array([r['word'].count(2) for r in rows])
    B=np.sin(((nx+ny)*np.pi+(nx-ny)*delta)/2)**2
    means=np.array([R['W'][i]@B for i in [0,6]])
    centers=np.array([R['groups'][i]['k']/R['groups'][i]['n'] for i in [0,6]])
    l,u=np.linalg.solve(np.c_[1-means,means],centers)
    if not 0<=l<=u<=1:raise ValueError('aggregate effect not physical')
    channels=[[np.eye(2)],[unitary(np.array([np.pi+delta,0,0]))],[unitary(np.array([np.pi-delta,0,0]))]]
    x=pack_channels(channels,[1,l,u,0]);return x
if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--cache',required=True);args=a.parse_args();rows=sources(args.cache);x=construct(rows)
    saved=json.loads((HERE/'results/aggregate-model.json').read_text());p=GeneralEvaluator(rows)(x)
    if np.max(abs(p-np.array(saved['probabilities'])))>1e-9:raise ValueError('aggregate reconstruction changed')
    print('Aggregate-only constructive recipe reproduced; saved dyadic factors remain certificate definition.')
