"""Finite-family adapted count bands and conservative nested-law enclosures."""
from math import log, expm1
import numpy as np
from common import (LAMBDAS, ALPHA, CDF_LABELS, PARTITION_LABELS, MAPS,
                    LEVELS, CELLS, grouped)


def count_bounds(k, u, n, labels):
    k,u=np.asarray(k),np.asarray(u)
    if not isinstance(n,int) or n<=0 or labels<=0 or np.any(k<0) or np.any(u<k) or np.any(u>n):
        raise ValueError('invalid count completion')
    b=log(2*len(LAMBDAS)*labels/ALPHA)
    lo=np.zeros(np.broadcast_shapes(k.shape,u.shape)); hi=np.full(lo.shape,float(n))
    # Scalar libm avoids NumPy SIMD-dispatch-dependent last-bit differences.
    for v in LAMBDAS:
        lo=np.maximum(lo,(v*k-b)/expm1(float(v)))
        hi=np.minimum(hi,(-v*u-b)/expm1(-float(v)))
    return np.stack([lo,hi],axis=-1)


def probabilities(k, unknown, n, epsilon, labels):
    if not np.isfinite(epsilon) or not 0<=epsilon<.25:
        raise ValueError('invalid assignment sensitivity')
    a=count_bounds(k,np.asarray(k)+np.asarray(unknown),n,labels)
    a[...,0]/=n*(.25+epsilon)
    a[...,1]=np.minimum(1.,a[...,1]/(n*(.25-epsilon)))
    return a


def compatible(a):
    return bool(np.all(a[...,0]<=a[...,1]+1e-14))


def delta(a):
    # arm axis is third from last: (..., arm, feature, endpoint)
    return np.stack([a[...,1,:,0]-a[...,0,:,1],a[...,1,:,1]-a[...,0,:,0]],axis=-1)


def absolute(a):
    return np.stack([np.maximum(0,np.maximum(a[...,0],-a[...,1])),
                     np.maximum(np.abs(a[...,0]),np.abs(a[...,1]))],axis=-1)


def cdf(counts, unknown, n, epsilon=0.):
    counts=np.asarray(counts); unknown=np.asarray(unknown)
    if counts.shape[-2:]!=(2,CELLS) or unknown.shape!=counts.shape[:-1]:
        raise ValueError('invalid CDF shapes')
    a=probabilities(np.cumsum(counts,axis=-1),unknown[...,None],n,epsilon,CDF_LABELS)
    a[...,0]=np.maximum.accumulate(a[...,0],axis=-1)
    a[...,1]=np.minimum.accumulate(a[...,1][...,::-1],axis=-1)[...,::-1]
    if not compatible(a): raise ValueError('incompatible CDF confidence/model intersection')
    band=delta(a); mag=absolute(band)
    grid=np.stack([mag[...,0].max(axis=-1),mag[...,1].max(axis=-1)],axis=-1)
    previous=np.concatenate([np.zeros_like(a[...,:1,0]),a[...,:-1,0]],axis=-1)
    bridge=np.maximum(a[...,1,:,1]-previous[...,0,:],a[...,0,:,1]-previous[...,1,:])
    full=np.stack([grid[...,0],np.maximum(grid[...,1],bridge.max(axis=-1))],axis=-1)
    return dict(band=band,arm_bands=a,k_grid=grid,k_full=full,
                incidence=band[...,-1,:])


def parents(side, level):
    child=MAPS[side][level]; parent=MAPS[side][level-1]
    return np.array([parent[np.flatnonzero(child==i)[0]] for i in range(child.max()+1)])

PARENTS={s:[None]+[parents(s,j) for j in range(1,len(LEVELS))] for s in MAPS}


def tighten_tree(nodes, root, side):
    """Exact interval propagation for each arm's additive partition tree.

    It is NOT a sharp joint optimization over missing-record assignments or TV.
    """
    for level in range(len(nodes)-1,0,-1):
        child=nodes[level]; parent=nodes[level-1]; mapping=PARENTS[side][level]
        for i in range(parent.shape[-2]):
            summed=child[...,mapping==i,:].sum(axis=-2)
            parent[...,i,0]=np.maximum(parent[...,i,0],summed[...,0])
            parent[...,i,1]=np.minimum(parent[...,i,1],summed[...,1])
    total=nodes[0].sum(axis=-2)
    root[...,0]=np.maximum(root[...,0],total[...,0]); root[...,1]=np.minimum(root[...,1],total[...,1])
    for level,child in enumerate(nodes):
        mapping=np.zeros(child.shape[-2],int) if level==0 else PARENTS[side][level]
        parent=root[...,None,:] if level==0 else nodes[level-1]
        for i in range(parent.shape[-2]):
            ix=np.flatnonzero(mapping==i); values=child[...,ix,:].copy(); sums=values.sum(axis=-2)
            child[...,ix,0]=np.maximum(values[...,0],parent[...,i,0,None]-sums[...,1,None]+values[...,1])
            child[...,ix,1]=np.minimum(values[...,1],parent[...,i,1,None]-sums[...,0,None]+values[...,0])
    if not compatible(root) or any(not compatible(v) for v in nodes):
        raise ValueError('incompatible partition confidence/model intersection')
    return nodes,root


def partitions(counts, unknown, n, side, epsilon=0.):
    counts=np.asarray(counts); unknown=np.asarray(unknown)
    if counts.shape[-2:]!=(2,CELLS) or unknown.shape!=counts.shape[:-1]:
        raise ValueError('invalid partition shapes')
    root=probabilities(counts.sum(axis=-1),unknown,n,epsilon,PARTITION_LABELS)
    nodes=[probabilities(grouped(counts,m),unknown[...,None],n,epsilon,PARTITION_LABELS) for m in MAPS[side]]
    nodes,root=tighten_tree(nodes,root,side)
    event_delta=delta(root[...,None,:])[...,0,:]
    no_event=-event_delta[...,::-1]
    tvs=[]
    for node in nodes:
        d=np.concatenate([no_event[...,None,:],delta(node)],axis=-2)
        mag=absolute(d)
        lo=np.maximum(mag[...,0].max(axis=-1),.5*mag[...,0].sum(axis=-1))
        hi=np.minimum(root[...,1].max(axis=-1),.5*mag[...,1].sum(axis=-1))
        tvs.append(np.stack([lo,hi],axis=-1))
    tvs=np.stack(tvs,axis=-2)
    tvs[...,0]=np.maximum.accumulate(tvs[...,0],axis=-1)
    tvs[...,1]=np.minimum.accumulate(tvs[...,1][...,::-1],axis=-1)[...,::-1]
    if not compatible(tvs): raise ValueError('incompatible nested TV bounds')
    gain=np.stack([np.maximum(0,tvs[...,0]-tvs[...,:1,1]),
                   np.minimum(1,tvs[...,1]-tvs[...,:1,0])],axis=-1)
    gain[...,0,:]=0
    full=np.stack([tvs[...,-1,0],root[...,1].max(axis=-1)],axis=-1)
    return dict(tv=tvs,gain=gain,event_mass=root,full_tv=full,incidence=event_delta)
