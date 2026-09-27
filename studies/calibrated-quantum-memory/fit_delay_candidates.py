#!/usr/bin/env python3
"""Optional bounded alternating likelihood search; not a null-exclusion solver.

The saved seed was itself fitted to these data. Inference uses pointwise
confidence-region inversion, never guessed fitted degrees of freedom.
The rational JSON output must pass the independent probability and statistical verifier.
"""
import argparse,json
from pathlib import Path
import identifiability as v
import numpy as np, cvxpy as cp
from scipy.optimize import least_squares
from scipy.spatial.transform import Rotation
from scipy.special import xlogy
def rationalize(values,x,settings,delays,mapping):
 D=10**12
 re=np.rint(values.real*D).astype(np.int64);im=np.rint(values.imag*D).astype(np.int64)
 for y in range(18):
  for b in range(2):re[y,b]=(re[y,b]+re[y,b].T)//2;im[y,b]=(im[y,b]-im[y,b].T)//2
 for axis in 'xyz':
  ref=next(i for i,k in enumerate(settings) if k[0]==axis)
  er=v.partial_trace_output(re[ref,0]);ei=v.partial_trace_output(im[ref,0])
  for y,s in enumerate(settings):
   if s[0]!=axis:continue
   for b in range(2):
    tr=er if b==0 else D*np.eye(2,dtype=int)-er;ti=ei if b==0 else -ei
    dr=tr-v.partial_trace_output(re[y,b]);di=ti-v.partial_trace_output(im[y,b])
    for i in range(2):
     for j in range(2):re[y,b,2*i+1,2*j+1]+=int(dr[i,j]);im[y,b,2*i+1,2*j+1]+=int(di[i,j])
 witness=dict(mapping=mapping,denominator=D,memory_weight=[1,10000],settings=settings,choi=[[dict(real=re[y,b].ravel().tolist(),imag=im[y,b].ravel().tolist()) for b in range(2)] for y in range(18)])
 v.check_instruments(witness)
 def ball(vec):
  vec=np.asarray(vec);radius=np.linalg.norm(vec)
  if radius>1-1e-8:vec=vec*(1-1e-8)/radius
  return np.rint(vec*D).astype(np.int64).tolist()
 original=np.array([[np.trace(v.STATES[name]@s).real for s in v.PAULIS] for name in v.STATES]);prep=x[54:78].reshape(6,4)
 preps=np.einsum('nij,nj->ni',Rotation.from_rotvec(prep[:,:3]).as_matrix(),original)*prep[:,3,None]
 measurement=x[78:].reshape(3,5);directions=np.einsum('zij,zj->zi',Rotation.from_rotvec(measurement[:,:3]).as_matrix(),np.eye(3))
 effects=[]
 for axis in range(3):
  bias=measurement[axis,4]-measurement[axis,3];vec=(1-measurement[axis,3:].sum())*directions[axis]
  if abs(bias)+np.linalg.norm(vec)>1-1e-8:vec=vec*(1-abs(bias)-1e-8)/np.linalg.norm(vec)
  effects.append([int(round(bias*D)),*np.rint(vec*D).astype(np.int64).tolist()])
 cayley=[]
 for rotvec in x[:54].reshape(18,3):
  theta=np.linalg.norm(rotvec);t=np.zeros(3) if theta==0 else rotvec*np.tan(theta/2)/theta
  cayley.append(np.rint(t*10**8).astype(np.int64).tolist())
 return dict(instrument=witness,spam_denominator=D,preparations={k:ball(vec) for k,vec in zip(v.STATES,preps)},effects=dict(zip('xyz',effects)),rotation_denominator=10**8,delays=delays,rotation_cayley=np.array(cayley).reshape(9,2,3).tolist())

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--source-dir',required=True,type=Path)
parser.add_argument('--mapping',choices=['regroup_first_1','literal_rows'],required=True)
parser.add_argument('--output',required=True,type=Path)
parser.add_argument('--rounds',type=int,default=8)
args=parser.parse_args()
mapping=args.mapping
counts,labels,delays=v.load_counts(args.source_dir,mapping)
n=counts.sum(-1);freq=counts/n[...,None];settings=sorted({','.join(k.split(',')[1:3]) for k in labels});axes='xyz'
r=np.array([[np.trace(v.STATES[k.split(',')[0]]@s).real for s in v.PAULIS] for k in labels]);z=np.array([axes.index(k[-1]) for k in labels]);y=np.array([settings.index(','.join(k.split(',')[1:3])) for k in labels]);sig=np.array([v.I,*v.PAULIS]);w=.0001
catalog=json.loads((v.HERE/'results/ibm-delay-candidates.json').read_text())
x=np.array(catalog['discovery_initial_parameters'],dtype=float)

names=list(v.STATES);initial=np.array([names.index(k.split(',')[0]) for k in labels])
lo=np.r_[np.full(54,-.5),np.tile([-.3,-.3,-.3,.8],6),np.tile([-.3,-.3,-.3,0,0],3)]
hi=np.r_[np.full(54,.5),np.tile([.3,.3,.3,1.],6),np.tile([.3,.3,.3,.1,.1],3)]
def pt(j):return cp.bmat([[cp.trace(j[2*i:2*i+2,2*k:2*k+2]) for k in range(2)] for i in range(2)])
def probes(x):
 rot=Rotation.from_rotvec(x[:54].reshape(18,3)).as_matrix().reshape(9,2,3,3)
 prep=x[54:78].reshape(6,4)
 orig=np.array([[np.trace(v.STATES[name]@s).real for s in v.PAULIS] for name in names])
 states=np.einsum('nij,nj->ni',Rotation.from_rotvec(prep[:,:3]).as_matrix(),orig)*prep[:,3,None]
 rr=np.einsum('dij,rj->dri',rot[:,0],states[initial])
 rho=(v.I+np.einsum('dri,iab->drab',rr,v.PAULIS))/2
 measurement=x[78:].reshape(3,5)
 directions=np.einsum('zij,zj->zi',Rotation.from_rotvec(measurement[:,:3]).as_matrix(),np.eye(3))
 directions=np.einsum('zj,djk->dzk',directions,rot[:,1])
 obs=(measurement[:,4]-measurement[:,3])[None,:,None,None]*v.I+(1-measurement[:,3:].sum(1))[None,:,None,None]*np.einsum('dzk,kab->dzab',directions,v.PAULIS)
 ee=np.stack([(v.I+s*obs[:,z])/2 for s in [1,-1]],axis=2)
 return rho,ee
js=[[cp.Variable((4,4),hermitian=True) for b in range(2)] for _ in settings]
cons=[];params=[];loss=[]
for yi,pair in enumerate(js):
 cons += [pt(sum(pair))==v.I]
 ref=next(i for i,k in enumerate(settings) if k[0]==settings[yi][0])
 ix=np.where(y==yi)[0]
 for b,j in enumerate(pair):
  cons += [j-w*cp.real(cp.trace(j))/2*v.IDENTITY_CHOI >> 1e-8*np.eye(4),pt(j-js[ref][b])==0]
  param=cp.Parameter((9*len(ix)*2,16),complex=True);params.append((param,ix))
  pred=cp.real(param@cp.vec(j,order='C'))
  observed=counts[:,ix,2*b:2*b+2].reshape(-1)
  exposures=np.broadcast_to(n[:,ix,None],(9,len(ix),2)).reshape(-1)
  loss.append(cp.sum(cp.kl_div(observed,cp.multiply(exposures,pred))))
prob=cp.Problem(cp.Minimize(sum(loss)),cons)
def predict(x,J):
 rho,eff=probes(x)
 out=np.einsum('drij,r b i a j c->drbac',rho,J[y].reshape(324,2,2,2,2,2))
 return np.einsum('drbac,drkca->drbk',out,eff).real.reshape(9,324,4)
def resid(x,J):
 p=np.maximum(predict(x,J),1e-12)
 d=2*(xlogy(freq,freq/p)-freq+p)*n[...,None]
 return (np.sign(p-freq)*np.sqrt(np.maximum(d,0))).ravel()
for iteration in range(args.rounds):
 rho,eff=probes(x)
 for param,ix in params:
  mat=np.einsum('drji,drkab->drkiajb',rho[:,ix],eff[:,ix]).reshape(9*len(ix)*2,4,4).transpose(0,2,1).reshape(-1,16)
  param.value=mat
 prob.solve(ignore_dpp=True,solver='CLARABEL',tol_gap_abs=1e-6,tol_feas=1e-9,tol_gap_rel=1e-8,max_iter=150)
 J=np.array([[j.value for j in pair] for pair in js]);print('SDP',iteration,prob.status,2*prob.value,flush=True)
 fit=least_squares(lambda a:resid(a,J),x,bounds=(lo,hi),max_nfev=60,ftol=1e-8,xtol=1e-8,gtol=1e-6)
 x=fit.x;p=predict(x,J)
 print('rot',iteration,2*fit.cost,'pearson',np.sum((counts-n[...,None]*p)**2/(n[...,None]*p)),flush=True)

candidate=rationalize(J,x,settings,delays,mapping)
args.output.write_text(json.dumps(dict(candidates=[candidate]),indent=2)+'\n')
print('Candidate only: no optimizer status certifies fit, physicality or exclusion.')
