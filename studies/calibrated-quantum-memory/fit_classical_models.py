#!/usr/bin/env python3
"""Optional unrestricted classical discovery; never impose a quantum counterpart.

Alternate convex CP/TP blocks. Local failure is not a class exclusion.
The saved model is a data-fitted starting point, not independent calibration.
"""
import argparse,json,copy
import numpy as np
import cvxpy as cp
import identifiability as v
import classical_models as cm
from pathlib import Path
D=10**12

def rounded(values):
 r=np.rint(values.real*D).astype(object);i=np.rint(values.imag*D).astype(object)
 r=np.vectorize(int)(r).astype(object);i=np.vectorize(int)(i).astype(object)
 for j in range(len(r)):r[j]=(r[j]+r[j].T)//2;i[j]=(i[j]-i[j].T)//2
 return r,i

def record(r,i,den):return dict(denominator=int(den),choi=[dict(real=[int(x) for x in a.ravel()],imag=[int(x) for x in b.ravel()]) for a,b in zip(r,i)])
def mix_and_check(r,i,instrument=True):
 branches=len(r) if instrument else 1
 for scale in [10**10,10**9,10**8,10**7,10**6,10**5,10**4]:
  rr=2*branches*(scale-1)*r+D*np.eye(4,dtype=object)[None]
  ii=2*branches*(scale-1)*i;den=D*scale*2*branches
  group=record(rr,ii,den)
  try:cm.check_group(group,instrument)
  except ValueError:continue
  group['depolarizing_repair_weight']=[1,scale]
  return group
 raise ValueError('cannot rationalize')
def channel_group(values,instrument=True):
 r,i=rounded(values)
 groups=[list(range(len(r)))] if instrument else [[j] for j in range(len(r))]
 for ix in groups:
  dr=D*np.eye(2,dtype=object)-sum(v.partial_trace_output(r[j]) for j in ix)
  di=-sum(v.partial_trace_output(i[j]) for j in ix)
  j=ix[0]
  for a in range(2):
   for b in range(2):r[j,2*a+1,2*b+1]+=dr[a,b];i[j,2*a+1,2*b+1]+=di[a,b]
 return mix_and_check(r,i,instrument)
def instruments(B,settings):
 r,i=rounded(B.reshape(-1,4,4));r=r.reshape(18,2,4,4);i=i.reshape(18,2,4,4)
 for axis in 'xyz':
  ref=next(j for j,s in enumerate(settings) if s[0]==axis);er=v.partial_trace_output(r[ref,0]);ei=v.partial_trace_output(i[ref,0])
  for j,s in enumerate(settings):
   if s[0]!=axis:continue
   for b in range(2):
    dr=(er if b==0 else D*np.eye(2,dtype=object)-er)-v.partial_trace_output(r[j,b]);di=(ei if b==0 else -ei)-v.partial_trace_output(i[j,b])
    for a in range(2):
     for c in range(2):r[j,b,2*a+1,2*c+1]+=dr[a,c];i[j,b,2*a+1,2*c+1]+=di[a,c]
 # Use one repair weight/denominator for every setting.
 for scale in [10**10,10**9,10**8,10**7,10**6,10**5,10**4]:
  groups=[record(4*(scale-1)*a+D*np.eye(4,dtype=object)[None],4*(scale-1)*b,D*scale*4) for a,b in zip(r,i)]
  try:
   for g in groups:cm.check_group(g)
  except ValueError:continue
  for g in groups:g['depolarizing_repair_weight']=[1,scale]
  return groups
 raise ValueError('instrument rationalization')
def spam(rho,E):
 states={};effects={}
 for name,R in zip(v.STATES,rho):
  vec=np.array([np.trace(R@s).real for s in v.PAULIS]);norm=np.linalg.norm(vec)
  if norm>1-1e-8:vec*= (1-1e-8)/norm
  states[name]=np.rint(vec*D).astype(np.int64).tolist()
 for axis,pair in zip('xyz',E):
  obs=2*pair[0]-v.I;bias=np.trace(obs).real/2;vec=np.array([np.trace(obs@s).real/2 for s in v.PAULIS]);norm=np.linalg.norm(vec)
  if abs(bias)+norm>1-1e-8:vec*=(1-abs(bias)-1e-8)/norm
  effects[axis]=[int(round(bias*D)),*np.rint(vec*D).astype(np.int64).tolist()]
 return states,effects

def apply(J,X):return np.einsum('ij,iajb->ab',X,J.reshape(2,2,2,2))
def adj(J,F):return np.einsum('ab,jbia->ij',F,J.reshape(2,2,2,2))
def pt(J):return cp.bmat([[cp.trace(J[2*i:2*i+2,2*j:2*j+2]) for j in range(2)] for i in range(2)])
def solve(problem):
    problem.solve(solver='CLARABEL',tol_gap_abs=1e-9,tol_feas=1e-9,tol_gap_rel=1e-8,max_iter=150)
    if problem.value is None:raise RuntimeError('No candidate; not a null exclusion')
def loss(obs,N,p):return cp.sum(cp.multiply(N/8000,cp.kl_div(obs/N,p)))
def floating(group):
    mats,den=cm.group_matrices(group)
    return np.array([(np.asarray(r,float)+1j*np.asarray(i,float))/den for r,i in mats])

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--source-dir',type=Path,required=True)
p.add_argument('--mapping',choices=['regroup_first_1','literal_rows'],required=True)
p.add_argument('--initial',type=Path,default=v.HERE/'results/ibm-classical-candidates.json')
p.add_argument('--output',type=Path,required=True)
p.add_argument('--rounds',type=int,default=1)
p.add_argument('--fixed-return',action='store_true')
p.add_argument('--fit-spam',action='store_true')
args=p.parse_args()
model=next(m for m in json.loads(args.initial.read_text())['models'] if m['mapping']==args.mapping)
cm.check_model(model)
counts,labels,delays=v.load_counts(args.source_dir,args.mapping);n=counts.sum(-1)
settings=model['settings'];yi=np.array([settings.index(','.join(k.split(',')[1:3])) for k in labels])
ai=np.array([list(v.STATES).index(k.split(',')[0]) for k in labels]);zi=np.array(['xyz'.index(k[-1]) for k in labels])
B=np.array([floating(g) for g in model['instruments']]);A=np.array([floating(g['before']) for g in model['processes']]);C=np.array([floating(g['after']) for g in model['processes']]);S=model['spam_denominator']
rho=np.array([(v.I+sum(k*s/S for k,s in zip(model['preparations'][a],v.PAULIS)))/2 for a in v.STATES])
E=np.array([[(v.I+sign*(value[0]/S*v.I+sum(k*s/S for k,s in zip(value[1:],v.PAULIS))))/2 for sign in (1,-1)] for value in [model['effects'][axis] for axis in 'xyz']])
for iteration in range(args.rounds):
    js=[[cp.Variable((4,4),hermitian=True) for b in range(2)] for y in settings];cons=[];losses=[]
    for y,pair in enumerate(js):
        cons += [j >> 0 for j in pair]+[pt(sum(pair))==v.I]
        ref=next(i for i,k in enumerate(settings) if k[0]==settings[y][0])
        cons += [pt(pair[b]-js[ref][b])==0 for b in range(2)]
        ix=np.where(yi==y)[0]
        for b,J in enumerate(pair):
            matrix=[]
            for d in range(9):
                for r in ix:
                    for c in range(2):matrix.append(sum(np.kron(apply(a,rho[ai[r]]).T,adj(ret,E[zi[r],c])) for a,ret in zip(A[d],C[d])).T.ravel())
            pred=cp.real(np.array(matrix)@cp.vec(J,order='C'))
            losses.append(loss(counts[:,ix,2*b:2*b+2].ravel(),np.repeat(n[:,ix],2).ravel(),pred))
    problem=cp.Problem(cp.Minimize(sum(losses)),cons);solve(problem);B=np.array([[j.value for j in pair] for pair in js])
    for stage in ('before','after'):
        if stage=='after' and args.fixed_return:continue
        for d in range(9):
            js=[cp.Variable((4,4),hermitian=True) for _ in A[d]];cons=[j >> 0 for j in js]
            cons += [pt(sum(js))==v.I] if stage=='before' else [pt(j)==v.I for j in js]
            pred=0
            for e,J in enumerate(js):
                matrix=[]
                for r in range(324):
                    for b in range(2):
                        for c in range(2):
                            if stage=='before':state=rho[ai[r]];effect=adj(B[yi[r],b],adj(C[d,e],E[zi[r],c]))
                            else:state=apply(B[yi[r],b],apply(A[d,e],rho[ai[r]]));effect=E[zi[r],c]
                            matrix.append(np.kron(state.T,effect).T.ravel())
                pred+=cp.real(np.array(matrix)@cp.vec(J,order='C'))
            problem=cp.Problem(cp.Minimize(loss(counts[d].ravel(),np.repeat(n[d],4),pred)),cons);solve(problem)
            if stage=='before':A[d]=np.array([j.value for j in js])
            else:C[d]=np.array([j.value for j in js])
    if args.fit_spam:
        for a in range(6):
            R=cp.Variable((2,2),hermitian=True);matrix=[];ix=np.where(ai==a)[0]
            for d in range(9):
                for r in ix:
                    for b in range(2):
                        for c in range(2):matrix.append(sum(adj(pre,adj(B[yi[r],b],adj(ret,E[zi[r],c]))) for pre,ret in zip(A[d],C[d])).T.ravel())
            pred=cp.real(np.array(matrix)@cp.vec(R,order='C'))
            problem=cp.Problem(cp.Minimize(loss(counts[:,ix].ravel(),np.repeat(n[:,ix],4).ravel(),pred)),[R >> 0,cp.trace(R)==1]);solve(problem);rho[a]=R.value
        for z in range(3):
            F=cp.Variable((2,2),hermitian=True);matrix=[];offset=[];ix=np.where(zi==z)[0]
            for d in range(9):
                for r in ix:
                    for b in range(2):
                        out=sum(apply(ret,apply(B[yi[r],b],apply(pre,rho[ai[r]]))) for pre,ret in zip(A[d],C[d]))
                        for c in range(2):matrix.append((1-2*c)*out.T.ravel());offset.append(c*np.trace(out).real)
            pred=cp.real(np.array(matrix)@cp.vec(F,order='C'))+offset
            problem=cp.Problem(cp.Minimize(loss(counts[:,ix].ravel(),np.repeat(n[:,ix],4).ravel(),pred)),[F >> 0,v.I-F >> 0]);solve(problem);E[z]=[F.value,v.I-F.value]
    print('completed classical-only round',iteration+1,flush=True)
states,effects=spam(rho,E)
out=copy.deepcopy(model)
out.pop('counterpart_certificate',None);out.pop('counterpart_post_fit_mixture',None)
out.update(instruments=instruments(B,settings),spam_denominator=D,preparations=states,effects=effects,
           processes=[dict(before=channel_group(a),after=channel_group(c,False)) for a,c in zip(A,C)],
           discovery=dict(rounds=args.rounds,initial_model_is_data_fitted=True,fixed_return=args.fixed_return,fit_spam=args.fit_spam,quantum_constraint=False))
cm.check_model(out)
args.output.write_text(json.dumps(dict(models=[out]),indent=2)+'\n')
print('Rational physical candidate written; run classical_audit.py for complete statistical acceptance.')
