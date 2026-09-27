"""Nonseparable diagonal-hopping follow-up, with date-separated numerical checks."""
import argparse,json,time
from itertools import combinations,permutations
from pathlib import Path
import extended_dynamics as e
import numpy as np
from scipy.linalg import eigh
from numpy.polynomial.hermite import hermgauss
from scipy.optimize import minimize
from scipy.special import xlogy

def channel_ensemble(p,times,order=7,half=10):
 jx,jy,j2x,j2y,fx,fy,kx,ky,off,sig,eta,jd=p
 c=np.arange(-half,half+1);n=len(c);A=np.eye(n,k=1)+np.eye(n,k=-1);B=np.eye(n,k=2)+np.eye(n,k=-2);I=np.eye(n)
 hx=np.diag(fx*c+kx*c*c)-jx*A-j2x*B;hy=np.diag(fy*c+ky*c*c)-jy*A-j2y*B
 H=np.kron(hy,I)+np.kron(I,hx)-jd*np.kron(A,A);en,V=eigh(H)
 z,w=hermgauss(order);w/=np.sqrt(np.pi);t=(np.array(times)[:,None]+off)*(1+sig*np.sqrt(2)*z[None,:])/1000
 inp=np.array([(half+j)*n+half for j in [-1,0,1]])
 U=np.einsum('oe,tze,je->tzjo',V,np.exp(-2j*np.pi*t[:,:,None]*en),V[inp],optimize=True)
 return np.sqrt(eta)*U,w,n


def prediction(p,records,order=7,half=10):
 times=sorted(set(r['time_ms'] for r in records));U,w,n=channel_ensemble(p,times,order,half)
 qs=np.einsum('tzjo,z->tjo',abs(U)**2,w).reshape(len(times),3,n,n)
 result=[]
 for r in records:
  y,x=r['y_start']+half,r['x_start']+half;W=r['width'];q=qs[times.index(r['time_ms']),r['input_y']+1,y:y+W,x:x+W].ravel();assert len(q)==W*W
  result.append(np.r_[q,1-q.sum()])
 return result


def july2_bunch(parameters):
 U,w,n=channel_ensemble(parameters,[2.45,4.65],31,16);perms=list(permutations(range(3)));out=[]
 for i,(W,ys,xs) in enumerate([(12,-6,-5),(18,-9,-8)]):
  pats=np.array([[(y+ys+16)*n+x+xs+16 for x in xx] for y in range(W) for xx in combinations(range(W),3)])
  amplitudes=np.stack([U[i,:,perm[0],pats[:,0]]*U[i,:,perm[1],pats[:,1]]*U[i,:,perm[2],pats[:,2]] for perm in perms],axis=-1)
  # Shape: patterns x quadrature points x assignment permutations.
  distinguishable=np.sum(abs(amplitudes)**2,axis=(0,2));boson=np.sum(abs(amplitudes.sum(axis=2))**2,axis=0)
  clusters=[]
  for singleton in range(3):
   groups=[[k for k,p in enumerate(perms) if p[position]==singleton] for position in range(3)]
   clusters.append(sum(np.sum(abs(amplitudes[:,:,g].sum(axis=2))**2,axis=0) for g in groups))
  out.append(dict(time_ms=[2.45,4.65][i],distinguishable=float(w@distinguishable),boson=float(w@boson),
       cluster_partition_upper_at_this_point=float(w@np.max(clusters+[distinguishable],axis=0))))
 return out


def fit(records,start):
 scale=np.r_[e.SCALE,5];bounds=list(e.BOUNDS)+[(-15,15)];bounds[8]=(max(-1.5,-min(r['time_ms'] for r in records)+.01),.2)
 def loss(v):return -sum(np.sum(xlogy(r['counts'],np.maximum(q,1e-200))) for r,q in zip(records,prediction(v*scale,records)))
 result=minimize(loss,np.r_[start,0.]/scale,bounds=[(a/s,b/s) for (a,b),s in zip(bounds,scale)],method='L-BFGS-B',options={'maxiter':100,'ftol':1e-10,'gtol':1e-4,'maxls':20})
 return dict(parameters=(result.x*scale).tolist(),nll=float(result.fun),solver_success=bool(result.success))


if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');p.add_argument('--fit',action='store_true');a=p.parse_args()
 params=json.loads((e.HERE/'results/diagonal-dynamics-parameters.json').read_text())
 groups={'July2':json.loads((e.HERE/'results/july2-dynamics-records.json').read_text()),'July7':e.july7_records()}
 if a.fit:
  initial=json.loads((e.HERE/'results/extended-dynamics-parameters.json').read_text())
  params={date:fit(rs,initial[date]['parameters']) for date,rs in groups.items()}
  (e.HERE/'results/diagonal-dynamics-parameters.json').write_text(json.dumps(params,sort_keys=True,indent=2)+'\n')
 result={date:e.evaluate(rs,params[date]['parameters'],predictor=prediction,orders=(21,31),halves=(14,16)) for date,rs in groups.items()}
 result['July2']['many_body_bunch_predictions']=july2_bunch(params['July2']['parameters'])
 out=e.HERE/'results/diagonal-dynamics.json'
 if a.check:e.close(result,json.loads(out.read_text()))
 else:out.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
 print(json.dumps(result,indent=2))
