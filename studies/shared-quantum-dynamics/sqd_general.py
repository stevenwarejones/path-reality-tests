"""General stationary qubit CPTP gates via normalized positive process matrices.

Unlike the earlier unital depolarizing model, this representation covers every
qubit CPTP map. Cholesky rank deficiency is allowed; a positive input marginal is
required for normalization. Every already-TP map has that marginal equal to I.
"""
import numpy as np
import ctypes,hashlib,subprocess,tempfile
from pathlib import Path
from sqd_sources import HERE
from scipy.optimize import least_squares
from sqd_models import PAULI,IDEAL,Evaluator,kraus as restricted_kraus
from sqd_inference import residual

# Column-by-column lower-triangular Cholesky: positive real diagonal plus complex
# off-diagonal. Scaling changes optimizer coordinates, not the represented set.
SLOTS=[(i,j,part) for j in range(4) for i in range(j,4) for part in ([0] if i==j else [0,1])]
SCALE=np.array([1. if i==j==0 else (.01 if i==j else .001) for i,j,part in SLOTS])


def unpack(v):
    L=np.zeros((4,4),complex)
    for a,(i,j,part),scale in zip(v,SLOTS,SCALE):L[i,j]+=a*scale*(1j if part else 1)
    return L


def pack(L):return np.array([(L[i,j].imag if part else L[i,j].real)/scale for (i,j,part),scale in zip(SLOTS,SCALE)])


def bounds():
    lo=np.array([0. if i==j else -1./scale for (i,j,part),scale in zip(SLOTS,SCALE)])
    hi=1./SCALE
    return np.r_[np.tile(lo,3),0.,0.,0.,0.],np.r_[np.tile(hi,3),1.,1.,1.,np.pi]


def unitary(vec):
    t=np.linalg.norm(vec)
    return np.eye(2,dtype=complex) if t==0 else np.cos(t/2)*PAULI[0]-1j*np.sin(t/2)*np.einsum('a,aij->ij',vec/t,PAULI[1:])


def channel_kraus(v,g):
    L=unpack(v)
    raw=np.einsum('ij,iab->jab',L,PAULI)
    S=sum(K.conj().T@K for K in raw)
    val,V=np.linalg.eigh(S)
    if val.min()<1e-14:raise ValueError('singular input marginal')
    inv=(V*(val**-.5))@V.conj().T
    U=unitary(IDEAL[g])
    return np.array([K@inv@U for K in raw])


def operations(x):
    x=np.asarray(x,float);lo,hi=bounds()
    if x.shape!=(52,) or np.any(~np.isfinite(x)) or np.any(x<lo-1e-8) or np.any(x>hi+1e-8):raise ValueError('invalid general-channel parameters')
    gates=np.zeros((2,3,5,5));channels=[]
    for g in range(3):
        kk=channel_kraus(x[16*g:16*(g+1)],g);channels.append(kk)
        M=np.zeros((5,5));M[0,0]=1;M[4,4]=1
        # Bloch affine coordinates; trace preservation is also checked separately.
        for j in range(4):
            image=sum(K@PAULI[j]@K.conj().T for K in kk)
            for i in range(1,4):M[i,j]=np.trace(PAULI[i]@image).real/2
        gates[:,g]=M
    contrast,l,u,theta=x[48:]
    initial=np.array([1.,0.,0.,contrast,0.])
    effect=np.array([(l+u)/2,(l-u)*np.sin(theta)/2,0.,(l-u)*np.cos(theta)/2,.5])
    return gates,initial,effect,channels


def from_restricted(x):
    out=[]
    for g in range(3):
        U=unitary(IDEAL[g]);kk=restricted_kraus(x,'qubit',0,g)
        coeff=np.array([[np.trace(sigma@K@U.conj().T)/2 for sigma in PAULI] for K in kk]).T
        chi=coeff@coeff.conj().T
        L=np.linalg.cholesky(chi+1e-15*np.eye(4))
        out.extend(pack(L))
    out.extend([x[14],x[12],x[13],0.])
    return np.array(out)


class GeneralEvaluator(Evaluator):
    def __init__(self,rows):
        self.order=sorted(range(len(rows)),key=lambda i:tuple(rows[i]['word']))
        self.inverse=np.argsort(self.order)
        ordered=[rows[i] for i in self.order]
        super().__init__(ordered)
        previous=();reuse=[]
        for row in ordered:
            w=row['word'];common=0
            for a,b in zip(previous,w):
                if a!=b:break
                common+=1
            reuse.append(common);previous=w
        self.reuse=np.array(reuse,dtype=np.int64)
        src=HERE/'propagate_general.cpp';sha=hashlib.sha256(src.read_bytes()).hexdigest()[:16]
        directory=Path(tempfile.gettempdir())/('sqd-general-'+sha);directory.mkdir(exist_ok=True);libpath=directory/'propagate.so'
        if not libpath.exists():subprocess.run(['g++','-O3','-std=c++11','-shared','-fPIC',str(src),'-o',str(libpath)],check=True)
        self.fast=ctypes.CDLL(str(libpath));pd=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS');pi=np.ctypeslib.ndpointer(dtype=np.int64,flags='C_CONTIGUOUS');pw=np.ctypeslib.ndpointer(dtype=np.int8,flags='C_CONTIGUOUS')
        self.fast.qubit_prefix.argtypes=[pi,pw,pi,ctypes.c_int64,pd,pd,pd,pd];self.fast.qubit_prefix.restype=None

    def __call__(self,x):
        gates,initial,effect,_=operations(x);out=np.empty(self.count)
        self.fast.qubit_prefix(self.offsets,self.words,self.reuse,self.count,gates,initial,effect,out)
        out=out[self.inverse]
        if np.any(out< -2e-8) or np.any(out>1+2e-8):raise ValueError('nonphysical general probabilities')
        return np.clip(out,0,1)


def certificate(x):
    _,_,_,channels=operations(x);tp=0.;cp=1.
    for kk in channels:
        tp=max(tp,float(np.max(abs(sum(K.conj().T@K for K in kk)-np.eye(2)))))
        J=sum(np.outer(K.reshape(-1,order='F'),K.reshape(-1,order='F').conj()) for K in kk)
        cp=min(cp,float(np.linalg.eigvalsh(J).min()))
    return {'trace_preservation_max_error':tp,'choi_min_eigenvalue':cp}


def density_probability(word,x):
    _,_,_,channels=operations(x);c,l,u,theta=x[48:]
    rho=(PAULI[0]+c*PAULI[3])/2
    E=(l+u)*PAULI[0]/2+(l-u)*(np.sin(theta)*PAULI[1]+np.cos(theta)*PAULI[3])/2
    for g in word:rho=sum(K@rho@K.conj().T for K in channels[g])
    return float(np.trace(E@rho).real)


def fit(rows,initial,starts=2,max_nfev=180,seed=260931):
    E=GeneralEvaluator(rows);k=np.array([r['k'] for r in rows]);n=np.array([r['n'] for r in rows]);lo,hi=bounds();rng=np.random.default_rng(seed);best=None
    for start in range(starts):
        x=np.array(initial).copy()
        if start:
            for g in range(3):
                x[g*16:(g+1)*16]+=rng.normal(0,.03,16)
        x=np.clip(x,lo+1e-10,hi-1e-10)
        result=least_squares(lambda v:residual(k,n,E(v)),x,bounds=(lo,hi),max_nfev=max_nfev,
            jac='3-point',diff_step=1e-3,ftol=2e-7,xtol=2e-7,gtol=2e-5,x_scale='jac')
        out={'x':result.x.tolist(),'deviance':float(2*result.cost),'success':bool(result.success),'nfev':int(result.nfev),'optimality':float(result.optimality),'physical_certificate':certificate(result.x)}
        if best is None or out['deviance']<best['deviance']:best=out
        print('general start',start,'deviance',out['deviance'],'nfev',out['nfev'],'success',out['success'],flush=True)
    return best
