"""CPTP shared operations in a fixed computational basis, plus independent witnesses."""
import ctypes
import hashlib
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from scipy.spatial.transform import Rotation

HERE = Path(__file__).resolve().parent
PAULI = np.array([[[1,0],[0,1]], [[0,1],[1,0]], [[0,-1j],[1j,0]], [[1,0],[0,-1]]], complex)
IDEAL = np.array([[0,0,0], [np.pi/2,0,0], [0,np.pi/2,0]])


def parameters(model):
    # delta rotation(9), depolarizing exponent(3), plus response to 0/1(2), prep contrast(1)
    x = np.r_[np.zeros(9), np.full(3, .8), .003, .997, .999]
    lo = np.r_[np.full(9, -10.), np.zeros(3), 0., .9, .9]
    hi = np.r_[np.full(9, 10.), np.full(3, 10.), .1, 1., 1.]
    if model == 'leakage':
        x = np.r_[x, .1, .1, .5]
        lo = np.r_[lo, 0., 0., 0.]
        hi = np.r_[hi, 10., 10., 1.]
    elif model in ('alternating', 'quasistatic'):
        x = np.r_[x, 1.]
        lo = np.r_[lo, 0.]
        hi = np.r_[hi, 5.]
    elif model != 'qubit':
        raise ValueError('unknown model')
    return x, lo, hi


def operations(x, model):
    x = np.asarray(x, float)
    _, lo, hi = parameters(model)
    if x.shape != lo.shape or np.any(~np.isfinite(x)) or np.any(x < lo-1e-10) or np.any(x > hi+1e-10):
        raise ValueError('parameters outside declared physical model')
    delta = x[:9].reshape(3,3)*.001
    eta = np.exp(-x[9:12]*.0001)
    l, u, contrast = x[12:15]
    leak, ret, z = (x[15]*.0001, x[16]*.0001, x[17]) if model == 'leakage' else (0.,0.,.5)
    theta = x[15]*.01 if model in ('alternating','quasistatic') else 0.
    mats = np.zeros((2,3,5,5))
    rotations = np.zeros((2,3,3,3))
    for b, sign in enumerate((1.,-1.)):
        for g in range(3):
            vec = IDEAL[g] + delta[g]
            if g > 0:
                vec[g-1] += sign*theta
            R = Rotation.from_rotvec(vec).as_matrix()
            rotations[b,g] = R
            M = mats[b,g]
            M[0,0], M[0,4] = 1-leak, ret
            M[1:4,1:4] = (1-leak)*eta[g]*R
            M[4,0], M[4,4] = leak, 1-ret
    initial = np.array([1.,0.,0.,contrast,0.])
    effect = np.array([(l+u)/2,0.,0.,(l-u)/2,z])
    return mats, initial, effect, rotations


def reference_probability(word, x, model):
    gates, initial, effect, _ = operations(x, model)
    p = 0.
    for branch in range(2):
        state = initial.copy()
        for t,g in enumerate(word):
            phase = branch ^ (t%2) if model == 'alternating' else branch
            state = gates[phase,g] @ state
        p += effect @ state / 2
    return float(p)


class Evaluator:
    def __init__(self, rows):
        # Build locally from reviewed source. Never execute downloaded source files.
        src = HERE/'propagate.cpp'
        sha = hashlib.sha256(src.read_bytes()).hexdigest()[:16]
        build = Path(tempfile.gettempdir()) / ('sqd-native-'+sha)
        build.mkdir(exist_ok=True)
        libpath = build/'propagate.so'
        if not libpath.exists():
            subprocess.run(['g++','-O3','-std=c++11','-shared','-fPIC',str(src),'-o',str(libpath)],check=True)
        self.lib = ctypes.CDLL(str(libpath))
        ptrd = np.ctypeslib.ndpointer(dtype=np.float64, flags='C_CONTIGUOUS')
        self.lib.propagate.argtypes = [np.ctypeslib.ndpointer(dtype=np.int64,flags='C_CONTIGUOUS'),np.ctypeslib.ndpointer(dtype=np.int8,flags='C_CONTIGUOUS'),ctypes.c_int64,ptrd,ptrd,ptrd,ctypes.c_int,ptrd]
        self.lib.propagate.restype = None
        self.offsets = np.r_[0,np.cumsum([len(r['word']) for r in rows])].astype(np.int64)
        self.words = np.fromiter((g for r in rows for g in r['word']),np.int8)
        self.count = len(rows)

    def __call__(self, x, model):
        gates, initial, effect, _ = operations(x,model)
        out = np.empty(self.count)
        self.lib.propagate(self.offsets,self.words,self.count,gates,initial,effect,int(model=='alternating'),out)
        if np.any(out < -1e-9) or np.any(out > 1+1e-9):
            raise ValueError('invalid propagated probability')
        return np.clip(out,0,1)


def kraus(x, model, phase, gate):
    """Independent complex Hilbert-space realization, not the Bloch propagator."""
    _,_,_,rots = operations(x,model)
    # Build unitary directly from rotation vector, not from the real transfer matrix.
    vec = IDEAL[gate] + np.asarray(x[:9]).reshape(3,3)[gate]*.001
    if model in ('alternating','quasistatic') and gate:
        vec[gate-1] += (1 if phase==0 else -1)*x[15]*.01
    a = np.linalg.norm(vec)
    U = np.eye(2,dtype=complex) if a == 0 else np.cos(a/2)*PAULI[0] - 1j*np.sin(a/2)*np.einsum('a,aij->ij',vec/a,PAULI[1:])
    eta = np.exp(-x[9+gate]*.0001)
    probs = [(1+3*eta)/4]+[(1-eta)/4]*3
    comp = [np.sqrt(p)*sigma@U for p,sigma in zip(probs,PAULI)]
    if model != 'leakage':
        return comp
    lam,ret = x[15]*.0001,x[16]*.0001
    out=[]
    for K in comp:
        A=np.zeros((3,3),complex);A[:2,:2]=np.sqrt(1-lam)*K;out.append(A)
    for i in range(2):
        A=np.zeros((3,3),complex);A[2,i]=np.sqrt(lam);out.append(A)
        A=np.zeros((3,3),complex);A[i,2]=np.sqrt(ret/2);out.append(A)
    A=np.zeros((3,3),complex);A[2,2]=np.sqrt(1-ret);out.append(A)
    return out


def density_probability(word,x,model):
    """Slow independent physical certificate including sign/reset convention."""
    d = 3 if model=='leakage' else 2
    rho=np.zeros((d,d),complex);rho[:2,:2]=(PAULI[0]+x[14]*PAULI[3])/2
    E=np.diag([x[12],x[13]]+([x[17]] if d==3 else []))
    channels = [[kraus(x,model,b,g) for g in range(3)] for b in range(2)]
    ans=0.
    for branch in range(2):
        r=rho.copy()
        for t,g in enumerate(word):
            b=branch^(t%2) if model=='alternating' else branch
            r=sum(K@r@K.conj().T for K in channels[b][g])
        ans+=np.trace(E@r).real/2
    return float(ans)


def certificate(x,model):
    residual=0.;mineig=1.
    for b in range(2):
        for g in range(3):
            kk=kraus(x,model,b,g);d=kk[0].shape[0]
            residual=max(residual,float(np.max(np.abs(sum(K.conj().T@K for K in kk)-np.eye(d)))))
            # Unnormalized Choi, column vectorization; marginal identity.
            C=sum(np.outer(K.reshape(-1,order='F'),K.reshape(-1,order='F').conj()) for K in kk)
            mineig=min(mineig,float(np.linalg.eigvalsh(C).min()))
    return {'trace_preservation_max_error':residual,'choi_min_eigenvalue':mineig}


def dormant_flag_probability(word,x):
    """Two classical memory states with 2x2 subnormalized quantum blocks.

    Active flag applies the qubit channel, loses trace lam to inactive I/2.
    Inactive flag returns trace ret to active I/2. Same gates at every occurrence.
    """
    lam,ret,z=x[15]*.0001,x[16]*.0001,x[17]
    active=(PAULI[0]+x[14]*PAULI[3])/2
    inactive=np.zeros((2,2),complex)
    E=np.diag(x[12:14])
    base=np.array(x[:15]); channels=[kraus(base,'qubit',0,g) for g in range(3)]
    for g in word:
        a=(1-lam)*sum(K@active@K.conj().T for K in channels[g])+ret*np.trace(inactive)*PAULI[0]/2
        b=lam*np.trace(active)*PAULI[0]/2+(1-ret)*inactive
        active,inactive=a,b
    return float((np.trace(E@active)+z*np.trace(inactive)).real)
