"""Exact all-word similarity, interval physicality, and an EB-status obstruction.

The trusted Bell-pair task is an additional operation, not part of the archived
binary gate-set experiment. No gate parameter is called an observed resource.
"""
from fractions import Fraction as F
import numpy as np
from mpmath import iv
from dynamics_common import *
from interval_physics import zero,add,scale,mul,adj,trace,real
from sqd_general import SLOTS
from sqd_models import PAULI

iv.dps=60

def matmul(A,B):
    return [[sum(A[i][k]*B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]
def dagger(A):return [[iv.mpc(A[j][i].real,-A[j][i].imag) for j in range(len(A))] for i in range(len(A[0]))]
def eye(n):return [[iv.mpc(int(i==j)) for j in range(n)] for i in range(n)]
def power(A,n):
    if not isinstance(n,int) or n<0 or n>20000:raise ValueError('unsupported target length')
    out=eye(len(A))
    while n:
        if n&1:out=matmul(A,out)
        n//=2
        if n:A=matmul(A,A)
    return out

def exact_gates(x):
    """Exact dyadic parameters, rational factor scales and nominal quarter turns.

Uses the validated PR17 interval arithmetic helpers, returning the matrices
needed for a new resource task instead of only its binary-word error bound.
    """
    P=[[[iv.mpc(complex(v).real,complex(v).imag) for v in row] for row in A] for A in PAULI]
    U=[P[0]]+[scale(1/iv.sqrt(2),add(P[0],scale(iv.mpc(0,-1),P[g]))) for g in [1,2]]
    out=[]
    for g in range(3):
        L=[[iv.mpc(0) for _ in range(4)] for _ in range(4)]
        for value,(i,j,part) in zip(x[g*16:(g+1)*16],SLOTS):
            factor=F(1) if i==j==0 else F(1,100) if i==j else F(1,1000)
            z=iv.mpf(float(value))*factor.numerator/factor.denominator
            L[i][j]+=iv.mpc(0,z) if part else z
        raw=[]
        for j in range(4):
            A=zero()
            for i in range(4):A=add(A,scale(L[i][j],P[i]))
            raw.append(A)
        S=zero()
        for A in raw:S=add(S,mul(adj(A),A))
        det=real(S[0][0]*S[1][1]-S[0][1]*S[1][0]);tr=real(trace(S))
        if float(det.a)<=0 or float(tr.a)<=0:raise ValueError('singular exact Kraus marginal')
        d=iv.sqrt(det);C=add(S,scale(d,P[0]));dc=real(C[0][0]*C[1][1]-C[0][1]*C[1][0])
        inv=scale(iv.sqrt(tr+2*d)/dc,[[C[1][1],-C[0][1]],[-C[1][0],C[0][0]]])
        ks=[mul(mul(A,inv),U[g]) for A in raw]
        M=eye(4)
        for i in range(1,4):
            for j in range(4):
                B=zero()
                for K in ks:B=add(B,mul(mul(K,P[j]),adj(K)))
                M[i][j]=real(trace(mul(P[i],B)))/2
        out.append(M)
    return out,P

def similarity(M,s):
    S=[iv.mpf(1)]+[iv.mpf(s.numerator)/s.denominator]*3
    return [[M[i][j]*S[i]/S[j] for j in range(4)] for i in range(4)]

def choi(M,P):
    # J=sum_j sigma_j^T tensor Phi(sigma_j)/2, input-first convention.
    J=[[iv.mpc(0) for _ in range(4)] for _ in range(4)]
    for j in range(4):
        B=zero()
        for i in range(4):B=add(B,scale(M[i][j],P[i]))
        for a in range(2):
            for b in range(2):
                for c in range(2):
                    for d in range(2):J[2*a+b][2*c+d]+=P[j][c][a]*B[b][d]/2
    return J

def partial_transpose(J):return [[J[2*(i//2)+j%2][2*(j//2)+i%2] for j in range(4)] for i in range(4)]

def positive_definite(A,margin=F(0)):
    """Interval LDL*: certify A-margin*I is strictly positive, no eigenvalue tolerance."""
    n=len(A);L=eye(n);D=[]
    for j in range(n):
        d=(A[j][j]-iv.mpf(margin.numerator)/margin.denominator-sum(L[j][k]*D[k]*iv.mpc(L[j][k].real,-L[j][k].imag) for k in range(j))).real
        if float(d.a)<=0:raise ValueError('positive definiteness not certified')
        D.append(d)
        for i in range(j+1,n):L[i][j]=(A[i][j]-sum(L[i][k]*D[k]*iv.mpc(L[j][k].real,-L[j][k].imag) for k in range(j)))/d
    return True

def midpoint(A):return np.array([[complex(float(v.real.mid),float(v.imag.mid)) for v in row] for row in A])
def rayleigh(A,v):
    z=[[iv.mpc(complex(a).real,complex(a).imag)] for a in v]
    return (matmul(matmul(dagger(z),A),z)[0][0]/matmul(dagger(z),z)[0][0]).real

def verify_resource(rayleigh_vector=None):
    x=baseline_models()['joint']['x'];gates,P=exact_gates(x);n=10749;scales=[F(1),F(1009,1000)];out={}
    c,l,u,theta=[iv.mpf(float(v)) for v in x[48:]];a=(l+u)/2;b=(u-l)/2
    for s in scales:
        ss=iv.mpf(s.numerator)/s.denominator
        if float((c*ss).b)>=1 or float((a-b/ss).a)<0 or float((1-a-b/ss).a)<0:raise ValueError('similarity produces nonphysical SPAM')
        G=[similarity(g,s) for g in gates]
        for g in G:positive_definite(choi(g,P),F(1,10**6))
        J=choi(power(G[2],n),P);PT=partial_transpose(J)
        # J/2 is the output for an independently trusted Bell-pair input.
        densityPT=[[v/2 for v in row] for row in PT]
        if s==1:
            positive_definite(densityPT,F(3,10**6))
            status={'status':'entanglement breaking','partial_transpose_eigenvalue_lower':3e-6,'method':'interval LDL of partially transposed normalized Choi minus 3e-6 I'}
        else:
            if rayleigh_vector is None:
                _,V=np.linalg.eigh(midpoint(densityPT));v=V[:,0]
            else:v=np.array([complex(*z) for z in rayleigh_vector])
            q=rayleigh(densityPT,v)
            if float(q.b)>=-3e-6:raise ValueError('negative partial transpose not certified')
            status={'status':'not entanglement breaking','negative_rayleigh_upper':-3e-6,'rayleigh_vector':[[float(z.real),float(z.imag)] for z in v],'method':'directed-interval Rayleigh quotient of normalized Choi partial transpose'}
        out[str(s)]={'scale':[s.numerator,s.denominator],'state_bloch_length_interval':[float(np.nextafter(float((c*ss).a),-np.inf)),float(np.nextafter(float((c*ss).b),np.inf))],'native_choi_eigenvalue_lower':1e-6,**status}
    # Deliberately broad rational outward bounds, independently checked below.
    delta=(iv.mpf(9)/1000)*c/2
    if not (float(delta.a)>.0044564 and float(delta.b)<.0044565):raise ValueError('calibration contrast changed')
    return {'target_expression':'(Gy)^10749','target_length':n,'word_sha256':target_hash((2,)*n),
        'models':out,'exact_word_equivalence':'S_s G_w S_s^-1, S_s rho, E S_s^-1 for every word; all original binary probabilities equal the joint-compatible base point',
        'trusted_Z_preparation_probability_gap':[.0044564,.0044565],
        'scope':'nonidentification of EB status for this fixed block from arbitrary archived-style terminal words, not impossibility of interval tightening or a measured memory resource'}

if __name__=='__main__':
    result=verify_resource();save(HERE/'results/resource-certificate.json',result);print({k:v['status'] for k,v in result['models'].items()})
