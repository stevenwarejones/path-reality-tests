"""Directed-interval gate construction and rational propagation error bound.

Physicality is structural (normalized Kraus factors), not inferred from small
negative floating eigenvalues. See theory.md for the accumulated error bound.
"""
from fractions import Fraction as F
import math
import numpy as np
from mpmath import iv
from gain_core import operations
from sqd_general import SLOTS

iv.dps=60

def zero():return [[iv.mpc(0),iv.mpc(0)],[iv.mpc(0),iv.mpc(0)]]
def add(A,B):return [[A[i][j]+B[i][j] for j in range(2)] for i in range(2)]
def scale(a,A):return [[a*A[i][j] for j in range(2)] for i in range(2)]
def mul(A,B):return [[sum(A[i][k]*B[k][j] for k in range(2)) for j in range(2)] for i in range(2)]
def adj(A):return [[iv.mpc(A[j][i].real, -A[j][i].imag) for j in range(2)] for i in range(2)]
def trace(A):return A[0][0]+A[1][1]
def real(x):return x.real

def upper_abs(v):
    # Convert an interval upper magnitude outward to binary64, then to an exact
    # rational grid. nextafter protects the conversion from downward rounding.
    a=max(abs(float(v.a)),abs(float(v.b)))
    v=F.from_float(float(np.nextafter(a,np.inf)))
    return F((v.numerator*10**30+v.denominator-1)//v.denominator,10**30)


def certify(x,max_length):
    x=np.asarray(x,float)
    if len(x)!=52 or not np.all(np.isfinite(x)):raise ValueError('invalid factor parameters')
    if not (0<=x[48]<=1 and 0<=x[49]<=1 and 0<=x[50]<=1 and 0<=x[51]<=np.pi):raise ValueError('nonphysical SPAM')
    P=[[[1,0],[0,1]],[[0,1],[1,0]],[[0,-1j],[1j,0]],[[1,0],[0,-1]]]
    P=[[[iv.mpc(complex(v).real,complex(v).imag) for v in row] for row in A] for A in P]
    U=[P[0],scale(1/iv.sqrt(2),add(P[0],scale(iv.mpc(0,-1),P[1]))),scale(1/iv.sqrt(2),add(P[0],scale(iv.mpc(0,-1),P[2])))]
    gates,initial,effect,_=operations(x);epsilon=F(0);s_margins=[]
    for g in range(3):
        L=[[iv.mpc(0) for j in range(4)] for i in range(4)]
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
        if float(det.a)<=0 or float(tr.a)<=0:raise ValueError('input marginal not certified positive')
        s_margins.append({'determinant_lower':float(np.nextafter(float(det.a),-np.inf)),'trace_lower':float(np.nextafter(float(tr.a),-np.inf))})
        d=iv.sqrt(det);C=add(S,scale(d,P[0]));dc=real(C[0][0]*C[1][1]-C[0][1]*C[1][0])
        inv=scale(iv.sqrt(tr+2*d)/dc,[[C[1][1],-C[0][1]],[-C[1][0],C[0][0]]])
        ks=[mul(mul(A,inv),U[g]) for A in raw]
        for j in range(4):
            B=zero()
            for K in ks:B=add(B,mul(mul(K,P[j]),adj(K)))
            for i in range(1,4):
                exact=real(trace(mul(P[i],B)))/2
                epsilon=max(epsilon,upper_abs(exact-iv.mpf(float(gates[0,g,i,j]))))
    l,u,theta=[iv.mpf(float(v)) for v in x[49:52]];a=(l+u)/2;b=(l-u)/2
    exact_effect=[a,b*iv.sin(theta),iv.mpf(0),b*iv.cos(theta)]
    eps_effect=max(upper_abs(v-iv.mpf(float(w))) for v,w in zip(exact_effect,effect[:4]))
    unit=F(1,2**53);gamma=7*unit/(1-7*unit)
    D=epsilon+gamma*(1+epsilon)+F(32,2**1022)
    # sqrt(3) < 7/4, so a=3 sqrt(3) D <=21 D/4, b<=7 D.
    A=F(21,4)*D;B=7*D
    if max_length*A>=1:raise ValueError('propagation bound unavailable')
    state_error=max_length*B/(1-max_length*A)
    D_effect=eps_effect+gamma*(1+eps_effect)+F(32,2**1022)
    probability_error=state_error/2+(4+3*state_error)*D_effect
    # Save a rational outward envelope for all subsequent exact comparisons.
    denominator=10**20;numerator=(probability_error.numerator*denominator+probability_error.denominator-1)//probability_error.denominator
    return {'method':'60-digit directed intervals for the exact normalized Kraus model; rational binary64 error accumulation using CPTP trace-distance contraction',
        'interval_digits':60,'maximum_word_length':int(max_length),'input_marginals':s_margins,
        'gate_entry_error_upper':float(epsilon),'effect_entry_error_upper':float(eps_effect),
        'probability_error_numerator':int(numerator),'probability_error_denominator':denominator,
        'probability_error_upper':numerator/denominator,'physicality':'positive input marginal; normalized Kraus maps exactly CPTP; physical state and binary effect eigenvalues'}
