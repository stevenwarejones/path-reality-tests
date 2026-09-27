"""Exact block-leakage similarity with gate dependence and nonunital dynamics.

Computational subspace is stipulated. This is an operational nonidentification
construction, not an estimate of a measured population.
"""
import numpy as np
from sqd_models import PAULI

TAU=np.eye(2)/2

def apply(kk,rho):return sum(K@rho@K.conj().T for K in kk)

def reset_kraus(tau):
    vals,V=np.linalg.eigh(tau)
    if vals.min() < -1e-12 or abs(np.trace(tau)-1)>1e-12:raise ValueError('nonphysical return state')
    return [np.sqrt(max(v,0))*np.outer(V[:,i],np.eye(2)[j]) for i,v in enumerate(vals) for j in range(2)]

def validate(kk,loss,ret,tau):
    if not 0<=loss<=1 or not 0<=ret<=1:raise ValueError('nonphysical transfer rate')
    if np.max(abs(sum(K.conj().T@K for K in kk)-np.eye(2)))>1e-10:raise ValueError('not trace preserving')
    reset_kraus(tau)

def matrix(kk,loss,ret,tau=TAU):
    validate(kk,loss,ret,tau)
    T=np.array([[np.trace(P@apply(kk,Q)).real/2 for Q in PAULI] for P in PAULI])
    M=np.zeros((5,5));M[:4,:4]=(1-loss)*T
    M[:4,4]=ret*np.array([np.trace(P@tau).real for P in PAULI]);M[4,0]=loss;M[4,4]=1-ret
    return M

def qutrit_kraus(kk,loss,ret,tau=TAU):
    validate(kk,loss,ret,tau);out=[]
    for K in kk:
        A=np.zeros((3,3),complex);A[:2,:2]=np.sqrt(1-loss)*K;out.append(A)
    for i in range(2):
        A=np.zeros((3,3),complex);A[2,i]=np.sqrt(loss);out.append(A)
    vals,V=np.linalg.eigh(tau)
    for i,v in enumerate(vals):
        A=np.zeros((3,3),complex);A[:2,2]=np.sqrt(ret*max(v,0))*V[:,i];out.append(A)
    A=np.zeros((3,3),complex);A[2,2]=np.sqrt(1-ret);out.append(A)
    return out

def transform(kk,loss,ret,kappa,tau=TAU):
    if not 0<kappa<=1:raise ValueError('kappa outside proved interval')
    validate(kk,loss,ret,tau)
    lp=kappa*loss;ep=ret+(1-kappa)*loss
    if ep<=0 or lp>=1:raise ValueError('degenerate transformation; use continuous limit separately')
    kp=[np.sqrt((1-loss)/(1-lp))*K for K in kk]+[np.sqrt((1-kappa)*loss/(1-lp))*K for K in reset_kraus(TAU)]
    tp=TAU+ret/(kappa*ep)*(tau-TAU)+(1-kappa)*(1-loss)/(kappa*ep)*(TAU-apply(kk,TAU))
    validate(kp,lp,ep,tp)
    return kp,lp,ep,tp

def effect_transform(E,z,kappa):
    a=np.trace(E@TAU).real;zp=a+(z-a)/kappa
    if not 0<=zp<=1:raise ValueError('transformed leakage effect outside [0,1]')
    return zp

def lower_kappa(kk,loss,ret,E,z):
    d=np.array([np.trace(P@apply(kk,TAU)).real for P in PAULI[1:]])
    c=(1-loss)*np.linalg.norm(d);b=ret+loss+c
    physical=0. if c==0 else 2*c/(b+np.sqrt(max(0,b*b-4*loss*c)))
    a=np.trace(E@TAU).real
    readout=(a-z)/a if z<a else (z-a)/(1-a) if z>a else 0.
    rates=max(0.,1-(1-ret)/loss) if loss else 0.
    return max(float(physical),float(readout),float(rates))

def similarity(kappa):
    T=np.eye(5);T[0,4]=1-kappa;T[4,4]=kappa;return T

def amplitude_damping(rate):
    return [np.diag([1.,np.sqrt(1-rate)]),np.array([[0,np.sqrt(rate)],[0,0]])]


def memory_probability(word,channels,loss,ret,returns,rho,E,z):
    """Independent two-state classical-flag/qubit instrument realization.

    The inactive branch deliberately discards its qubit state. Its trace carries
    the flag population; no gate counter or circuit label is stored.
    """
    active=np.array(rho,dtype=complex);inactive=np.zeros((2,2),complex)
    for g in word:
        validate(channels[g],loss[g],ret[g],returns[g])
        mass_a=np.trace(active).real;mass_i=np.trace(inactive).real
        active_next=(1-loss[g])*apply(channels[g],active)+ret[g]*mass_i*returns[g]
        inactive_next=(loss[g]*mass_a+(1-ret[g])*mass_i)*TAU
        active,inactive=active_next,inactive_next
    return float(np.trace(E@active).real+z*np.trace(inactive).real)


def lift_qubit(kk,loss,ret):
    """Exact observable lift of an interior qubit channel to block leakage.

    Q(rho_C,L)=rho_C+L I/2 intertwines the lifted map with the input
    channel. The leakage effect must equal Tr(E I/2); return polarization is
    chosen explicitly. No inference about a device population is made.
    """
    if not 0<=loss<1 or not 0<ret<=1:raise ValueError('invalid lift rates')
    J=sum(np.outer(K.reshape(-1,order='F'),K.reshape(-1,order='F').conj()) for K in kk)
    B=(J-loss*np.kron(np.eye(2),TAU))/(1-loss)
    vals,V=np.linalg.eigh(B)
    if vals.min() < -1e-12:raise ValueError('loss exceeds CP replacement budget')
    base=[(np.sqrt(max(v,0))*V[:,i]).reshape((2,2),order='F') for i,v in enumerate(vals)]
    translation=apply(kk,TAU)-TAU
    # TP makes this exactly traceless; remove floating-point trace error before
    # dividing by a small return rate.
    translation-=np.trace(translation)*TAU
    return_state=TAU+translation/ret
    validate(base,loss,ret,return_state)
    return base,return_state
