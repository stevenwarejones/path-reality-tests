"""Exact physicality and rational probabilities for delay-dependent candidates."""
import numpy as np
import identifiability as physical


def rotation_matrix(numerators,denominator):
    if type(denominator) is not int or denominator<=0 or len(numerators)!=3 or any(type(k) is not int for k in numerators):
        raise ValueError('invalid Cayley rotation')
    a=np.array(numerators,dtype=object);d=denominator
    norm=sum(k*k for k in a);q=d*d+norm
    x,y,z=a;cross=np.array([[0,-z,y],[z,0,-x],[-y,x,0]],dtype=object)
    R=(d*d-norm)*np.eye(3,dtype=object)+2*np.outer(a,a)+2*d*cross
    if not np.array_equal(R.T@R,q*q*np.eye(3,dtype=object)):
        raise ValueError('rotation identity failed')
    det=R[0,0]*(R[1,1]*R[2,2]-R[1,2]*R[2,1])-R[0,1]*(R[1,0]*R[2,2]-R[1,2]*R[2,0])+R[0,2]*(R[1,0]*R[2,1]-R[1,1]*R[2,0])
    if det!=q**3:raise ValueError("rotation orientation failed")
    return R,q


def check_candidate(candidate):
    boundary=physical.check_instruments(candidate['instrument'])
    D=candidate['spam_denominator']
    if type(D) is not int or D<=0 or set(candidate['preparations'])!=set(physical.STATES) or set(candidate['effects'])!=set('xyz'):
        raise ValueError('invalid SPAM grid')
    for v in candidate['preparations'].values():
        if len(v)!=3 or any(type(k) is not int for k in v) or sum(k*k for k in v)>D*D:
            raise ValueError('nonphysical preparation')
    for effect in candidate['effects'].values():
        if len(effect)!=4 or any(type(k) is not int for k in effect):raise ValueError('invalid effect')
        bias,*vec=effect
        if abs(bias)>D or sum(k*k for k in vec)>(D-abs(bias))**2:
            raise ValueError('nonphysical effect')
    if len(candidate['delays'])!=9 or len(set(candidate['delays']))!=9 or len(candidate['rotation_cayley'])!=9:
        raise ValueError('invalid delay grid')
    for pair in candidate['rotation_cayley']:
        if len(pair)!=2:raise ValueError('missing pre/post process')
        for vec in pair:rotation_matrix(vec,candidate['rotation_denominator'])
    return boundary


def probabilities(candidate,labels):
    """Exact Pauli transfer contraction, with all four outcomes retained."""
    check_candidate(candidate)
    witness=candidate['instrument'];D=witness['denominator'];S=candidate['spam_denominator']
    sigmas=[physical.I,*physical.PAULIS]
    transfers=[]
    for pair in witness['choi']:
        ts=[]
        for block in pair:
            re,im=physical.arrays(block);T=np.zeros((4,4),dtype=object)
            for mu,out in enumerate(sigmas):
                for nu,inp in enumerate(sigmas):
                    coefficient=np.kron(inp.T,out).T
                    cr=np.rint(coefficient.real).astype(int);ci=np.rint(coefficient.imag).astype(int)
                    T[mu,nu]=int(np.sum(re*cr-im*ci))
            ts.append(T)
        transfers.append(ts)
    nums=[];dens=[]
    for pair in candidate['rotation_cayley']:
        before,db=rotation_matrix(pair[0],candidate['rotation_denominator'])
        after,da=rotation_matrix(pair[1],candidate['rotation_denominator'])
        denominator=4*D*S*S*db*da;rows=[]
        for label in labels:
            a,m,p,z=label.split(',');y=witness['settings'].index(m+','+p)
            inp=np.r_[S*db,before@np.array(candidate['preparations'][a],dtype=object)]
            bias,*vector=candidate['effects'][z]
            pull=np.array(vector,dtype=object)@after
            row=[]
            for b in range(2):
                output=transfers[y][b]@inp
                for sign in (1,-1):
                    effect=np.r_[(S+sign*bias)*da,sign*pull]
                    row.append(int(effect@output))
            if min(row)<=0 or sum(row)!=denominator:raise ValueError('invalid probabilities')
            rows.append(row)
        nums.append(rows);dens.append([denominator]*len(labels))
    return np.array(nums,dtype=object),np.array(dens,dtype=object)
