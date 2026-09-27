"""Independent exact verification of shared-instrument classical-memory models."""
import numpy as np
import identifiability as physical


def group_matrices(group):
    D=group['denominator']
    if type(D) is not int or D<=0:raise ValueError('invalid Choi denominator')
    out=[]
    for record in group['choi']:
        if any(type(x) is not int for key in ('real','imag') for x in record[key]):raise ValueError('noninteger Choi coefficient')
        r,i=physical.arrays(record);physical.ldl_positive(physical.realify(r,i));out.append((r,i))
    if not out:raise ValueError('empty CP group')
    return out,D


def check_group(group,instrument=True):
    matrices,D=group_matrices(group)
    traces=[(physical.partial_trace_output(r),physical.partial_trace_output(i)) for r,i in matrices]
    check=[(sum(r for r,i in traces),sum(i for r,i in traces))] if instrument else traces
    for r,i in check:
        if not np.array_equal(r,D*np.eye(2,dtype=object)) or np.any(i):raise ValueError('not trace preserving')
    return matrices,D


def check_model(model):
    settings=model['settings']
    if len(settings)!=18 or set(settings)!={m+','+p for m in 'xyz' for p in physical.STATES}:raise ValueError('setting grid')
    if len(model['instruments'])!=18:raise ValueError('missing instruments')
    first={}
    for y,group in zip(settings,model['instruments']):
        matrices,D=check_group(group)
        if len(matrices)!=2:raise ValueError('both recorded flags required')
        for b,(r,i) in enumerate(matrices):
            effect=(physical.partial_trace_output(r),physical.partial_trace_output(i),D);key=(y[0],b)
            if key in first:
                oldr,oldi,oldD=first[key]
                if not np.array_equal(effect[0]*oldD,oldr*D) or not np.array_equal(effect[1]*oldD,oldi*D):raise ValueError('unshared first effect')
            first[key]=effect
    S=model['spam_denominator']
    if type(S) is not int or S<=0 or set(model['preparations'])!=set(physical.STATES) or set(model['effects'])!=set('xyz'):raise ValueError('SPAM grid')
    for state in model['preparations'].values():
        if len(state)!=3 or any(type(x) is not int for x in state) or sum(x*x for x in state)>S*S:raise ValueError('nonphysical preparation')
    for effect in model['effects'].values():
        if len(effect)!=4 or any(type(x) is not int for x in effect):raise ValueError('invalid effect')
        bias,*vec=effect
        if abs(bias)>S or sum(x*x for x in vec)>(S-abs(bias))**2:raise ValueError('nonphysical effect')
    if len(model['processes'])!=9 or len(set(model['delays']))!=9 or len(model['delays'])!=9:raise ValueError('delay grid')
    for process in model['processes']:
        before,_=check_group(process['before']);after,_=check_group(process['after'],False)
        if len(before)!=len(after):raise ValueError('missing classical record branch')
    return True


def transfer(matrix):
    re,im=matrix;sigmas=[physical.I,*physical.PAULIS];out=np.zeros((4,4),dtype=object)
    for mu,a in enumerate(sigmas):
        for nu,b in enumerate(sigmas):
            c=np.kron(b.T,a).T
            out[mu,nu]=int(np.sum(re*np.rint(c.real).astype(int)-im*np.rint(c.imag).astype(int)))
    return out  # denominator 2*Choi denominator


def probabilities(model,labels):
    check_model(model);S=model['spam_denominator']
    middle=[]
    denominators={group['denominator'] for group in model['instruments']}
    if len(denominators)!=1:raise ValueError('unsupported mixed instrument denominators')
    DB=denominators.pop()
    for group in model['instruments']:middle.append([transfer(m) for m in group_matrices(group)[0]])
    nums=[];dens=[]
    for process in model['processes']:
        before,DA=group_matrices(process['before']);after,DC=group_matrices(process['after'])
        pre=[transfer(m) for m in before];post=[transfer(m) for m in after]
        denominator=16*S*S*DA*DB*DC;rows=[]
        for label in labels:
            a,m,p,z=label.split(',');y=model['settings'].index(m+','+p)
            state=np.array([S,*model['preparations'][a]],dtype=object)
            bias,*vec=model['effects'][z];row=[]
            for b in range(2):
                output=sum(C@middle[y][b]@A@state for A,C in zip(pre,post))
                for sign in (1,-1):row.append(int(np.array([S+sign*bias,*[sign*x for x in vec]],dtype=object)@output))
            if min(row)<0 or sum(row)!=denominator:raise ValueError('nonphysical probability table')
            rows.append(row)
        nums.append(rows);dens.append([denominator]*len(labels))
    return np.array(nums,dtype=object),np.array(dens,dtype=object)


def mix_instruments(model,weight):
    """Post-fit physical probe direction; never required by classical discovery."""
    import copy
    from fractions import Fraction as F
    t=F(*weight)
    if not 0<=t<1:raise ValueError('invalid probe mixture')
    out=copy.deepcopy(model)
    for group in out['instruments']:
        D=group['denominator'];records=[]
        for r,i in group_matrices(group)[0]:
            re=4*(t.denominator-t.numerator)*r+t.numerator*D*np.eye(4,dtype=object)
            im=4*(t.denominator-t.numerator)*i
            records.append(dict(real=[int(x) for x in re.ravel()],imag=[int(x) for x in im.ravel()]))
        group['denominator']=4*t.denominator*D;group['choi']=records
    out['post_fit_instrument_mixture']=weight
    check_model(out)
    return out


def quantum_instrument(model,weight=None):
    """Assess the exact counterpart only after an independent classical fit."""
    from fractions import Fraction as F
    witness=dict(settings=model['settings'],denominator=model['instruments'][0]['denominator'],
                 choi=[g['choi'] for g in model['instruments']],memory_weight=[1,10**40])
    if len({g['denominator'] for g in model['instruments']})!=1:raise ValueError('mixed denominators')
    endpoint=physical.check_instruments(witness)
    w=F(*weight) if weight is not None else F(99,100)*endpoint
    witness['memory_weight']=[w.numerator,w.denominator]
    physical.check_instruments(witness)
    return witness,endpoint


def process_choi(process):
    """Exact classical and SWAP-return numerators, common denominator 2 DA DC."""
    before,DA=group_matrices(process['before']);after,DC=group_matrices(process['after'])
    wr=np.zeros((16,16),dtype=object);wi=wr.copy();jr=np.zeros((4,4),dtype=object);ji=jr.copy()
    for (ar,ai),(cr,ci) in zip(before,after):
        wr+=2*(np.kron(ar,cr)-np.kron(ai,ci));wi+=2*(np.kron(ar,ci)+np.kron(ai,cr))
        for a in range(2):
            for aa in range(2):
                for c in range(2):
                    for cc in range(2):
                        for r in range(2):
                            for s in range(2):
                                av,bv=ar[2*a+r,2*aa+s],ai[2*a+r,2*aa+s]
                                cv,dv=cr[2*r+c,2*s+cc],ci[2*r+c,2*s+cc]
                                jr[2*a+c,2*aa+cc]+=av*cv-bv*dv
                                ji[2*a+c,2*aa+cc]+=av*dv+bv*cv
    sr=np.zeros((16,16),dtype=object);si=sr.copy()
    for a in range(2):
        for aa in range(2):
            for b in range(2):
                for o in range(2):
                    for c in range(2):
                        for cc in range(2):
                            row,col=8*a+4*b+2*o+c,8*aa+4*b+2*o+cc
                            sr[row,col]=jr[2*a+c,2*aa+cc];si[row,col]=ji[2*a+c,2*aa+cc]
    denominator=2*DA*DC
    if np.trace(wr)!=4*denominator or np.trace(sr)!=4*denominator:raise ValueError('process normalization')
    return wr,wi,sr,si,denominator


def npt_expectation(model,certificate):
    """Exact negative expectation, not a floating eigenvalue acceptance test."""
    from fractions import Fraction as F
    w=F(*certificate['memory_weight']);quantum_instrument(model,certificate['memory_weight'])
    index=certificate['delay_index'];wr,wi,sr,si,D=process_choi(model['processes'][index])
    qr=(w.denominator-w.numerator)*wr+w.numerator*sr
    qi=(w.denominator-w.numerator)*wi+w.numerator*si
    transpose=lambda a:a.reshape([2]*8).transpose(4,5,2,3,0,1,6,7).reshape(16,16)
    qr,qi=transpose(qr),transpose(qi)
    re=np.array(certificate['vector_real'],dtype=object);im=np.array(certificate['vector_imag'],dtype=object)
    if re.shape!=(16,) or im.shape!=(16,) or any(type(x) is not int for x in [*re,*im]):raise ValueError('invalid vector')
    norm=int(re@re+im@im)
    if norm==0:raise ValueError('zero vector')
    value=F(int(re@qr@re+im@qr@im-re@qi@im+im@qi@re),4*D*w.denominator*norm)
    if value>=0:raise ValueError('no certified quantum process resource')
    return value


def discover_npt_vector(model):
    """Numerical vector discovery followed by independent exact verification."""
    from fractions import Fraction as F
    witness,endpoint=quantum_instrument(model);w=F(*witness['memory_weight']);best=None
    for d,process in enumerate(model['processes']):
        wr,wi,sr,si,D=process_choi(process)
        qr=(w.denominator-w.numerator)*wr+w.numerator*sr;qi=(w.denominator-w.numerator)*wi+w.numerator*si
        Q=(np.asarray(qr,float)+1j*np.asarray(qi,float))/float(4*D*w.denominator)
        Q=Q.reshape([2]*8).transpose(4,5,2,3,0,1,6,7).reshape(16,16)
        values,vectors=np.linalg.eigh(Q);vec=vectors[:,0]
        certificate=dict(memory_weight=witness['memory_weight'],delay_index=d,
                         vector_real=np.rint(vec.real*10**10).astype(np.int64).tolist(),
                         vector_imag=np.rint(vec.imag*10**10).astype(np.int64).tolist())
        try:value=npt_expectation(model,certificate)
        except ValueError:continue
        if best is None or value<best[0]:best=(value,certificate)
    if best is None:raise ValueError('No NPT certificate found; this is not a classicality proof')
    return best[1]


def isolated_probe(model,certificate):
    from fractions import Fraction as F
    import math
    w=F(*certificate['memory_weight']);best=None
    for setting,group in zip(model['settings'],model['instruments']):
        matrices,D=group_matrices(group);transfers=[transfer(m) for m in matrices]
        for a,rho in physical.STATES.items():
            state=np.array([1,*[int(round(np.trace(rho@s).real)) for s in physical.PAULIS]],dtype=object)
            for z in range(3):
                differences=[]
                for b,((re,_),T) in enumerate(zip(matrices,transfers)):
                    q=F(int(np.trace(re)),2*D)
                    for sign in (1,-1):
                        effect=np.array([1,*[sign*int(i==z) for i in range(3)]],dtype=object)
                        p=F(int(effect@T@state),4*D);identity=F(1+sign*int(state[z+1]),2)
                        differences.append(p-(p-w*q*identity)/(1-w))
                tv=sum(abs(d) for d in differences)/2
                if best is None or tv>best[0]:best=(tv,a+','+setting+','+'xyz'[z],differences)
    tv,label,differences=best
    return dict(label=label,event=[physical.audit.OUTCOMES[i] for i,d in enumerate(differences) if d>0],
                total_variation=float(tv),total_variation_exact=[tv.numerator,tv.denominator],
                half_diamond_upper=float(w),independent_measured_calibration_available=False,
                ideal_external_probes_assumed=True,
                shots_per_candidate_for_two_95_percent_hoeffding_radii=math.ceil(2*math.log(80)/float(tv)**2)+1)
