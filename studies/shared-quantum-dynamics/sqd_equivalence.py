"""Exact observable-coordinate fiber for uniform incoherent leakage/return.

No optimizer is used to derive the compatibility interval. Values applied to
fitted invariants remain conditional examples, not data confidence intervals.
"""
import numpy as np
from sqd_models import operations,parameters,reference_probability


def invariants(x):
    lam,ret,z=x[15]*1e-4,x[16]*1e-4,x[17]
    s=lam+ret
    if not 0<s<1:raise ValueError('nondegenerate 0 < loss+return < 1 required')
    a=(x[12]+x[13])/2
    B=(a-z)*lam/s
    return {'a':a,'A':a-B,'B':B,'C':(x[12]-x[13])*x[14]/2,
            's':s,'gamma':(np.asarray(x[9:12])*1e-4-np.log1p(-lam)).tolist()}


def fiber_interval(x,declared_box=True):
    v=invariants(x);a,B,s=v['a'],v['B'],v['s']
    lower=B*s/a if B>=0 else -B*s/(1-a)
    upper=min(s,-np.expm1(-min(v['gamma'])))
    if declared_box:
        # Fit family also bounds return <= .001 and damping exponents <= .001.
        lower=max(lower,s-.001,max(0,1-np.exp(.001-max(v['gamma']))))
        upper=min(upper,.001)
    return max(0,float(lower)),float(upper)


def fiber_point(x,lam):
    v=invariants(x);lo,hi=fiber_interval(x)
    if not lo-1e-13<=lam<=hi+1e-13 or lam<=0:raise ValueError('outside nonzero physical fiber')
    out=np.array(x,float)
    out[9:12]=(np.asarray(v['gamma'])+np.log1p(-lam))/1e-4
    out[15]=lam/1e-4;out[16]=(v['s']-lam)/1e-4
    out[17]=v['a']-v['B']*v['s']/lam
    _,low,high=parameters('leakage')
    if np.any(out<low-1e-9) or np.any(out>high+1e-9):raise ValueError('nonphysical transformed point')
    return np.clip(out,low,high)


def formula_probability(word,x):
    v=invariants(x);_,_,_,rots=operations(x,'leakage');r=np.array([0.,0.,1.])
    for g in word:r=rots[0,g]@r
    return float(v['A']+v['B']*(1-v['s'])**len(word)+v['C']*np.exp(-sum(v['gamma'][g] for g in word))*r[2])


def leakage_population(length,x):
    v=invariants(x)
    return float((x[15]*1e-4/v['s'])*(-np.expm1(length*np.log1p(-v['s']))))


def monitored_loss(population,length,s):
    if length<1 or not 0<s<1:raise ValueError('degenerate population probe')
    return float(population*s/(-np.expm1(length*np.log1p(-s))))


def clock_probability(word,x,wait_after=None):
    """External gate-slot clock. A wait advances clock but is identity on qubit.

    No recorded timing is inferred. Agreement without waits is algebraic under
    the declared one-slot-per-primitive schedule, not a claim about hardware clocks.
    """
    gates,initial,effect,_=operations(x,'alternating');total=0.
    for branch in range(2):
        state=initial.copy();tick=branch
        for t,g in enumerate(word):
            state=gates[tick%2,g]@state;tick+=1
            if wait_after==t:tick+=1
        total+=effect@state/2
    return float(total)


def synthetic_certificate():
    # Interior, intentionally larger than the measured near-null examples.
    x,_,_=parameters('leakage');x[9:12]=[4.,5.,6.];x[15:]=[1.,2.,.85]
    lo,hi=fiber_interval(x);a=fiber_point(x,lo+.1*(hi-lo));b=fiber_point(x,lo+.9*(hi-lo))
    words=[(),(1,),(1,2),(2,1),(1,1,2,0)*19,(1,2,0)*333]
    maxerr=max(abs(reference_probability(w,a,'leakage')-reference_probability(w,b,'leakage')) for w in words)
    y,_,_=parameters('alternating');y[15]=3.
    w=(1,1)
    no_wait=clock_probability(w,y);wait=clock_probability(w,y,0)
    return {'leakage_parameters_a':a.tolist(),'leakage_parameters_b':b.tolist(),
        'sharp_conditional_loss_interval':[lo,hi], 'terminal_probability_max_error':maxerr,
        'one_gate_monitor_a':leakage_population(1,a),'one_gate_monitor_b':leakage_population(1,b),
        'monitor_loss_reconstruction_error':abs(monitored_loss(leakage_population(7,a),7,invariants(a)['s'])-a[15]*1e-4),
        'clock_pair':{'parameters':y.tolist(),'word':[1,1],'memory_probability':reference_probability(w,y,'alternating'),
                      'clock_no_wait':no_wait,'clock_identity_wait_after_first_gate':wait,'separation':abs(no_wait-wait)},
        'status':'analytic identities; floating-point checks are witnesses, not proof substitutes'}
