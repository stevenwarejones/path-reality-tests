#!/usr/bin/env python3
"""Prospective finite-pointer study. Synthetic probabilities, no apparatus data.

All statistical guarantees are conditional on the ontological representation
premises and characterized calibration controls stated in protocol.md.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import json
import math
from pathlib import Path
import numpy as np
from scipy.stats import beta as beta_dist

ROOT = Path(__file__).resolve().parent
SEED = 260926
AUDIT_CONTEXTS = 44
AUDIT_COORDINATES = 88


def _unit(name, value):
    if not math.isfinite(float(value)) or not 0 <= value <= 1:
        raise ValueError(f"{name} must be in [0,1]")


def rational_instrument(t=F(1, 3), u=F(1, 3), v=F(1, 3)):
    """Rational unit-circle coordinates; chosen reference t=1/3 gives (a,b)=(3/5,4/5).
    Source (x,y), success (z,-w) with independently tunable rational coordinates.
    """
    for name, value in [('t', t), ('u', u), ('v', v)]:
        _unit(name, value)
    a, b = 2*t/(1+t*t), (1-t*t)/(1+t*t)
    x, y = (1-u*u)/(1+u*u), 2*u/(1+u*u)
    z, w = (1-v*v)/(1+v*v), 2*v/(1+v*v)
    if a > b:
        raise ValueError("requires a <= b (negative outcome biased toward Q)")
    table = [[(b*x*z-a*y*w)**2, (b*x*w+a*y*z)**2],
             [(a*x*z-b*y*w)**2, (a*x*w+b*y*z)**2]]
    assert a*a+b*b == x*x+y*y == z*z+w*w == 1
    assert sum(map(sum, table)) == 1
    q, d, f = b*b, (a-b)**2/2, (x*z-y*w)**2
    return {'a': a, 'b': b, 'x': x, 'y': y, 'z': z, 'w': w,
            'table': table, 'q': q, 'd': d, 'f': f,
            'gap': table[0][0]-q*f-d*(1-f)}


def lossy_table(model, probe_eff=1., final_eff=1., flip=0.):
    """3x3 table: pointer (-,+,lost), final (success,failure,no detection).
    Independent classical erasures after a complete instrument and a fixed
    symmetric final-bit flip. This is a specified noise channel, not fair sampling.
    """
    for key, val in [('probe_eff', probe_eff), ('final_eff', final_eff), ('flip', flip)]:
        _unit(key, val)
    raw = np.asarray(model['table'], float)
    noise = np.array([[1-flip, flip], [flip, 1-flip]])
    final = raw @ noise
    p = np.zeros((3, 3))
    p[:2, :2] = probe_eff*final_eff*final
    p[:2, 2] = probe_eff*(1-final_eff)*raw.sum(axis=1)
    p[2, :2] = (1-probe_eff)*final_eff*final.sum(axis=0)
    p[2, 2] = (1-probe_eff)*(1-final_eff)
    assert np.all(p >= -1e-15) and abs(p.sum()-1) < 1e-12
    f = final_eff*((1-2*flip)*float(model['f'])+flip)
    q = probe_eff*float(model['q'])
    d = float(model['d'])
    return p, (float(p[0, 0]), f, q, d)


def witness(a, f, q, d):
    return a-min(q, q*f+d*(1-f))


def cp(k, n, error):
    if not isinstance(k, (int, np.integer)) or not isinstance(n, (int, np.integer)):
        raise ValueError("counts must be integers")
    if n < 0 or k < 0 or k > n or not 0 < error < 1:
        raise ValueError("invalid count or error probability")
    if n == 0:
        return (0., 1.)
    lo = 0. if k == 0 else float(beta_dist.ppf(error/2, k, n-k+1))
    hi = 1. if k == n else float(beta_dist.ppf(1-error/2, k+1, n-k))
    return lo, hi


def decision(counts, alpha=0.01, allowance=0.):
    """4 separate eligible-trial Bernoulli counts (a,f,q,d); Bonferroni CP.
    No independence between confidence intervals is needed for the union bound.
    IID Bernoulli trials conditional on context counts are needed for each CP.
    """
    if len(counts) != 4 or not 0 < alpha < 1 or allowance < 0:
        raise ValueError("four counts, alpha in (0,1), and nonnegative allowance required")
    intervals = [cp(k, n, alpha/4) for k, n in counts]
    return interval_decision(intervals, allowance)


def interval_decision(intervals, allowance=0.):
    (al, _), (fl, fu), (_, qu), (_, du) = intervals
    ceiling = min(qu, max(qu*fl+du*(1-fl), qu*fu+du*(1-fu)))
    lower_gap = al-ceiling-allowance
    return {'reject': lower_gap > 0, 'lower_gap': lower_gap, 'intervals': intervals}


def hoeffding_decision(counts, alpha=.005, allowance=0.):
    if len(counts) != 4 or not 0 < alpha < 1 or allowance < 0:
        raise ValueError("invalid decision parameters")
    intervals=[]
    for k,n in counts:
        cp(k,n,alpha/4)  # validate integer trial counts and domain
        r=math.sqrt(math.log(8/alpha)/(2*n)) if n else 1.
        p=k/n if n else 0.
        intervals.append((max(0.,p-r),min(1.,p+r)))
    return interval_decision(intervals,allowance)


def certified_budget(gap, alpha=.005, beta=.1, tolerance=.01, cal_alpha=.005):
    """Conservative fixed-count IID Hoeffding power certificate, including audit cost.
    With four simultaneous estimates, |error_i| <= r_beta with probability 1-beta.
    Hoeffding decision intervals add r_alpha. The witness is 1-Lipschitz in each
    of its four coordinates on [0,1]^4. Thus gap > 4*(r_alpha+r_beta) suffices.
    Audit intervals cover 88 coordinates; passing tolerance cannot prove exact
    operational equivalence or bound ontic variation. Their trial cost is counted.
    """
    if gap <= 0 or not 0 < tolerance < 1:
        raise ValueError("positive gap and audit tolerance required")
    for x in [alpha, beta, cal_alpha]:
        if not 0 < x < 1:
            raise ValueError("risk budgets must be in (0,1)")
    c = math.sqrt(math.log(8/alpha)/2)+math.sqrt(math.log(8/beta)/2)
    n = math.floor((4*c/gap)**2)+1
    ncal = math.ceil(math.log(2*AUDIT_COORDINATES/cal_alpha)/(2*tolerance*tolerance))
    return {'per_main_context': n, 'main_contexts': 4, 'per_audit_context': ncal,
            'audit_contexts': AUDIT_CONTEXTS, 'eligible_trials': 4*n+AUDIT_CONTEXTS*ncal,
            'alpha_main': alpha, 'alpha_calibration': cal_alpha, 'beta_main': beta,
            'audit_radius': tolerance,
            'guarantee': 'conditional Hoeffding power; calibration equivalences remain premises'}


def random_context_budget(gap, contexts=48, alpha=.005, beta=.1):
    """Fixed total N, independent uniform context selection; random counts retained.
    Chernoff: Pr[min n_s < N/(2K)] <= K exp(-N/(8K)). Allocate beta/2
    to counts and beta/2 to outcome estimation. No stopping or dropped trials.
    """
    if gap <= 0 or not isinstance(contexts, int) or contexts < 4+AUDIT_CONTEXTS:
        raise ValueError("positive gap and at least 48 integral contexts required")
    budget = certified_budget(gap, alpha=alpha, beta=beta/2)
    n = max(budget['per_main_context'], budget['per_audit_context'])
    total = max(2*contexts*n, math.ceil(8*contexts*math.log(2*contexts/beta)))
    return {'eligible_trials': total, 'contexts': contexts,
            'minimum_count_event': total/(2*contexts), 'log_counts_failure_bound':
            math.log(contexts)-total/(8*contexts), 'beta_total': beta}


def exact_countermodels():
    """Independent exact verification of complete finite kernels and their moments."""
    def evaluate(mu, kernel, response):
        n=len(mu)
        assert len(kernel)==len(response)==n
        assert all(0<=r<=1 for r in response)
        assert all(len(row)==2 and all(len(branch)==n for branch in row) for row in kernel)
        assert sum(mu) == 1 and min(mu) >= 0
        for row in kernel:
            assert sum(map(sum, row)) == 1 and min(sum(row, [])) >= 0
        table = [[sum(mu[l]*kernel[l][m][j]*(response[j] if f == 0 else 1-response[j])
                      for l in range(n) for j in range(n)) for f in range(2)] for m in range(2)]
        assert sum(map(sum, table)) == 1
        return table, sum(mu[l]*response[l] for l in range(n))
    h=F(1,2)
    f0=F(49,625); q0=F(16,25); d0=F(1,50)
    identity_probe=[[[h,0],[h,0]], [[0,h],[0,h]]]
    null_table,null_f=evaluate([1,0],identity_probe,[f0,f0])
    assert null_f==f0 and null_table[0][0]==F(49,1250)
    assert h<=q0 and (1-d0)+d0==1 and null_table[0][0]<=q0*f0+d0*(1-f0)
    # Ordinary disturbance: uniform preparation, always reset to the pointer bit.
    # Every final marginal is unchanged for this preparation, but outcome branches change.
    k=[[[h,0],[0,h]], [[h,0],[0,h]]]
    table,f=evaluate([h,h], k, [1,0])
    assert table == [[h,0],[0,h]] and f==h
    # Stable identity channel but context-dependent pointer correlated with final property.
    k2=[[[1,0],[0,0]], [[0,0],[0,1]]]
    t2,f2=evaluate([h,h],k2,[1,0])
    assert t2==table and f2==f
    # Postselection inflation under a legitimate q=1/2,d=1/100 null.
    d=F(1,100); q=h
    k3=[[[q-d,d],[1-q,0]], [[0,q],[0,1-q]]]
    t3,f3=evaluate([1,0],k3,[0,1])
    assert t3[0][0]==d and sum(t3[m][0] for m in range(2))==d and f3==0
    # Exact reference-data rivals: drop exactly one representation premise.
    target=rational_instrument()['table']
    cap_kernel=[[[F(1081,1225),0],[F(144,1225),0]],
                [[d0,F(49,100)],[0,F(49,100)]]]
    cap_table,cap_f=evaluate([f0,1-f0],cap_kernel,[1,0])
    assert cap_table==target and cap_f==f0
    for l in range(2):
        for j in range(2):
            assert sum(cap_kernel[l][m][j] for m in range(2)) == (1-d0)*(l==j)+d0*(j==0)
    assert sum(cap_kernel[0][0])>q0
    disturbance_kernel=[[[0,target[0][0],target[0][1]],
                         [0,target[1][0],target[1][1]]],
                        [[0,h,0],[0,h,0]], [[0,0,h],[0,0,h]]]
    disturbance_table,disturbance_f=evaluate([1,0,0],disturbance_kernel,[f0,1,0])
    assert disturbance_table==target and disturbance_f==f0
    max_negative=max(sum(row[0]) for row in disturbance_kernel)
    assert max_negative==F(337,625) and max_negative<=q0
    # L has zero probability of staying L, contradicting the diagonal lower bound.
    assert sum(disturbance_kernel[0][m][0] for m in range(2)) < 1-d0
    return {'drop_cap_quantum_joint': [[str(x) for x in row] for row in cap_table],
            'drop_cap_bypass': str(cap_f),
            'drop_cap_negative_response': str(sum(cap_kernel[0][0])),
            'drop_disturbance_quantum_joint': [[str(x) for x in row] for row in disturbance_table],
            'drop_disturbance_bypass': str(disturbance_f),
            'drop_disturbance_max_negative_response': str(max_negative),
            'reference_parameter_null_joint':str(null_table[0][0]),
            'reference_parameter_null_bypass':str(null_f),
            'undisturbed_marginal_invasive_joint': [[str(x) for x in row] for row in table],
            'marginal_disturbance': '0', 'invasive_gap_if_d_misidentified_as_zero': '1/4',
            'contextual_pointer_same_joint': [[str(x) for x in row] for row in t2],
            'legitimate_null_postselected_negative_fraction': '1',
            'legitimate_null_joint_negative_success': str(d),
            'legitimate_null_bypass': str(f3)}


def matrix_check(model):
    """Independent complex-matrix Born/channel checks, including imaginary coherence."""
    a,b,x,y,z,w=[float(model[k]) for k in ['a','b','x','y','z','w']]
    ks=[np.diag([b,a]).astype(complex), np.diag([a,b]).astype(complex)]
    psi=np.array([x,y],complex); rho=np.outer(psi,psi.conj())
    fs=[np.array([z,-w],complex),np.array([w,z],complex)]
    result=np.array([[np.vdot(v,k@rho@k.conj().T@v).real for v in fs] for k in ks])
    np.testing.assert_allclose(result, np.array(model['table'],float),atol=1e-14,rtol=0)
    zz=np.diag([1,-1])
    # Basis of ALL complex 2x2 matrices, not just real density matrices.
    for i in range(2):
        for j in range(2):
            X=np.zeros((2,2),complex); X[i,j]=1+2j
            np.testing.assert_allclose(sum(k@X@k.conj().T for k in ks),
                (1-float(model['d']))*X+float(model['d'])*zz@X@zz,atol=1e-14,rtol=0)
    return float(np.max(np.abs(result-np.array(model['table'],float))))


def calibration_tables(model, probe_eff=.95):
    """All 44 ideal audit contexts with complete outcomes, not fitted data.
    Final calibration readout is characterized/complete; prospective main loss
    is separate. No finite statistical fit is asserted to prove these identities.
    """
    _unit('probe_eff',probe_eff)
    a,b=float(model['a']),float(model['b'])
    q,d=float(model['q']),float(model['d']); pm=2*q-1
    I=np.eye(2,dtype=complex); Z=np.diag([1,-1]).astype(complex)
    X=np.array([[0,1],[1,0]],complex); Y=np.array([[0,-1j],[1j,0]],complex)
    states={'Q':np.array([1,0],complex),'P':np.array([0,1],complex),
            'X+':np.array([1,1],complex)/np.sqrt(2),
            'Y+':np.array([1,1j],complex)/np.sqrt(2)}
    ks=[np.diag([b,a]),np.diag([a,b])]
    tables={}; residuals=[]
    for name,psi in states.items():
        rho=np.outer(psi,psi.conj())
        pointer=[float(np.trace(k@rho@k.conj().T).real) for k in ks]
        simulator=[pm*float(rho[0,0].real)+(1-pm)/2,
                   pm*float(rho[1,1].real)+(1-pm)/2]
        for label,p in [('probe',pointer),('simulator',simulator)]:
            tables[f'pointer/{name}/{label}']=[probe_eff*p[0],probe_eff*p[1],1-probe_eff]
        residuals.extend(np.array(pointer)-np.array(simulator))
        outputs={'M':sum(k@rho@k.conj().T for k in ks),'I':rho,'Z':Z@rho@Z}
        for axis,A in [('X',X),('Y',Y),('Z',Z)]:
            for op,out in outputs.items():
                plus=float(np.trace((I+A)@out/2).real)
                tables[f'channel/{name}/{axis}/{op}']=[plus,1-plus,0.]
            residuals.extend(np.array(tables[f'channel/{name}/{axis}/M'])-
                ((1-d)*np.array(tables[f'channel/{name}/{axis}/I'])+
                 d*np.array(tables[f'channel/{name}/{axis}/Z'])))
    assert len(tables)==AUDIT_CONTEXTS
    for row in tables.values():
        assert abs(sum(row)-1)<1e-12 and min(row)>-1e-14
    # Exact informational completeness in (trace,X,Y,Z) coordinates.
    from itertools import permutations
    rows=[[1,0,0,1],[1,0,0,-1],[1,1,0,0],[1,0,1,0]]
    determinant=sum((-1)**sum(p[i]>p[j] for i in range(4) for j in range(i+1,4))*
                    math.prod(rows[i][p[i]] for i in range(4)) for p in permutations(range(4)))
    assert abs(determinant)==2
    return {'contexts':tables,'max_equivalence_residual':float(max(abs(x) for x in residuals)),
            'preparation_coordinate_determinant':determinant,
            'status':'ideal synthetic controls; not a finite-data equivalence certificate'}


def simulate_power(params, n, repetitions=2000, seed=SEED):
    """CP rule Monte Carlo; simulation uncertainty is separate from analytic guarantee."""
    rng=np.random.default_rng(seed)
    ks=rng.binomial(n,params,size=(repetitions,4))
    rejects=sum(decision([(int(k),n) for k in row],alpha=.005)['reject'] for row in ks)
    lo,hi=cp(rejects,repetitions,.01)
    return {'seed':seed,'repetitions':repetitions,'rejections':rejects,
            'power':rejects/repetitions,'mc_99pct_interval':[lo,hi]}


def report():
    ref=rational_instrument(); matrix_check(ref)
    noise_scenarios=[]
    for eta,beta,flip in [(1.,1.,0.),(.9,.95,0.),(.8,.9,.01),(.5,.8,.02)]:
        table, params=lossy_table(ref,beta,eta,flip)
        gap=witness(*params)
        noise_scenarios.append({'final_eff':eta,'probe_eff':beta,'final_flip':flip,
            'table':table.tolist(),'a_f_q_d':list(params),'gap':gap,
            'budget':certified_budget(gap) if gap>0 else None})
    # Finite rational grid search, explicitly not a global optimum.
    grid=[]
    for t in [F(i,40) for i in range(2,17)]:
        for u in [F(i,20) for i in range(2,11)]:
            for v in [F(i,20) for i in range(2,11)]:
                model=rational_instrument(t,u,v)
                _,p=lossy_table(model,.95,.9,.01)
                gap=witness(*p)
                if gap>0:
                    grid.append((certified_budget(gap)['eligible_trials'],t,u,v,gap))
    best=min(grid)
    budget=noise_scenarios[0]['budget']
    _,params=lossy_table(ref)
    return {'schema':1,'status':'prospective synthetic; no apparatus evidence',
        'reference':{k:([[str(x) for x in r] for r in v] if k=='table' else str(v))
                     for k,v in ref.items()}, 'noise_scenarios':noise_scenarios,
        'grid_search':{'total_candidates':15*9*9,'positive_candidates':len(grid),
            'best_parameters':[str(x) for x in best[1:4]],'gap':best[4],
            'eligible_trials':best[0],'claim':'minimum over this declared grid only'},
        'reference_random_allocation':random_context_budget(float(ref['gap'])),
        'cp_power':simulate_power(params, budget['per_main_context']),
        'countermodels':exact_countermodels(), 'calibration':calibration_tables(ref)}


def compare(expected,actual,path='$'):
    """Exact keys/discrete values; explicit absolute/relative float roundoff tolerance."""
    if type(expected) is not type(actual): raise AssertionError(f'{path}: type mismatch')
    if isinstance(expected,dict):
        if expected.keys()!=actual.keys(): raise AssertionError(f'{path}: keys mismatch')
        for k in expected: compare(expected[k],actual[k],f'{path}.{k}')
    elif isinstance(expected,list):
        if len(expected)!=len(actual): raise AssertionError(f'{path}: length mismatch')
        for i,(x,y) in enumerate(zip(expected,actual)): compare(x,y,f'{path}[{i}]')
    elif isinstance(expected,float):
        if not math.isclose(expected,actual,rel_tol=1e-11,abs_tol=1e-13):
            raise AssertionError(f'{path}: {expected} != {actual}')
    elif expected!=actual: raise AssertionError(f'{path}: {expected} != {actual}')


if __name__=='__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('--check',action='store_true')
    args=parser.parse_args(); result=report(); target=ROOT/'results.json'
    if args.check: compare(json.loads(target.read_text()),result); print('Path contextuality snapshot verified')
    else: target.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
