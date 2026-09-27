#!/usr/bin/env python3
"""Exact archive-compatible pair-control evasion and third-order control certificate."""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations,product
from math import log
from pathlib import Path
import argparse,json
import full_models as fm
import robustness as rb
from certificate import check_endpoint
HERE=Path(__file__).resolve().parent


def pair(T,d):
    q=[[F(fm.norm(r[j]),d*d) for r in T[108:168]] for j in [0,1]]
    return sum(q[0])*sum(q[1])-sum(a*b for a,b in zip(*q))


def independent_bunch(T,d):
    """Actual three-distinct-cell event; rows 8..19 of the reference array."""
    total=F(0)
    for y in range(8,20):
        a,b,c=[[F(fm.norm(r[j]),d*d) for r in T[12*y:12*(y+1)]] for j in range(3)]
        A,B,C=map(sum,[a,b,c])
        total+=A*B*C-C*sum(x*z for x,z in zip(a,b))-B*sum(x*z for x,z in zip(a,c))-A*sum(x*z for x,z in zip(b,c))+2*sum(x*z*t for x,z,t in zip(a,b,c))
    return total


def gaussian_bound(U,gamma,kappa):
    """kappa = time * spectral bandwidth * Gaussian scale SD / 2 (hbar=1)."""
    return U+12*gamma*kappa*kappa+20*kappa**3


def selected_mixture(read,models,weights):
    assert len(models)==len(weights) and sum(weights)==1 and all(w>=0 for w in weights)
    Ts=[fm.transport(m)[:2] for m in models]
    q=[[sum(w*F(m['base_amplitude_numerators'][y][x]**2,m['amplitude_denominator']**2) for w,m in zip(weights,models)) for x in range(12)] for y in range(28)]
    rb.single_probabilities(q,read('cell-certificate.json'))
    counts=read('counts.json');pats=fm.histogram(read('parity-histogram.json'),counts);reg=read('robust-region.json')
    iv={int(k):tuple(map(F,v)) for k,v in reg['pattern_intervals'].items()}
    comps=[(w,T[96:240],d,'cluster') for w,(T,d) in zip(weights,Ts)]
    result=rb.pattern_check(comps,pats,iv,12,F(reg['bunch_lower']))
    hits=rb.row_hits(comps,pats,reg);assert all(h['compatible'] for h in hits)
    means=[[sum(w*rb.column_probabilities(T,d,j)[i] for w,(T,d) in zip(weights,Ts)) for i in range(337)] for j in range(3)]
    drift=max(sum(rb.tv(rb.column_probabilities(T,d,j),means[j]) for j in range(3)) for T,d in Ts)
    pair_pr=sum(w*pair(T,d) for w,(T,d) in zip(weights,Ts))
    third=sum(w*independent_bunch(T,d) for w,(T,d) in zip(weights,Ts))
    result.update(mean_singleton_region_pass=True,row_hit_constraints_pass=True,drift_budget_sum_tv=float(drift),pair_control_exact=str(pair_pr),labelled_triple_control=float(third))
    # Descriptive fitted-point statistic only: neither a p-value nor a calibrated fit test.
    ps=[fm.ParityModel(T,d,kind) for w,T,d,kind in comps];N=sum(pats.values())
    dev=2*sum(n*log(n/(N*float(sum(w*F(p.parity(pat),p.den) for w,p in zip(weights,ps))))) for pat,n in pats.items())
    result['many_pattern_deviance_descriptive']=dev
    return result


def pair_moment_obstruction():
    """Loss-dilated deterministic isometries: complete 1/2-body laws agree, 3-body not."""
    laws=[[z for z in product([0,1],repeat=3) if sum(z)%2==parity] for parity in [0,1]]
    for law in laws:
        for z in law:
            # Three detected cells in one row, and three distinct unobserved loss modes.
            T=[[int(o==j+3*z[j]) for j in range(3)] for o in range(6)]
            assert all(sum(T[o][j]*T[o][k] for o in range(6))==int(j==k) for j in range(3) for k in range(3))
    def distribution(law,inputs):
        return Counter(tuple(j for j in inputs if not z[j]) for z in law)
    for n in [1,2]:
        for inputs in combinations(range(3),n):assert distribution(laws[0],inputs)==distribution(laws[1],inputs)
    bunch=[F(distribution(law,range(3))[(0,1,2)],4) for law in laws]
    assert bunch==[F(1,4),F(0)]
    return {'all_singleton_and_pair_distributions_equal':True,'triple_bunch_probabilities':list(map(str,bunch)),
            'scope':'Structural counterexample, not a witness for the archive confidence region.'}


def verify(read):
    rb.selected_verify(read) # Recheck all retained confidence endpoints and the joint quantum model.
    Tq,dq,_=fm.transport(read('full-boson.json'));target=pair(Tq,dq)
    above=[read('pair-above-'+s+'.json') for s in ['minus','plus']]
    below=[read('pair-below-'+s+'.json') for s in ['minus','plus']]
    pa=sum(pair(*fm.transport(m)[:2]) for m in above)/2;pb=sum(pair(*fm.transport(m)[:2]) for m in below)/2
    assert pb<target<pa
    w=(target-pb)/(pa-pb);weights=[w/2,w/2,(1-w)/2,(1-w)/2]
    result=selected_mixture(read,above+below,weights)
    assert F(result['pair_control_exact'])==target
    result['weights_exact']=list(map(str,weights));result['quantum_pair_control_exact']=str(target)
    result['pair_control_separation_gap_exact']='0'
    # Both endpoints are feasible, hence the whole interval between their predictions is feasible.
    endpoints=[selected_mixture(read,models,[F(1,2)]*2) for models in [below,above]]
    result['certified_pair_prediction_inner_interval']=[float(pb),float(pa)]
    result['endpoint_region_checks']=[r['all_pattern_constraints_pass'] and r['row_hit_constraints_pass'] for r in endpoints]
    L=F(read('robust-region.json')['bunch_lower']);qd=independent_bunch(Tq,dq)
    # A future fixed-size control, not measured counts. Conservative rational CP endpoint.
    upper=F('0.0115');check_endpoint(3000,19,upper,F(1,320),'upper');assert 2*upper<L
    assert 2*qd<L
    q=[[F(v*v,read('full-boson.json')['amplitude_denominator']**2) for v in row] for row in read('full-boson.json')['base_amplitude_numerators']]
    rb.single_probabilities(q,read('cell-certificate.json'))
    assert qd<upper
    c=read('cell-certificate.json');gamma=max(F(p[1]) for p in c['row_intervals']);kappa=F('0.00925')
    bound=2*gaussian_bound(F(c['cell_upper']),gamma,kappa);assert bound<L
    return {'pair_evading_cluster':result,'pair_moment_obstruction':pair_moment_obstruction(),
        'labelled_triple_control':{'quantum_example_probability':float(qd),'exclusion_threshold':float(L/2),
          'universal_cluster_factor':2,'prospective_trials':3000,'prospective_max_successes':19,
          'certified_control_upper':str(upper),'one_sided_tail':'1/320','joint_two_event_error_upper':'1/160',
          'status':'Prospective conditional protocol; no new acquisition has occurred.'},
        'prospective_source_deletions':{
          'control_removed':'The exact pair-evading C2 mixture passes the complete existing selected region.',
          'original_many_body_source_removed':'The full-boson.json channel with cluster labels passes singleton calibration and the prospective control upper bound.',
          'joint_quantum_feasible':True,'singleton_source_needed_for_new_factor_two_test':False,
          'scope':'Conditional on the declared future upper region [0,0.0115]; not an observed two-source result.'},
        'gaussian_scale_envelope':{'excluded_through_kappa':str(kappa),'cluster_upper':float(bound),
          'scope':'Conditional bandwidth/noise product bound, not a measured apparatus noise limit.'}}


def main():
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args()
    read=lambda n:json.loads((HERE/'results'/n).read_text());r=verify(read)
    # Round descriptive floating logarithms only; every acceptance comparison above is exact.
    r['pair_evading_cluster']['many_pattern_deviance_descriptive']=round(r['pair_evading_cluster']['many_pattern_deviance_descriptive'],6)
    text=json.dumps(r,sort_keys=True,indent=2)+'\n';out=HERE/'results/control-design.json'
    if a.check:assert text==out.read_text()
    else:out.write_text(text)
    print(text)
if __name__=='__main__':main()
