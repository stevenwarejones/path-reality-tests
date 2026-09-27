#!/usr/bin/env python3
"""Numerical fixed-size IID binomial power; exact acceptance examples are separate."""
from functools import lru_cache
from pathlib import Path
from fractions import Fraction as F
import argparse,json,math
import numpy as np
from scipy.stats import beta,binom
from fresh_protocol import coefficients
from dynamics import close
HERE=Path(__file__).resolve().parent
TAIL=1/320


@lru_cache(None)
def endpoints(n):
    k=np.arange(n+1);lo=np.zeros(n+1);hi=np.ones(n+1)
    lo[1:]=beta.ppf(TAIL,k[1:],n-k[1:]+1);hi[:-1]=beta.isf(TAIL,k[:-1]+1,n-k[:-1])
    return lo,hi


def critical_control_counts(n,m,factor=2.,offset=0.):
    lo=endpoints(n)[0];hi=endpoints(m)[1]
    return np.searchsorted(hi,(lo-offset)/factor,side='left')-1


def power(n,m,p,q,factor=2.,offset=0.):
    assert n>0 and m>0 and 0<=p<=1 and 0<=q<=1 and factor>0 and offset>=0
    critical=critical_control_counts(n,m,factor,offset)
    return float(np.dot(binom.pmf(np.arange(n+1),n,p),binom.cdf(critical,m,q)))


def zero_failure_cost(bound,tail):
    # Exact rational verification after logarithms propose a sample count.
    b,t=F(bound),F(tail);assert 0<b<1 and 0<t<1
    n=math.ceil(math.log(float(t))/math.log1p(-float(b)))
    while (1-b)**n>t:n+=1
    while n>0 and (1-b)**(n-1)<=t:n-=1
    return n


def evaluate():
    alts=json.loads((HERE/'results/fresh-protocol.json').read_text())['alternatives']
    nominal=next(r for r in alts if F(r['amplitude_scale'])==1);p,q=nominal['original'],nominal['control']
    costs=[dict(original_trials=n,control_trials=n,total_science_trials=2*n,power=power(n,n,p,q)) for n in [1500,3000,5000,7500,10000]]
    variation=[dict(amplitude_scale=r['amplitude_scale'],original_probability=r['original'],control_probability=r['control'],power_10000=power(5000,5000,r['original'],r['control'])) for r in alts]
    sensitivities=[]
    for a,delta,eO,eC in [('0','0','0','0'),('0.1','0','0','0'),('0.2','0','0','0'),('0','0.001','0','0'),('0','0.003','0','0'),('0.1','0.001','0','0'),('0','0','0.001','0.001')]:
        A,b,lam=coefficients(a,delta,eO,eC)
        # Gershgorin upper envelope also bounds control probability by (1+r)D.
        # Adverse response shifts are a conservative alternative-rate envelope.
        p_min=max(0,p-float(eO));q_max=min(1,(2-float(lam))*(q+float(delta))+float(eC))
        sensitivities.append(dict(overlap_amplitude=a,transport_event_tv=delta,original_error=eO,control_error=eC,
            factor=float(A),offset=float(b),fixed_rate_power_10000=power(5000,5000,p,q,float(A),float(b)),
            adverse_rate_original=p_min,adverse_rate_control=q_max,
            adverse_rate_power_10000=power(5000,5000,p_min,q_max,float(A),float(b))))
    pmin=min(r['original'] for r in alts);qmax=max(r['control'] for r in alts)
    # Monotonicity supplies a numerical conservative bound for the declared finite convex hull.
    hull=power(5000,5000,pmin,qmax)
    return dict(iid_fixed_size_rule='Reject if lower_CP(original) > factor * upper_CP(control) + offset; each tail is 1/320.',
        numerical_status='SciPy binomial integration, not an interval-certified power calculation.',
        ideal_nominal_costs=costs,certified_alternative_examples=variation,finite_hull_conservative_power_10000=hull,
        finite_hull_conservative_power_15000=power(7500,7500,pmin,qmax),
        mismatch_sensitivity=sensitivities,
        historical_3000_control_only_power=float(binom.cdf(19,3000,q)),
        calibration_cost_example=dict(per_assay_error_target='1/1000',per_assay_tail='1/1000',zero_failure_trials=zero_failure_cost('1/1000','1/1000'),
            number_of_assays=3,status='Cost illustration for validated binary error assays only; not a transport-TV or coherent-state certification.'),
        total_cost_formula='ceil(N_original/yield_original)+ceil(N_control/yield_control)+N_validation. Validation cost and usable yields are unestablished.',
        uniform_over_all_compatible_quantum_models=False)


def main():
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();result=evaluate();out=HERE/'results/fresh-power.json'
    if a.check:close(result,json.loads(out.read_text()))
    else:out.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
