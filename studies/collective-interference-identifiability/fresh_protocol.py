#!/usr/bin/env python3
"""Exact physical alternatives and mismatch coefficients for a fresh two-arm protocol."""
from fractions import Fraction as F
from pathlib import Path
import argparse,copy,json
import full_models as fm
import robustness as rb
from control_design import independent_bunch
from certificate import check_endpoint
HERE=Path(__file__).resolve().parent


def coefficients(overlap_amplitude,transport_tv,original_error,control_error):
    """Uniform product-label overlap amplitude bound, plus separately calibrated errors.

transport_tv bounds the change in the fully distinguishable third-order event,
not merely TV between averaged singleton distributions. Errors are event-level
trace-distance/readout/contamination allowances for each preparation arm.
"""
    a,delta,eO,eC=map(F,[overlap_amplitude,transport_tv,original_error,control_error])
    assert 0<=a<1 and all(0<=v<=1 for v in [delta,eO,eC])
    row_sum=3*a*a+2*a**3;lam=1-row_sum;assert lam>0
    return F(2)/lam,2*delta+eO+2*eC/lam,lam


def scaled_model(model,scale):
    m=copy.deepcopy(model);m['amplitude_denominator']*=scale.denominator
    m['base_amplitude_numerators']=[[v*scale.numerator for v in row] for row in m['base_amplitude_numerators']]
    return m


def axial_contact_ratios():
    """Diagonal harmonic-oscillator contact integrals, relative to I_00.

    These are not bounds on off-diagonal motional transitions or real-apparatus errors.
    Polynomials are densities relative to the n=0 Gaussian, in powers of x^2.
    """
    densities=[[F(1)],[F(0),F(2)],[F(1,2),F(-2),F(2)]]
    moments=[F(1)]
    for k in range(1,5):moments.append(moments[-1]*F(2*k-1,4))
    ratios=[[sum(a*b*moments[i+j] for i,a in enumerate(p) for j,b in enumerate(q)) for q in densities] for p in densities]
    assert ratios==[[F(1),F(1,2),F(3,8)],[F(1,2),F(3,4),F(7,16)],[F(3,8),F(7,16),F(41,64)]]
    return [[str(v) for v in row] for row in ratios]


def alternatives(read):
    cert=read('cell-certificate.json');pats=fm.histogram(read('parity-histogram.json'),read('counts.json'));reg=read('robust-region.json')
    iv={int(k):tuple(map(F,v)) for k,v in reg['pattern_intervals'].items()};result=[]
    for numerator in [985,990,995,1000,1005,1010]:
        scale=F(numerator,1000);m=scaled_model(read('full-boson.json'),scale);T,d,_=fm.transport(m)
        q=[[F(v*v,m['amplitude_denominator']**2) for v in row] for row in m['base_amplitude_numerators']]
        rb.single_probabilities(q,cert);comps=[(F(1),T[96:240],d,'boson')]
        check=rb.pattern_check(comps,pats,iv,12,F(reg['bunch_lower']));assert all(h['compatible'] for h in rb.row_hits(comps,pats,reg))
        # Three surviving clicks scale with the sixth power of the amplitude multiplier.
        p=fm.ParityModel(T[96:240],d,'boson')
        from itertools import combinations
        original=F(sum(p.occupation(tuple(12*y+x for x in xs)) for y in range(12) for xs in combinations(range(12),3)),p.den)
        result.append(dict(amplitude_scale=str(scale),original_exact=str(original),control_exact=str(independent_bunch(T,d)),
                           original=float(original),control=float(independent_bunch(T,d)),full_selected_region_pass=True))
    return result


def verify(read):
    values=alternatives(read)
    # Example observed count pair only, never treated as acquired data.
    lo=F('0.0201');hi=F('0.0094')
    check_endpoint(5000,130,lo,F(1,320),'lower');check_endpoint(5000,23,hi,F(1,320),'upper');assert lo>2*hi
    A,b,lam=coefficients(F(1,10),F(1,1000),0,0)
    assert lam==F(121,125) and A==F(250,121) and b==F(1,500)
    return dict(alternatives=values,fresh_example=dict(original_trials=5000,original_events=130,control_trials=5000,control_events=23,
        lower=str(lo),upper=str(hi),one_sided_tail='1/320',ideal_rejects=True,status='Hypothetical count pair, not acquired data.'),
        mismatch_example=dict(overlap_amplitude='1/10',transport_event_tv='1/1000',factor=str(A),offset=str(b),gram_lower=str(lam)),
        axial_contact_ratios=axial_contact_ratios(),
        power_scope='Finite certified bosonic examples and their convex hull, not all archive-compatible quantum laws.')


def main():
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();read=lambda n:json.loads((HERE/'results'/n).read_text())
    result=verify(read);text=json.dumps(result,sort_keys=True,indent=2)+'\n';out=HERE/'results/fresh-protocol.json'
    if a.check:assert out.read_text()==text
    else:out.write_text(text)
    print(text)
if __name__=='__main__':main()
