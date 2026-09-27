#!/usr/bin/env python3
"""Exact rational / interval physicality and event-probability verification."""
from fractions import Fraction as F
from itertools import combinations,permutations
from math import isqrt
from pathlib import Path
import argparse,json
from certificate import check_endpoint

HERE=Path(__file__).resolve().parent
SCALE=10**12

def sqrt_interval(q):
    assert q>=0
    k=isqrt(q.numerator*SCALE*SCALE//q.denominator)
    a,b=F(k,SCALE),F(k+1,SCALE)
    assert a*a<=q<=b*b
    return a,b

def mul(z,w):return (z[0]*w[0]-z[1]*w[1], z[0]*w[1]+z[1]*w[0])

def phase(t):return ((1-t*t)/(1+t*t),2*t/(1+t*t))

def verify(w,counts,cert):
    n=w['single_shots'];assert n==counts['single_shots']
    p={}
    for y,x,c in w['nonzero_cell_counts']:
        assert 9<=y<=20 and 8<=x<=19 and isinstance(c,int) and c>0
        assert (y,x) not in p
        p[y,x]=F(c,n)
    assert [sum(p.get((y,x),F(0)) for x in range(28))*n for y in range(9,21)]==counts['reference_row_counts']
    ts=list(map(F,w['phase_tangents']));assert len(ts)==56
    phases=[[ (F(1),F(0)),phase(ts[2*y]),phase(ts[2*y+1])] for y in range(28)]
    assert all(a*a+b*b==1 for row in phases for a,b in row)
    shifts=[0,1,-1]
    def prob(y,x,j):return p.get((y-shifts[j],x),F(0))
    diag=[sum(prob(y,x,j) for y in range(28) for x in range(28)) for j in range(3)]
    off={}
    for j,k in combinations(range(3),2):
        real=imag=F(0);error=F(0)
        for y in range(28):
            z=mul((phases[y][j][0],-phases[y][j][1]),phases[y][k])
            for x in range(28):
                l,u=sqrt_interval(prob(y,x,j)*prob(y,x,k))
                mid=(l+u)/2
                real+=mid*z[0];imag+=mid*z[1]
                error+=(u-l)/2*(abs(z[0])+abs(z[1]))
        off[j,k]=abs(real)+abs(imag)+error
    slack=min(1-diag[j]-sum(off[tuple(sorted((j,k)))] for k in range(3) if k!=j) for j in range(3))
    assert slack>0, 'Not a certified contraction'
    # Each detected row has a common phase factor across all six assignments.
    # Hence the permanent is a sum of positive square roots times a unit phase.
    boson_lo=boson_hi=dist=F(0)
    for y in counts['many_row_indices']:
        cols=[x for x in range(28) if any(prob(y,x,j)>0 for j in range(3))]
        for xs in combinations(cols,3):
            lo=hi=F(0)
            for perm in permutations(range(3)):
                prod=F(1)
                for r,j in enumerate(perm):prod*=prob(y,xs[r],j)
                dist+=prod
                l,u=sqrt_interval(prod);lo+=l;hi+=u
            boson_lo+=lo*lo;boson_hi+=hi*hi
    assert boson_lo>F(cert['bunch_lower'])
    # Also check an upper confidence endpoint for the model's bunching event.
    # This is an inner example for row calibration + bunching, not a fit to all many-body images.
    assert boson_hi < F(1,25) # .04; verified CP upper endpoint below
    check_endpoint(counts['many_shots'],counts['bunch_shots'],F(1,25),F(1,160),'upper')
    return {'contraction_slack_lower':float(slack),'bosonic_bunch_interval':[float(boson_lo),float(boson_hi)],
            'distinguishable_bunch_probability':float(dist),
            'singleton_calibration_exact':True,'joint_summary_feasible':True,
            'complete_many_body_record_feasibility_established':False}

def main():
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');args=p.parse_args()
    read=lambda f:json.loads((HERE/'results'/f).read_text())
    result=verify(read('physical-witness.json'),read('counts.json'),read('certificate.json'))
    text=json.dumps(result,indent=2,sort_keys=True)+'\n';path=HERE/'results/witness-check.json'
    if args.check:assert path.read_text()==text
    else:path.write_text(text)
    print(text)

if __name__=='__main__':main()
