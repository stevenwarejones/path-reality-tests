#!/usr/bin/env python3
"""Conditional binomial power; does not pretend to resample calibration."""
import argparse,json
from pathlib import Path
from scipy.stats import binom
HERE=Path(__file__).resolve().parent

def calculate():
    c=json.loads((HERE/'results/counts.json').read_text())
    b=json.loads((HERE/'results/bounds.json').read_text())
    w=json.loads((HERE/'results/witness-check.json').read_text())
    n=c['many_shots'];p0=b['cluster_two_ceiling']
    k=next(k for k in range(n+1) if binom.sf(k-1,n,p0)<=1/160)
    point=c['bunch_shots']/n
    inner=sum(w['bosonic_bunch_interval'])/2
    return {'conditional_on_frozen_calibration':True,'trials':n,'critical_count':k,
       'null_tail':round(float(binom.sf(k-1,n,p0)),12),
       'power_at_observed_frequency':round(float(binom.sf(k-1,n,point)),12),
       'power_at_physical_inner_example':round(float(binom.sf(k-1,n,inner)),12),
       'power_at_probability_0_04':round(float(binom.sf(k-1,n,.04)),12),
       'synthetic_common_drift':{'singleton_row_probabilities':[.5,.5],
           'product_of_averages_D':.25,'C2_ceiling_if_sharing_were_valid':.5,
           'actual_distinguishable_bunch_probability':1.,'violates_fixed_channel_assumption':True}}

def main():
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args()
    text=json.dumps(calculate(),indent=2,sort_keys=True)+'\n';out=HERE/'results/diagnostics.json'
    if a.check:assert out.read_text()==text
    else:out.write_text(text)
    print(text)

if __name__=='__main__':main()
