#!/usr/bin/env python3
"""Reproducible shared-scale sensitivity boundary; synthetic, no source data."""
import argparse
import hashlib
import json
from pathlib import Path
from joint import certify_shared_gap
from model import required_trials,calibration_trials
from design import compare

ROOT=Path(__file__).resolve().parent


def refined_certificate(sites,*,refinements=6,**kwargs):
    """Retain only successful certificates; failed searches are not upper bounds.

    Refine between the best certified target and the preceding unsuccessful
    target (or the search cap). Each saved partition belongs to the returned gap.
    """
    certificate=None; lower=0.; upper=.02; attempts=[]
    def check(target):
        trial=certify_shared_gap(sites,target,**kwargs)
        attempts.append(dict(target=target,certified=trial['certified'],cells=trial['cells']))
        return trial
    for target in (.02,.01,.005,.001,.0001):
        trial=check(target)
        if trial['certified']:
            certificate=trial; lower=target; break
        upper=target
    if certificate is not None:
        for _ in range(refinements):
            if upper<=lower: break
            target=(lower+upper)/2
            trial=check(target)
            if trial['certified']: certificate=trial; lower=target
            else: upper=target
    return certificate,attempts,upper


def generate():
    rows=[]; alpha_cal=.025; alpha_main=.025; beta=.1
    for offset in (0.,.005):
        for sites in (96,128,192,256,512):
            certificate,attempts,search_upper=refined_certificate(sites,
                phase_offset=offset,calibration_phase_radius=offset,max_cells=4096)
            gap=None if certificate is None else certificate['certified_gap']
            n=None if gap is None else required_trials(gap,72,alpha_main,beta)
            nc,cal=calibration_trials(.001,alpha_cal)
            rows.append(dict(sites=sites,phase_offset=offset,fractional_scale=.001,
                efficiency=.8,visibility=.9,calibration_radius=.001,calibration_phase_radius=offset,
                certified_gap=gap,gap_search_attempts=attempts,search_upper_target=search_upper,
                trials_per_setting=n,total_main_trials=None if n is None else 24*n,
                calibration_trials_per_stratum=nc,total_calibration_trials=cal,
                total_eligible_trials=None if n is None else cal+24*n,
                interval_cells=None if certificate is None else certificate['cells'],
                decision_partition=None if certificate is None else certificate['decision_partition'],
                disposition='certified-shared-gap' if certificate else 'inconclusive-cell-search'))
    return dict(schema=1,kind='synthetic-shared-nuisance-design',data_inputs=[],
        alpha_total=alpha_cal+alpha_main,alpha_calibration=alpha_cal,alpha_main=alpha_main,
        beta_main=beta,guaranteed_unconditional_power_lower=1-alpha_cal-beta,
        search='coarse seed then six bisection refinements; failures are inconclusive, not upper bounds on separation',
        calibration_reference='calibration phase radius equals the main phase-offset radius; prospective envelopes include center drift plus interval bias',
        source_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
          for p in ('model.py','decision.py','joint.py','joint_design.py','protocol.md')},rows=rows)


def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--check',action='store_true'); args=parser.parse_args()
    result=generate(); path=ROOT/'results/shared-sensitivity.json'
    if args.check:
        compare(json.loads(path.read_text()),result)
        print(f"Verified {len(result['rows'])} shared-nuisance designs")
    else: path.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')


if __name__=='__main__': main()
