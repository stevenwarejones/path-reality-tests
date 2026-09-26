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


def generate():
    rows=[]; alpha_cal=.025; alpha_main=.025; beta=.1
    for offset in (0.,.005):
        for sites in (96,128,192,256,512):
            certificate=None
            for target in (.02,.01,.005,.001,.0001):
                trial=certify_shared_gap(sites,target,phase_offset=offset,max_cells=4096)
                if trial['certified']:
                    certificate=trial; break
            gap=None if certificate is None else certificate['certified_gap']
            n=None if gap is None else required_trials(gap,72,alpha_main,beta)
            nc,cal=calibration_trials(.001,alpha_cal)
            rows.append(dict(sites=sites,phase_offset=offset,fractional_scale=.001,
                efficiency=.8,visibility=.9,calibration_radius=.001,calibration_phase_radius=0.,
                certified_gap=gap,trials_per_setting=n,total_main_trials=None if n is None else 24*n,
                calibration_trials_per_stratum=nc,total_calibration_trials=cal,
                total_eligible_trials=None if n is None else cal+24*n,
                interval_cells=None if certificate is None else certificate['cells'],
                decision_partition=None if certificate is None else certificate['decision_partition'],
                disposition='certified-shared-gap' if certificate else 'inconclusive-cell-search'))
    return dict(schema=1,kind='synthetic-shared-nuisance-design',data_inputs=[],
        alpha_total=alpha_cal+alpha_main,alpha_calibration=alpha_cal,alpha_main=alpha_main,
        beta_main=beta,guaranteed_unconditional_power_lower=1-alpha_cal-beta,
        search='descending finite target menu; lower certificates, not optimal gaps',
        calibration_reference='external exact phase calibration in this targeted comparison; baseline map propagates phase radius',
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
