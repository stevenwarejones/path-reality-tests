#!/usr/bin/env python3
"""Deterministic synthetic sensitivity map; --check never rewrites snapshots."""
import argparse
import hashlib
import json
from pathlib import Path
from math import pi
import numpy as np
from model import Ring, phase_gap, outcome_box, box_gap, required_trials, calibration_trials, calibration_phase_bias

ROOT=Path(__file__).resolve().parent


def generate():
    ring=Ring()
    alpha_calibration=.025; alpha_main=.025; beta_main=.1
    # Every candidate uses the SAME physical acquisition menu and calibration.
    menu=[(j,t,q) for j in (1,2,3) for t in (.5,2.,8.,32.) for q in (0.,pi/2)]
    rows=[]
    for eta,visibility in ((1.,1.),(.8,.9),(.25,.5),(0.,1.),(.8,0.)):
        for phase_offset in (0.,.005,.05):
            for fractional_scale in (0.,.001):
                for cal_radius in (.001,.005):
                    for sites in (8,12,16,24,32,48,64,96,128,192,256):
                        gaps=[]
                        for j,t,q in menu:
                            delta=phase_gap(ring,sites,0,j,t)
                            # Common scale = t*hbar/(m*L^2). Null and alternative
                            # boxes may each vary over its shared calibrated set.
                            ec=float(t*ring.continuum(j)/ring.hbar)
                            ea=float(t*ring.lattice(j,sites)/ring.hbar)
                            # Power envelopes allow calibration centers to move by r:
                            # intervals around centers add another r (2r and 4r).
                            nuisance=dict(eta=eta,eta_radius=2*cal_radius,
                                contrast=eta*visibility,contrast_radius=4*cal_radius+2*calibration_phase_bias(phase_offset))
                            continuum=outcome_box(-q,phase_offset+fractional_scale*abs(ec),**nuisance)
                            lattice=outcome_box(delta-q,phase_offset+fractional_scale*abs(ea),**nuisance)
                            gaps.append(box_gap(continuum,lattice))
                        best=int(np.argmax(gaps)); gap=gaps[best]
                        # One simultaneous confidence box excludes any number of
                        # models: no extra factor for number of candidate N values.
                        n=required_trials(gap,3*len(menu),alpha=alpha_main,beta=beta_main)
                        per_cal,total_cal=calibration_trials(cal_radius,error=alpha_calibration)
                        eligible=None if n is None else len(menu)*n+total_cal
                        rows.append(dict(sites=sites,spacing=ring.length/sites,
                            efficiency=eta,visibility=visibility,phase_offset=phase_offset,
                            fractional_scale=fractional_scale,calibration_radius=cal_radius,
                            calibration_phase_radius=phase_offset,
                            certified_coordinate_gap=gap, best_setting=list(menu[best]),
                            trials_per_setting=n,main_strata=len(menu),
                            calibration_trials_per_stratum=per_cal,calibration_strata=3,
                            total_calibration_trials=total_cal,total_main_trials=None if n is None else len(menu)*n,
                            total_eligible_trials=eligible,
                            excludes_with_100000000_eligible=(eligible is not None and eligible<=100000000),
                            total_source_attempts=None,
                            decision='conditional-power-guarantee' if n else 'no-certified-gap'))
    hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in ('model.py','decision.py','joint.py','design.py','protocol.md','model.md')}
    return dict(schema=1,kind='synthetic-prospective-design',data_inputs=[],source_sha256=hashes,
        alpha_total=alpha_calibration+alpha_main,alpha_calibration=alpha_calibration,alpha_main=alpha_main,beta_main=beta_main,
        guaranteed_unconditional_power_lower=1-alpha_calibration-beta_main,
        units='dimensionless L=2pi, hbar=m=1; no apparatus performance is inferred',
        source_attempts_status='requires external preparation-efficiency and scale/phase calibration certificates',
        settings=[list(s) for s in menu],rows=rows)


def compare(expected,actual,path='$'):
    if type(expected) is not type(actual):
        raise AssertionError(f'{path}: type mismatch')
    if isinstance(actual,dict):
        if expected.keys()!=actual.keys():
            raise AssertionError(f'{path}: keys differ')
        for k in actual: compare(expected[k],actual[k],f'{path}.{k}')
    elif isinstance(actual,list):
        if len(expected)!=len(actual): raise AssertionError(f'{path}: length differs')
        for i,(e,a) in enumerate(zip(expected,actual)): compare(e,a,f'{path}[{i}]')
    elif isinstance(actual,float):
        if not np.isclose(expected,actual,rtol=2e-12,atol=2e-14):
            raise AssertionError(f'{path}: {expected} != {actual}')
    elif expected!=actual:
        raise AssertionError(f'{path}: {expected!r} != {actual!r}')


def svg(result):
    rows=[r for r in result['rows'] if r['efficiency']==.8 and r['visibility']==.9
          and r['fractional_scale']==0 and r['calibration_radius']==.001]
    parts=['<svg xmlns="http://www.w3.org/2000/svg" width="820" height="390" viewBox="0 0 820 390">',
        '<rect width="820" height="390" fill="white"/>',
        '<g font-family="sans-serif" font-size="13" fill="#172635">',
        '<text x="40" y="26" font-size="18">Synthetic guaranteed coordinate separation</text>',
        '<text x="40" y="48">24 shared settings; efficiency 0.8; visibility 0.9; calibration radius 0.001</text>',
        '<path d="M65 70 V320 H770" stroke="#172635" fill="none"/>']
    for value in (0,.2,.4,.6,.8):
        y=320-280*value
        parts.append(f'<text x="28" y="{y+4:.3f}">{value:.1f}</text>')
    ns=sorted(set(r['sites'] for r in rows))
    for i,n in enumerate(ns):
        parts.append(f'<text x="{65+i*69:.3f}" y="342">{n}</text>')
    for offset,color in ((0.,'#086e9d'),(.005,'#a74b00'),(.05,'#75439c')):
        selected=[r for r in rows if r['phase_offset']==offset]
        points=' '.join(f"{65+ns.index(r['sites'])*69:.3f},{320-280*r['certified_coordinate_gap']:.3f}" for r in selected)
        parts.append(f'<polyline points="{points}" stroke="{color}" stroke-width="2" fill="none"/>')
        parts.append(f'<text x="{80+list((0.,.005,.05)).index(offset)*235}" y="374" fill="{color}">Phase radius {offset:g} rad</text>')
    parts.append('<text x="660" y="358">Integer site count N</text></g></svg>\n')
    return '\n'.join(parts)


def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--check',action='store_true')
    args=parser.parse_args(); result=generate()
    output=ROOT/'results'/'sensitivity.json'; figure=ROOT/'results'/'sensitivity.svg'
    if args.check:
        compare(json.loads(output.read_text()),result)
        if figure.read_text()!=svg(result): raise AssertionError('sensitivity.svg differs')
        print(f"Verified {len(result['rows'])} synthetic scenarios, hashes and figure")
    else:
        output.parent.mkdir(exist_ok=True)
        header={k:v for k,v in result.items() if k!='rows'}
        prefix=json.dumps(header,indent=2,allow_nan=False)[:-2]
        output.write_text(prefix+',\n  \"rows\": [\n'+',\n'.join('    '+json.dumps(row,allow_nan=False) for row in result['rows'])+'\n  ]\n}\n')
        figure.write_text(svg(result))


if __name__=='__main__': main()
