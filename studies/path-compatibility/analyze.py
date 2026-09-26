#!/usr/bin/env python3
"""Reproduce exact full-table compatibility and conditional sensitivity."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
from fractions import Fraction as F
import boundary

ROOT=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('contextuality_design',ROOT.parent/'path-contextuality/design.py')
comparison=importlib.util.module_from_spec(spec);spec.loader.exec_module(comparison)


def stringify(value):
    if isinstance(value,F): return str(value)
    if isinstance(value,dict): return {k:stringify(v) for k,v in value.items()}
    if isinstance(value,(tuple,list)): return [stringify(v) for v in value]
    return value


def analyze():
    caps=[boundary.G,F(16,25),boundary.CROSS,F(1081,1225),F(1)]
    models=[boundary.attaining_model(q) for q in caps]
    sensitivity=[]
    for eps in [F(0),F(1,1000),F(1,100),F(1,20)]:
        sensitivity.append({'b_probability_error':eps,'bypass_interval':[boundary.F0,boundary.F0],
                            'preparation_tv':eps,'readout_probability_error':eps,
                            'positive_success_floor':boundary.positive_floor(F(16,25),F(1,50),
                                  (boundary.F0,boundary.F0),eps,eps,eps),
                            'reference_observed_positive_success':boundary.B})
    source_hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                   for p in sorted(ROOT.glob('*.py'))}
    return stringify({'scope':'fixed ideal reference-table compatibility; no experimental evidence',
            'source_hashes':source_hashes,
            'boundary':{'reference_joint':boundary.TARGET,'bypass':boundary.F0,
                        'minimum_cap':boundary.G,'branch_crossing':boundary.CROSS,
                        'cap_at_minimum_disturbance':boundary.Q_MAX,
                        'reference_cap_minimum_disturbance':boundary.required_disturbance(F(16,25)),
                        'reference_cap_negative_only_necessary_disturbance':F(13,320),
                        'attaining_models':models,'conditional_sensitivity':sensitivity}})


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--check',action='store_true');args=ap.parse_args()
    result=analyze();path=ROOT/'results.json'
    if args.check:
        comparison.compare(json.loads(path.read_text()),result)
        print('Sharp compatibility verified')
    else:
        path.write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__': main()
