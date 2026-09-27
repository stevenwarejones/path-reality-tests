#!/usr/bin/env python3
"""Optional local proposals; control_design.py independently certifies acceptance.

Different solver versions may give different feasible points. Reproduction of
committed exact certificates does not rerun a floating-point optimization.
"""
import json
from build_robustness import drift_models,HERE
read=lambda n:json.loads((HERE/'results'/n).read_text())
if __name__=='__main__':
    for label,target,seeds in [('above',.1432359739972339,8),('below',.14320,2)]:
        for sign,model in zip(['minus','plus'],drift_models(read,target,seeds)):
            (HERE/'results'/f'pair-{label}-{sign}.json').write_text(json.dumps(model,sort_keys=True,indent=2)+'\n')
    print('Proposals saved; run control_design.py to certify exact physical feasibility.')
