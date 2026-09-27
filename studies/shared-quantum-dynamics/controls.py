"""Retrospective nuisance-refitted ordinary adversaries at actual exposures."""
import argparse
import json
import time
import numpy as np
from sqd_sources import HERE,sources,split
from sqd_models import Evaluator,parameters
from sqd_inference import fit


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--cache',required=True);a=ap.parse_args()
    rows=sources(a.cache);rows=[r for r,m in zip(rows,split(rows)) if m]
    E=Evaluator(rows);n=np.array([r['n'] for r in rows]);fits=json.loads((HERE/'results/fits.json').read_text())
    base=np.array(fits['qubit_joint']['x']);rng=np.random.default_rng(260928)
    out={'seed':260928,'null_repetitions':39,'repetitions_per_control':20,'training_rows':len(rows),'training_shots':int(n.sum()),
         'null_scores':[],'controls':{},'scope':'conditional independent binomial shots; local nuisance refits; no mechanism classifier and no probability deformation'}
    path=HERE/'results/controls.json'
    def save():path.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    for rep in range(39):
        counts=rng.binomial(n,E(base,'qubit'));f=fit(rows,'qubit',starts=1,counts=counts,initial=base)
        out['null_scores'].append(f['deviance']);save()
        if rep%10==0:print('null',rep,flush=True)
    threshold=sorted(out['null_scores'])[37]
    out['empirical_95_threshold']=threshold
    cases={}
    cases['leakage']=('leakage',np.r_[base,[3.,2.,.9]])
    cases['alternating_memory']=('alternating',np.r_[base,2.])
    cases['quasistatic_drift']=('quasistatic',np.r_[base,.5])
    imperfect=base.copy();imperfect[14]=.93;cases['imperfect_reset']=('qubit',imperfect)
    readout=base.copy();readout[12:14]=[.035,.96];cases['readout_bias']=('qubit',readout)
    # Weakest non-gauge direction from the earlier declared generic Jacobian is
    # not a detectable alternative along the exact leakage fiber: power is zero
    # above test size for ANY count statistic, proved separately.
    for name,(model,x) in cases.items():
        p=E(x,model);entry={'model':model,'injected_parameters':x.tolist(),'qubit_scores':[], 'generating_class_scores':[], 'generating_class_parameters':[],'qubit_success':[],'generating_class_success':[]}
        out['controls'][name]=entry
        for rep in range(20):
            counts=rng.binomial(n,p)
            f=fit(rows,'qubit',starts=1,counts=counts,initial=base)
            right=f if model=='qubit' else fit(rows,model,starts=1,counts=counts,initial=x)
            entry['qubit_scores'].append(f['deviance']);entry['generating_class_scores'].append(right['deviance']);entry['generating_class_parameters'].append(right['x'])
            entry['qubit_success'].append(f['success']);entry['generating_class_success'].append(right['success'])
            save()
        entry['narrow_qubit_rejections']=int(np.sum(np.array(entry['qubit_scores'])>threshold))
        entry['generating_class_above_qubit_threshold']=int(np.sum(np.array(entry['generating_class_scores'])>threshold))
        save();print(name,'narrow rejections',entry['narrow_qubit_rejections'],'of 20',flush=True)
    out['weakest_structural_direction']={'mechanism':'exact leakage/readout/depolarization fiber','terminal_total_variation':0,'power_upper_bound':'test size, for every exposure under matched sampling law','status':'analytic, not a Monte Carlo estimate'}
    save()


if __name__=='__main__':main()
