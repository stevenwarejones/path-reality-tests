"""Pilot-based nuisance-adapted parametric calibration; exploratory, not a uniform test."""
import argparse,json
import numpy as np
from sqd_sources import HERE,sources,split
from sqd_models import Evaluator
from sqd_inference import fit


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--cache',required=True);a=ap.parse_args();rows=sources(a.cache);rows=[r for r,m in zip(rows,split(rows)) if m];E=Evaluator(rows);n=np.array([r['n'] for r in rows]);rng=np.random.default_rng(260930)
    controls=json.loads((HERE/'results/controls.json').read_text());fits=json.loads((HERE/'results/fits.json').read_text());base=np.array(fits['qubit_joint']['x'])
    out={'seed':260930,'calibration_replicates':39,'scope':'separate synthetic pilot at actual exposures, fitted qubit null, all bootstrap samples nuisance-refitted; empirical conditional calibration, NOT uniform coverage or global class exclusion','cases':{}}
    path=HERE/'results/calibration.json'
    for name,c in controls['controls'].items():
        p=E(c['injected_parameters'],c['model']);pilot=rng.binomial(n,p);f=fit(rows,'qubit',starts=1,counts=pilot,initial=base);null=np.array(f['x']);pn=E(null,'qubit');scores=[]
        for rep in range(39):
            counts=rng.binomial(n,pn);g=fit(rows,'qubit',starts=1,counts=counts,initial=null);scores.append(g['deviance'])
        threshold=sorted(scores)[37];out['cases'][name]={'pilot_null_parameters':null.tolist(),'pilot_fit_success':f['success'],'null_scores':scores,'threshold':threshold,
            'qubit_exceedances':int(np.sum(np.array(c['qubit_scores'])>threshold)),'trials':len(c['qubit_scores'])}
        path.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(name,'adapted exceedances',out['cases'][name]['qubit_exceedances'],flush=True)


if __name__=='__main__':main()
