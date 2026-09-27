"""Hypothetical correlated blocks, explicitly not reconstructed acquisition jobs."""
import argparse,json
import numpy as np
from sqd_sources import HERE,sources,split
from sqd_models import Evaluator,operations
from sqd_inference import fit


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--cache',required=True);a=ap.parse_args()
    rows=sources(a.cache);rows=[r for r,m in zip(rows,split(rows)) if m];E=Evaluator(rows);n=np.array([r['n'] for r in rows])
    fits=json.loads((HERE/'results/fits.json').read_text());base=np.array(fits['qubit_joint']['x']);x=np.r_[base,1.]
    gates,initial,effect,_=operations(x,'quasistatic');curves=[]
    for b in range(2):
        same=np.repeat(gates[b:b+1],2,axis=0);p=np.empty(len(rows));E.lib.propagate(E.offsets,E.words,E.count,same,initial,effect,0,p);curves.append(p)
    rng=np.random.default_rng(260929);out={'seed':260929,'hypothetical_block_size_rows':64,'status':'synthetic dependence stress; file row adjacency is not verified chronology; no measured block calibration', 'qubit_scores':[],'marginal_drift_scores':[],'success':[]}
    for rep in range(20):
        signs=np.repeat(rng.integers(0,2,size=(len(rows)+63)//64),64)[:len(rows)]
        p=np.where(signs,curves[1],curves[0]);counts=rng.binomial(n,p)
        q=fit(rows,'qubit',starts=1,counts=counts,initial=base)
        d=fit(rows,'quasistatic',starts=1,counts=counts,initial=x)
        out['qubit_scores'].append(q['deviance']);out['marginal_drift_scores'].append(d['deviance']);out['success'].append([q['success'],d['success']])
        (HERE/'results/dependence.json').write_text(json.dumps(out,indent=2)+'\n')
    print('Correlated ordinary-drift control complete',flush=True)


if __name__=='__main__':main()
