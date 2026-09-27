"""External-source refit and independent source-to-prediction reproduction."""
import argparse
import json
from pathlib import Path
import time
import numpy as np
from sqd_sources import HERE,sources,split,digest
from sqd_models import Evaluator,parameters,certificate,reference_probability,density_probability,dormant_flag_probability
from sqd_inference import fit,deviance,diagnostic,context_certificate


def write(path,x):
    path.write_text(json.dumps(x,indent=2,sort_keys=True,allow_nan=False)+'\n')


def analyze(rows,out):
    train=split(rows);families=np.array([r['family'] for r in rows]);all_eval=Evaluator(rows)
    k=np.array([r['k'] for r in rows]);n=np.array([r['n'] for r in rows]);fits={}
    for model in ('qubit','leakage','alternating','quasistatic'):
        for source in ('GST','RB','joint'):
            mask=train if source=='joint' else train&(families==source)
            rr=[r for r,m in zip(rows,mask) if m];t=time.time()
            result=fit(rr,model)
            x=result['x'];p=all_eval(x,model)
            result['model']=model;result['source']=source
            result['physical_certificate']=certificate(x,model)
            result['observability']=diagnostic(Evaluator(rr),x,model,n[mask])
            result['scores']={}
            for fam in ('GST','RB'):
                for part,tt in [('train',train),('validation',~train)]:
                    m=tt&(families==fam)
                    result['scores'][fam+'_'+part]={'rows':int(m.sum()),'shots':int(n[m].sum()),
                        'deviance':float(deviance(k[m],n[m],p[m]).sum()),
                        'rms_probability_error':float(np.sqrt(np.mean((k[m]/n[m]-p[m])**2)))}
            result['probabilities']=p.tolist();fits[model+'_'+source]=result
            write(out/'fits.json',fits)
            print(model,source,'seconds',round(time.time()-t,1),'deviance',round(result['deviance'],2),'success',result['success'],flush=True)
    return fits


def audit(rows):
    ans={'families':{}}
    train=split(rows)
    for fam in ('GST','RB'):
        rr=[r for r in rows if r['family']==fam]
        ns=sorted(set(r['n'] for r in rr))
        ans['families'][fam]={'rows':len(rr),'unique_words':len(set(r['word_sha256'] for r in rr)),
            'total_shots':sum(r['n'] for r in rr),'denominators':ns,'length_range':[min(r['length'] for r in rr),max(r['length'] for r in rr)],
            'gate_occurrences':{str(g):sum(r['word'].count(g) for r in rr) for g in range(3)},
            'count_order_sha256':digest(json.dumps([[r['expression'],r['k'],r['n']] for r in rr],separators=(',',':')).encode())}
    ans['context_certificate']=context_certificate(rows)
    return ans


def verify(rows,fits):
    E=Evaluator(rows)
    for label,fit0 in fits.items():
        x,model=fit0['x'],fit0['model']
        p=E(x,model)
        err=float(np.max(np.abs(p-np.array(fit0['probabilities']))))
        if err>2e-9:raise ValueError('source prediction mismatch '+label)
        # Include short and long words, both families; independent density/Kraus.
        for idx in (0,1,5,60,len(rows)//2,4657,len(rows)-1):
            a=reference_probability(rows[idx]['word'],x,model)
            b=density_probability(rows[idx]['word'],x,model)
            if max(abs(p[idx]-a),abs(p[idx]-b))>2e-9:raise ValueError('physical propagation mismatch')
        if model=='leakage':
            for idx in (0,60,4657,len(rows)-1):
                if abs(dormant_flag_probability(rows[idx]['word'],x)-p[idx])>2e-9:raise ValueError('flag equivalence mismatch')
    bank_path=HERE/'results/removal.json'
    if bank_path.exists():
        bank=json.loads(bank_path.read_text())['bank']
        for label,b in bank.items():
            actual=E(b['x'],b['model'])
            if np.max(np.abs(actual-np.array(b['probabilities'])))>2e-9:
                raise ValueError('source bank prediction mismatch '+label)
    print('All source and witness-bank predictions and independent Kraus checks reproduce',flush=True)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--cache',required=True);ap.add_argument('--refit',action='store_true');ap.add_argument('--check',action='store_true');ap.add_argument('--output-dir',type=Path,default=HERE/'results');args=ap.parse_args()
    rows=sources(args.cache);out=args.output_dir;out.mkdir(parents=True,exist_ok=True)
    if args.refit:
        fits=analyze(rows,out)
        write(out/'observations.json',[{k:r[k] for k in ('family','row','word_sha256','k','n','length')} for r in rows])
        write(out/'audit.json',audit(rows))
    else:fits=json.loads((out/'fits.json').read_text())
    verify(rows,fits)
    if args.check:
        actual=[{k:r[k] for k in ('family','row','word_sha256','k','n','length')} for r in rows]
        if actual!=json.loads((out/'observations.json').read_text()):raise ValueError('count reduction mismatch')
        if json.loads(json.dumps(audit(rows)))!=json.loads((out/'audit.json').read_text()):raise ValueError('audit mismatch')


if __name__=='__main__':main()
