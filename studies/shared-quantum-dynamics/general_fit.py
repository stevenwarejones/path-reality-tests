"""General CPTP escalation: retained split, plus explicitly in-sample diagnosis."""
import argparse,json
import numpy as np
from sqd_sources import HERE,sources,split
from sqd_general import GeneralEvaluator,from_restricted,fit,certificate,density_probability
from sqd_inference import deviance,simultaneous_intervals
from format_artifacts import formatted


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--cache',required=True);ap.add_argument('--check',action='store_true');ap.add_argument('--all-lengths',action='store_true');a=ap.parse_args();rows=sources(a.cache);train=split(rows);family=np.array([r['family'] for r in rows]);old=json.loads((HERE/'results/fits.json').read_text());E=GeneralEvaluator(rows);k=np.array([r['k'] for r in rows]);n=np.array([r['n'] for r in rows]);path=HERE/'results/general-fits.json'
    results=json.loads(path.read_text()) if path.exists() else {}
    if a.check:
        for key,r in results.items():
            p=E(r['x'])
            if np.max(abs(p-r['probabilities']))>2e-8:raise ValueError('general prediction mismatch')
            for i in [0,10,4657,len(rows)-1]:
                if abs(density_probability(rows[i]['word'],r['x'])-p[i])>2e-8:raise ValueError('independent density mismatch')
        print('General CPTP source predictions reproduce');return
    for source in (['all_joint'] if a.all_lengths else ['GST','RB','joint']):
        if source in results:continue
        mask=np.ones(len(rows),bool) if source=='all_joint' else train if source=='joint' else train&(family==source)
        initial=np.array(results['joint']['x']) if source=='all_joint' else from_restricted(old['qubit_'+source]['x'])
        r=fit([row for row,m in zip(rows,mask) if m],initial);p=E(r['x']);r['probabilities']=p.tolist();r['scores']={}
        for f in ['GST','RB']:
            for part,tt in [('train',train),('validation',~train)]:
                mm=tt&(family==f);lo,hi=simultaneous_intervals(k[mm],n[mm]);r['scores'][f+'_'+part]={'rows':int(mm.sum()),'deviance':float(deviance(k[mm],n[mm],p[mm]).sum()),'cell_violations':int(np.sum((p[mm]<lo)|(p[mm]>hi))), 'rms_probability_error':float(np.sqrt(np.mean((k[mm]/n[mm]-p[mm])**2)))}
        r['status']='in-sample diagnostic, no held-out claim' if source=='all_joint' else 'frozen training-only general CPTP fit; local feasible example'
        results[source]=r;path.write_text(formatted(results)+'\n');print(source,r['scores'],flush=True)


if __name__=='__main__':main()
