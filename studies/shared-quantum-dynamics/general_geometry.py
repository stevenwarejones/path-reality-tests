"""Observable Jacobian diagnostic with nuisance and gate-set gauges retained."""
import argparse,json
import numpy as np
from sqd_sources import HERE,sources,split
from sqd_general import GeneralEvaluator,bounds
from format_artifacts import formatted


def calculate(rows):
    train=split(rows);rr=[r for r,t in zip(rows,train) if t];families=np.array([r['family'] for r in rr]);E=GeneralEvaluator(rr);x=np.array(json.loads((HERE/'results/general-fits.json').read_text())['joint']['x']);lo,hi=bounds();p=E(x);cols=[]
    for j in range(52):
        a=x.copy();b=x.copy();a[j]=max(lo[j],x[j]-1e-4);b[j]=min(hi[j],x[j]+1e-4);cols.append((E(b)-E(a))/(b[j]-a[j]))
    J=np.array(cols).T;W=J*np.sqrt(np.array([r['n'] for r in rr])/np.maximum(p*(1-p),1e-10))[:,None];out={}
    for source in ['GST','RB','joint']:
        mask=np.ones(len(rr),bool) if source=='joint' else families==source;s=np.linalg.svd(W[mask],compute_uv=False)
        out[source]={'singular_values':s.tolist(),'rank_relative_1e-6':int(np.sum(s>s[0]*1e-6)),'rank_relative_1e-4':int(np.sum(s>s[0]*1e-4))}
    return {'status':'local observable derivative diagnostic at joint training fit; scaled redundant coordinates; not global identification',
        'coordinate_count':52,'generic_observable_dimension_upper':31,
        'dimension_count':'3 channels * 12 + preparation 3 + effect 4 - TP similarity gauge 12 = 31 at regular interior points. Factor coordinates retain additional normalization and unitary gauge redundancies.',
        'step':1e-4,'sources':out}

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--cache',required=True);a.add_argument('--check',action='store_true');args=a.parse_args();r=calculate(sources(args.cache));path=HERE/'results/general-geometry.json'
    if args.check:
        old=json.loads(path.read_text())
        for source in r['sources']:np.testing.assert_allclose(r['sources'][source]['singular_values'],old['sources'][source]['singular_values'],rtol=2e-3,atol=1e-5)
    else:path.write_text(formatted(r)+'\n')
    print({s:(v['rank_relative_1e-6'],v['rank_relative_1e-4']) for s,v in r['sources'].items()})
