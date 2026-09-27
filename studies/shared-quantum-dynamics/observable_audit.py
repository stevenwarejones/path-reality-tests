"""A constructive failure of one interval-inversion relaxation, not a no-go."""
import argparse,json,itertools,math
import numpy as np
from sqd_sources import HERE,sources
from sqd_inference import simultaneous_intervals
from format_artifacts import formatted


def calculate(rows):
    # Standard nominally complete preparation/measurement fiducials, fixed by
    # design rather than outcome: identity, X90, Y90, X180.
    fid=[(),(1,),(2,),(1,1)];lookup={tuple(r['word']):r for r in rows if r['family']=='GST'}
    rr=[lookup[a+b] for b in fid for a in fid];k=np.array([r['k'] for r in rr]);n=np.array([r['n'] for r in rr])
    # Same conservative primary source-removal cell budget, no spending reset.
    lo,hi=simultaneous_intervals(k,n,alpha=.025*len(rr)/3785);lo=lo.reshape(4,4);hi=hi.reshape(4,4)
    point=k.reshape(4,4)/n.reshape(4,4)
    # Find two realizations in this box with opposite exact determinant signs.
    # Equal circuit words share one variable, including repeated fiducial words.
    unique=sorted(set(tuple(r['word']) for r in rr));index={w:i for i,w in enumerate(unique)}
    positions=np.array([index[tuple(r['word'])] for r in rr]).reshape(4,4)
    scale=10**8;low=[];high=[]
    for w in unique:
        r=lookup[w];l,h=simultaneous_intervals(np.array([r['k']]),np.array([r['n']]),alpha=.025/3785)
        low.append(math.ceil((float(l[0])+1e-10)*scale));high.append(math.floor((float(h[0])-1e-10)*scale))
    def determinant(A):
        return sum((-1)**sum(p[i]>p[j] for i in range(4) for j in range(i+1,4))*math.prod(int(A[i,p[i]]) for i in range(4)) for p in itertools.permutations(range(4)))
    first=np.rint(point*scale).astype(np.int64);sign=determinant(first);second=None;rng=np.random.default_rng(191)
    for _ in range(10000):
        v=np.where(rng.integers(0,2,len(unique)),high,low);candidate=v[positions]
        if determinant(candidate)*sign<0:second=candidate;break
    feasible=second is not None
    return {'status':'failure of chosen observable-basis interval relaxation only; not impossibility of stronger certification',
        'fiducials':[list(w) for w in fid],'rows':rr,'lower':lo.tolist(),'upper':hi.tolist(),
        'point_singular_values':np.linalg.svd(point,compute_uv=False).tolist(),
        'singular_witness_found':feasible,'determinant_segment_denominator':scale,'segment_first_numerators':first.tolist(),'segment_second_numerators':second.tolist() if feasible else None,'first_exact_determinant':str(sign),'second_exact_determinant':str(determinant(second)) if feasible else None,
        'note':'opposite exact determinant signs certify a singular H on the segment inside the entrywise count box; this relaxation cannot uniformly invert H. Witness need not obey a common CPTP realization or all other data.'}

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--cache',required=True);a.add_argument('--check',action='store_true');args=a.parse_args();result=calculate(sources(args.cache));path=HERE/'results/observable-audit.json'
    if args.check:
        from certificates import close
        close(result,json.loads(path.read_text()))
    else:path.write_text(formatted(result)+'\n')
    print(result['singular_witness_found'],result['point_singular_values'])
