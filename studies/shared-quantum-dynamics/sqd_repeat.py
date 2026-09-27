"""Optimizer-independent, all-stationary-qubit outer bound for a doubled word."""
import argparse,json,math
import numpy as np
from scipy.stats import beta
from sqd_sources import HERE,sources,split
from format_artifacts import formatted

BINS=[(2,16),(17,128),(129,256),(257,512),(513,1000)]

def repeat_upper(p0,p1):
    if not 0<=p0<=1 or not 0<=p1<=1:raise ValueError('probability outside [0,1]')
    return min(1.,(2*np.sqrt(p1)+np.sqrt(p0))**2)

GRID=10**10

def exact_upper(k,n,tail_denominator):
    """Outward rational binomial upper bound; integer arithmetic certifies tail."""
    if k==n:return GRID
    a=math.ceil(float(beta.ppf(1-1/tail_denominator,k+1,n-k))*GRID)
    def accepted(v):
        return tail_denominator*sum(math.comb(n,j)*v**j*(GRID-v)**(n-j) for j in range(k+1))<=GRID**n
    while not accepted(a):a+=1
    while a>0 and accepted(a-1):a-=1
    return a

def rational_repeat(a0,a1):
    # ceil(sqrt(a / GRID) * GRID), followed by an outward squared sum.
    roots=[]
    for a in [a0,a1]:
        b=math.isqrt(a*GRID)
        roots.append(b+(b*b<a*GRID))
    numerator=(roots[0]+2*roots[1])**2
    return min(GRID,(numerator+GRID-1)//GRID)

def calculate(rows):
    empty=[r for r in rows if r['family']=='GST' and r['length']==0]
    if len(empty)!=1:raise ValueError('expected one empty GST word')
    selected=[]
    for low,high in BINS:
        candidates=[r for r,t in zip(rows,split(rows)) if t and r['family']=='RB' and low<=r['length']<=high]
        selected.append(min(candidates,key=lambda r:(r['word_sha256'],r['row'])))
    # One-sided Bonferroni, six fixed rows; independent shots within each row.
    # No independence between rows is needed by the union bound.
    tail_denominator=20*(len(selected)+1);tail=1/tail_denominator
    upper=lambda r:exact_upper(r['k'],r['n'],tail_denominator)
    a0=upper(empty[0]);u0=a0/GRID;records=[];known={tuple(r['word']) for r in rows}
    for r in selected:
        records.append({**r,'doubled_length':2*r['length'],'target_already_recorded':tuple(r['word'])*2 in known,'rb_upper_numerator':upper(r),'rb_upper':upper(r)/GRID,
            'joint_upper_numerator':rational_repeat(a0,upper(r)), 'joint_upper':rational_repeat(a0,upper(r))/GRID,'GST_removed_relaxation_upper':1.,'RB_removed_relaxation_upper':1.})
    return {'status':'certified full-class outer bounds, conditional on stationary qubit CPTP and row-binomial sampling',
        'confidence':.95,'one_sided_tail':tail,'empty':empty[0],'empty_upper':u0,'empty_upper_numerator':a0,'rational_denominator':GRID,'tail_denominator':tail_denominator,'targets':records,
        'rounding':'rational outward endpoints; binomial tails checked by exact integer arithmetic','dependence_unrestricted_upper':1.,
        'source_removal_scope':'removing a source makes this relaxation vacuous; not a proof that its full feasible model class attains one'}

def main():
    a=argparse.ArgumentParser();a.add_argument('--cache',required=True);a.add_argument('--check',action='store_true');args=a.parse_args()
    result=calculate(sources(args.cache));path=HERE/'results/repeat-bounds.json';serialized=formatted(result)+'\n'
    if args.check:
        if json.loads(path.read_text())!=json.loads(serialized):raise ValueError('repeat certificate mismatch')
    else:path.write_text(serialized)
    print([(r['length'],r['k'],r['n'],r['joint_upper']) for r in result['targets']])

if __name__=='__main__':main()
