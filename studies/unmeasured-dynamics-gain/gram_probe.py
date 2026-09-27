"""Exploratory SDP outer relaxation; solver objectives are NOT certificates."""
import argparse,time
from dynamics_common import *
from dynamics_verify import TARGET

def graph(rows,limit,source):
    ids=[i for i,r in enumerate(rows) if r['length']<=limit and (source=='joint' or r['family']==source)]
    words={TARGET[:j] for j in range(len(TARGET)+1)}
    for i in ids:
        w=rows[i]['word'];words.update(w[:j] for j in range(len(w)+1))
    words=sorted(words,key=lambda w:(len(w),w));index={w:i for i,w in enumerate(words)}
    pairs=[[(index[w],index[w+(g,)]) for w in words if w+(g,) in index] for g in range(3)]
    return ids,words,index,pairs

def contraction_matrix(Q,pairs,t):
    before=[i for i,j in pairs];after=[j for i,j in pairs]
    return Q[np.ix_(before,before)]-Q[np.ix_(after,after)]+Q[after,t,None]+Q[None,t,after]-Q[t,t]

def solve(rows,R,limit,source):
    import cvxpy as cp
    ids,words,index,pairs=graph(rows,limit,source);n=len(words);e=n;N=n+4
    if n>160:raise ValueError('exploratory graph size exceeds declared resource cap')
    Q=cp.Variable((N,N),symmetric=True);a=cp.Variable()
    cons=[Q>>0,cp.diag(Q)[:n]<=1,cp.diag(Q)[n+1:]<=1,Q[e,e]<=.25,Q[e,e]<=a/2,Q[e,e]<=(1-a)/2,a>=0,a<=1]
    p=a+Q[:n,e];cons +=[p>=0,p<=1]
    for i in ids:
        q=p[index[rows[i]['word']]];cons.extend([q>=R['lo'][i],q<=R['hi'][i]])
    for g,pair in enumerate(pairs):
        before=[i for i,j in pair];after=[j for i,j in pair];t=n+1+g
        M=Q[np.ix_(before,before)]-Q[np.ix_(after,after)]+cp.reshape(Q[after,t],(len(after),1),order='F')+cp.reshape(Q[t,after],(1,len(after)),order='F')-Q[t,t]
        cons.append(M>>0)
        for j in after:cons.append(4*Q[t,t]-4*Q[t,j]+Q[j,j]<=1)
    prob=cp.Problem(cp.Minimize(p[index[TARGET]]),cons);prob.solve(solver='SCS',eps=1e-6,max_iters=10000)
    if prob.status not in ['optimal','optimal_inaccurate']:raise ValueError('probe solver did not return a usable discovery result')
    return {'source':source,'max_word_length':limit,'rows_used':len(ids),'prefix_states':n,'solver_status':prob.status,'numerical_objective':float(prob.value),'readout_bias':float(a.value),
        'scope':'Uncertified numerical objective of an outer relaxation; not a physical witness, dual bound, or full-source exclusion'}

def run(cache):
    import cvxpy,scs
    rows=sources(cache);R=region(rows)
    return {'versions':{'cvxpy':cvxpy.__version__,'scs':scs.__version__},'results':[solve(rows,R,L,s) for L,s in [(4,'GST'),(5,'GST'),(4,'joint')]],
        'missing_constraints':'No rank-three Gram restriction, full effect positivity, complete positivity of native Choi matrices, or full measured long-word/group constraints',
        'containment':'Every full physical stationary qubit model maps into this relaxation; converse not asserted'}
if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--cache',required=True);a.add_argument('--check',action='store_true');args=a.parse_args();out=run(args.cache);path=HERE/'results/gram-probes.json'
    if args.check:
        old=json.loads(path.read_text())
        for a,b in zip(old['results'],out['results']):
            if (a['source'],a['rows_used'],a['prefix_states'])!=(b['source'],b['rows_used'],b['prefix_states']) or abs(a['numerical_objective']-b['numerical_objective'])>1e-4:raise ValueError('exploratory probe changed')
    else:save(path,out)
    print(out)
