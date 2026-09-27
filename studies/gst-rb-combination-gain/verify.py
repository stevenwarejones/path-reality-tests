"""Exact count-region verification; optional independent full-source propagation."""
import argparse
import re
from fractions import Fraction as F
import numpy as np
from gain_core import HERE,BASE,GRID,region,save,sources,operations,GeneralEvaluator
from interval_physics import certify
from sqd_models import PAULI
import json

ERROR=F(1,10**9)

def fraction_record(x):
    return {'numerator':x.numerator,'denominator':x.denominator,'decimal':float(x)}

def exact_audit(rows,model,R):
    p=[F.from_float(float(v)) for v in model['probabilities']]
    if len(p)!=len(rows) or any(v<0 or v>1 for v in p):raise ValueError('invalid probability artifact')
    physical=certify(model['x'],max(r['length'] for r in rows))
    bound=F(physical['probability_error_numerator'],physical['probability_error_denominator'])
    if bound>ERROR:raise ValueError('roundoff certificate exceeds common envelope')
    means={};slacks={f:[] for f in ['GST','RB']}
    for i,(r,v) in enumerate(zip(rows,p)):
        slacks[r['family']].extend([v-F(R['lo_num'][i],GRID)-ERROR,F(R['hi_num'][i],GRID)-v-ERROR])
    for g in R['groups']:
        ids=[i for i,r in enumerate(rows) if r['family']==g['family'] and g['min_length']<=r['length']<=g['max_length']]
        mean=sum(rows[i]['n']*p[i] for i in ids)/g['n']
        slacks[g['family']].extend([mean-F(g['lower_numerator'],g['lower_denominator'])-ERROR,F(g['upper_numerator'],g['upper_denominator'])-mean-ERROR])
        if g['min_length']==0 and g['max_length']==8198:means[g['family']]=mean
    q=(1+means['RB']-means['GST'])/2
    retained=['GST','RB'] if model['retained']=='joint' else [model['retained']]
    for f in retained:
        if min(slacks[f])<=F(2,10**8):raise ValueError('retained constraints not strictly certified: '+f)
    gap=q-ERROR-R['joint_upper_q_exact']
    if model['retained']!='joint' and gap<=F(2,10**8):raise ValueError('source-removal separation not certified')
    return {'retained':model['retained'],'q_interval':[float(np.nextafter(float(q-ERROR),-np.inf)),float(np.nextafter(float(q+ERROR),np.inf))],
        'q_lower_exact':fraction_record(q-ERROR),'q_upper_exact':fraction_record(q+ERROR),'separation_lower_exact':fraction_record(gap),
        'mean_intervals':{f:[float(np.nextafter(float(v-ERROR),-np.inf)),float(np.nextafter(float(v+ERROR),np.inf))] for f,v in means.items()},
        'minimum_certified_slack':{f:float(min(v)) for f,v in slacks.items()},
        'violated_inequalities':{f:sum(s<0 for s in v) for f,v in slacks.items()},
        'checked_two_sided_constraints':{f:len(v)//2 for f,v in slacks.items()},
        'physicality':'exact normalized Kraus construction; directed-interval input marginals positive',
        'numerical_probability_envelope':fraction_record(ERROR)}

def independent_probabilities(rows,x):
    """Complex density superoperators and compressed expression powers, no Bloch kernel."""
    _,_,_,kk=operations(x)
    gates=[sum(np.kron(K.conj(),K) for K in ks) for ks in kk]
    c,l,u,theta=x[48:]
    rho=((PAULI[0]+c*PAULI[3])/2).reshape(-1,order='F')
    effect=((l+u)*PAULI[0]/2+(l-u)*(np.sin(theta)*PAULI[1]+np.cos(theta)*PAULI[3])/2).T.reshape(-1,order='F')
    cache={};out=[];alphabet={'Gi':0,'Gx':1,'Gy':2}
    for row in rows:
        state=rho.copy()
        if row['expression']!='{}':
            for token in re.findall(r'\((?:G[ixy])+\)\^[1-9][0-9]*|G[ixy]',row['expression']):
                if token not in cache:
                    if token.startswith('('):body,exponent=token[1:].split(')^');power=int(exponent)
                    else:body=token;power=1
                    M=np.eye(4,dtype=complex)
                    for name in re.findall(r'G[ixy]',body):M=gates[alphabet[name]]@M
                    cache[token]=np.linalg.matrix_power(M,power)
                state=cache[token]@state
        out.append(float((effect@state).real))
    return np.array(out)

def run(cache=None,download=False):
    rows=json.loads((BASE/'results/observations.json').read_text())
    if isinstance(rows,dict):rows=rows['rows']
    if cache:
        raw=sources(cache,download)
        for a,b in zip(rows,raw):
            if any(a[k]!=b[k] for k in ['family','row','k','n','length','word_sha256']):raise ValueError('observation/source mismatch')
        if len(raw)!=len(rows):raise ValueError('source row count mismatch')
        rows=raw
    R=region(rows);models=json.loads((HERE/'results/witnesses.json').read_text())
    result={'claim':'strict full-class source-removal gain for a recorded-context mixture',
        'joint_upper_exact':fraction_record(R['joint_upper_q_exact']),
        'confidence':.95,'cell_budget':6661,'group_budget':12,'groups':R['groups'],'models':{}}
    evaluator=GeneralEvaluator(rows) if cache else None
    for name,model in models.items():
        if cache:
            p=evaluator(model['x']);saved=np.array(model['probabilities']);ind=independent_probabilities(rows,model['x'])
            # Add any cross-platform artifact difference to the kernel envelope.
            cert=certify(model['x'],max(r['length'] for r in rows))
            artifact_difference=max(abs(F.from_float(float(a))-F.from_float(float(b))) for a,b in zip(p,saved))
            kernel_bound=F(cert['probability_error_numerator'],cert['probability_error_denominator'])
            if artifact_difference+kernel_bound>ERROR:raise ValueError('saved propagation exceeds common envelope')
            if np.max(abs(ind-p))>float(ERROR):raise ValueError('independent density propagation mismatch')
            print(name,'full-source max independent discrepancy',float(np.max(abs(ind-p))),flush=True)
        result['models'][name]=exact_audit(rows,model,R)
    return result

def main():
    a=argparse.ArgumentParser();a.add_argument('--cache');a.add_argument('--download',action='store_true');a.add_argument('--check',action='store_true');args=a.parse_args()
    if args.download and not args.cache:a.error('--download requires --cache')
    result=run(args.cache,args.download);path=HERE/'results/certificate.json'
    if args.check:
        if json.loads(path.read_text())!=result:raise ValueError('certificate artifact mismatch')
    else:save(path,result)
    print('certified joint upper:',result['joint_upper_exact']['decimal'])
    for name,r in result['models'].items():print(name,r['q_interval'],r['minimum_certified_slack'])

if __name__=='__main__':main()
