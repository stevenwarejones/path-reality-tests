"""Synthetic generalized-equivalence controls on the actual recorded geometry."""
import argparse,json
import numpy as np
from sqd_sources import HERE,sources
from sqd_models import PAULI,IDEAL,Evaluator
from sqd_general import unitary
from sqd_transfer import amplitude_damping,matrix,qutrit_kraus,transform,effect_transform,lower_kappa,similarity,TAU
from format_artifacts import formatted


def setup():
    loss=[.001,.002,.0015];ret=[.01,.012,.008];damping=[.0001,.0002,.0003]
    E=np.diag([.02,.98]);z=.4;kappa=.6
    old=[];new=[];locked=[];kk=[];kp=[];tp=[];lower=[]
    for g in range(3):
        ks=[K@unitary(IDEAL[g]) for K in amplitude_damping(damping[g])]
        ks2,l2,e2,t2=transform(ks,loss[g],ret[g],kappa)
        old.append(matrix(ks,loss[g],ret[g]));new.append(matrix(ks2,l2,e2,t2));locked.append(matrix(ks2,l2,e2,TAU))
        kk.append(qutrit_kraus(ks,loss[g],ret[g]));kp.append(qutrit_kraus(ks2,l2,e2,t2));tp.append(t2)
        lower.append(lower_kappa(ks,loss[g],ret[g],E,z))
    state=np.array([1.,0,0,1.,0]);effect=np.array([.5,0,0,-.48,z]);effect2=effect.copy();effect2[4]=effect_transform(E,z,kappa)
    return np.array(old),np.array(new),np.array(locked),state,effect,effect2,kk,kp,tp,lower


def probabilities(rows,gates,state,effect):
    ev=Evaluator(rows);p=np.empty(len(rows));ev.lib.propagate(ev.offsets,ev.words,len(rows),np.array([gates,gates]),state,effect,0,p);return p


def calculate(rows):
    old,new,locked,state,E,E2,kk,kp,tp,lower=setup();T=similarity(.6)
    p=probabilities(rows,old,state,E);q=probabilities(rows,new,state,E2);r=probabilities(rows,locked,state,E2)
    examples=[]
    for family in ['GST','RB']:
        ids=[i for i,row in enumerate(rows) if row['family']==family];i=max(ids,key=lambda j:abs(p[j]-r[j]));row=rows[i]
        examples.append({**row,'original_probability':float(p[i]),'free_return_probability':float(q[i]),'locked_return_probability':float(r[i]),'absolute_break':float(abs(p[i]-r[i]))})
    detectability={}
    # This is a rigorous simple-pair power ceiling under independent binomial
    # exposures, not a claim that either synthetic model fits the acquisition.
    for family in ['GST','RB','joint']:
        mask=np.array([family=='joint' or row['family']==family for row in rows]);nn=np.array([row['n'] for row in rows])[mask];pp=p[mask];qq=r[mask]
        kl=float(np.sum(nn*(pp*np.log(pp/qq)+(1-pp)*np.log((1-pp)/(1-qq)))))
        detectability[family]={'independent_binomial_KL_original_vs_locked':kl,'level_005_power_upper_by_pinsker':min(1.,.05+np.sqrt(max(0,kl)/2))}
    tp_error=max(float(np.max(abs(sum(K.conj().T@K for K in ks)-np.eye(3)))) for ks in kk+kp)
    cp_min=min(float(np.linalg.eigvalsh(sum(np.outer(K.reshape(-1),K.reshape(-1).conj()) for K in ks)).min()) for ks in kk+kp)
    return {'status':'synthetic ordinary controls evaluated on measured sequence geometry; not fits to measured probabilities',
        'locked_return_detectability':detectability,'kappa':.6,'loss':[.001,.002,.0015],'return':[.01,.012,.008],'damping':[.0001,.0002,.0003],
        'common_kappa_lower':max(lower),'transformed_return_state_min_eigenvalue':min(float(np.linalg.eigvalsh(t).min()) for t in tp),
        'trace_preservation_max_error':tp_error,'choi_min_eigenvalue':cp_min,
        'similarity_max_error':float(max(np.max(abs(B-T@A@np.linalg.inv(T))) for A,B in zip(old,new))),
        'all_recorded_words_max_error':float(np.max(abs(p-q))),'locked_return_break_examples':examples,
        'monitor_one_gate_original':.001,'monitor_one_gate_transformed':.0006}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--cache',required=True);ap.add_argument('--check',action='store_true');a=ap.parse_args()
    result=calculate(sources(a.cache));path=HERE/'results/transfer-audit.json'
    if a.check:
        from certificates import close
        close(result,json.loads(path.read_text()))
    else:path.write_text(formatted(result)+'\n')
    print({k:v for k,v in result.items() if k!='locked_return_break_examples'})

if __name__=='__main__':main()
