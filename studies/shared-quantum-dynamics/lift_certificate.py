"""Exact leakage/memory lifts attached to the compatible all-length qubit fit."""
import argparse,json,hashlib
import numpy as np
from sqd_sources import HERE,sources
from sqd_general import operations,GeneralEvaluator
from sqd_models import PAULI,Evaluator
from sqd_transfer import lift_qubit,matrix,qutrit_kraus,TAU,apply
from transfer_audit import probabilities
from format_artifacts import formatted


def calculate(rows=None):
    fit=json.loads((HERE/'results/general-fits.json').read_text())['all_joint'];gg,initial,effect,channels=operations(fit['x']);comp=json.loads((HERE/'results/general-compatibility.json').read_text());assert comp['compatible']
    budgets=[];minimum_returns=[]
    for kk in channels:
        J=sum(np.outer(K.reshape(-1,order='F'),K.reshape(-1,order='F').conj()) for K in kk)
        budgets.append(float(2*np.linalg.eigvalsh(J).min()))
        minimum_returns.append(float(np.linalg.norm([np.trace(P@apply(kk,TAU)).real for P in PAULI[1:]])))
    # Fixed interior fractions; no measured residual is optimized here.
    loss=np.array(budgets)*.8;ret=np.maximum(2*np.array(minimum_returns),1e-5)
    if np.any(ret>1):raise ValueError('chosen interior return multiplier unavailable')
    Q=np.eye(5)[:4];Q[0,4]=1.;models=[]
    for fraction in [0.,.25,1.]:
        matrices=[];cp=[];tp=[];rt=[]
        for g,kk in enumerate(channels):
            base,tau=lift_qubit(kk,float(loss[g]*fraction),float(ret[g]));G=matrix(base,float(loss[g]*fraction),float(ret[g]),tau);matrices.append(G)
            ks=qutrit_kraus(base,float(loss[g]*fraction),float(ret[g]),tau)
            J=sum(np.outer(K.reshape(-1,order='F'),K.reshape(-1,order='F').conj()) for K in ks)
            cp.append(float(np.linalg.eigvalsh(J).min()));tp.append(float(np.max(abs(sum(K.conj().T@K for K in ks)-np.eye(3)))));rt.append(float(np.linalg.eigvalsh(tau).min()))
        effect2=effect.copy();effect2[4]=effect[0]
        # General model has zero leakage initially. Its Bloch coordinates and
        # arbitrary binary effect are shared exactly by the lifted model.
        err=max(float(np.max(abs(Q@G-gg[0,g,:4,:4]@Q))) for g,G in enumerate(matrices))
        populations={}
        for g in range(3):
            lam=loss[g]*fraction;eta=ret[g]
            populations[str(g)]={str(L):float(lam/(lam+eta)*(1-(1-lam-eta)**L)) for L in [1,1000,8192]}
        item={'loss':(loss*fraction).tolist(),'return':ret.tolist(),'fraction_of_selected_loss':fraction,
            'intertwining_max_error':err,'choi_min_eigenvalue':min(cp),'trace_preservation_max_error':max(tp),
            'return_state_min_eigenvalue':min(rt),'repeated_gate_leakage_populations':populations,
            'leakage_effect':float(effect2[4]),'same_all_word_probabilities':True}
        if rows is not None:
            actual=probabilities(rows,np.array(matrices),initial,effect2)
            item['full_source_max_probability_error']=float(np.max(abs(actual-np.array(fit['probabilities']))))
        models.append(item)
    observations=json.loads((HERE/'results/observations.json').read_text());key=hashlib.sha256(bytes([1])*8192).hexdigest()
    index=next(i for i,r in enumerate(observations) if r['word_sha256']==key)
    witness={**observations[index],'word':'Gx repeated 8192 times','common_probability':fit['probabilities'][index],
        'conditional_leakage_populations':[m['repeated_gate_leakage_populations']['1']['8192'] for m in models]}
    return {'recorded_example':witness,'status':'exact model-class construction, numerically evaluated at the full-acquisition compatible qubit point; populations are conditional examples, not device estimates',
        'computational_subspace':'stipulated block, using the saved qubit representation; channel budget is gauge-conditional',
        'CP_loss_upper_per_gate':budgets,'positivity_return_lower_per_gate':minimum_returns,'models':models,
        'same_joint_region_minimum_slack':comp['minimum_all_constraint_slack'],
        'missing_observation':'calibrated subspace population monitor after one gate; unknown terminal binary words cannot distinguish these models'}


if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--cache');a.add_argument('--check',action='store_true');args=a.parse_args();r=calculate(sources(args.cache) if args.cache else None);path=HERE/'results/lift-certificate.json'
    if args.check:
        from certificates import close
        old=json.loads(path.read_text())
        if not args.cache:
            for m in old['models']:m.pop('full_source_max_probability_error',None)
        close(r,old)
    else:path.write_text(formatted(r)+'\n')
    print({'CP_loss_upper':r['CP_loss_upper_per_gate'],'return_lower':r['positivity_return_lower_per_gate'],'models':r['models']})
