"""Audit new external portfolios and reproduce complementary-observable feasibility."""
import argparse,csv,hashlib,io,json,re,zipfile
from pathlib import Path
import numpy as np
from scipy.sparse import load_npz
from br_sources import verify
from br_combinations import models,feasibility
# Avoid importing another study's generic analyze module in tests.
HERE=Path(__file__).resolve().parent


def clean(x):
    if isinstance(x,dict):return {k:clean(v) for k,v in x.items()}
    if isinstance(x,(list,tuple,np.ndarray)):return [clean(v) for v in x]
    if isinstance(x,(float,np.floating)):return float(format(float(x),'.12g'))
    if isinstance(x,np.integer):return int(x)
    return x


def portfolio_manifest():return json.loads((HERE/'portfolio-manifest.json').read_text())


def ibm_archive(path):
    with zipfile.ZipFile(path) as z:
        names=[n for n in z.namelist() if not n.startswith('__MACOSX/')]
        meta=json.loads(z.read(next(n for n in names if n.endswith('/data.json'))))
        jobs=sorted([n for n in names if re.search(r'/wyniki_testy_\d+\.csv$',n)],
                    key=lambda n:int(re.search(r'_(\d+)\.csv$',n).group(1)))
        sums={}; job_probs={}; controls=[];rows_total=0;depths_present=True
        for member in jobs:
            rows=list(csv.DictReader(io.StringIO(z.read(member).decode('utf-8-sig'))))
            job={};control=[]
            for r in rows:
                k0,k1=int(r['0']),int(r['1']); n=k0+k1
                if min(k0,k1)<0 or n!=meta['shots']:raise ValueError('Bad counts/shots')
                rows_total+=1
                if r['theta']=='TEST':control.append([k0,k1]);continue
                if 'SX number' not in r:
                    depths_present=False;continue
                depth=int(float(r['SX number']));angle=float(r['theta'])
                key=(depth,angle)
                for dest in [sums,job]:
                    dest.setdefault(key,np.zeros(2,dtype=np.int64));dest[key]+=[k0,k1]
            if len(control)!=2*meta['tests circuits']:raise ValueError('Control count mismatch')
            controls.append(control)
            for key,v in job.items():job_probs.setdefault(key,[]).append(v[0]/v.sum())
        if len(jobs)!=meta['jobs']:raise ValueError('Job metadata mismatch')
        out={'backend':meta['backend'],'start_time':meta['start time'],'jobs':len(jobs),
             'shots_per_circuit':meta['shots'],'circuit_rows':rows_total,
             'evidence_type':'simulation' if 'simulator' in meta['backend'] else 'observed',
             'depth_column_present':depths_present,'control_count_total':np.array(controls).sum(axis=(0,1)).tolist(),
             'control_preparation_success':[
                 float(np.array(controls)[:,0::2,0].sum()/np.array(controls)[:,0::2,:].sum()),
                 float(np.array(controls)[:,1::2,1].sum()/np.array(controls)[:,1::2,:].sum())],
             'depths':{}}
        for depth in sorted({k[0] for k in sums}):
            keys=sorted([k for k in sums if k[0]==depth]); angles=np.array([k[1] for k in keys])
            counts=np.array([sums[k] for k in keys]);prob=counts[:,0]/counts.sum(1)
            X=np.column_stack([np.ones(len(angles)),np.cos(angles),np.sin(angles)])
            fit=np.linalg.lstsq(X,prob,rcond=None)[0];res=prob-X@fit
            out['depths'][str(depth)]={'angles_rad':angles,'aggregate_counts_0_1':counts,
                'harmonic_fit_offset_cos_sin':fit,'residual_max_abs':float(abs(res).max()),
                'residual_rms':float(np.sqrt(np.mean(res*res))),
                'job_level_se':np.array([np.std(job_probs[k],ddof=1)/np.sqrt(len(job_probs[k])) for k in keys])}
        return clean(out)


def photon_archive(path,label):
    with zipfile.ZipFile(path) as z:
        def ld(name):return np.load(io.BytesIO(z.read(label+'/'+name)),allow_pickle=False)
        if label=='8-bin TMD':pn,rn='P_8bin.npy','POVM_8bin.npy'
        else:pn,rn='P_EMG_4+.npy','POVMs_EMG_4+_smoothing_1e-6_D40.npy'
        mu=ld('mpn.npy');p=ld(pn);response=ld(rn); audits=[]
        for name in sorted(n for n in z.namelist() if n.endswith('.npz')):
            index=int(re.search(r'_(\d+)\.npz$',name).group(1)); a=load_npz(io.BytesIO(z.read(name)))
            if label=='8-bin TMD':
                if a.shape[1]!=8 or not np.all(a.data==1):raise ValueError('Invalid click encoding')
                outcomes=np.asarray(a.sum(1)).ravel().astype(int)
            else:
                if a.shape[1]!=1 or not np.equal(a.data,a.data.astype(int)).all():raise ValueError('Invalid event encoding')
                outcomes=np.asarray(a.toarray()).ravel().astype(int)
            if outcomes.min()<0 or outcomes.max()>=p.shape[1]:raise ValueError('Outcome outside range')
            counts=np.bincount(outcomes,minlength=p.shape[1]);freq=counts/counts.sum()
            audits.append({'input_index':index,'mean_photon_label':float(mu[index]),'triggers':len(outcomes),
                          'max_frequency_discrepancy':float(abs(freq-p[index]).max())})
        return clean({'evidence_type':'observed_events_and_derived_tomography','mean_photon_labels':mu,
            'outcome_matrix_shape':list(p.shape),'response_matrix_shape':list(response.shape),
            'response_row_sum_max_error':float(abs(response.sum(1)-1).max()),
            'response_min':float(response.min()),'response_max':float(response.max()),
            'outcome_normalization_error':float(abs(p.sum(1)-1).max()),'event_reconciliation':audits,
            'same_data_calibration_is_independent':False,
            'use_in_physical_joint_fit':False,
            'reason':'POVM is inferred using assumed coherent-state probability model from these probes; not independent calibration of a Born-rule deformation or the waveguide detector'})


def audit(cache):
    m=portfolio_manifest()
    for e in m['files']:verify(cache/e['cache_name'],e)
    return {'source':m,'ibm':{n:ibm_archive(cache/'7538941'/f'{n}.zip') for n in ['lima1-5','lima1-5s','lagos_bench']},
            'photon':{n:photon_archive(cache/'18483031'/f'{n}.zip',n) for n in ['8-bin TMD','PNR SNSPD']}}


def compute(summary):
    hardware=summary['ibm']['lima1-5']
    counts=np.array(hardware['depths']['1']['aggregate_counts_0_1'])
    shots=int(np.median(counts.sum(1)))
    fit=hardware['depths']['1']['harmonic_fit_offset_cos_sin']
    contrast=2*float(np.hypot(fit[1],fit[2]))
    result={'qubit_equivalence':[models(t) for t in [-.1,.01,.1]],
            'qubit_prospective':feasibility(shots,contrast)}
    # Explicit global gauge for photon calibration: a common stochastic latent
    # transformation T can always be absorbed by independent detector responses.
    p=np.array([[.8,.15,.05],[.2,.5,.3]])
    t=np.array([[1,0,0],[.1,.9,0],[0,.1,.9]])
    r1=np.array([[.99,.01],[.4,.6],[.1,.9]])
    r2=np.array([[.95,.04,.01],[.3,.5,.2],[.05,.25,.7]])
    errs=[float(abs((p@t)@r-p@(t@r)).max()) for r in [r1,r2]]
    result['photon_gauge_control']={'evidence_type':'exact_algebra_numerically_checked_synthetic',
        'separate_max_errors':errs,'joint_max_error':max(errs),
        'claim':'For any row-stochastic P,T,R_d, (P T)R_d=P(T R_d). Separate detectors do not break this gauge; T is an ordinary preparation/response change, not a genuine new probability theory.',
        'remove_detector_1_still_ambiguous':True,'remove_detector_2_still_ambiguous':True}
    return clean(result)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--cache',type=Path);p.add_argument('--check',action='store_true');a=p.parse_args()
    summary=audit(a.cache) if a.cache else json.loads((HERE/'results/portfolio-audit.json').read_text())
    if summary['source']!=portfolio_manifest():raise ValueError('Portfolio manifest mismatch')
    outputs={'portfolio-audit.json':summary,'combination-feasibility.json':compute(summary)}
    # Import this file's sibling explicitly: generic module names exist elsewhere.
    import importlib.util
    spec=importlib.util.spec_from_file_location('born_analysis_helpers',HERE/'analyze.py');helper=importlib.util.module_from_spec(spec);spec.loader.exec_module(helper)
    for name,data in outputs.items():
        if a.check:helper.compare(data,json.loads((HERE/'results'/name).read_text()),name)
        else:(HERE/'results'/name).write_text(helper.dumps(data))
    from br_combination_report import render
    report=render(summary,outputs['combination-feasibility.json'])
    target=HERE/'results/combination-report.md'
    if a.check:
        if report!=target.read_text():raise ValueError('Combination report differs')
    else:target.write_text(report)
    print('Verified' if a.check else 'Generated','portfolio audit and two-candidate feasibility')
if __name__=='__main__':main()
