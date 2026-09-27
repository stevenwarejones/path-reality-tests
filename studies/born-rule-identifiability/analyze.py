"""End-to-end source reproduction and offline conditional identifiability checks."""
import argparse
import json
import math
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from scipy.stats import norm, binomtest
from br_core import (statistics, gram_from_lower, imaginary_limit, complete_gram, table,
    phase_shift, phase_envelope, cancellation_witness, phase_error_bound, budget,
    hac_covariance, leakage_response_bound)
from br_sources import manifest
from br_io import read, retained, EXCLUSIONS
HERE = Path(__file__).resolve().parent


def clean(x):
    if isinstance(x, dict): return {k: clean(v) for k,v in x.items()}
    if isinstance(x, (list, tuple, np.ndarray)): return [clean(v) for v in x]
    if isinstance(x, (float, np.floating)):
        if not np.isfinite(x): raise ValueError('Nonfinite result')
        return float(format(float(x), '.12g'))
    if isinstance(x, np.integer): return int(x)
    return x


def dumps(x):
    return json.dumps(clean(x), indent=2, sort_keys=True, allow_nan=False)+'\n'


def compare(a, e, path='root'):
    if isinstance(a, dict):
        if not isinstance(e, dict) or a.keys()!=e.keys(): raise ValueError(path+': keys differ')
        for k in a: compare(a[k], e[k], path+'.'+k)
    elif isinstance(a, list):
        if not isinstance(e, list) or len(a)!=len(e): raise ValueError(path+': lengths differ')
        for i,(x,y) in enumerate(zip(a,e)): compare(x,y,f'{path}[{i}]')
    elif isinstance(a,float) and isinstance(e,(int,float)):
        if not math.isclose(a,e,rel_tol=2e-9,abs_tol=2e-11): raise ValueError(f'{path}: {a} != {e}')
    elif a!=e: raise ValueError(f'{path}: {a!r} != {e!r}')


def source_reproduction(cache):
    out={'evidence_type':'observed','units':'detector volts','source':manifest(),'datasets':{},
         'uncertainty_claim':'Descriptive Bartlett HAC sensitivity, not coverage or physics p-values'}
    for e in manifest()['files']:
        y,inv=read(cache,e)
        if e['role']=='contrast':
            out['contrast_inventory']=inv
            continue
        t=int(e['name'][:2]); yy,ids=retained(y,t); s=statistics(yy)
        values=np.column_stack([s['peres'],s['epsilon_V'],s['denominator_V']]); means=values.mean(0)
        hac={}
        for lag in [0,5,10,20,40]:
            cov=hac_covariance(values,ids,lag)
            gradient=np.array([0,1/means[2],-means[1]/means[2]**2])
            hac[str(lag)]={'F_se':math.sqrt(max(0,cov[0,0])),
                'epsilon_se_V':math.sqrt(max(0,cov[1,1])),
                'kappa_delta_se':math.sqrt(max(0,float(gradient@cov@gradient)))}
        sensitivity={}
        for label,base,cuts in [('no_exclusions',1,False),('zero_based_exclusions',0,True)]:
            alt,_=retained(y,t,base,cuts); ss=statistics(alt)
            sensitivity[label]={'cycles':len(alt),'F':float(ss['peres'].mean()),
                'kappa':float(ss['epsilon_V'].mean()/ss['denominator_V'].mean())}
        ci=ids-ids.mean()
        out['datasets'][str(t)]={'inventory':inv,'retained_cycles':len(yy),'excluded_one_based_ids':EXCLUSIONS[t],
            'mean_voltage_by_mask':yy.mean(0),'cycle_covariance_V2':np.cov(yy,rowvar=False,ddof=1),
            'mean_F':means[0],'epsilon_V':means[1],'kappa_ratio_of_means':means[1]/means[2],
            'mean_cycle_kappa_not_paper_definition':float(np.mean(s['epsilon_V']/s['denominator_V'])),
            'mean_correlations_AB_AC_BC':s['correlations_AB_AC_BC'].mean(0),
            'hac_sensitivity':hac,'exclusion_sensitivity':sensitivity,
            'F_epsilon_denominator_slopes_per_cycle':ci@(values-values.mean(0))/(ci@ci),
            'temperature_status':'absent_in_file' if t==23 else 'housing_voltage_recorded_one_missing',
            'published_summary':{'evidence_type':'published_summary','F':.9553 if t==23 else .9683,
                'F_se':.0004 if t==23 else .0009,'kappa':-.00140 if t==23 else .0011,
                'epsilon_nW':-2.58 if t==23 else 1.5},
            'status':'Kappa agrees at reported precision; small F discrepancy unresolved; voltage gain absent'}
    return clean(out)


def identification(rep):
    out={'evidence_type':'observed_summary_conditional_construction','datasets':{},'joint':[],
         'claim':'Mean-table ambiguity only; no fit of complete time-series distributions',
         'alternative':'Additive all-open intensity theta*trace(G), not a complete probability theory'}
    intervals={}
    for t,d in rep['datasets'].items():
        y=np.array(d['mean_voltage_by_mask']); g=gram_from_lower(y); scale=float(np.trace(g))
        u=.5*imaginary_limit(g); gc=complete_gram(g,u)
        eps=float(y[7]-table(gc,background=y[0])[7]); grid=np.linspace(-.5,.5,2001); roots=[]
        for l,r in zip(grid[:-1],grid[1:]):
            if (phase_shift(g,u,l)-eps)*(phase_shift(g,u,r)-eps)<=0:
                roots.append(brentq(lambda z:phase_shift(g,u,z)-eps,l,r,xtol=1e-14))
        if not roots: raise ValueError('No ordinary phase witness in declared domain')
        delta=min(roots,key=abs)
        ordinary=table(gc,delta=delta,background=y[0]); alt=table(gc,theta=eps,background=y[0])
        regions=[]
        for radius in [0.,.01,.03,.1]:
            lo,hi=phase_envelope(g,radius); interval=[(eps-hi)/scale,(eps-lo)/scale]
            regions.append({'phase_radius_rad':radius,'theta_interval':interval,'contains_zero':interval[0]<=0<=interval[1]})
            intervals.setdefault(radius,[]).append(interval)
        cov=np.array(d['cycle_covariance_V2'])/d['retained_cycles']
        j=np.zeros((8,2));j[7]=[scale,2*u]; w=np.linalg.solve(np.linalg.cholesky(cov),j)
        proj=w[:,0]-w[:,1]*(w[:,1]@w[:,0])/(w[:,1]@w[:,1])
        out['datasets'][t]={'gram_real':g,'gram_eigenvalues':np.linalg.eigvalsh(gc),
            'imaginary_AB_V':u,'imaginary_AB_limit_V':imaginary_limit(g),'ordinary_phase_rad':delta,
            'alternative_theta':eps/scale,'intensity_scale_V':scale,
            'max_table_residual_V':float(max(abs(ordinary-y).max(),abs(alt-y).max())),
            'calibration_regions':regions,'local_whitened_rank':int(np.linalg.matrix_rank(w,tol=1e-9)),
            'projected_theta_relative_norm':float(np.linalg.norm(proj)/np.linalg.norm(w[:,0])),
            'added_witness_ordinary_V':cancellation_witness(gc,delta),
            'added_witness_alternative_V':cancellation_witness(gc,0.,eps)}
    for radius,iv in intervals.items():
        lo,hi=max(a[0] for a in iv),min(a[1] for a in iv)
        out['joint'].append({'phase_radius_rad':radius,'theta_interval':[lo,hi],
            'nonempty':lo<=hi,'contains_zero':lo<=0<=hi,
            'evidence_type':'conditional_feasibility_not_confidence_interval',
            'shared':'theta only; independent coherence and phase nuisance per temperature'})
    return clean(out)


def prospective(rep):
    rows=[]
    for t,d in rep['datasets'].items():
        g=gram_from_lower(np.array(d['mean_voltage_by_mask'])); signal=.002*float(np.trace(g))
        for error in [0.,.001,.005]:
            for r in [0.,.0001]:
                bound=phase_error_bound(g,error)+4*r; gap=signal-2*bound
                rows.append({'temperature_reference':t,'theta_target':.002,'signal_V':signal,
                    'phase_error_rad':error,'per_reading_systematic_V':r,'witness_bound_V':bound,
                    'separation_gap_V':gap,'budget':budget(gap)})
    return clean({'evidence_type':'prospective_design','alpha':.01,'power_lower_bound':.9,
        'voltage_range_V':1.,'assumptions':'Independent bounded trials and external deterministic calibration; no achieved power','rows':rows})


def simulations(rep):
    rng=np.random.default_rng(260926); nrep=4000
    g=np.array([[.18,-.10,.12],[-.10,.36,-.08],[.12,-.08,.16]])
    gc=complete_gram(g,.5*imaginary_limit(g)); sigma=.003;n=200;se=sigma*math.sqrt(3.5/n)
    z=norm.ppf(1-.01/2);bound=leakage_response_bound(gc,.014,.002,1.02);rows=[]
    cases=[('ordinary_phase',.03,0.,0j,0.,0.),('ordinary_large_phase',-.2,0.,0j,0.,0.),
        ('injected_intensity',.03,.002,0j,0.,0.),
        ('leakage_nonlinearity_unprotected',.03,0.,.014j,.002,0.),
        ('leakage_nonlinearity_protected',.03,0.,.014j,.002,bound),
        ('injected_with_conservative_envelope',.03,.002,.014j,.002,bound)]
    cases=[(*case,np.full(5,sigma),n) for case in cases]
    for temperature,data in rep['datasets'].items():
        sd=np.sqrt(np.diag(data['cycle_covariance_V2']))[[7,7,1,6,0]]
        cases.append((f'archive_scale_{temperature}_injection',.03,abs(data['epsilon_V']),0j,0.,0.,sd,data['retained_cycles']))
    for label,delta,theta,leak,quad,b,noise_sd,mean_n in cases:
        se=float(np.sqrt(np.sum(np.array([.5,.5,-1,-1,1])**2*noise_sd**2)/mean_n))
        rejected=0
        for _ in range(nrep):
            scale=1+rng.uniform(-.02,.02);phase=delta+rng.uniform(-.02,.02)
            y=table(gc*scale,phase,theta,background=.01,leakage=leak,quadratic=quad)
            yp=table(gc*scale,phase,theta,phase=math.pi,background=.01,leakage=leak,quadratic=quad)
            raw=np.array([y[7],yp[7],y[1],y[6],y[0]])+rng.normal(0,noise_sd/math.sqrt(mean_n),5)
            rejected+=abs(raw@np.array([.5,.5,-1,-1,1]))>b+z*se
        ci=binomtest(rejected,nrep).proportion_ci(.95,method='exact')
        rows.append({'scenario':label,'rejected':int(rejected),'replications':nrep,
            'binomial_95_interval':[ci.low,ci.high],'systematic_bound_V':b,
            'readings_per_setting':mean_n,'setting_sd_V':noise_sd,'injected_intensity_V':theta})
    return clean({'evidence_type':'simulation','seed':260926,'rows':rows,
        'noise':'known Gaussian sample means with common block power/phase drift',
        'decision':'two-sided alpha=.01 normal threshold plus analytic synthetic calibration envelope',
        'scope':'Archive-scale scenarios use measured cycle SDs, retained cycle counts and absolute I3 as injection scale, assuming independent setting noise and equal all-open SD at the unmeasured antipode. Not validated temporal coverage or achieved sensitivity.'})


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--cache',type=Path)
    p.add_argument('--output-dir',type=Path,default=HERE/'results');p.add_argument('--check',action='store_true');a=p.parse_args()
    rep=source_reproduction(a.cache) if a.cache else json.loads((HERE/'results/reproduction.json').read_text())
    if rep['source']!=manifest(): raise ValueError('Manifest mismatch')
    out={'reproduction.json':rep,'identifiability.json':identification(rep),
         'prospective.json':prospective(rep),'simulations.json':simulations(rep)}
    from br_report import render
    report=render(out)
    if a.check:
        for name,data in out.items(): compare(data,json.loads((HERE/'results'/name).read_text()),name)
        if report!=(HERE/'results/report.md').read_text(): raise ValueError('Report differs')
        print('Verified full source reproduction' if a.cache else 'Verified offline derived diagnostics; raw check requires --cache')
    else:
        a.output_dir.mkdir(parents=True,exist_ok=True)
        for name,data in out.items(): (a.output_dir/name).write_text(dumps(data))
        (a.output_dir/'report.md').write_text(report)
        print('Wrote',a.output_dir)
if __name__=='__main__': main()
