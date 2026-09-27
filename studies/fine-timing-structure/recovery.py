#!/usr/bin/env python3
"""Paired empirical-background recovery at frozen resolutions; not hardware power."""
import argparse
import json
import numpy as np
from scipy.stats import beta
from common import HERE,SIDES,LEVELS,MAPS,CELLS,PHASE_CELLS,PEAKS,SEED,REPETITIONS,grouped,json_text,artifact_json
from inference import cdf,partitions


def rate(k,n=REPETITIONS):
    return dict(count=int(k),repetitions=n,rate=k/n,interval95=[
        0. if k==0 else float(beta.ppf(.025,k,n-k+1)),
        1. if k==n else float(beta.ppf(.975,k+1,n-k))])


def transform(p,side,mode,amplitude,x):
    out=p.copy()
    if mode=='within':
        core=MAPS[side][0]==1
        mass=p[core].sum(); out[core]*=1-amplitude
        out[PHASE_CELLS+PEAKS[side]+x]+=amplitude*mass
    elif mode=='translation' and x:
        out[:]=0
        cell=np.arange(CELLS); shifted=cell//PHASE_CELLS*PHASE_CELLS+np.minimum(PHASE_CELLS-1,cell%PHASE_CELLS+int(amplitude))
        np.add.at(out,shifted,p)
    elif mode=='rate' and x: out*=1+amplitude
    if out.sum()>1 or np.any(out<0): raise ValueError('invalid synthetic probability')
    return out


def truths(p0,p1,side):
    d=p1-p0; incidence=float(d.sum())
    tv=[.5*(abs(incidence)+float(np.abs(grouped(d,m)).sum())) for m in MAPS[side]]
    return dict(k_grid=float(np.abs(np.cumsum(d)).max()),incidence=incidence,tv=tv,
                gain=[v-tv[0] for v in tv])


def evaluate(counts,u,n,side,epsilon=0.):
    a=cdf(counts,u,n,epsilon); b=partitions(counts,u,n,side,epsilon)
    results=dict(cdf=rate(int((a['k_grid'][:,0]>0).sum())),
        tv={level:rate(int((b['tv'][:,j,0]>0).sum())) for j,level in enumerate(LEVELS)},
        gain={level:rate(int((b['gain'][:,j,0]>0).sum())) for j,level in enumerate(LEVELS)},
        median_k_upper=float(np.median(a['k_grid'][:,1])),
        median_tv_upper=np.median(b['tv'][:,:,1],axis=0).tolist())
    return results


def run():
    rng=np.random.default_rng(SEED); cases=[]
    for side in SIDES:
        source=json.loads((HERE/f'{side}-counts.json').read_text())
        blockcounts=np.array(source['pooled_block_counts']); exposures=np.array(source['pooled_block_exposures'])
        unknown=np.array(source['unknown']).sum(axis=0)
        for y in (0,1):
            laws=blockcounts[:,y]/exposures[:,y,None]
            for scale in (.25,1.,4.):
                n=int(round(source['interior_rows']*scale))
                sizes=np.rint(exposures.sum(axis=1)*scale).astype('i8')
                weights=sizes/n
                for mode,amps in [('within',(0.,.1,.25,.5,1.)),('translation',(1.,2.)),('rate',(.1,.5))]:
                    for amplitude in amps:
                        counts=np.zeros((REPETITIONS,2,CELLS),dtype='i8')
                        average=np.zeros((2,CELLS))
                        for block,p in enumerate(laws):
                            arm=np.array([transform(p,side,mode,amplitude,x) for x in (0,1)])
                            average+=weights[block]*arm
                            probabilities=np.r_[.25*arm.ravel(),1-.25*arm.sum()]
                            counts+=rng.multinomial(int(sizes[block]),probabilities,size=REPETITIONS)[:,:-1].reshape(REPETITIONS,2,CELLS)
                        truth=truths(*average,side)
                        boundary_mass=float(sum(w*p[np.arange(CELLS)%PHASE_CELLS > PHASE_CELLS-1-int(amplitude)].sum() for w,p in zip(weights,laws))) if mode=='translation' else 0.
                        for envelope in ('no_missing','event_supported'):
                            u=np.zeros((REPETITIONS,2),dtype='i8') if envelope=='no_missing' else np.broadcast_to(np.rint(unknown[:,y,1]*scale).astype('i8'),(REPETITIONS,2))
                            cases.append(dict(receiver=side,local_setting=y,scale=scale,n=n,mode=mode,
                                amplitude=amplitude,envelope=envelope,truth=truth,saturated_boundary_mass=boundary_mass,
                                uncertainty=u[0].tolist(),outcomes=evaluate(counts,u,n,side)))
    controls=memory_controls(rng)
    return dict(schema_version=1,seed=SEED,repetitions=REPETITIONS,cases=cases,controls=controls,
        model='Independent multinomial trusted counts within 16 pooled chronological blocks; fixed archive completion budgets; not a temporal resampling model',
        monte_carlo_precision='Worst-case SE 0.025; exact pointwise 95% binomial intervals; no tail guarantee inferred from simulation')


def memory_controls(rng):
    """Explicit row models, separate from archive histogram power simulations."""
    n=200_000; side='alice'; allcounts={k:np.zeros((REPETITIONS,2,CELLS),dtype='i8') for k in ('fresh_memory','persistent_memory','misalignment','selective_missing')}
    missing=np.zeros((REPETITIONS,2),dtype='i8')
    for rep in range(REPETITIONS):
        x=rng.integers(0,2,n+1); local=rng.integers(0,2,n)==0
        event=rng.random(n)<.02
        phase=PHASE_CELLS+PEAKS[side]+x[:-1]
        def tally(labels,coords,keep):
            return np.bincount(labels[keep]*CELLS+coords[keep],minlength=2*CELLS).reshape(2,CELLS)
        allcounts['fresh_memory'][rep]=tally(x[1:],phase,local & event)
        allcounts['misalignment'][rep]=tally(x[:-1],phase,local & event)
        flips=rng.random(n)<.1
        persistent=np.r_[rng.integers(0,2),np.zeros(n,dtype='i8')]
        persistent[1:]=(persistent[0]+np.cumsum(flips))%2
        persistent_phase=PHASE_CELLS+PEAKS[side]+persistent[:-1]
        allcounts['persistent_memory'][rep]=tally(persistent[1:],persistent_phase,local & event)
        # Current-outcome-dependent removal is allowed by pathwise count completion.
        removed=local & event & (x[1:]==x[:-1])
        allcounts['selective_missing'][rep]=tally(x[1:],phase,local & event & ~removed)
        missing[rep]=np.bincount(x[1:][removed],minlength=2)
    result=[]
    for name,counts in allcounts.items():
        for epsilon in ((0.,.2) if name=='persistent_memory' else (0.,)):
            u=missing if name=='selective_missing' else np.zeros((REPETITIONS,2),dtype='i8')
            result.append(dict(model=name,n=n,epsilon=epsilon,outcomes=evaluate(counts,u,n,side,epsilon)))
    # Unsafe deletion is a deliberate assumption violation, not an inferential option.
    result.append(dict(model='selective_missing_wrongly_ignored',n=n,epsilon=0.,outcomes=evaluate(allcounts['selective_missing'],np.zeros_like(missing),n,side)))
    return result


def report(r):
    lines=['# Paired recovery and assumption stress tests','',
        f'{r["repetitions"]} repetitions per case; seed {r["seed"]}. Exact pointwise 95% binomial intervals are in recovery-results.json.',
        'Worst-case Monte Carlo SE is 2.5 percentage points. Zero/400 has a 95% upper',
        'limit of about 0.92%; it is not proof of the analytic 1% family guarantee.', '',
        'The empirical background pools remote labels within 16 chronological blocks',
        'and local settings. Events are sampled independently within blocks. Drift is',
        'retained only through block rates/shapes; arbitrary temporal memory is NOT',
        'recreated. Trusted sample size and fixed uncertain-event candidate budgets',
        'match the archive (scaled with N). The no-missing control uses the same trusted',
        'counts but removes the completion penalty; it does not add missing events.',
        'The remaining N minus trusted rows are no-selected-event padding in this model;',
        'truth values use trusted block sizes divided by full N. Fixed completion budgets',
        'are conservative sensitivity penalties, not a fitted physical missingness law.',
        'Truth values describe the specified pooled law, not an estimated real effect.', '',
        'Every resolution sees the identical simulated realization. Within-core moves',
        'preserve coarse mass exactly. Translation shifts phase cells within each pulse',
        'and saturates at 161; event-rate changes preserve conditional timing shape.',
        'The mass that would cross that upper boundary is recorded per translation case.', '',
        '## Sparse archive-size within-core redistributions','',
        'Event-supported completion, original N. Counts are recoveries /400; J is',
        'positive resolved TV gain beyond coarse. All sizes, amplitudes and resolutions',
        'are retained in the machine-readable result.','',
        '| Receiver/y | Moved core fraction | True fine TV (ppm) | Coarse | CDF | width8 | width2 | width1 | Positive J width1 |',
        '|---|---:|---:|---:|---:|---:|---:|---:|---:|']
    for v in r['cases']:
        if v['mode']=='within' and v['scale']==1 and v['envelope']=='event_supported':
            o=v['outcomes']; k=lambda key:o['tv'][key]['count']
            lines.append(f'| {v["receiver"]}/{v["local_setting"]} | {v["amplitude"]} | {v["truth"]["tv"][-1]*1e6:.1f} | {k("coarse")} | {o["cdf"]["count"]} | {k("width8")} | {k("width2")} | {k("width1")} | {o["gain"]["width1"]["count"]} |')
    lines+=['','## Explicit row-model controls','',
        'N=200,000, event probability .02, random local setting; previous remote bit',
        'selects one of the two core cells. This is a denser control, not archive power.',
        'Fresh assignments are independent of history. Persistent assignments stay',
        'with probability .9; their joint per-history propensity needs epsilon=.2.',
        'Misalignment deliberately substitutes the previous causal bit as current.',
        'Selective missingness removes events whose current and previous bits match;',
        'completion budgets restore a valid enclosure without predictable missingness.','',
        '| Model | Epsilon | CDF /400 | Coarse /400 | width1 /400 |',
        '|---|---:|---:|---:|---:|']
    for v in r['controls']:
        o=v['outcomes']; lines.append(f'| {v["model"]} | {v["epsilon"]} | {o["cdf"]["count"]} | {o["tv"]["coarse"]["count"]} | {o["tv"]["width1"]["count"]} |')
    lines+=['','Rejecting after misalignment, ignoring selective missingness, or falsely',
        'assuming fresh bits is an assumption failure, not evidence of unusual physics.',
        'The ideal-quantizer witness remains exactly invisible to every digital method;',
        'it is verified algebraically in tests, not assigned a fictitious power curve.','',
        '![Paired recovery](figures/recovery.svg)','']
    return '\n'.join(lines)


def main():
    p=argparse.ArgumentParser(); p.add_argument('--check',action='store_true'); a=p.parse_args()
    r=run()
    for name,text in [('recovery-results.json',artifact_json(r)),('recovery-report.md',report(r))]:
        if a.check:
            if (HERE/name).read_text()!=text: raise SystemExit(name+' differs')
        else: (HERE/name).write_text(text)
    print('fine recovery '+('verified' if a.check else 'written'))
if __name__=='__main__': main()
