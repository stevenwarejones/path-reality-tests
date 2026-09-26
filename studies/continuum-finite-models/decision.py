"""Executable fixed-allocation candidate exclusion using calibration/main counts.

Input counts include failures. External phase/scale/transport certificates are
explicit parameters, never estimated from a separately refitted null model.
"""
import argparse
import json
from math import isfinite, pi
from pathlib import Path
import numpy as np
from model import Ring, check_sites, phase_gap, cosine_interval, radius, calibration_phase_bias

DEFAULT_SETTINGS=tuple((j,t,q) for j in (1,2,3) for t in (.5,2.,8.,32.) for q in (0.,pi/2))


def _count(value):
    if isinstance(value,bool) or not isinstance(value,(int,np.integer)) or value<0:
        raise ValueError('counts must be nonnegative integers')
    return int(value)


def _binomial(record):
    if len(record)!=2:
        raise ValueError('calibration records are [successes,total]')
    success,total=map(_count,record)
    if total<1 or success>total:
        raise ValueError('calibration needs 0 <= successes <= positive total')
    return success,total


def calibration_box(calibration,alpha=.025,efficiency_transport=0.,contrast_transport=0.):
    """Intersect simultaneous Bernoulli intervals with 0 <= w <= eta <= 1.

    Returns an enclosing rectangle of the physical intersection. Empirical
    centers need not themselves be physical; no point-estimate clipping occurs.
    An empty intersection yields an inconclusive calibration result.
    """
    if set(calibration)!={'detection','phase0','phase_pi'}:
        raise ValueError('three declared calibration strata required')
    if any(not isfinite(x) or x<0 for x in (efficiency_transport,contrast_transport)):
        raise ValueError('finite nonnegative transport certificates required')
    data={key:_binomial(value) for key,value in calibration.items()}
    means={key:x/n for key,(x,n) in data.items()}
    radii={key:radius(n,3,alpha) for key,(_,n) in data.items()}
    e=means['detection']; er=radii['detection']+efficiency_transport
    w=means['phase0']-means['phase_pi']
    wr=radii['phase0']+radii['phase_pi']+contrast_transport
    el,eh=max(0.,e-er),min(1.,e+er)
    wl,wh=max(0.,w-wr),min(eh,w+wr)
    el=max(el,wl)
    empty=el>eh or wl>wh
    return dict(efficiency=[el,eh],contrast=[wl,wh],empty=empty,
                estimates=means,radii=radii,confidence_error=alpha)


def calibrated_outcome_box(phase,phase_radius,cal,contamination=0.):
    if cal['empty']:
        raise ValueError('empty calibration set')
    if not isfinite(phase_radius) or phase_radius<0 or not 0<=contamination<=1:
        raise ValueError('invalid external certificate')
    el,eh=cal['efficiency']; wl,wh=cal['contrast']
    cl,ch=cosine_interval(phase-phase_radius,phase+phase_radius)
    products=[w*c for w in (wl,wh) for c in (cl,ch)]
    lower=np.array([(el+min(products))/2,(el-max(products))/2,1-eh])-contamination
    upper=np.array([(eh+max(products))/2,(eh-min(products))/2,1-el])+contamination
    return np.maximum(0,lower),np.minimum(1,upper)


def outside_enlarged_box(frequency,box,statistical_radius):
    """Strict decision: touching an enlarged endpoint does not reject."""
    p=np.asarray(frequency,float); lo,hi=map(np.asarray,box)
    if p.shape!=(3,) or lo.shape!=(3,) or hi.shape!=(3,):
        raise ValueError('three complete outcome coordinates required')
    if not isfinite(statistical_radius) or statistical_radius<0:
        raise ValueError('nonnegative finite statistical radius required')
    return (p<lo-statistical_radius)|(p>hi+statistical_radius)


def exclude_candidates(calibration,main_counts,candidates,*,ring=None,settings=DEFAULT_SETTINGS,
                       phase_offset=0.,fractional_scale=0.,contamination=0.,
                       efficiency_transport=0.,contrast_transport=0.,
                       alpha_calibration=.025,alpha_main=.025,calibration_phase_radius=0.):
    """Return all rejected integer N with a simultaneous coordinate witness.

    main_counts has one [plus,minus,failure] row per PREDECLARED setting.
    Counts may differ across strata, but allocations must be fixed externally;
    the function does not authorize optional stopping at the observed counts.
    """
    ring=Ring() if ring is None else ring
    if not settings or len(main_counts)!=len(settings):
        raise ValueError('one complete count row per setting required')
    if any(not isfinite(x) or x<0 for x in (phase_offset,fractional_scale,contamination,calibration_phase_radius)) or contamination>1:
        raise ValueError('invalid external phase/scale/contamination certificate')
    rows=[]
    for row in main_counts:
        if len(row)!=3:
            raise ValueError('plus, minus and failure counts required')
        counts=list(map(_count,row)); total=sum(counts)
        if total<1:
            raise ValueError('every predeclared stratum requires positive allocation')
        rows.append((counts,total))
    candidates=list(candidates)
    if len(set(candidates))!=len(candidates):
        raise ValueError('duplicate candidate site counts')
    for n in candidates:
        check_sites(n)
        for j,t,q in settings:
            if not isfinite(t) or not isfinite(q):
                raise ValueError('finite controls required')
            phase_gap(ring,n,0,j,t) # Validates integer modes and strict alias cutoff.
    phase_bias=calibration_phase_bias(calibration_phase_radius)
    cal=calibration_box(calibration,alpha_calibration,efficiency_transport,contrast_transport+phase_bias)
    radii=[radius(total,3*len(settings),alpha_main) for _,total in rows]
    if cal['empty']:
        return dict(status='inconclusive-calibration',rejected=[],retained=candidates,
                    calibration=cal,models=[],alpha_total=alpha_calibration+alpha_main)
    model_results=[]
    for n in candidates:
        witnesses=[]; boxes=[]
        for index,((j,t,q),(counts,total),stat_radius) in enumerate(zip(settings,rows,radii)):
            delta=phase_gap(ring,n,0,j,t)
            phase_uncertainty=phase_offset+fractional_scale*abs(float(t*ring.lattice(j,n)/ring.hbar))
            box=calibrated_outcome_box(delta-q,phase_uncertainty,cal,contamination)
            observed=np.asarray(counts)/total
            rejected=outside_enlarged_box(observed,box,stat_radius)
            boxes.append(dict(lower=box[0].tolist(),upper=box[1].tolist(),radius=stat_radius))
            for outcome in np.flatnonzero(rejected):
                witnesses.append(dict(setting=index,outcome=('plus','minus','failure')[outcome],
                                      frequency=float(observed[outcome]),
                                      lower=float(box[0][outcome]),upper=float(box[1][outcome]),
                                      statistical_radius=stat_radius))
        model_results.append(dict(sites=n,rejected=bool(witnesses),witnesses=witnesses,boxes=boxes))
    return dict(status='evaluated',rejected=[r['sites'] for r in model_results if r['rejected']],
                retained=[r['sites'] for r in model_results if not r['rejected']],calibration=cal,
                models=model_results,alpha_total=alpha_calibration+alpha_main,
                external_certificates=dict(phase_offset=phase_offset,fractional_scale=fractional_scale,
                    contamination=contamination,efficiency_transport=efficiency_transport,
                    contrast_transport=contrast_transport,calibration_phase_radius=calibration_phase_radius),
                scope='fixed independent allocations; conditional on external certificates')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('counts',type=Path,help='JSON with calibration, main_counts, candidates and optional certificates')
    parser.add_argument('--output',type=Path)
    parser.add_argument('--shared',action='store_true',help='enforce one shared phase/scale across settings')
    args=parser.parse_args(); data=json.loads(args.counts.read_text())
    if args.shared:
        from joint import exclude_shared_candidates
        partitions=data.get('partitions')
        partitions=None if partitions is None else {int(n):cuts for n,cuts in partitions.items()}
        result=exclude_shared_candidates(data['calibration'],data['main_counts'],data['candidates'],
            partitions=partitions,**data.get('certificates',{}))
    else:
        result=exclude_candidates(data['calibration'],data['main_counts'],data['candidates'],
                                  **data.get('certificates',{}))
    rendered=json.dumps(result,indent=2,allow_nan=False)+'\n'
    if args.output: args.output.write_text(rendered)
    else: print(rendered,end='')


if __name__=='__main__': main()
