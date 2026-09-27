"""Shared-nuisance interval test, without fitting a different scale per setting.

Every discarded parameter cell has an enclosing interval witness. Exhausting
the cell budget is INCONCLUSIVE, never evidence of overlap or rejection.
"""
from math import isfinite
from itertools import product
import numpy as np
from model import Ring, phase_gap, calibration_phase_bias
from decision import DEFAULT_SETTINGS, exclude_candidates, calibrated_outcome_box


def _phases(ring,sites,settings):
    ec=np.array([t*float(ring.continuum(j))/ring.hbar for j,t,q in settings])
    ea=np.array([t*float(ring.lattice(j,sites))/ring.hbar for j,t,q in settings])
    q=np.array([q for j,t,q in settings])
    delta=np.array([phase_gap(ring,sites,0,j,t) for j,t,q in settings])
    return ec,ea,q,delta


def _boxes(center,coeff,cell,cal,contamination=0.):
    """Two shared parameters: additive phase and fractional kinetic scale."""
    mid=np.mean(cell,axis=1); rad=(cell[:,1]-cell[:,0])/2
    phases=center+mid[0]+coeff*mid[1]
    radii=rad[0]+np.abs(coeff)*rad[1]
    # Pad arithmetic endpoints outward; ambiguous sub-picoprobability gaps are
    # never certified by this floating implementation.
    pairs=[calibrated_outcome_box(float(p),float(r)+1e-13,cal,contamination)
           for p,r in zip(phases,radii)]
    return np.array([p[0] for p in pairs])-1e-13,np.array([p[1] for p in pairs])+1e-13


def _cover(initial,enclosure,weights,threshold,max_cells):
    """Prove every cell has a separating coordinate; retain unresolved cells."""
    if isinstance(max_cells,bool) or not isinstance(max_cells,int) or max_cells<1:
        raise ValueError('positive integer cell budget required')
    pending=[initial]; examined=0; witnesses=[]
    cuts=[set(map(float,row)) for row in initial]
    while pending and examined<max_cells:
        cell=pending.pop(); examined+=1
        gap,index=enclosure(cell)
        if gap>threshold+1e-12:
            witnesses.append((float(gap),int(index)))
            continue
        widths=cell[:,1]-cell[:,0]
        axis=int(np.argmax(widths*weights))
        mid=(cell[axis,0]+cell[axis,1])/2
        if not cell[axis,0]<mid<cell[axis,1]:
            return dict(certified=False,cells=examined,reason='unresolved-point',witnesses=witnesses,cuts=[sorted(c) for c in cuts])
        cuts[axis].add(float(mid))
        low=cell.copy(); high=cell.copy(); low[axis,1]=mid; high[axis,0]=mid
        pending.extend((high,low))
    return dict(certified=not pending,cells=examined,
                reason='covered' if not pending else 'cell-budget',witnesses=witnesses,cuts=[sorted(c) for c in cuts])


def exclude_shared_candidates(calibration,main_counts,candidates,*,max_cells=4096,partitions=None,**kwargs):
    """Executable shared-scale/phase refinement of the SAME simultaneous test.

    Rejection requires covering the entire calibrated shared nuisance set. No
    extra alpha is spent on cell subdivision, settings or candidate choices.
    Complete failures and strict confidence thresholds are preserved.
    """
    baseline=exclude_candidates(calibration,main_counts,candidates,**kwargs)
    if baseline['status']!='evaluated': return baseline
    ring=kwargs.get('ring') or Ring(); settings=kwargs.get('settings',DEFAULT_SETTINGS)
    offset=kwargs.get('phase_offset',0.); scale=kwargs.get('fractional_scale',0.)
    contam=kwargs.get('contamination',0.)
    observed=np.array([np.array(row)/sum(row) for row in main_counts])
    initial=np.array([[-offset,offset],[-scale,scale]])
    for model in baseline['models']:
        if model['rejected']:
            model['shared_test']=dict(certified=True,cells=0,reason='coordinate-witness')
            continue
        ec,ea,q,delta=_phases(ring,model['sites'],settings)
        radii=np.array([b['radius'] for b in model['boxes']])[:,None]
        def enclosure(cell):
            lo,hi=_boxes(delta-q,ea,cell,baseline['calibration'],contam)
            gaps=np.maximum(observed-radii-hi,lo-observed-radii)
            index=int(np.argmax(gaps))
            return float(gaps.flat[index]),index
        if partitions is not None and model['sites'] in partitions:
            cuts=partitions[model['sites']]
            if len(cuts)!=2: raise ValueError('two shared parameter partitions required')
            axes=[]
            for values,bounds in zip(cuts,initial):
                values=np.asarray(values,float)
                if len(values)<1 or not np.isfinite(values).all() or values[0]!=bounds[0] or values[-1]!=bounds[1] or (len(values)>1 and not np.all(np.diff(values)>0)):
                    raise ValueError('partition must exactly cover calibrated parameter bounds')
                axes.append([(values[0],values[0])] if len(values)==1 else list(zip(values[:-1],values[1:])))
            checks=[enclosure(np.array(cell)) for cell in product(*axes)]
            result=dict(certified=all(g>1e-12 for g,i in checks),cells=len(checks),
                        reason='predeclared-partition',witnesses=[(float(g),int(i)) for g,i in checks])
        else:
            result=_cover(initial,enclosure,np.array([1.,max(np.max(np.abs(ea)),1.)]),0.,max_cells)
        model['shared_test']=result; model['rejected']=result['certified']
    baseline['rejected']=[m['sites'] for m in baseline['models'] if m['rejected']]
    baseline['retained']=[m['sites'] for m in baseline['models'] if not m['rejected']]
    baseline['scope']+='; one common phase offset and kinetic scale across settings'
    return baseline


def certify_shared_gap(sites,target_gap,*,ring=None,settings=DEFAULT_SETTINGS,
                       eta=.8,visibility=.9,calibration_radius=.001,
                       phase_offset=0.,fractional_scale=.001,max_cells=4096,calibration_phase_radius=0.):
    """Prospective uniform gap, comparing BOTH shared nuisance families.

    Four-dimensional interval coverage: (offset_c,scale_c,offset_a,scale_a).
    The calibration envelope includes center drift as in design.py. A positive
    certificate supplies power using the existing radius/budget arithmetic.
    The observed-count shared test uses the same cells and may refine further.
    """
    if not isfinite(target_gap) or target_gap<0: raise ValueError('nonnegative finite gap required')
    if not 0<=eta<=1 or not 0<=visibility<=1: raise ValueError('physical efficiency and visibility required')
    if any(not isfinite(x) or x<0 for x in (calibration_radius,phase_offset,fractional_scale,calibration_phase_radius)):
        raise ValueError('nonnegative finite certificates required')
    ring=Ring() if ring is None else ring
    ec,ea,q,delta=_phases(ring,sites,settings)
    # Rectangular enlargement of the common calibration region is conservative.
    contrast_radius=4*calibration_radius+2*calibration_phase_bias(calibration_phase_radius)
    cal=dict(efficiency=[max(0,eta-2*calibration_radius),min(1,eta+2*calibration_radius)],
             contrast=[max(0,eta*visibility-contrast_radius),min(1,eta*visibility+contrast_radius)],empty=False)
    initial=np.array([[-phase_offset,phase_offset],[-fractional_scale,fractional_scale]]*2)
    def enclosure(cell):
        lc,hc=_boxes(-q,ec,cell[:2],cal)
        la,ha=_boxes(delta-q,ea,cell[2:],cal)
        gaps=np.maximum(lc-ha,la-hc); index=int(np.argmax(gaps))
        return float(gaps.flat[index]),index
    weights=np.array([1,max(np.max(np.abs(ec)),1),1,max(np.max(np.abs(ea)),1)])
    result=_cover(initial,enclosure,weights,target_gap,max_cells)
    result['certified_gap']=target_gap if result['certified'] else None
    result['decision_partition']=result.pop('cuts')[2:] if result['certified'] else None
    return result
