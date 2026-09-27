#!/usr/bin/env python3
"""Finite-sample multinomial aggregate confidence-region audit.

No chi-square approximation, fitted degrees of freedom, or bootstrap is used.
Invert a Cantelli test at each complete probability point. See joint-statistics.md.
"""
import argparse
from decimal import Context, Decimal, ROUND_CEILING, ROUND_FLOOR
import json
from pathlib import Path
import numpy as np
import identifiability as physical
import delay_candidates


def pearson_region(counts, numerators, denominators, alpha=Decimal('0.025')):
    """Directed bounds on Pearson X and its exact multinomial variance.

    Independent rows, each with four nonnegative probabilities and a fixed trial
    count. A complete point can be chosen after observing the data: acceptance
    means membership in an inverted confidence region, not a fitted-null test.
    """
    counts=np.asarray(counts)
    nums=np.broadcast_to(np.asarray(numerators,dtype=object),counts.shape)
    dens=np.broadcast_to(np.asarray(denominators,dtype=object),counts.shape[:-1])
    if counts.shape[-1]!=4 or counts.dtype.kind not in 'iu' or np.any(counts<0):
        raise ValueError('expected four nonnegative integer counts per row')
    down=Context(prec=60,rounding=ROUND_FLOOR)
    up=Context(prec=60,rounding=ROUND_CEILING)
    xl=xu=vl=vu=Decimal(0)
    mean_int=0;impossible=False
    for row,pnums,den in zip(counts.reshape(-1,4),nums.reshape(-1,4),dens.reshape(-1)):
        if type(den) not in (int,np.int64):raise ValueError("noninteger denominator")
        den=int(den);N=int(row.sum())
        if N<=0 or den<=0 or any(type(p) not in (int,np.int64) for p in pnums):
            raise ValueError('invalid rational probability or exposure')
        pnums=[int(p) for p in pnums]
        if min(pnums)<0 or sum(pnums)!=den:
            raise ValueError('probabilities must be nonnegative and normalized')
        support=sum(p>0 for p in pnums);mean_int+=support-1
        leading=2*(support-1);correction=support*support+leading
        reciprocal_l=reciprocal_u=Decimal(0)
        for k,p in zip(row,pnums):
            if p==0:
                impossible=impossible or int(k)>0
                continue
            a=Decimal((int(k)*den-N*p)**2);b=Decimal(N*p*den)
            xl=down.add(xl,down.divide(a,b));xu=up.add(xu,up.divide(a,b))
            reciprocal_l=down.add(reciprocal_l,down.divide(Decimal(den),Decimal(p)))
            reciprocal_u=up.add(reciprocal_u,up.divide(Decimal(den),Decimal(p)))
        row_vl=down.add(Decimal(leading),down.divide(down.subtract(reciprocal_l,Decimal(correction)),Decimal(N)))
        row_vu=up.add(Decimal(leading),up.divide(up.subtract(reciprocal_u,Decimal(correction)),Decimal(N)))
        vl=down.add(vl,row_vl);vu=up.add(vu,row_vu)
    mean=Decimal(mean_int)
    if vl<0:raise ValueError('unsupported negative variance bound')
    def tail_bound(stat,var,ctx,other):
        if stat<=mean:return Decimal(1)
        delta=other.subtract(stat,mean)
        square=other.multiply(delta,delta)
        return ctx.divide(var,other.add(var,square))
    # Cantelli is increasing in V and decreasing in X above the null mean.
    if impossible:
        xl=xu=Decimal('Infinity');lower=upper=Decimal(0)
    else:
        lower=tail_bound(xu,vl,down,up)
        upper=tail_bound(xl,vu,up,down)
    decision='inside' if lower>alpha else 'excluded' if upper<=alpha else 'rounding_unresolved'
    return dict(rows=int(counts.size//4),outcomes_per_row=4,aggregate_component_alpha=str(alpha),
                mean_exact=int(mean),pearson_lower=str(xl),pearson_upper=str(xu),
                variance_lower=str(vl),variance_upper=str(vu),
                cantelli_tail_lower=str(lower),cantelli_tail_upper=str(upper),
                aggregate_region_membership=decision,
                classical_memory_class_excluded=False,
                fitted_degrees_of_freedom_used=False)


def binomial_cell_membership(n,k,pnum,pden,cells):
    """Certify the 97.5% simultaneous CP box, tail threshold 1/(80M)."""
    if not (n>0 and 0<=k<=n and 0<=pnum<=pden and pden>0 and cells>0):
        raise ValueError('invalid binomial cell')
    if pnum==0:return 'inside' if k==0 else 'excluded'
    if pnum==pden:return 'inside' if k==n else 'excluded'
    down=Context(prec=60,rounding=ROUND_FLOOR);up=Context(prec=60,rounding=ROUND_CEILING)
    tl=down.divide(Decimal(1),Decimal(80*cells));tu=up.divide(Decimal(1),Decimal(80*cells))
    masses=[]
    for ctx in (down,up):
        p=ctx.divide(Decimal(pnum),Decimal(pden));q=ctx.divide(Decimal(pden-pnum),Decimal(pden))
        masses.append(ctx.multiply(Decimal(physical.choose(n,k)),ctx.multiply(physical.floor_power(ctx,p,k),physical.floor_power(ctx,q,n-k))))
    decisions=[]
    for direction in (-1,1):
        if (k==0 and direction==1) or (k==n and direction==-1):continue
        low,high=masses;term_l,term_u=masses;j=k
        while low<tu and 0<=j+direction<=n:
            top,bottom=((n-j)*pnum,(j+1)*(pden-pnum)) if direction==1 else (j*(pden-pnum),(n-j+1)*pnum)
            rl=down.divide(Decimal(top),Decimal(bottom));ru=up.divide(Decimal(top),Decimal(bottom))
            # Future ratios decrease in either outward direction.
            if ru<1:
                remainder=up.divide(up.multiply(term_u,ru),down.subtract(Decimal(1),ru))
                if up.add(high,remainder)<tl:return 'excluded'
            term_l=down.multiply(term_l,rl);term_u=up.multiply(term_u,ru)
            low=down.add(low,term_l);high=up.add(high,term_u);j+=direction
        if low>=tu:decisions.append('inside')
        elif high<tl:return 'excluded'
        else:decisions.append('rounding_unresolved')
    return 'inside' if all(d=='inside' for d in decisions) else 'rounding_unresolved'


def cell_region(counts,numerators,denominators):
    counts=np.asarray(counts);nums=np.broadcast_to(np.asarray(numerators,dtype=object),counts.shape)
    dens=np.broadcast_to(np.asarray(denominators,dtype=object),counts.shape[:-1])
    tally=dict(inside=0,excluded=0,rounding_unresolved=0);failures=[]
    for row,(obs,ps,den) in enumerate(zip(counts.reshape(-1,4),nums.reshape(-1,4),dens.reshape(-1))):
        N=int(obs.sum())
        for outcome,(k,p) in enumerate(zip(obs,ps)):
            decision=binomial_cell_membership(N,int(k),int(p),int(den),counts.size)
            tally[decision]+=1
            if decision!='inside' and len(failures)<10:failures.append(dict(row=row,outcome=outcome,decision=decision))
    membership='excluded' if tally['excluded'] else 'rounding_unresolved' if tally['rounding_unresolved'] else 'inside'
    return dict(cell_region_membership=membership,cell_component_alpha='0.025',
                cell_tail_denominator=80*counts.size,cell_decisions=tally,first_failed_cells=failures)


def analyze(source_dir,candidate_file=physical.HERE/'results/ibm-delay-candidates.json'):
    catalog=json.loads((physical.HERE/'results/ibm-instrument-witnesses.json').read_text())
    reports=[]
    for witness in catalog['witnesses']:
        physical.check_instruments(witness)
        counts,labels,_=physical.load_counts(source_dir,witness['mapping'])
        nums,den=physical.probability_numerators(witness,labels)
        reports.append(dict(mapping=witness['mapping'],point='previous_shared_identity_process_pair',
                            **pearson_region(counts,nums,den),**cell_region(counts,nums,den)))
    candidates=json.loads(candidate_file.read_text())
    for candidate in candidates['candidates']:
        witness=candidate['instrument']
        boundary=delay_candidates.check_candidate(candidate)
        counts,labels,delays=physical.load_counts(source_dir,witness['mapping'])
        if delays!=candidate['delays']:raise ValueError('delay order mismatch')
        nums,dens=delay_candidates.probabilities(candidate,labels)
        p=np.asarray(nums/dens[...,None],float)
        expected=counts.sum(-1)[...,None]*p
        deviance=float(2*np.where(counts>0,counts*np.log(np.maximum(counts,1)/expected),0).sum())
        reports.append(dict(mapping=witness['mapping'],point='delay_dependent_unitaries_shared_fitted_spam_pair',
                            exact_physical_classical_and_quantum_pair=True,
                            quantum_npt_expectation=str(-physical.F(*witness['memory_weight'])/8),
                            fixed_instrument_weight_endpoint=float(boundary),
                            multinomial_deviance_display=round(deviance,6),
                            **pearson_region(counts,nums,dens),**cell_region(counts,nums,dens)))
    for report in reports:
        decisions=[report['aggregate_region_membership'],report['cell_region_membership']]
        report['joint_region_membership']='excluded' if 'excluded' in decisions else 'inside' if decisions==['inside','inside'] else 'rounding_unresolved'
    return dict(kind='finite_sample_joint_point_region_audit',familywise_alpha=.05,
                confidence_region='Intersection of 97.5% simultaneous cell intervals and 97.5% aggregate Cantelli region; both components evaluated.',
                assumptions='Independent fixed-exposure multinomial rows; conditional on record mapping and stable within-row probabilities.',
                retrospective_fitting='Confidence-region inversion; no exclusion of the composite fitted model family.',
                previous_cell_interval_certificate_unchanged=True,reports=reports)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-dir',required=True,type=Path)
    parser.add_argument('--check',action='store_true')
    parser.add_argument('--candidate-file',type=Path,default=physical.HERE/'results/ibm-delay-candidates.json')
    parser.add_argument('--output',type=Path,default=physical.HERE/'results/ibm-joint-statistics.json')
    args=parser.parse_args()
    output=json.dumps(analyze(args.source_dir,args.candidate_file),indent=2,sort_keys=True)+'\n'
    target=args.output
    if args.check:
        if target.read_text()!=output:raise SystemExit('joint statistics differ')
        print('Finite-sample aggregate region audit reproduced')
    else:target.write_text(output);print(target)


if __name__=='__main__':main()
