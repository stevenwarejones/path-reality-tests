"""Exact prospective binomial decision and explicitly hypothetical sensitivity."""
import argparse,json,math
from fractions import Fraction as F
from dynamics_common import HERE,save,GRID
from dynamics_verify import lower_bound

def cdf(n,k,p):return sum(F(math.comb(n,j))*p**j*(1-p)**(n-j) for j in range(k+1))
def run():
    c=json.loads((HERE/'results/certificate.json').read_text());L=F(c['joint_lower']['numerator'],c['joint_lower']['denominator'])
    z=c['models']['RB_only']['target_upper'];p=F(z['numerator'],z['denominator']);decision=None
    for n in range(1,201):
        for k in range(n+1):
            size=cdf(n,k,L);power=cdf(n,k,p)
            if size<=F(1,20) and power>=F(19,20):
                decision={'shots':n,'reject_if_plus_count_at_most':k,'conditional_size_upper':float(F(math.ceil(size*10**12),10**12)),'power_at_RB_witness_lower':float(F(math.floor(power*10**12),10**12)), 'probability_rounding_grid':10**12};break
        if decision:break
    if decision is None:raise ValueError('no decision in declared search')
    inp=c['inputs'];a=inp['0']['upper_numerator'];b=inp['4662']['upper_numerator'];f=inp['10']['lower_numerator']
    sensitivity=[]
    for delta in [F(0),F(1,10000),F(1,1000),F(1,100)]:
        # Two gates in observed II, two in YY, four in the target. Empty word unchanged.
        num=int(2*delta*GRID)
        bound=max(F(0),lower_bound(a,b+num,f)-6*delta)
        sensitivity.append({'half_diamond_per_occurrence':float(delta),'target_lower':float(bound)})
    return {'prospective':decision,'prospective_scope':'Independent future shots, conditional on training-region membership. Training coverage 95%; future conditional test size 5%; union bound gives at most 10% unconditional error, not 5%. No experiment conducted.',
        'hypothetical_context_variation':sensitivity,'context_scope':'Fixed circuit-context channels within half-diamond delta of one stationary qubit gate set; fixed shared SPAM and independent shots. Delta is not inferred from this acquisition.',
        'composition_removed_target_interval':[0,1],
        'within_row_dependence':'Unrestricted dependence invalidates the binomial/Hoeffding confidence claim. Algebraic implications of region membership remain; no 95% empirical certificate is asserted in that enlargement.'}
if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--check',action='store_true');args=a.parse_args();out=run();p=HERE/'results/controls.json'
    if args.check:
        if json.loads(p.read_text())!=out:raise ValueError('controls artifact mismatch')
    else:save(p,out)
    print(out['prospective']);print(out['hypothetical_context_variation'])
