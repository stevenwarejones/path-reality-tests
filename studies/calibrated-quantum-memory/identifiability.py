#!/usr/bin/env python3
"""Independently certify physical, observationally equal IBM memory models.

Discovery uses an optional SDP; this verifier uses rational matrices and
one-sided Decimal bounds on binomial tails, not optimizer status/tolerances.
"""
import argparse
from decimal import Context, Decimal, ROUND_CEILING, ROUND_FLOOR
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
from scipy.stats import beta

HERE = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location('nmn_audit', HERE/'audit.py')
audit = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(audit)
I = np.eye(2, dtype=complex)
PAULIS = np.array([[[0,1],[1,0]], [[0,-1j],[1j,0]], [[1,0],[0,-1]]])
STATES = {a+s:(I+sign*PAULIS[i])/2 for i,a in enumerate('xyz')
          for s,sign in [('p',1),('m',-1)]}
IDENTITY_CHOI = np.outer([1,0,0,1], [1,0,0,1]).astype(int)


def ldl_positive(matrix):
    """Exact LDL^T of a real symmetric rational matrix; reject non-PD input."""
    n = len(matrix)
    if any(matrix[i][j] != matrix[j][i] for i in range(n) for j in range(n)):
        raise ValueError('not symmetric')
    lower = [[F(int(i == j)) for j in range(n)] for i in range(n)]
    diagonal = []
    for i in range(n):
        d = F(matrix[i][i]) - sum(lower[i][k]**2*diagonal[k] for k in range(i))
        if d <= 0:
            raise ValueError('matrix not strictly positive definite')
        diagonal.append(d)
        for j in range(i+1,n):
            lower[j][i] = (F(matrix[j][i]) -
                           sum(lower[j][k]*lower[i][k]*diagonal[k] for k in range(i)))/d
    return lower, diagonal


def realify(real, imag):
    if not np.array_equal(real,real.T) or not np.array_equal(imag,-imag.T):
        raise ValueError('non-Hermitian Choi matrix')
    return np.block([[real,-imag],[imag,real]]).astype(object).tolist()


def inverse_quadratic(lower, diagonal, vector):
    y = []
    for i,v in enumerate(vector):
        y.append(F(v)-sum(lower[i][k]*y[k] for k in range(i)))
    return sum(y[i]**2/diagonal[i] for i in range(len(y)))


def partial_trace_output(matrix):
    return np.array([[sum(matrix[2*i+k,2*j+k] for k in range(2))
                      for j in range(2)] for i in range(2)],dtype=object)


def arrays(record):
    return (np.array(record['real'],dtype=object).reshape(4,4),
            np.array(record['imag'],dtype=object).reshape(4,4))


def check_instruments(witness):
    denominator = witness['denominator']
    if type(denominator) is not int or denominator <= 0:
        raise ValueError('invalid denominator')
    w = F(*witness['memory_weight'])
    if not 0 < w < 1:
        raise ValueError('invalid quantum mixture weight')
    expected={m+','+p for m in 'xyz' for p in audit.STATES}
    if set(witness['settings']) != expected:
        raise ValueError('invalid setting grid')
    effects = {}
    boundaries = []
    for setting, pair in zip(witness['settings'],witness['choi']):
        if len(pair) != 2:
            raise ValueError('missing retained flag')
        traces = []
        quantum_traces = []
        for b,record in enumerate(pair):
            if any(type(x) is not int for key in ('real','imag') for x in record[key]):
                raise ValueError('certificate coefficients must be integers')
            re,im = arrays(record)
            lower,diag = ldl_positive(realify(re,im))
            trace = sum(re[i,i] for i in range(4))
            # q=Tr(J)/2D, v=vec(I). The CP endpoint is exact by a rank-one update.
            quadratic = inverse_quadratic(lower,diag,[1,0,0,1,0,0,0,0])
            boundaries.append(F(2,1)/(trace*quadratic))
            qre = 2*w.denominator*re - w.numerator*trace*IDENTITY_CHOI.astype(object)
            qim = 2*w.denominator*im
            ldl_positive(realify(qre,qim))
            quantum_traces.append((partial_trace_output(qre),partial_trace_output(qim)))
            er,ei = partial_trace_output(re),partial_trace_output(im)
            key = (setting.split(',')[0],b)
            if key in effects and (not np.array_equal(er,effects[key][0]) or
                                   not np.array_equal(ei,effects[key][1])):
                raise ValueError('first measurement changes with later rotation')
            effects[key] = (er,ei)
            traces.append((er,ei))
            # Exact equality J_class = (1-w) J_quantum + w*q J_identity.
            if not np.array_equal(qre+w.numerator*trace*IDENTITY_CHOI.astype(object),
                                  2*w.denominator*re):
                raise AssertionError('affine instrument identity')
        if not np.array_equal(traces[0][0]+traces[1][0],denominator*np.eye(2,dtype=object)):
            raise ValueError('instrument sum is not trace preserving')
        qden=2*(w.denominator-w.numerator)*denominator
        if not np.array_equal(quantum_traces[0][0]+quantum_traces[1][0],qden*np.eye(2,dtype=object)) or np.any(quantum_traces[0][1]+quantum_traces[1][1]):
            raise ValueError('quantum instrument is not trace preserving')
        if np.any(traces[0][1]+traces[1][1]):
            raise ValueError('imaginary trace-preservation residual')
    if len(witness['settings']) != 18 or len(witness['choi']) != 18:
        raise ValueError('incomplete setting grid')
    return min(boundaries)


def probability_numerators(witness, labels):
    result = []
    for label in labels:
        a,m,p,z = label.split(',')
        yi = witness['settings'].index(m+','+p)
        row = []
        for b in range(2):
            re,im = arrays(witness['choi'][yi][b])
            for c in range(2):
                k = 4*np.kron(STATES[a].T,(I+(1-2*c)*PAULIS['xyz'.index(z)])/2)
                kr = np.rint(k.real).astype(int).T
                ki = np.rint(k.imag).astype(int).T
                row.append(int(np.sum(re*kr-im*ki)))
        result.append(row)
    nums = np.array(result,dtype=np.int64)
    den = 4*witness['denominator']
    if np.any(nums<0) or np.any(nums>den) or np.any(nums.sum(axis=1)!=den):
        raise ValueError('nonphysical joint probabilities')
    return nums,den


def floor_power(context, base, exponent):
    result = Decimal(1)
    while exponent:
        if exponent & 1:
            result = context.multiply(result,base)
        exponent >>= 1
        if exponent:
            base = context.multiply(base,base)
    return result


@lru_cache(maxsize=20000)
def choose(n,k):
    return math.comb(n,k)


def binomial_membership(n,k,pnum,pden,cells):
    """Prove both required tails >= alpha/(2*cells), alpha=1/20.

    Every mass, ratio, product and sum is rounded downward. The threshold is
    rounded upward. A partial tail sum is a lower bound on the complete tail.
    No floating-point inverse CDF is used to accept a model.
    """
    if not (0<=k<=n and 0<pnum<pden and cells>0):
        raise ValueError('unsupported or invalid binomial certificate')
    down = Context(prec=50,rounding=ROUND_FLOOR)
    up = Context(prec=50,rounding=ROUND_CEILING)
    threshold = up.divide(Decimal(1),Decimal(40*cells))
    pp = down.divide(Decimal(pnum),Decimal(pden))
    qq = down.divide(Decimal(pden-pnum),Decimal(pden))
    mass = down.multiply(Decimal(choose(n,k)),
                         down.multiply(floor_power(down,pp,k),floor_power(down,qq,n-k)))
    checked = []
    for direction in (-1,1):
        # At endpoints CP imposes only one tail.
        if (k==0 and direction==1) or (k==n and direction==-1):
            continue
        total,term,j = mass,mass,k
        while total < threshold and 0 <= j+direction <= n:
            if direction == 1:
                top,bottom = (n-j)*pnum,(j+1)*(pden-pnum)
            else:
                top,bottom = j*(pden-pnum),(n-j+1)*pnum
            ratio = down.divide(Decimal(top),Decimal(bottom))
            term = down.multiply(term,ratio)
            total = down.add(total,term)
            j += direction
        if total < threshold:
            raise ValueError(f'binomial region failure: N={n}, k={k}, p={pnum}/{pden}')
        checked.append(float(down.divide(total,threshold)))
    return min(checked)


def load_counts(source_dir, mapping):
    if mapping not in ('regroup_first_1','literal_rows'):
        raise ValueError('unsupported record mapping')
    manifest = json.loads((HERE/'sources.json').read_text())
    entry = next(e for e in manifest['nmn_files'] if e['name']=='NMN_tomog_rerun.json')
    runs = json.loads(audit.verify_file(source_dir/entry['name'],entry))
    labels = sorted(audit.LABELS)
    counts = []
    for key in sorted(runs):
        rows = audit.regroup(runs[key]) if mapping=='regroup_first_1' else runs[key]
        audit.validate_run(rows)
        counts.append([[rows[label][o] for o in audit.OUTCOMES] for label in labels])
    return np.array(counts,dtype=np.int64),labels,sorted(runs)


def isolated_probe(witness, labels):
    nums,den = probability_numerators(witness,labels)
    w = F(*witness['memory_weight'])
    best = None
    for r,label in enumerate(labels):
        a,m,p,z = label.split(',');yi=witness['settings'].index(m+','+p)
        differences = []
        for b in range(2):
            re,_=arrays(witness['choi'][yi][b]);q=F(int(np.trace(re)),2*witness['denominator'])
            for c in range(2):
                identity_probability=F(int(round(2*np.trace(STATES[a]@(I+(1-2*c)*PAULIS['xyz'.index(z)])/2).real)),2)
                classical=F(int(nums[r,2*b+c]),den)
                quantum=(classical-w*q*identity_probability)/(1-w)
                differences.append(classical-quantum)
        tv=sum(abs(x) for x in differences)/2
        if best is None or tv>best[0]: best=(tv,label,differences)
    tv,label,differences=best
    event=[audit.OUTCOMES[i] for i,d in enumerate(differences) if d>0]
    # Two separately estimated event probabilities. Each Hoeffding radius < gap/2.
    shots=math.ceil(2*math.log(4/.05)/float(tv)**2)+1
    return dict(label=label,positive_difference_event=event,
                total_variation_exact=[tv.numerator,tv.denominator],
                total_variation=float(tv),half_diamond_upper=float(w),
                shots_per_candidate_for_two_radius_separation_at_95_percent=shots)



def assignment_ceiling_diagnostic(counts,labels,delays):
    # Fixed Table-4 readout: any first-bit-0 joint event has p <= .976.
    restricted=counts[:,:,:2]
    n=counts.sum(axis=-1)
    frequencies=restricted/n[...,None]
    d,r,o=np.unravel_index(frequencies.argmax(),frequencies.shape)
    k,total=int(counts[d,r,o]),int(n[d,r])
    up=Context(prec=50,rounding=ROUND_CEILING)
    down=Context(prec=50,rounding=ROUND_FLOOR)
    p,q=up.divide(Decimal(122),Decimal(125)),up.divide(Decimal(3),Decimal(125))
    mass=up.multiply(Decimal(choose(total,k)),up.multiply(floor_power(up,p,k),floor_power(up,q,total-k)))
    ratio=up.divide(Decimal((total-k)*122),Decimal((k+1)*3))
    if ratio>=1:raise ValueError('geometric tail bound not applicable')
    upper=up.divide(mass,down.subtract(Decimal(1),ratio))
    threshold=down.divide(Decimal(1),Decimal(40*counts.size))
    return dict(label=labels[r],delay=delays[d],outcome=audit.OUTCOMES[o],count=k,shots=total,
                assumed_exact_assignment_ceiling=.976,binomial_upper_tail_bound=str(upper),
                rejects_fixed_forward_assignment_point=upper<threshold,
                interpretation='Conditional forward-model consistency failure; not a quantum-memory exclusion.')


def analyze(source_dir, witness_file=HERE/'results/ibm-instrument-witnesses.json'):
    catalog=json.loads(witness_file.read_text())
    reports=[]
    for witness in catalog['witnesses']:
        boundary=check_instruments(witness)
        counts,labels,delays=load_counts(source_dir,witness['mapping'])
        nums,den=probability_numerators(witness,labels)
        cells=counts.size
        probabilities=nums/den
        minimum_tail_ratio=math.inf
        n=counts.sum(axis=-1)
        for d in range(len(delays)):
            for r in range(len(labels)):
                for o in range(4):
                    ratio=binomial_membership(int(n[d,r]),int(counts[d,r,o]),int(nums[r,o]),den,cells)
                    minimum_tail_ratio=min(minimum_tail_ratio,ratio)
        tail=1/(40*cells)
        lo=np.where(counts==0,0,beta.ppf(tail,np.maximum(counts,1),n[...,None]-counts+1))
        hi=np.where(counts==n[...,None],1,beta.ppf(1-tail,counts+1,np.maximum(n[...,None]-counts,1)))
        slack=np.minimum(probabilities[None]-lo,hi-probabilities[None])
        per_run=[dict(delay=d,minimum_interval_slack=round(float(slack[i].min()),12)) for i,d in enumerate(delays)]
        split={}
        for name,axis_test in [('retain_z_only',lambda m:m=='z'),('retain_x_y_only',lambda m:m!='z')]:
            ix=[i for i,k in enumerate(labels) if axis_test(k.split(',')[1])]
            split[name]=dict(cells=int(counts[:,ix].size),confidence_budget_cells=cells,
                            same_complete_physical_models_fit=True,minimum_interval_slack=round(float(slack[:,ix].min()),12))
        expected=n[...,None]*probabilities[None]
        deviance=float(2*np.where(counts>0,counts*np.log(np.maximum(counts,1)/expected),0).sum())
        empirical_shared=counts.sum(axis=0)/n.sum(axis=0)[:,None]
        shared_deviance=float(2*np.where(counts>0,counts*np.log(np.maximum(counts,1)/np.maximum(n[...,None]*empirical_shared[None],1e-300)),0).sum())
        reports.append(dict(multinomial_deviance_to_saturated=round(deviance,6),
                            unconstrained_shared_row_deviance=round(shared_deviance,6),
                            saturated_multinomial_dimension=int(3*counts.shape[0]*counts.shape[1]),
                            aggregate_goodness_of_fit_not_certified=True,
                            exact_table4_diagnostic=assignment_ceiling_diagnostic(counts,labels,delays),
                            mapping=witness['mapping'],outcome='physical_classical_and_quantum_models_certified',
                            cells=cells,shots_per_row_min=int(n.min()),shots_per_row_max=int(n.max()),
                            memory_weight=witness['memory_weight'],exact_classical_cp_tp=True,exact_quantum_instrument_cp_tp=True,
                            shared_first_effects_across_repreparations=True,
                            all_binomial_tails_certified=True,minimum_certified_tail_to_threshold_ratio=round(minimum_tail_ratio,10),
                            interval_slack_is_display_only=True,per_run=per_run,source_deletion=split,
                            exact_probability_identity=True,
                            fixed_classical_instrument_memory_weight_boundary=[boundary.numerator,boundary.denominator],
                            boundary_decimal=float(boundary),
                            quantum_partial_transpose_witness_expectation=str(-F(*witness['memory_weight'])/8),
                            isolated_instrument_probe=isolated_probe(witness,labels),
                            probability_numerators_sha256=hashlib.sha256(json.dumps(nums.tolist(),separators=(',',':')).encode()).hexdigest()))
    return dict(kind='measured_conditional_full_record_identifiability_certificate',alpha=.05,
                empirical_quantum_memory_excluded_or_certified=False,
                mapping_historical_provenance_verified=False,
                first_bit_zero_branch='Same certified models by opposite-repreparation relabeling; all constraints permute.',
                interpretation='A shared-instrument classical point remains; no positive joint exclusion over the mapping union.',
                reports=reports)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source-dir',required=True,type=Path)
    p.add_argument('--check',action='store_true')
    p.add_argument('--witness-file',type=Path,default=HERE/'results/ibm-instrument-witnesses.json')
    p.add_argument('--output',type=Path,default=HERE/'results/ibm-identifiability.json')
    args=p.parse_args()
    result=json.dumps(analyze(args.source_dir,args.witness_file),indent=2,sort_keys=True)+'\n'
    target=args.output
    if args.check:
        if target.read_text()!=result: raise SystemExit('identifiability result differs')
        print('Exact physicality, every binomial constraint, equal flagged records and discriminating probe verified')
    else:
        target.write_text(result);print(target)


if __name__=='__main__':main()
