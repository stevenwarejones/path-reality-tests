#!/usr/bin/env python3
"""Exact verifier: no optimizer, NumPy, or floating-point feasibility decisions."""
from fractions import Fraction as F
from pathlib import Path
import argparse
import json
import math

HERE = Path(__file__).resolve().parent

def cdf_numerator(n, k, p):
    """Integer numerator of Bin(n,p) CDF, denominator p.denominator**n."""
    a, d = p.numerator, p.denominator
    b = d - a
    if k < 0:
        return 0, d**n
    if a == 0:
        return 1, 1
    if b == 0:
        return int(k >= n), 1
    term = b**n
    total = term
    for j in range(min(k, n)):
        term, rem = divmod(term * (n-j) * a, (j+1) * b)
        assert rem == 0
        total += term
    return total, d**n

def check_endpoint(n, k, p, tail, side):
    assert 0 <= k <= n and 0 <= p <= 1 and 0 < tail < 1
    if side == 'lower':
        if k == 0:
            assert p == 0
            return
        num, den = cdf_numerator(n, k-1, p)
        assert (den-num)*tail.denominator <= den*tail.numerator
    elif side == 'upper':
        if k == n:
            assert p == 1
            return
        num, den = cdf_numerator(n, k, p)
        assert num*tail.denominator <= den*tail.numerator
    else:
        raise ValueError(side)

def constraints(lo, hi, loss_interval):
    """q = rows 7..20, residual loss; z_i=q_i q_(i+1), t_i=z_i q_(i+2)."""
    assert len(lo) == len(hi) == 15
    rows = [{0:F(1), 1:F(1), 14:F(1)}, {0:F(-1), 1:F(-1), 14:F(-1)}]
    rhs = [loss_interval[1], -loss_interval[0]]
    bounds = list(zip(lo, hi))
    def mc(i, j, k, li, ui, lj, uj):
        for a,b,c,d in [(lj,li,-1,li*lj),(uj,ui,-1,ui*uj),
                        (-lj,-ui,1,-ui*lj),(-uj,-li,1,-li*uj)]:
            rows.append({i:F(a), j:F(b), k:F(c)})
            rhs.append(F(d))
    for i in range(12):
        zl, zu = lo[i]*lo[i+1], hi[i]*hi[i+1]
        bounds.append((zl, zu))
        mc(i,i+1,15+i,lo[i],hi[i],lo[i+1],hi[i+1])
        mc(15+i,i+2,27+i,zl,zu,lo[i+2],hi[i+2])
    bounds.extend((lo[i]*lo[i+1]*lo[i+2], hi[i]*hi[i+1]*hi[i+2]) for i in range(12))
    return rows, rhs, bounds

def dual_bound(lo, hi, loss, multipliers, mu):
    rows, rhs, bounds = constraints(lo,hi,loss)
    assert len(multipliers) == len(rows)
    assert all(v >= 0 for v in multipliers)
    residual = [F(0)]*27 + [F(1)]*12
    bound = mu
    for j in range(15):
        residual[j] -= mu
    for lam, row, b in zip(multipliers, rows, rhs):
        bound += lam*b
        for j,v in row.items():
            residual[j] -= lam*v
    # Residual repair makes arbitrary rounded duals rigorous, even if not feasible.
    for r,(l,u) in zip(residual,bounds):
        bound += r*(u if r >= 0 else l)
    return bound

def verify(cert, counts):
    # Fix all allocation choices here, not in untrusted certificate metadata.
    assert cert['schema'] == 1
    assert counts['selected_n'] == 3
    assert len(counts['reference_row_counts']) == 12
    assert counts['reference_row_indices'] == list(range(9,21))
    assert counts['many_row_indices'] == list(range(8,20))
    assert counts['copies'] == {'31':1, '32':-1}
    cats = counts['reference_row_counts'] + [counts['single_shots']-sum(counts['reference_row_counts'])]
    assert counts['single_shots'] == 931 and counts['many_shots'] == 2999
    assert counts['bunch_shots'] == 93
    intervals = [[F(v) for v in pair] for pair in cert['category_intervals']]
    assert len(intervals) == len(cats)
    for k,(l,u) in zip(cats,intervals):
        assert l <= u
        check_endpoint(counts['single_shots'], k, l, F(1,4160), 'lower')
        check_endpoint(counts['single_shots'], k, u, F(1,4160), 'upper')
    lower = F(cert['bunch_lower'])
    check_endpoint(counts['many_shots'],counts['bunch_shots'],lower,F(1,160),'lower')
    loss = intervals[-1]
    lo = [F(0),F(0)] + [x[0] for x in intervals[:-1]] + [F(0)]
    hi = [loss[1],loss[1]] + [x[1] for x in intervals[:-1]] + [loss[1]]
    leaves = 0
    def visit(node,l,u):
        nonlocal leaves
        if 'split' in node:
            j = node['split']; mid = F(node['at'])
            assert isinstance(j,int) and 0 <= j < 15 and l[j] < mid < u[j]
            left_u = u.copy();left_u[j]=mid
            right_l = l.copy();right_l[j]=mid
            return max(visit(node['left'],l,left_u),visit(node['right'],right_l,u))
        leaves += 1
        nr = len(constraints(l,u,loss)[0])
        lam = [F(0)]*nr
        for key,v in node['lambda'].items():
            assert 0 <= int(key) < nr
            lam[int(key)] = F(v)
        return dual_bound(l,u,loss,lam,F(node['mu']))
    upper = visit(cert['tree'],lo,hi)
    reported = F(cert['polynomial_upper'])
    assert upper <= reported
    margin = lower-2*reported
    assert margin > 0
    return {'leaf_count':leaves, 'bunch_lower':float(lower),
            'distinguishable_row_upper':float(reported),
            'cluster_two_ceiling':float(2*reported),
            'separation_margin':float(margin),
            'sum_transfer_tv_budget':float(margin/2),
            'non_cluster_weight_lower':float(margin/(4*reported)),
            'confidence_familywise':0.95}

def main():
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');args=p.parse_args()
    counts=json.loads((HERE/'results/counts.json').read_text())
    cert=json.loads((HERE/'results/certificate.json').read_text())
    result=verify(cert,counts)
    path=HERE/'results/bounds.json'
    text=json.dumps(result,indent=2,sort_keys=True)+'\n'
    if args.check:
        assert path.read_text()==text, 'Stale bounds.json'
    else:
        path.write_text(text)
    print(json.dumps(result,indent=2,sort_keys=True))

if __name__=='__main__':
    main()
