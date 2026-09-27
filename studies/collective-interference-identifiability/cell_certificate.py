#!/usr/bin/env python3
"""Exact 2D uncertainty-region and repaired-dual verification, without SciPy."""
from fractions import Fraction as F
from pathlib import Path
import argparse,json
from certificate import check_endpoint, constraints as row_constraints

HERE=Path(__file__).resolve().parent

def product_rows(rows,rhs,i,j,k,li,ui,lj,uj):
    for a,b,c,d in [(lj,li,-1,li*lj),(uj,ui,-1,ui*uj),
                    (-lj,-ui,1,-ui*lj),(-uj,-li,1,-li*uj)]:
        rows.append({i:F(a),j:F(b),k:F(c)});rhs.append(F(d))

def conditional_problem(lo,hi):
    """3 conditional x-distributions; objective is 1 - pair collisions + 2 triple collisions."""
    assert len(lo)==len(hi)==3 and all(len(v)==12 for v in lo+hi)
    bounds=[(lo[j][x],hi[j][x]) for j in range(3) for x in range(12)]
    pairs=[(0,1),(0,2),(1,2)]
    bounds += [(lo[a][x]*lo[b][x],hi[a][x]*hi[b][x]) for a,b in pairs for x in range(12)]
    bounds += [(lo[0][x]*lo[1][x]*lo[2][x],hi[0][x]*hi[1][x]*hi[2][x]) for x in range(12)]
    rows=[];rhs=[]
    for j,(a,b) in enumerate(pairs):
        for x in range(12):
            product_rows(rows,rhs,12*a+x,12*b+x,36+12*j+x,lo[a][x],hi[a][x],lo[b][x],hi[b][x])
    for x in range(12):
        product_rows(rows,rhs,36+x,24+x,72+x,*bounds[36+x],lo[2][x],hi[2][x])
    eq=[{12*j+x:F(1) for x in range(12)} for j in range(3)]
    c=[F(0)]*36+[F(-1)]*36+[F(2)]*12
    return rows,rhs,bounds,eq,c,F(1)

def row_problem(lo,hi,loss,weights):
    rows,rhs,bounds=row_constraints(lo,hi,loss)
    return rows,rhs,bounds,[{j:F(1) for j in range(15)}],[F(0)]*27+weights,F(0)

def repaired_bound(problem,dual):
    rows,rhs,bounds,eq,c,constant=problem
    mu=list(map(F,dual['mu']));assert len(mu)==len(eq)
    lam=[F(0)]*len(rows)
    for i,v in dual['lambda'].items():
        assert 0<=int(i)<len(rows);lam[int(i)]=F(v)
    assert all(v>=0 for v in lam)
    residual=c.copy();bound=constant+sum(mu)
    for m,row in zip(mu,eq):
        for j,v in row.items():residual[j]-=m*v
    for m,row,b in zip(lam,rows,rhs):
        bound+=m*b
        for j,v in row.items():residual[j]-=m*v
    for r,(l,u) in zip(residual,bounds):bound+=r*(u if r>=0 else l)
    return bound

def tree_bound(tree,lo,hi,loss,weights):
    if 'split' not in tree:return repaired_bound(row_problem(lo,hi,loss,weights),tree),1
    j=tree['split'];m=F(tree['at'])
    assert isinstance(j,int) and 0<=j<15 and lo[j]<m<hi[j]
    lu=hi.copy();lu[j]=m;rl=lo.copy();rl[j]=m
    a,na=tree_bound(tree['left'],lo,lu,loss,weights)
    b,nb=tree_bound(tree['right'],rl,hi,loss,weights)
    return max(a,b),na+nb

def verify(cert,counts,witness):
    assert cert['schema']==1
    C=[[0]*12 for _ in range(12)]
    for y,x,n in witness['nonzero_cell_counts']:
        assert 9<=y<=20 and 8<=x<=19 and n>0
        assert C[y-9][x-8]==0
        C[y-9][x-8]=n
    row_counts=list(map(sum,C));assert row_counts==counts['reference_row_counts']
    intervals=[[F(v) for v in p] for p in cert['row_intervals']]
    cats=row_counts+[counts['single_shots']-sum(row_counts)]
    assert len(intervals)==13
    for k,(l,u) in zip(cats,intervals):
        assert l<=u
        check_endpoint(counts['single_shots'],k,l,F(1,8320),'lower')
        check_endpoint(counts['single_shots'],k,u,F(1,8320),'upper')
    conditionals=cert['conditional_intervals'];assert len(conditionals)==12
    cints=[]
    for y,row in enumerate(conditionals):
        assert len(row)==12
        pairs=[[F(v) for v in p] for p in row]
        for x,(l,u) in enumerate(pairs):
            assert l<=u
            check_endpoint(row_counts[y],C[y][x],l,F(1,92160),'lower')
            check_endpoint(row_counts[y],C[y][x],u,F(1,92160),'upper')
        cints.append(pairs)
    weights=[]
    assert len(cert['conditional_duals'])==12
    for y,dual in zip(range(8,20),cert['conditional_duals']):
        ls=[];us=[]
        for z in [y-1,y,y+1]:
            pairs=cints[z-9] if 9<=z<=20 else [[F(0),F(1)]]*12
            ls.append([p[0] for p in pairs]);us.append([p[1] for p in pairs])
        ub=min(F(1),repaired_bound(conditional_problem(ls,us),dual))
        reported=F(dual['upper']);assert 0<=reported<=1 and ub<=reported
        weights.append(reported)
    loss=intervals[-1]
    lo=[F(0),F(0)]+[p[0] for p in intervals[:-1]]+[F(0)]
    hi=[loss[1],loss[1]]+[p[1] for p in intervals[:-1]]+[loss[1]]
    old,old_n=tree_bound(cert['row_tree'],lo,hi,loss,[F(1)]*12)
    new,new_n=tree_bound(cert['cell_tree'],lo,hi,loss,weights)
    row_upper=F(cert['row_upper']);cell_upper=F(cert['cell_upper'])
    assert old<=row_upper and new<=cell_upper and cell_upper<row_upper
    lower=F(cert['bunch_lower'])
    check_endpoint(counts['many_shots'],counts['bunch_shots'],lower,F(1,160),'lower')
    assert lower>2*cell_upper
    # Exact empirical singleton distribution belongs to every full-cell constraint.
    for y in range(12):
        q=F(row_counts[y],counts['single_shots'])
        assert intervals[y][0]<=q<=intervals[y][1]
        if row_counts[y]:
            for x in range(12):assert cints[y][x][0]<=F(C[y][x],row_counts[y])<=cints[y][x][1]
    margin=lower-2*cell_upper
    return {'confidence_familywise':0.95,'bunch_lower':float(lower),
        'historical_row_upper':0.01166,'same_region_row_upper':float(row_upper),
        'cell_event_upper':float(cell_upper),'cluster_ceiling':float(2*cell_upper),
        'margin':float(margin),'row_only_margin_same_region':float(lower-2*row_upper),
        'sum_full_cell_transfer_tv_sufficient':float(margin/2),
        'non_cluster_weight_lower':float(margin/(4*cell_upper)),
        'conditional_distinct_upper_by_row':list(map(float,weights)),
        'row_leaves':old_n,'cell_leaves':new_n,
        'full_singleton_region_witness_compatible':True,
        'scope':'Full singleton region plus scalar bunching; see full_models.py for the expanded region'}

def main():
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args()
    read=lambda n:json.loads((HERE/'results'/n).read_text())
    result=verify(read('cell-certificate.json'),read('counts.json'),read('physical-witness.json'))
    text=json.dumps(result,indent=2,sort_keys=True)+'\n';out=HERE/'results/cell-bounds.json'
    if a.check:assert out.read_text()==text
    else:out.write_text(text)
    print(text)

if __name__=='__main__':main()
