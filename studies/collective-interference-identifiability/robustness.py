#!/usr/bin/env python3
"""Exact nuisance bounds, shared-drift model, row-hit and matched-control checks."""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
from math import comb,lcm
import argparse,json
import full_models as fm
from certificate import check_endpoint
from cell_certificate import verify as verify_cells
HERE=Path(__file__).resolve().parent

def preparation_bound(U,beta,gamma,d):return U+beta*d+gamma*d*d/3+d**3/27

def drift_bound(U,gamma,d):return U+2*gamma*d*d/3+4*d**3/27

def check_chernoff(n,k,p,tail,side):
    """Exact integer Chernoff tail envelope; avoids summing long binomial tails."""
    assert n>0 and 0<=k<=n and 0<=p<=1 and 0<tail<1
    if side=='lower':
        if k==0:assert p==0;return
        assert p<=F(k,n)
    else:
        assert side=='upper'
        if k==n:assert p==1;return
        assert p>=F(k,n)
    a,d=p.numerator,p.denominator
    num=a**k*(d-a)**(n-k)*n**n
    den=d**n*k**k*(n-k)**(n-k)
    assert num*tail.denominator<=den*tail.numerator

def single_probabilities(q,cert):
    assert len(q)==28 and all(len(row)==12 and all(p>=0 for p in row) for row in q)
    for y in range(9,21):
        mass=sum(q[y]);lo,hi=map(F,cert['row_intervals'][y-9]);assert lo<=mass<=hi,('calibration row',y)
        if mass:
            for x in range(12):
                lo,hi=map(F,cert['conditional_intervals'][y-9][x]);assert lo<=q[y][x]/mass<=hi,('calibration cell',y,x)
    lo,hi=map(F,cert['row_intervals'][-1]);assert lo<=1-sum(sum(row) for row in q[9:21])<=hi

def generic_transport(model):
    A=model['amplitude_denominator'];R=model['phase_radius'];amp=model['amplitudes'];ph=model['phases'];W=model['width']
    assert isinstance(A,int) and A>0 and isinstance(R,int) and R>0 and W>0
    assert len(amp)==len(ph) and all(len(row)==W and all(len(v)==3 for v in row) for row in amp)
    assert all(isinstance(a,int) and a>=0 for row in amp for v in row for a in v)
    assert all(len(row)==3 and all(len(z)==2 and all(isinstance(v,int) for v in z) and fm.norm(z)==R*R for z in row) for row in ph)
    T=[[(amp[y][x][j]*ph[y][j][0],amp[y][x][j]*ph[y][j][1]) for j in range(3)] for y in range(len(amp)) for x in range(W)]
    d=A*R;fm.positive_definite(fm.gram_complement(T,d));return T,d

def column_probabilities(T,d,j):return [F(fm.norm(row[j]),d*d) for row in T]+[1-sum(F(fm.norm(row[j]),d*d) for row in T)]
def tv(a,b):assert len(a)==len(b);return sum(abs(x-y) for x,y in zip(a,b))/2

def pattern_check(components,patterns,intervals,W,lower):
    """All parity outcomes for a convex mixture of certified physical channels."""
    M=W*W;assert all(len(T)==M and kind in ['cluster','boson'] for w,T,d,kind in components)
    assert sum(w for w,T,d,kind in components)==1 and all(w>=0 for w,T,d,kind in components)
    ps=[fm.ParityModel(T,d,kind) for w,T,d,kind in components]
    den=lcm(*(p.den*w.denominator for p,(w,*_) in zip(ps,components)))
    factors=[den//(p.den*w.denominator)*w.numerator for p,(w,*_) in zip(ps,components)]
    def numerator(pat):return sum(a*p.parity(pat) for a,p in zip(factors,ps))
    bunch=sum(a*sum(p.occupation(tuple(W*y+x for x in xs)) for y in range(W) for xs in combinations(range(W),3)) for a,p in zip(factors,ps))
    assert F(bunch,den)>=lower
    for n in range(3):
        for pat in combinations(range(M),n):
            pr=F(numerator(pat),den);lo,hi=intervals[patterns.get(pat,0)];assert lo<=pr<=hi,('parity',pat)
    for pat,k in patterns.items():
        if len(pat)==3:
            pr=F(numerator(pat),den);lo,hi=intervals[k];assert lo<=pr<=hi,('observed triple',pat)
    # An analytic bound checks every unobserved triple, including bins absent from the histogram.
    upper=F(0);probs=[];orders=[]
    for w,T,d,kind in components:
        q=[[fm.norm(z) for z in row] for row in T];probs.append(q);factor=6 if kind=='boson' else 2;orders.append(factor)
        maxima=sorted((max(row) for row in q),reverse=True)
        upper+=w*F(6*factor*maxima[0]*maxima[1]*maxima[2],d**6)
    refined=0
    if upper>intervals[0][1]:
        hi=intervals[0][1];maxnum=0
        for pat in combinations(range(M),3):
            if pat in patterns:continue
            bound=0
            for coeff,q,factor in zip(factors,probs,orders):
                a,b,c=[q[i] for i in pat]
                bound+=coeff*6*factor*sum(a[s[0]]*b[s[1]]*c[s[2]] for s in fm.PERMS)
            if bound*hi.denominator>den*hi.numerator:
                bound=numerator(pat);refined+=1
                assert bound*hi.denominator<=den*hi.numerator,('unobserved triple',pat)
            maxnum=max(maxnum,bound)
        upper=F(maxnum,den)
    return {'bunch_probability':float(F(bunch,den)),'all_pattern_constraints_pass':True,
        'unobserved_triple_upper':float(upper),'refined_unobserved_triples':refined}

def row_hits(components,patterns,region):
    result=[];N=sum(patterns.values())
    assert len(region['row_hit_intervals'])==12
    for y in range(12):
        k=sum(n for pat,n in patterns.items() if any(o//12==y for o in pat))
        lo,hi=map(F,region['row_hit_intervals'][y]);check_chernoff(N,k,lo,F(1,15360),'lower');check_chernoff(N,k,hi,F(1,15360),'upper')
        pr=sum(w*(1-F(fm.ParityModel(T[12*y:12*(y+1)],d,kind).parity(()),6*d**6)) for w,T,d,kind in components)
        result.append({'row':y+8,'count':k,'lower':float(lo),'upper':float(hi),'prediction':float(pr),'compatible':lo<=pr<=hi})
    return result

def selected_verify(read):
    cert=read('cell-certificate.json');counts=read('counts.json');verify_cells(cert,counts,read('physical-witness.json'))
    patterns=fm.histogram(read('parity-histogram.json'),counts);reg=read('robust-region.json');K=fm.K
    assert reg['schema']==1 and reg['alphabet_size']==K
    iv={int(k):tuple(map(F,p)) for k,p in reg['pattern_intervals'].items()};assert set(iv)==set(patterns.values())|{0}
    for k,(lo,hi) in iv.items():
        check_endpoint(counts['many_shots'],k,lo,F(1,1280*K),'lower');check_endpoint(counts['many_shots'],k,hi,F(1,1280*K),'upper')
    lower=F(reg['bunch_lower']);check_endpoint(counts['many_shots'],counts['bunch_shots'],lower,F(1,320),'lower')
    models=[read('drift-minus.json'),read('drift-plus.json')];Ts=[fm.transport(m)[:2] for m in models]
    q=[[sum(F(m['base_amplitude_numerators'][y][x]**2,m['amplitude_denominator']**2) for m in models)/2 for x in range(12)] for y in range(28)]
    single_probabilities(q,cert)
    means=[[sum(column_probabilities(T,d,j)[i] for T,d in Ts)/2 for i in range(337)] for j in range(3)]
    drift=max(sum(tv(column_probabilities(T,d,j),means[j]) for j in range(3)) for T,d in Ts)
    comps=[(F(1,2),T[96:240],d,'cluster') for T,d in Ts]
    result=pattern_check(comps,patterns,iv,12,lower);hits=row_hits(comps,patterns,reg);assert all(r['compatible'] for r in hits)
    result.update({'shared_law':'Two equiprobable channels, identical acquisition law in both preparation families.',
        'mean_singleton_region_pass':True,'row_hit_constraints_pass':True,'drift_budget_sum_tv':float(drift),'row_hits':hits})
    # An event-specific polynomial envelope, uniform over the declared nuisance classes.
    u={y+9:F(p[1]) for y,p in enumerate(cert['row_intervals'][:12])};u[7]=u[8]=F(cert['row_intervals'][-1][1])
    beta=max(u[y+a]*u[y+b] for y in range(8,20) for a,b in [(-1,0),(-1,1),(0,1)])
    gamma=max(u.values());U=F(cert['cell_upper'])
    prep=F('0.00554');fluct=F('0.0397')
    assert 2*preparation_bound(U,beta,gamma,prep)<lower
    assert 2*drift_bound(U,gamma,fluct)<lower
    result['certified_bounds']={'beta':float(beta),'gamma':float(gamma),'fixed_preparation_change_excluded_through_sum_tv':float(prep),
        'shared_zero_mean_drift_excluded_through_sum_tv':float(fluct),'drift_critical_upper':float(drift),
        'fixed_preparation_change_upper':'Not established for the augmented row-hit region.'}
    # Existing quantum example and old higher-only focusing explanation on the expanded region.
    rows={};control={}
    for kind in ['boson','cluster']:
        T,d,_=fm.transport(read('full-'+kind+'.json'));cs=[(F(1),T[96:240],d,kind)]
        rows[kind]=row_hits(cs,patterns,reg)
        if kind=='boson':assert all(r['compatible'] for r in rows[kind]);pattern_check(cs,patterns,iv,12,lower)
    result['old_focusing_rejected_rows']=[r['row'] for r in rows['cluster'] if not r['compatible']]
    result['quantum_row_hits']=rows['boson']
    # Proposed matched control: orthogonally labelled inputs 0 and 1, two distinct hits in rows 9..13.
    def pair_control(T,d):
        subset=T[9*12:14*12];a=[F(fm.norm(t[0]),d*d) for t in subset];b=[F(fm.norm(t[1]),d*d) for t in subset]
        return sum(a)*sum(b)-sum(x*y for x,y in zip(a,b))
    Tq,dq,_=fm.transport(read('full-boson.json'))
    result['proposed_distinguishable_pair_control']={'fixed_quantum_channel':float(pair_control(Tq,dq)),
        'shared_drift_cluster':float(sum(pair_control(T,d) for T,d in Ts)/2),
        'status':'Prediction for a new matched acquisition; not an existing control at this setting.'}
    return result

def matched_verify(case,model):
    W=case['width'];M=W*W;N=case['many_shots'];K=sum(comb(M,k) for k in range(4))
    T,d=generic_transport(model);assert len(T)==M and model['width']==W
    row_uppers=[];cell_uppers=[];empirical=[]
    for j,record in enumerate(case['singletons']):
        n=record['shots'];cts=record['counts'];assert len(cts)==M+1 and sum(cts)==n
        iv={int(k):tuple(map(F,p)) for k,p in record['intervals'].items()}
        assert set(iv)==set(cts)
        for k,(lo,hi) in iv.items():check_endpoint(n,k,lo,F(1,7680*(M+1)),'lower');check_endpoint(n,k,hi,F(1,7680*(M+1)),'upper')
        probs=column_probabilities(T,d,j)
        for pr,k in zip(probs,cts):assert iv[k][0]<=pr<=iv[k][1]
        rc=[sum(cts[y*W:(y+1)*W]) for y in range(W)]+[cts[-1]]
        rp=[sum(probs[y*W:(y+1)*W]) for y in range(W)]+[probs[-1]]
        assert len(record['row_intervals'])==W+1
        ri=[tuple(map(F,p)) for p in record['row_intervals']]
        for k,pr,(l,u) in zip(rc,rp,ri):
            check_endpoint(n,k,l,F(1,7680*(W+1)),'lower');check_endpoint(n,k,u,F(1,7680*(W+1)),'upper');assert l<=pr<=u,('matched row',case['time_ms'],j,float(pr),float(l),float(u))
        row_uppers.append([u for l,u in ri[:-1]])
        cell_uppers.append([iv[k][1] for k in cts[:M]])
        empirical.append([F(k,n) for k in cts[:M]])
    pats={tuple(p):k for p,k in case['patterns']};assert len(pats)==len(case['patterns']) and sum(pats.values())==N
    assert all(tuple(sorted(set(p)))==p and len(p)<=3 and all(0<=v<M for v in p) and k>0 for p,k in pats.items())
    B=sum(k for p,k in pats.items() if len(p)==3 and len({i//W for i in p})==1);assert B==case['bunch_count']
    lo=F(case['bunch_lower']);check_endpoint(N,B,lo,F(1,1280),'lower')
    iv={int(k):tuple(map(F,p)) for k,p in case['pattern_intervals'].items()};assert set(iv)==set(pats.values())|{0}
    for k,(l,u) in iv.items():check_endpoint(N,k,l,F(1,2560*K),'lower');check_endpoint(N,k,u,F(1,2560*K),'upper')
    def event(prob):
        total=F(0)
        for y in range(W):
            a,b,c=[p[y*W:(y+1)*W] for p in prob];A,B,C=map(sum,[a,b,c])
            total+=A*B*C-C*sum(x*z for x,z in zip(a,b))-B*sum(x*z for x,z in zip(a,c))-A*sum(x*z for x,z in zip(b,c))+2*sum(x*z*t for x,z,t in zip(a,b,c))
        return total
    upper=min(F(1),event(cell_uppers),sum(row_uppers[0][y]*row_uppers[1][y]*row_uppers[2][y] for y in range(W)))
    result=pattern_check([(F(1),T,d,'cluster')],pats,iv,W,lo)
    result.update({'empirical_distinguishable_event':float(event(empirical)),
        'certified_distinguishable_event_upper':float(upper),'certified_cluster_ceiling':float(min(F(1),2*upper))})
    result.update({'time_ms':case['time_ms'],'many_shots':N,'bunch_count':B,'bunch_lower':float(lo),'independent_singleton_shots':[r['shots'] for r in case['singletons']],
        'full_independent_singleton_region_pass':True,'translation_assumption':False})
    return result

def main():
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();read=lambda n:json.loads((HERE/'results'/n).read_text())
    result={'selected_shared_drift':selected_verify(read),'independently_calibrated_settings':[matched_verify(case,read(case['model_file'])) for case in read('matched-regions.json')['cases']]}
    text=json.dumps(result,sort_keys=True,indent=2)+'\n';out=HERE/'results/robustness-check.json'
    if a.check:assert out.read_text()==text
    else:out.write_text(text)
    print(text)
if __name__=='__main__':main()
