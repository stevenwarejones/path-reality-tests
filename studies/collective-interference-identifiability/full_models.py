#!/usr/bin/env python3
"""Exact physical and complete parity-histogram feasibility checks (stdlib only)."""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations,permutations
from math import comb,factorial
from pathlib import Path
import argparse,json
from certificate import check_endpoint
from cell_certificate import verify as verify_cells
HERE=Path(__file__).resolve().parent
PERMS=list(permutations(range(3)))
K=sum(comb(144,k) for k in range(4))
ZERO=(0,0)
def add(z,w):return z[0]+w[0],z[1]+w[1]
def mul(z,w):return z[0]*w[0]-z[1]*w[1],z[0]*w[1]+z[1]*w[0]
def conj(z):return z[0],-z[1]
def norm(z):return z[0]*z[0]+z[1]*z[1]
def gram_complement(T,d):
    G=[]
    for i in range(3):
        row=[]
        for j in range(3):
            z=ZERO
            for t in T:z=add(z,mul(conj(t[i]),t[j]))
            row.append((d*d*int(i==j)-z[0],-z[1]))
        G.append(row)
    return G

def positive_definite(G):
    assert all(G[i][j]==conj(G[j][i]) for i in range(3) for j in range(3))
    a,b,c=[G[i][i][0] for i in range(3)]
    minor=a*b-norm(G[0][1])
    det=a*b*c-a*norm(G[1][2])-b*norm(G[0][2])-c*norm(G[0][1])+2*mul(mul(G[0][1],G[1][2]),G[2][0])[0]
    assert a>0 and minor>0 and det>0,'Not a physical contraction'
    return [a,minor,det]

def transport(model):
    A=model['amplitude_denominator'];R=model['phase_radius'];d=A*R
    assert isinstance(A,int) and A>0 and isinstance(R,int) and R>0
    base=model['base_amplitude_numerators'];ph=model['row_phase_numerators']
    assert len(base)==len(ph)==28 and all(len(v)==12 for v in base)
    assert all(isinstance(v,int) and v>=0 for row in base for v in row)
    assert all(len(row)==3 and all(len(z)==2 and all(isinstance(v,int) for v in z) and norm(z)==R*R for z in row) for row in ph)
    assert all(v==0 for y in [0,27] for v in base[y]),'Boundary mass needs explicit translation accounting'
    T=[]
    for y in range(28):
        for x in range(12):
            row=[]
            for j,s in enumerate([0,1,-1]):
                a=base[y-s][x] if 0<=y-s<28 else 0
                row.append((a*ph[y][j][0],a*ph[y][j][1]))
            T.append(row)
    minors=positive_definite(gram_complement(T,d))
    return T,d,[F(v,d**(2*(i+1))) for i,v in enumerate(minors)]

class ParityModel:
    def __init__(self,T,d,kind):
        assert kind in ['boson','cluster','distinguishable']
        self.T=T;self.d=d;self.den=6*d**6
        G=gram_complement(T,d);positive_definite(G)
        def j(s,t):
            return kind=='boson' or s==t or (kind=='cluster' and tuple(1-v if v<2 else v for v in s)==t)
        self.L={}
        for m in range(4):
            L=[]
            for s in PERMS:
                for t in PERMS:
                    if not j(s,t):L.append(ZERO);continue
                    v=(1,0)
                    for k in range(m,3):v=mul(v,G[t[k]][s[k]])
                    L.append(v)
            self.L[m]=L
    def occupation(self,outs):
        """Integer numerator at common denominator 6*d^6; traces all loss modes."""
        m=len(outs);V=[]
        for s in PERMS:
            v=(1,0)
            for i,o in enumerate(outs):v=mul(v,self.T[o][s[i]])
            V.append(v)
        total=ZERO
        for i,v in enumerate(V):
            for j,w in enumerate(V):
                l=self.L[m][6*i+j]
                if l!=ZERO:total=add(total,mul(mul(v,conj(w)),l))
        assert total[1]==0 and total[0]>=0
        f=factorial(3-m)
        for n in Counter(outs).values():f*=factorial(n)
        assert 6%f==0
        return total[0]*(6//f)
    def parity(self,pat):
        if len(pat)==0:return self.occupation(())+sum(self.occupation((x,x)) for x in range(len(self.T)))
        if len(pat)==1:
            i=pat[0]
            return self.occupation((i,))+self.occupation((i,i,i))+sum(self.occupation((i,k,k)) for k in range(len(self.T)) if k!=i)
        return self.occupation(pat)
    def bunch(self):
        return sum(self.occupation(tuple(12*y+x for x in xs)) for y in range(12) for xs in combinations(range(12),3))

def histogram(h,counts):
    assert h['schema']==1 and h['mode_count']==144
    patterns={}
    for p,k in h['patterns']:
        p=tuple(p)
        assert p==tuple(sorted(set(p))) and len(p)<=3 and all(isinstance(v,int) and 0<=v<144 for v in p)
        assert isinstance(k,int) and k>0 and p not in patterns
        patterns[p]=k
    assert sum(patterns.values())==counts['many_shots']
    assert sum(k for p,k in patterns.items() if len(p)==3 and len({x//12 for x in p})==1)==counts['bunch_shots']
    return patterns

def single_region(model,cert):
    A=model['amplitude_denominator'];base=model['base_amplitude_numerators']
    q=[[F(v*v,A*A) for v in row] for row in base]
    for y in range(9,21):
        mass=sum(q[y]);l,u=map(F,cert['row_intervals'][y-9]);assert l<=mass<=u
        if mass:
            for x in range(12):
                l,u=map(F,cert['conditional_intervals'][y-9][x]);assert l<=q[y][x]/mass<=u
    l,u=map(F,cert['row_intervals'][-1]);assert l<=1-sum(sum(q[y]) for y in range(9,21))<=u

def verify(model,h,region,counts,cell_cert,check_histogram=True):
    patterns=histogram(h,counts);N=counts['many_shots'];T,d,minors=transport(model)
    assert region['schema']==1 and region['alphabet_size']==K
    intervals={int(k):tuple(map(F,p)) for k,p in region['pattern_intervals'].items()}
    assert set(intervals)==set(patterns.values())|{0}
    for k,(l,u) in intervals.items():
        assert l<=u
        check_endpoint(N,k,l,F(1,640*K),'lower');check_endpoint(N,k,u,F(1,640*K),'upper')
    lower=F(region['bunch_lower']);check_endpoint(N,counts['bunch_shots'],lower,F(1,320),'lower')
    assert lower>2*F(cell_cert['cell_upper'])
    obs=T[8*12:20*12];p=ParityModel(obs,d,model['kind']);pb=F(p.bunch(),p.den);assert pb>=lower
    if model['kind']=='boson':single_region(model,cell_cert)
    if check_histogram:
        # Every empty, one- and two-click outcome, observed or not.
        checked=set()
        for m in range(3):
            for pat in combinations(range(144),m):
                pr=F(p.parity(pat),p.den);l,u=intervals[patterns.get(pat,0)];assert l<=pr<=u,(pat,pr,l,u)
                checked.add(pat)
        for pat,k in patterns.items():
            if len(pat)==3:
                pr=F(p.parity(pat),p.den);l,u=intervals[k];assert l<=pr<=u,(pat,pr,l,u)
                checked.add(pat)
        # For any unobserved three-click set, p <= |H| sum_sigma prod p_sj.
        maxima=sorted((max(norm(z) for z in row) for row in obs),reverse=True)
        factor={'boson':6,'cluster':2,'distinguishable':1}[model['kind']]
        zero_upper=F(factor*6*maxima[0]*maxima[1]*maxima[2],d**6)
        refined=0
        if zero_upper>intervals[0][1]:
            probabilities=[[norm(z) for z in row] for row in obs]
            u=intervals[0][1];max_num=0
            for pat in combinations(range(144),3):
                if pat in patterns:continue
                rows=[probabilities[o] for o in pat]
                bound=6*factor*sum(rows[0][s[0]]*rows[1][s[1]]*rows[2][s[2]] for s in PERMS)
                if bound*u.denominator>p.den*u.numerator:
                    bound=p.occupation(pat);refined+=1
                    assert bound*u.denominator<=p.den*u.numerator,(pat,'unobserved triple')
                max_num=max(max_num,bound)
            zero_upper=F(max_num,p.den)
    return {'kind':model['kind'],'bunch_probability':float(pb),
        'physical_principal_minors':list(map(float,minors)),
        'full_parity_region_compatible':check_histogram,
        'full_singleton_region_compatible':model['kind']=='boson',
        'bunch_lower':float(lower),'cluster_ceiling':float(2*F(cell_cert['cell_upper'])),
        'margin':float(lower-2*F(cell_cert['cell_upper'])),
        'non_cluster_weight_lower':float((lower-2*F(cell_cert['cell_upper']))/(4*F(cell_cert['cell_upper']))),
        'sum_transfer_tv_sufficient':float((lower-2*F(cell_cert['cell_upper']))/2),
        'unobserved_triple_probability_upper':float(zero_upper) if check_histogram else None,
        'individually_refined_unobserved_triples':refined if check_histogram else None,
        'alphabet_size':K,'observed_patterns':len(patterns)}

def main():
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args()
    read=lambda n:json.loads((HERE/'results'/n).read_text())
    counts=read('counts.json');cert=read('cell-certificate.json')
    verify_cells(cert,counts,read('physical-witness.json'))
    result={kind:verify(read('full-'+kind+'.json'),read('parity-histogram.json'),read('full-region.json'),counts,cert) for kind in ['boson','cluster']}
    out=HERE/'results/full-model-check.json';text=json.dumps(result,sort_keys=True,indent=2)+'\n'
    if a.check:assert out.read_text()==text
    else:out.write_text(text)
    print(text)
if __name__=='__main__':main()
