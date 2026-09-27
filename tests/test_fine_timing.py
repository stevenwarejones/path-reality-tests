"""Independent row, decoder, probability and additive-tree oracles."""
import importlib.util
import io
import os
import subprocess
import itertools
import json
from math import exp,expm1,log
from pathlib import Path
import sys
import unittest
import numpy as np
from scipy.optimize import linprog
HERE=Path(__file__).resolve().parents[1]/'studies/fine-timing-structure'
# Study scripts are standalone; prevent their generic module names leaking to other suites.
saved={name:sys.modules.get(name) for name in ('common','inference','extract')}
sys.path.insert(0,str(HERE))
try:
    def load(name,file):
        spec=importlib.util.spec_from_file_location(name,HERE/file); m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m
    common=load('fine_common_test','common.py'); sys.modules['common']=common
    inference=load('fine_inference_test','inference.py'); sys.modules['inference']=inference
    extract=load('fine_extract_test','extract.py'); sys.modules['extract']=extract
    recovery=load('fine_recovery_test','recovery.py')
    witnesses=common.load('fine_baseline_witness',common.PREVIOUS/'witnesses.py')
finally:
    sys.path.pop(0)
    for name,old in saved.items():
        if old is None: sys.modules.pop(name,None)
        else: sys.modules[name]=old


class FineRepresentationTests(unittest.TestCase):
    def test_domain_and_nested_partitions(self):
        from math import ceil
        self.assertEqual((ceil(799.5*129102/800),ceil(800.5*129102/800)-1),(129022,129182))
        self.assertLess(129182/800,common.PHASE_CELLS)
        self.assertEqual(common.CDF_LABELS,3888); self.assertEqual(common.PARTITION_LABELS,6512)
        for side,maps in common.MAPS.items():
            for coarse,fine in zip(maps,maps[1:]):
                for i in np.unique(fine): self.assertEqual(len(np.unique(coarse[fine==i])),1)
            self.assertEqual(len(np.unique(maps[-1])),486)

    def test_boundaries_ties_first_order_multiplicity_and_tail(self):
        pk=90; phases=np.array([np.nextafter(90.,0.),90.,91.999,92.,161.47,3.,4.])
        events=dict(row=np.array([0,1,2,3,4,4,5]),pulse=np.array([32,33,34,32,34,32,35]),phase=phases)
        rows,cells,bits,phase,count,multi=extract.select(events,dict(pk=pk,bitoffset=28),0,6)
        self.assertEqual(count,6); self.assertEqual(multi,1)
        self.assertEqual(rows.tolist(),[0,1,2,3,4])
        self.assertEqual(bits[-1],6) # raw first, not smallest pulse or phase
        np.testing.assert_array_equal(phase,phases[:5])
        self.assertEqual(common.MAPS['alice'][0][cells].tolist(),[0,1,1,2,2])
        bad=dict(events); bad['phase']=phases.copy(); bad['phase'][0]=162
        with self.assertRaises(ValueError): extract.select(bad,dict(pk=90,bitoffset=28),0,6)

    def test_audited_decoder_chunk_order_and_censored_tail(self):
        config=dict(pk=90,radius=4,bitoffset=28); records=[]
        for i in range(5):
            t=i*129102; records.extend([(6,t,0),(2,t+1,0)])
            if i<4:
                for pulse,phase in [(32,0),(33,90),(33,90),(34,161)]:
                    records.append((0,t+int(round(pulse*129102/800+phase)),0))
        raw=np.array(records,dtype=common.decoder.DTYPE)
        results=[]
        for chunk in (1,7,1000):
            previous=None; base=0; got=[]; tails=0
            for a,next_tag,tail in common.decoder.closed_blocks(io.BytesIO(raw.tobytes()),chunk):
                if tail: tails+=1; continue
                settings,_,_,previous,_,events=common.decoder.decode_block(a,next_tag,previous,config,include_events=True)
                sy=np.flatnonzero(a['ch']==6); det=np.flatnonzero(a['ch']==0)
                for j,k in enumerate(det):
                    row=np.searchsorted(sy,k,side='right')-1
                    delay=int(a['ttag'][k])-int(a['ttag'][sy[row]])
                    self.assertEqual(events['pulse'][j],int(np.floor(delay/(129102/800))))
                    self.assertEqual(events['phase'][j],delay%(129102/800))
                    got.append((base+row,int(events['pulse'][j]),float(events['phase'][j])))
                base+=len(settings)
            self.assertEqual(tails,1); results.append(got)
        self.assertEqual(results[0],results[1]); self.assertEqual(results[0],results[2])

    def test_independent_row_oracle_and_chunking(self):
        rng=np.random.default_rng(38); n=301; rows=np.arange(n)
        local=rng.choice([0,1,2,3],n); remote=rng.choice([0,1,2,3],n)
        bad=rng.random(n)<.03; cells=rng.integers(-1,486,n); detector=rng.random(n)<.7
        cells[~detector]=-1
        shapes=[(2,2,2,486),(2,2,2),(2,2,2,2),(16,2,486),(16,2)]
        oracle=[np.zeros(s,dtype='i8') for s in shapes]
        for i in range(64,n-64):
            h=int(i>=n//2); block=(i-64)*16//(n-128)
            uncertain=bad[i] or local[i] not in (1,2) or remote[i] not in (1,2)
            if not uncertain:
                x,y=remote[i]-1,local[i]-1; oracle[1][h,x,y]+=1; oracle[4][block,y]+=1
                if cells[i]>=0: oracle[0][h,x,y,cells[i]]+=1; oracle[3][block,y,cells[i]]+=1
            else:
                for x,y in itertools.product((0,1),repeat=2):
                    if (remote[i] not in (1,2) or remote[i]==x+1) and (local[i] not in (1,2) or local[i]==y+1):
                        oracle[2][h,x,y,0]+=1; oracle[2][h,x,y,1]+=int(detector[i])
        for chunk in (1,17,1000):
            total=[np.zeros(s,dtype='i8') for s in shapes]
            for start in range(0,n,chunk):
                sl=slice(start,start+chunk)
                got=extract.summarize(rows[sl],local[sl],remote[sl],bad[sl],cells[sl],detector[sl],n)
                for a,b in zip(total,got): a+=b
            for a,b in zip(total,oracle): np.testing.assert_array_equal(a,b)

    def test_committed_fine_to_coarse_counts(self):
        for side in common.SIDES:
            fine=json.loads((HERE/f'{side}-counts.json').read_text())
            old=json.loads((common.PREVIOUS/f'{side}-counts.json').read_text())
            extract.compare_coarse(fine,old)
            self.assertEqual(np.array(fine['counts']).sum(),np.array(fine['pooled_block_counts']).sum())


class FineInferenceTests(unittest.TestCase):
    def test_exponential_process_with_adaptive_outcomes(self):
        # Exhaust all paths; each next Bernoulli probability depends on its past.
        for stake in (-2.,-.16,.16,2.):
            expectation=0.
            for bits in itertools.product((0,1),repeat=7):
                probability=1.; compensator=0
                for i,bit in enumerate(bits):
                    mu=(.25+.01*(1 if sum(bits[:i])%2 else -1))*(.1+.7*(sum(bits[:i])%2))
                    probability*=mu if bit else 1-mu; compensator+=mu
                expectation+=probability*exp(stake*sum(bits)-expm1(stake)*compensator)
            self.assertLessEqual(expectation,1+1e-12)

    def test_count_boundary_and_pathwise_completions(self):
        n=10000; b=log(2*10*common.CDF_LABELS/.005)
        bounds=inference.count_bounds(np.array([0,20,600]),np.array([5,27,608]),n,common.CDF_LABELS)
        for k,u,(lo,hi) in zip((0,20,600),(5,27,608),bounds):
            for v in common.LAMBDAS:
                self.assertLessEqual(v*k-expm1(v)*lo,b+1e-9)
                self.assertLessEqual(-v*u-expm1(-v)*hi,b+1e-9)
            for completed in range(k,u+1):
                complete=inference.count_bounds(completed,completed,n,common.CDF_LABELS)
                self.assertLessEqual(lo,complete[0]+1e-12); self.assertGreaterEqual(hi,complete[1]-1e-12)

    def test_tree_projections_against_linear_program(self):
        side='alice'; rng=np.random.default_rng(10)
        law=rng.dirichlet(np.ones(486))*.02
        nodes=[]
        for mapping in common.MAPS[side]:
            mass=common.grouped(law,mapping)
            nodes.append(np.stack([np.maximum(0,mass-rng.random(len(mass))*.0001),mass+rng.random(len(mass))*.0001],axis=-1))
        root=np.array([.019,.021]); original=[v.copy() for v in nodes]
        tightened,total=inference.tighten_tree(nodes,root.copy(),side)
        mat=[]; rhs=[]
        for mapping,node in zip(common.MAPS[side],original):
            for i,(lo,hi) in enumerate(node):
                indicator=(mapping==i).astype(float); mat.extend([indicator,-indicator]); rhs.extend([hi,-lo])
        mat.extend([np.ones(486),-np.ones(486)]); rhs.extend([root[1],-root[0]])
        for level,index in [(0,0),(1,4),(4,200)]:
            cost=(common.MAPS[side][level]==index).astype(float)
            l=linprog(cost,A_ub=mat,b_ub=rhs,bounds=(0,None),method='highs')
            u=linprog(-cost,A_ub=mat,b_ub=rhs,bounds=(0,None),method='highs')
            self.assertTrue(l.success and u.success)
            np.testing.assert_allclose(tightened[level][index],[l.fun,-u.fun],atol=1e-8)

    def test_incompatible_tree_flag(self):
        nodes=[np.zeros((int(m.max())+1,2)) for m in common.MAPS['alice']]
        for v in nodes: v[:,1]=1
        nodes[-1][0,0]=.8; nodes[-1][1,0]=.8
        with self.assertRaises(ValueError): inference.tighten_tree(nodes,np.array([0.,1.]),'alice')

    def test_equal_coarse_and_detectable_fine(self):
        n=100_000_000; counts=np.zeros((2,486),dtype='i8')
        counts[0,252]=25000; counts[1,253]=25000
        a=inference.cdf(counts,np.zeros(2,int),n); b=inference.partitions(counts,np.zeros(2,int),n,'alice')
        self.assertEqual(b['tv'][0,0],0); self.assertGreater(a['k_grid'][0],0)
        self.assertGreater(b['gain'][-1,0],0)
        self.assertEqual(witnesses.quantization()['recorded_tv'],'0')

    def test_missingness_assignment_monotonic_and_batch(self):
        c=np.zeros((2,486),dtype='i8'); c[0,252]=1000; c[1,253]=1000; n=10_000_000
        for fun in (lambda c,u,e:inference.cdf(c,u,n,e)['k_grid'],lambda c,u,e:inference.partitions(c,u,n,'bob',e)['tv']):
            tight=fun(c,np.zeros(2,int),0); wide=fun(c,np.full(2,200),.01)
            self.assertTrue(np.all(wide[...,0]<=tight[...,0]+1e-12)); self.assertTrue(np.all(wide[...,1]>=tight[...,1]-1e-12))
            np.testing.assert_allclose(fun(np.stack([c,c]),np.zeros((2,2),int),0),np.stack([tight,tight]))

    def test_mass_cap_and_monotone_full_grid_bridge(self):
        c=np.zeros((2,486),dtype='i8'); c[:,252]=400
        a=inference.cdf(c,np.zeros(2,int),1_000_000)
        self.assertGreaterEqual(a['k_full'][1],a['k_grid'][1])
        b=inference.partitions(c,np.zeros(2,int),1_000_000,'alice')
        self.assertLess(b['full_tv'][1],.01)
        self.assertTrue(np.all(np.diff(b['tv'][:,1])>=-1e-12))

    def test_transform_truth_and_tail_mass(self):
        p=np.full(486,1e-6)
        for side in common.SIDES:
            a,b=[recovery.transform(p,side,'within',.5,x) for x in (0,1)]
            np.testing.assert_allclose(common.grouped(a,common.MAPS[side][0]),common.grouped(b,common.MAPS[side][0]))
            self.assertAlmostEqual(recovery.truths(a,b,side)['gain'][-1],.5*p[common.MAPS[side][0]==1].sum())
            shifted=recovery.transform(p,side,'translation',2,1)
            self.assertAlmostEqual(shifted.sum(),p.sum())
            self.assertAlmostEqual(shifted[161],3e-6)
            scaled=recovery.transform(p,side,'rate',.5,1)
            np.testing.assert_allclose(scaled/scaled.sum(),p/p.sum())

    def test_subcell_difference_is_not_mistaken_for_grid_equality(self):
        # Opposite point masses inside the same grid cell have true full K=p.
        # Equal cell counts cannot detect that, but the bridge must allow it.
        n=100_000_000; c=np.zeros((2,486),dtype='i8'); c[:,252]=25000
        a=inference.cdf(c,np.zeros(2,int),n)
        self.assertLess(a['k_grid'][1],.001)
        self.assertGreaterEqual(a['k_full'][1],.001)
        self.assertEqual(a['k_grid'][0],0)

    def test_cpu_dispatch_does_not_change_serialized_bounds(self):
        # Reproduces the differing local/full-source runners without tolerating drift.
        code="""import json, numpy as np
from inference import cdf, partitions
c=np.zeros((2,486),dtype='i8'); c[0,252]=400; c[1,253]=521
u=np.array([17,23]); n=107109468
r=[cdf(c,u,n)['band'].tolist(),partitions(c,u,n,'alice')['tv'].tolist()]
print(json.dumps(r,sort_keys=True))
"""
        base=os.environ.copy(); base.pop('NPY_DISABLE_CPU_FEATURES',None)
        disabled=dict(base,NPY_DISABLE_CPU_FEATURES='AVX512F,AVX2,FMA3')
        outputs=[subprocess.check_output([sys.executable,'-c',code],cwd=HERE,env=e) for e in (base,disabled)]
        self.assertEqual(*outputs)

    def test_invalid_shapes_and_assumptions(self):
        with self.assertRaises(ValueError): inference.count_bounds(2,1,10,10)
        with self.assertRaises(ValueError): inference.cdf(np.zeros((2,10)),np.zeros(2),100)
        with self.assertRaises(ValueError): inference.probabilities(0,0,100,.25,10)

if __name__=='__main__': unittest.main()
