"""Independent decoder, missingness, score identity, and source reconciliation checks."""
import importlib.util
import io
import itertools
import json
from pathlib import Path
import sys
import unittest

import numpy as np

STUDY=Path(__file__).resolve().parents[1]/'studies/nist-bell-causal-audit'
sys.path.insert(0,str(STUDY))
# Other studies have modules with similarly generic names. Load this one explicitly.
def load(name,file):
    spec=importlib.util.spec_from_file_location(name,STUDY/file)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module
inf=load('nist_inference','inference.py')
rec=load('nist_reconstruct','reconstruct.py')
audit=load('nist_audit','audit.py')


class DecoderTests(unittest.TestCase):
    @staticmethod
    def synthetic_events():
        # Invented times: exact 161.375-bin laser period, 800 pulses per sync.
        # Two distinct eligible bits, plus a duplicate, in the second trial.
        a=np.array([(6,0,0),(2,1,0),(6,129100,0),(4,128001,0),
                    (0,129100+90+round(28*161.375),0),(0,129100+90+round(29*161.375),0),
                    (0,129100+90+round(29*161.375),0)],dtype=rec.DTYPE)
        return a,258200,{'pk':90,'radius':4,'bitoffset':28}

    def test_synthetic_multiclick_semantics(self):
        a,tag,config=self.synthetic_events()
        s,w,legacy,_,_=rec.decode_block(a,tag,None,config)
        self.assertEqual(s.tolist(),[1,2])
        self.assertEqual(w['nominal'].tolist(),[0,3])
        self.assertEqual(legacy.tolist(),[0,2])

    def test_chunk_boundaries_preserve_records_and_censored_tail(self):
        a,tag,config=self.synthetic_events();data=a.tobytes()
        tail=np.array([(6,tag,0)],dtype=rec.DTYPE).tobytes()
        for chunk in (1,2,3,7,100):
            blocks=list(rec.closed_blocks(io.BytesIO(data+tail),chunk))
            self.assertEqual(b''.join(a.tobytes() for a,_,_ in blocks),data+tail)
            self.assertTrue(blocks[-1][2]);self.assertEqual(len(blocks[-1][0]),1)
            words=[];previous=None
            for a,tag,censored in blocks:
                if not censored:
                    _,w,_,previous,_=rec.decode_block(a,tag,previous,config)
                    words.extend(w['nominal'].tolist())
            self.assertEqual(words,[0,3])

    def test_truncated_and_unhandled_sync_fail(self):
        with self.assertRaises(ValueError):list(rec.closed_blocks(io.BytesIO(b'x')))
        a=np.array([(6,0,0)],dtype=rec.DTYPE)
        with self.assertRaises(ValueError):rec.decode_block(a,10,None,{'pk':90,'radius':4,'bitoffset':28})

    def test_patch_validation_and_no_click_conservation(self):
        old=np.array([0,16,0],dtype='u2')
        got=audit.apply_patches(old,10,np.array([[11,16,48]]))
        self.assertEqual(got.tolist(),[0,48,0])
        with self.assertRaises(ValueError):audit.apply_patches(old,10,np.array([[11,32,48]]))
        t=audit.counts(np.array([1,1,2,3]),np.array([1,2,1,2]),[0,0,1,0],[0,1,0,0])
        self.assertEqual(t.sum(),4);self.assertEqual(t[1,1,0,0],1);self.assertEqual(t[3,2,0,0],1)


class InferenceTests(unittest.TestCase):
    def test_uniform_joint_score_identity_exhaustive(self):
        # All 16 deterministic receiver response functions over the four settings.
        for response in itertools.product((0,1),repeat=4):
            for y in (0,1):
                mean=sum(2*(2*x-1)*(2*response[2*x+v]-1)/4 for x in (0,1) for v in (0,1) if v==y)
                self.assertEqual(mean,response[2+y]-response[y])

    def test_selection_creates_gap_from_a_complete_no_influence_table(self):
        # Fair latent bit R is receiver outcome, independent of setting X.
        # Sender clicks iff R=X. Complete B is fair; conditioning on A=1 gives B=X.
        a=np.array([1,1,2,2]);b=np.ones(4,dtype=int)
        receiver=np.array([0,1,0,1]);sender=np.array([1,0,0,1])
        t=audit.counts(a,b,sender,receiver)
        self.assertEqual(audit.contrast(t,t,4,'A_to_B',0)['descriptive_gap'],0)
        selected=audit.selected_contrast(t,'A_to_B',0)
        self.assertEqual(selected['descriptive_selected_gap'],1)
        self.assertFalse(selected['causal_interpretation'])

    def test_setting_tv_range_bound(self):
        q=np.array([.4,.1,.2,.3]);tv=.5*abs(q-.25).sum()
        for response in itertools.product((0,1),repeat=4):
            for y in (0,1):
                values=np.array([2*(2*x-1)*(2*response[2*x+v]-1) if v==y else 0 for x in (0,1) for v in (0,1)])
                self.assertLessEqual(abs(np.dot(q-.25,values)),4*tv+1e-15)

    def test_unknown_scores_cover_every_completion(self):
        interval=inf.score_interval(2,2,10,alpha=.01,comparisons=396,setting_tv=0,leakage=0)
        for unknown in itertools.product((-2,0,2),repeat=2):
            mean=(2+sum(unknown))/10
            self.assertLessEqual(abs(mean-interval['score_mean']),interval['unknown_radius']+1e-15)

    def test_missing_outcome_identification_by_enumeration(self):
        lo,hi=inf.missing_outcome_difference(1,5,2,2,6,3)
        completions=[(2+b)/6-(1+a)/5 for a in range(3) for b in range(4)]
        self.assertEqual((lo,hi),(min(completions),max(completions)))

    def test_interval_error_allocation_telescopes(self):
        # Finite partial sums are exactly below their infinite total, one.
        from fractions import Fraction
        total=sum(Fraction(1,n*(n+1)) for n in range(1,101))
        self.assertEqual(total,Fraction(100,101))
        base=dict(total=0,unknown=0,n=100000,alpha=.01,comparisons=396,setting_tv=0,leakage=0)
        self.assertGreater(inf.score_interval(**base,start=100000)['sampling_radius'],
                           inf.score_interval(**base,start=0)['sampling_radius'])

    def test_claim_gate_refuses_missing_physical_premises(self):
        ledger=json.loads((STUDY/'sufficiency.json').read_text())
        status=audit.claim_status(ledger)
        self.assertFalse(status['causal_claim_enabled'])
        self.assertIn('history_conditional_joint_setting_tv',status['causal_blockers'])
        self.assertIn('spacelike_endpoint_calibration',status['spacelike_blockers'])

    def test_bias_leakage_and_multiplicity_widen_bounds(self):
        base=dict(total=2,unknown=1,n=1_000_000,alpha=.01,comparisons=4,setting_tv=0,leakage=0)
        a=inf.score_interval(**base)
        b=inf.score_interval(**(base|dict(comparisons=396,setting_tv=.001,leakage=.002)))
        self.assertLess(b['lower'],a['lower']);self.assertGreater(b['upper'],a['upper'])
        self.assertEqual(b['setting_bias_allowance'],.004)
        for update in ({'n':0},{'unknown':-1},{'setting_tv':float('nan')},{'alpha':1},{'total':float('inf')}):
            with self.assertRaises(ValueError):inf.score_interval(**(base|update))

    def test_alternating_effects_do_not_bound_mean_absolute_effect(self):
        # Each block has all four setting pairs, B=X then B=1-X.
        a=np.tile([1,1,2,2],2);b=np.tile([1,2,1,2],2)
        outcomes=np.array([0,0,1,1,1,1,0,0])
        t=audit.counts(a,b,np.zeros(8,dtype=int),outcomes)
        full=audit.contrast(t,t,8,'A_to_B',0)
        self.assertEqual(full['score_sum'],0)
        effects=[]
        for sl in (slice(0,4),slice(4,8)):
            u=audit.counts(a[sl],b[sl],np.zeros(4,dtype=int),outcomes[sl])
            effects.append(audit.contrast(u,u,4,'A_to_B',0)['score_sum']/4)
        self.assertEqual(effects,[1,-1])


class SnapshotTests(unittest.TestCase):
    def test_reconciled_run_and_conservative_unknowns(self):
        b=json.loads((STUDY/'results/counts.json').read_text())
        self.assertEqual(b['n'],107109596)
        result=audit.summarize(b)
        self.assertEqual(len(result['rows']),396)
        self.assertFalse(result['causal_claim_enabled']);self.assertFalse(result['spacelike_claim_enabled'])
        self.assertTrue(all(r['unknown_rows']>0 for r in result['rows'] if r['block']=='all'))
        self.assertEqual(result,json.loads((STUDY/'results/analysis.json').read_text()))



if __name__=='__main__':unittest.main()
