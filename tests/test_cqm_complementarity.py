"""Leakage, independent mark oracle and normalized generative-law checks."""
from pathlib import Path
import sys
import unittest

import numpy as np

STUDY=Path(__file__).resolve().parents[1]/'studies/correlated-qubit-mechanisms'
sys.path.insert(0,str(STUDY))
from cqm_complementarity import (FIRST, FEATURE_NAMES, local_marks, run_features,
                                 environment_mixture, joint_component, smooth_tables)
from cqm_interventions import counts, normalized_flags


class ComplementarityChecks(unittest.TestCase):
    def test_no_future_or_other_acquisition_in_history(self):
        states=[np.zeros(FIRST+100,dtype=np.int8) for _ in range(2)]
        states[0][FIRST-20:FIRST-10]=1
        states[1][FIRST-100:FIRST-40]=2
        before, marks, cell=run_features(states,.3,starts=[FIRST])
        for s in states:
            s[FIRST+1:]=1
        after, changed, other=run_features(states,.3,starts=[FIRST])
        np.testing.assert_array_equal(before,after)
        np.testing.assert_array_equal(cell,other)
        self.assertTrue((marks==0).all())
        self.assertTrue((changed>0).all())
        self.assertEqual(before.shape,(1,len(FEATURE_NAMES)))
        # No history before the declared 0.75 s window can influence a feature.
        for state in states:
            state[FIRST+10]=0
        states[0][:10]=2
        same,_,_=run_features(states,.3,starts=[FIRST+10])
        states[0][:10]=0
        check,_,_=run_features(states,.3,starts=[FIRST+10])
        self.assertEqual(len(same),1)
        np.testing.assert_array_equal(same,check)
        with self.assertRaises(ValueError):
            run_features(states,0,starts=[FIRST-1])

    def test_marks_against_direct_oracle(self):
        rng=np.random.default_rng(20)
        states=rng.integers(0,3,size=2000,dtype=np.int8)
        starts=np.arange(0,1900,33)
        expected=[]
        for start in starts:
            future=states[start+1:start+34]
            candidates=[i for i in range(31) if all(future[i:i+3]!=0)]
            expected.append(0 if not candidates else 1+2*(candidates[0]//11)+int(sum(future!=0)>=12))
        np.testing.assert_array_equal(local_marks(states,starts),expected)
        explicit=np.zeros(100,dtype=np.int8)
        explicit[31:34]=[1,2,1]
        self.assertEqual(local_marks(explicit,np.array([0]))[0],5)

    def test_generative_models_normalize_and_improve_training_fit(self):
        low=np.array([.94,.01,.01,.01,.01,.01,.01])
        high=np.array([.4,.1,.1,.1,.1,.1,.1])
        shared=.7*np.outer(low,low)+.3*np.outer(high,high)
        independent=np.outer(low,high)
        tables=smooth_tables(np.array([shared,independent])*100000)
        total=tables.sum((1,2))
        product=(tables.sum(2)/total[:,None])[:,:,None]*(tables.sum(1)/total[:,None])[:,None,:]
        score=lambda p: float((tables*np.log(p)).sum())
        for model in (environment_mixture,joint_component):
            predicted=model(tables)
            np.testing.assert_allclose(predicted.sum((1,2)),1,atol=1e-12)
            self.assertTrue(np.isfinite(predicted).all())
            self.assertTrue((predicted>0).all())
            self.assertGreaterEqual(score(predicted),score(product)-1e-4)
        self.assertLess(np.max(np.abs(environment_mixture(tables)[0]-tables[0]/total[0])),.003)

    def test_intervention_windows_never_cross_acquisitions(self):
        flags=np.zeros((2,2,4),dtype=bool)
        flags[:,1,0]=True  # Must not count as the future of the preceding acquisition.
        result=counts(flags)
        self.assertEqual(result['total']['paired_flags'],1)
        self.assertEqual(result['total']['next_1ms_paired_flags'],0)
        r=np.zeros((2,1,32768))
        r[:,0,1]=2
        params=[{'ground_mean':-1,'excited_mean':1}]*2
        a=normalized_flags(r,params,np.ones(2))
        b=normalized_flags(r/3,params,np.full(2,3))
        np.testing.assert_array_equal(a,b)
        self.assertEqual(a.sum(),2)


if __name__=='__main__':
    unittest.main()
