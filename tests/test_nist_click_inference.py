"""Independent algebra, completion and nonstationary assignment tests for v2."""
import importlib.util
import itertools
import json
import math
from pathlib import Path
import sys
import unittest
import numpy as np

STUDY=Path(__file__).resolve().parents[1]/'studies/nist-bell-causal-audit'
sys.path.insert(0,str(STUDY))
from click_inference import count_interval, effect_interval, LAMBDAS
from click_audit import possible_context, summarize


class ClickInferenceTests(unittest.TestCase):
    def test_click_score_identity_all_response_tables(self):
        for response in itertools.product((0,1),repeat=4):
            for y in (0,1):
                mean=sum((2*x-1)*response[2*x+y] for x in (0,1))
                self.assertEqual(mean,response[2+y]-response[y])
                # Adding arbitrary no-click trials contributes exactly zero.
                self.assertEqual(sum(4*(2*x-1)*0 for x in [0]*100+[1]),0)

    def test_binary_exponential_process_for_every_history_probability(self):
        # Independent direct conditional expectation: arbitrary history may pick p.
        for p in (0.,1e-4,.03,.4,1.):
            for v in LAMBDAS:
                for lam in (-v,v):
                    e=(1-p)*math.exp(-math.expm1(lam)*p)+p*math.exp(lam-math.expm1(lam)*p)
                    self.assertLessEqual(e,1+1e-14)

    def test_inversion_endpoints_against_exponential_threshold(self):
        n=100000;k=100
        b=count_interval(k,k,n)
        self.assertGreater(b['lower'],0);self.assertLess(b['upper'],n)
        lo=max(v*k-b['log_threshold'] for v in LAMBDAS)
        self.assertTrue(math.isfinite(lo))
        # At either selected boundary, the largest corresponding process equals 1/delta.
        self.assertAlmostEqual(max(v*k-math.expm1(v)*b['lower'] for v in LAMBDAS),b['log_threshold'])
        self.assertAlmostEqual(max(-v*k-math.expm1(-v)*b['upper'] for v in LAMBDAS),b['log_threshold'])

    def test_unknown_completion_envelope(self):
        envelope=count_interval(100,103,100000)
        for k in range(100,104):
            exact=count_interval(k,k,100000)
            self.assertLessEqual(envelope['lower'],exact['lower'])
            self.assertGreaterEqual(envelope['upper'],exact['upper'])
        h=np.zeros((4,4),dtype=int)
        h[1,1]=2;h[2,1]=5;h[3,1]=7;h[0,3]=11;h[1,2]=13
        self.assertEqual(possible_context(h,0,0),20)
        self.assertEqual(possible_context(h,1,0),23)

    def test_history_tv_inversion_with_drift_and_adversarial_correlation(self):
        eps=.03
        q=np.array([.25-eps,.25+eps,.25,.25])
        for probabilities in itertools.product((0.,.001,.8,1.),repeat=4):
            p=np.array(probabilities)
            # Each of four histories can align rare outcomes with assignment bias.
            m=float(np.dot(q,p));target=float(p.mean())
            self.assertLessEqual(m/(4*(.25+eps)),target+1e-14)
            self.assertGreaterEqual(m/(4*(.25-eps)),target-1e-14)
        ideal=effect_interval([(100,101),(105,106)],100000)
        biased=effect_interval([(100,101),(105,106)],100000,per_history_tv=eps,leakage=.001)
        self.assertLess(biased['lower'],ideal['lower']);self.assertGreater(biased['upper'],ideal['upper'])

    def test_sparse_counts_gain_precision_and_zero_counts_are_not_certainty(self):
        rare=count_interval(0,0,1000000)
        self.assertEqual(rare['lower'],0);self.assertGreater(rare['upper'],0)
        dense=count_interval(500000,500000,1000000)
        self.assertLess(rare['upper'],dense['upper']-dense['lower'])
        later=count_interval(0,0,1000000,start=1000)
        self.assertGreater(later['upper'],rare['upper'])
        for args in ((2,1,10),(-1,0,10),(0,11,10),(0,0,0),(0.,0,10)):
            with self.assertRaises(ValueError):count_interval(*args)
        with self.assertRaises(ValueError):effect_interval([(1,1),(1,1)],10,per_history_tv=.25)

    def test_snapshot_and_event_support_conservation(self):
        b=json.loads((STUDY/'results/counts.json').read_text())
        r=summarize(b)
        self.assertEqual(r,json.loads((STUDY/'results/click-analysis.json').read_text()))
        self.assertEqual(len(r['rows']),396)
        e=np.array(b['uncertain_event_contexts']);u=np.array(b['unknown_contexts'])
        self.assertEqual(int(u.sum()),19968)
        self.assertEqual(e.sum(axis=(0,2,3)).tolist(),[877,797])
        self.assertTrue(np.all(e<=u[:,None]))
        self.assertFalse(r['causal_claim_enabled']);self.assertFalse(r['spacelike_claim_enabled'])
        for row in r['rows']:
            self.assertFalse(row['event_supported_interval']['empty'])
            self.assertLessEqual(row['event_supported_interval']['absolute_upper'],row['unrestricted_interval']['absolute_upper'])


if __name__=='__main__':unittest.main()
