import importlib.util
from pathlib import Path
from fractions import Fraction as F
import unittest
import math
p=Path(__file__).resolve().parents[1]/'studies/path-contextuality/design.py'
s=importlib.util.spec_from_file_location('dual_design',p);d=importlib.util.module_from_spec(s);s.loader.exec_module(d)

class DualWitnessTests(unittest.TestCase):
    def test_exact_gaps_and_budget_penalty(self):
        gaps=d.dual_witnesses(F(1369,15625),F(144,15625),F(49,625),F(16,25),F(1,50))
        self.assertEqual(gaps,(F(297,15625),F(109,6250)))
        budget=d.dual_certified_budget(gaps)
        self.assertGreater(budget['eligible_trials'],d.certified_budget(gaps[0])['eligible_trials'])
        n=budget['per_main_context']
        radii=sum(math.sqrt(math.log(10/risk)/(2*n)) for risk in [.005,.1])
        self.assertGreater(float(max(gaps)),4*radii)
        self.assertEqual(budget['eligible_trials'],4*n+44*budget['per_audit_context'])

    def test_positive_only_rejection_and_negative_coefficient(self):
        # d=.1 destroys the negative witness but the positive witness survives.
        vals=[.087616,.009216,.0784,.64,.1]
        result=d.dual_interval_decision([(x,x) for x in vals])
        self.assertLess(result['negative_lower_gap'],0)
        self.assertGreater(result['positive_lower_gap'],0)
        self.assertTrue(result['reject'])
        result=d.dual_interval_decision([(0,0),(0,0),(.1,.8),(.8,.8),(.5,.5)])
        self.assertAlmostEqual(result['positive_lower_gap'],-.24)
        self.assertFalse(result['reject'])

    def test_null_intervals_never_reject(self):
        # Fair identity pointer: a=b=f/2; cap .5, disturbance 0.
        for f in [0,.1,.5,1]:
            vals=[f/2,f/2,f,.5,0]
            boxes=[(max(0,x-.01),min(1,x+.01)) for x in vals]
            self.assertFalse(d.dual_interval_decision(boxes)['reject'])

    def test_loss_complement_and_count_validation(self):
        table,old=d.lossy_table(d.rational_instrument(),.8,.9)
        b=table[1:,0].sum()
        self.assertGreater(b,table[1,0])
        self.assertAlmostEqual(old[0]+b,table[:,0].sum())
        with self.assertRaises(ValueError):d.dual_decision([(6,10),(5,10),(1,10),(1,10),(1,10)])
        self.assertFalse(d.dual_decision([(0,0)]*5)['reject'])
        self.assertFalse(d.dual_decision([(0,0)]*5,method='hoeffding')['reject'])

if __name__=='__main__':unittest.main()
