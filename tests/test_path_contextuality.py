"""Independent scientific/analysis regressions for the prospective finite test."""
import importlib.util
from pathlib import Path
from fractions import Fraction as F
import unittest
import numpy as np

P=Path(__file__).resolve().parents[1]/'studies/path-contextuality/design.py'
spec=importlib.util.spec_from_file_location('path_contextuality_design',P)
d=importlib.util.module_from_spec(spec); spec.loader.exec_module(d)

class ContextualityTests(unittest.TestCase):
    def test_exact_reference_instrument(self):
        m=d.rational_instrument()
        self.assertEqual(m['table'],[[F(1369,15625),F(7056,15625)],
                                     [F(144,15625),F(7056,15625)]])
        self.assertEqual(m['f'],F(49,625)); self.assertEqual(m['d'],F(1,50))
        self.assertEqual(m['gap'],F(297,15625))
        self.assertLess(d.matrix_check(m),1e-14)

    def test_exact_channel_and_born_check_on_distinct_states(self):
        for t,u,v in [(F(1,5),F(1,4),F(2,5)),(F(2,5),F(1,2),F(1,3)),(F(0),F(1),F(0))]:
            d.matrix_check(d.rational_instrument(t,u,v))

    def test_loss_full_table_and_marginals(self):
        m=d.rational_instrument()
        for eta in [0.,.7,1.]:
            for beta in [0.,.8,1.]:
                for flip in [0.,.1,1.]:
                    p,(a,f,q,dist)=d.lossy_table(m,beta,eta,flip)
                    self.assertAlmostEqual(p.sum(),1)
                    self.assertAlmostEqual(p[:2,:].sum(),beta)
                    self.assertAlmostEqual(p[:,:2].sum(),eta)
                    self.assertEqual(a,p[0,0])
        self.assertLessEqual(d.witness(*d.lossy_table(m,0,1)[1]),0)
        self.assertLessEqual(d.witness(*d.lossy_table(m,1,0)[1]),0)

    def test_cp_empty_and_invalid_counts(self):
        self.assertEqual(d.cp(0,0,.01),(0,1))
        self.assertFalse(d.decision([(0,0)]*4)['reject'])
        for k,n in [(-1,5),(6,5),(1.5,5),(1,-1)]:
            with self.assertRaises(ValueError): d.cp(k,n,.01)

    def test_interval_bound_handles_reversed_slope(self):
        # q < d: the worst case is the LOWER endpoint of f.
        result=d.interval_decision([(0.18,.2),(.1,.9),(.1,.2),(.7,.8)])
        self.assertFalse(result['reject'])
        self.assertAlmostEqual(result['lower_gap'],-.02)  # trivial cap q=.2 dominates

    def test_real_positive_and_nonviolating_cases(self):
        for eta,beta,flip,expected in [(1,1,0,True),(.5,.8,.02,False)]:
            _,params=d.lossy_table(d.rational_instrument(),beta,eta,flip)
            counts=[(round(p*10_000_000),10_000_000) for p in params]
            self.assertEqual(d.decision(counts)['reject'],expected)

    def test_budget_certificate_matches_implemented_rule(self):
        m=d.rational_instrument(); gap=float(m['gap']); b=d.certified_budget(gap)
        n=b['per_main_context']
        self.assertEqual(b['eligible_trials'],4*n+44*b['per_audit_context'])
        ra=np.sqrt(np.log(8/.005)/(2*n)); rb=np.sqrt(np.log(8/.1)/(2*n))
        self.assertLess(4*(ra+rb),gap)
        _,p=d.lossy_table(m)
        # Worst-direction estimates at the declared power-event boundary.
        counts=[(round(max(0,min(1,x+s*rb))*n),n) for x,s in zip(p,[-1,1,1,1])]
        self.assertTrue(d.hoeffding_decision(counts)['reject'])

    def test_calibration_includes_imaginary_direction_and_full_channel(self):
        report=d.calibration_tables(d.rational_instrument())
        self.assertEqual(len(report['contexts']),44)
        self.assertEqual(abs(report['preparation_coordinate_determinant']),2)
        self.assertLess(report['max_equivalence_residual'],1e-14)
        self.assertIn('channel/Y+/Y/M',report['contexts'])

    def test_random_allocation_covers_audit_when_main_signal_is_large(self):
        random=d.random_context_budget(.9)
        fixed=d.certified_budget(.9,beta=.05)
        self.assertGreaterEqual(random['minimum_count_event'],fixed['per_audit_context'])
        self.assertGreaterEqual(random['minimum_count_event'],fixed['per_main_context'])
        for contexts in [4,47,48.5]:
            with self.assertRaises(ValueError): d.random_context_budget(.9,contexts=contexts)

    def test_countermodels(self):
        r=d.exact_countermodels()
        expected=[[str(x) for x in row] for row in d.rational_instrument()['table']]
        for name in ['drop_cap','drop_disturbance']:
            self.assertEqual(r[name+'_quantum_joint'],expected)
            self.assertEqual(r[name+'_bypass'],'49/625')
        self.assertGreater(F(r['drop_cap_negative_response']),F(16,25))
        self.assertLessEqual(F(r['drop_disturbance_max_negative_response']),F(16,25))
        self.assertEqual(r['reference_parameter_null_bypass'],'49/625')
        self.assertEqual(r['reference_parameter_null_joint'],'49/1250')
        self.assertEqual(r['marginal_disturbance'],'0')
        self.assertEqual(r['legitimate_null_joint_negative_success'],'1/100')
        self.assertEqual(r['legitimate_null_postselected_negative_fraction'],'1')

    def test_not_all_large_weak_values_violate(self):
        # Orthogonal pre/post: weak-value denominator vanishes; joint rule stays defined.
        m=d.rational_instrument(F(1,3),F(1,3),F(1,2))
        self.assertEqual(m['f'],0)
        self.assertLessEqual(m['gap'],0)

    def test_snapshot_does_not_tolerate_discrete_drift(self):
        with self.assertRaisesRegex(AssertionError,r'\$\.n'):
            d.compare({'n':100},{'n':101})
        with self.assertRaisesRegex(AssertionError,r'\$\.reject'):
            d.compare({'reject':True},{'reject':False})
        d.compare({'p':.2},{'p':.2+1e-14})

if __name__=='__main__': unittest.main()
