"""Failure checks and independent identities for the full-class gain certificate."""
import copy
from fractions import Fraction as F
import importlib.util
import json
import math
from pathlib import Path
import sys
import unittest
import numpy as np

STUDY=Path(__file__).resolve().parents[1]/'studies/gst-rb-combination-gain'
sys.path.insert(0,str(STUDY))
# A unique module name avoids the repository's other generic verify modules.
spec=importlib.util.spec_from_file_location('gain_verify',STUDY/'verify.py')
v=importlib.util.module_from_spec(spec);spec.loader.exec_module(v)
from gain_core import BASE,GRID,region,exp_lower,radius_numerator,GeneralEvaluator
from interval_physics import certify
from sqd_sources import parse_word
from sqd_general import density_probability

class CombinationGain(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rows=json.loads((BASE/'results/observations.json').read_text())
        cls.R=region(cls.rows)
        cls.models=json.loads((STUDY/'results/witnesses.json').read_text())

    def test_all_retained_counts_and_nonempty_joint(self):
        for model in self.models.values():
            result=v.exact_audit(self.rows,model,self.R)
            expected={'GST':4663,'RB':2010}
            self.assertEqual(result['checked_two_sided_constraints'],expected)
            if model['retained']!='joint':self.assertGreater(result['separation_lower_exact']['decimal'],.004)

    def test_joint_point_is_not_a_separation_witness(self):
        model=copy.deepcopy(self.models['joint']);model['retained']='GST'
        with self.assertRaisesRegex(ValueError,'separation'):v.exact_audit(self.rows,model,self.R)

    def test_witness_fails_when_deleted_source_restored(self):
        for name in ['GST_only','RB_only']:
            model=copy.deepcopy(self.models[name]);model['retained']='joint'
            with self.assertRaisesRegex(ValueError,'retained constraints'):v.exact_audit(self.rows,model,self.R)

    def test_bad_probability_artifact_fails(self):
        model=copy.deepcopy(self.models['GST_only']);model['probabilities'][0]=1.0
        with self.assertRaisesRegex(ValueError,'retained constraints'):v.exact_audit(self.rows,model,self.R)

    def test_physicality_rejects_invalid_spam_and_singular_marginal(self):
        x=np.array(self.models['joint']['x']);x[49]=-1e-5
        with self.assertRaisesRegex(ValueError,'SPAM'):certify(x,8198)
        x=np.array(self.models['joint']['x']);x[:16]=0
        with self.assertRaises(ValueError):certify(x,8198)

    def test_exact_binomial_tail_endpoints(self):
        # Distinct count pairs across both denominators, including boundaries.
        checked=set()
        for i,r in enumerate(self.rows):
            key=(r['k'],r['n'])
            if key in checked:continue
            checked.add(key)
            for k,a in [(r['k'],self.R['hi_num'][i]),(r['n']-r['k'],GRID-self.R['lo_num'][i])]:
                if k==r['n']:self.assertEqual(a,GRID);continue
                n=r['n'];tail=sum(math.comb(n,j)*a**j*(GRID-a)**(n-j) for j in range(k+1))
                self.assertLessEqual(80*6661*tail,GRID**n)

    def test_group_radii_and_coverage_budget(self):
        self.assertEqual(F(2*6661,80*6661)+F(12,480),F(1,20))
        for g in self.R['groups']:
            a=g['radius_numerator'];n=g['n']
            self.assertGreaterEqual(exp_lower(F(2*n*a*a,GRID*GRID)),960)
        self.assertEqual(len(self.R['groups']),12)

    def test_bound_is_exact_sum_of_two_constraints(self):
        g,r=self.R['groups'][0],self.R['groups'][6]
        expected=(1+F(r['upper_numerator'],r['upper_denominator'])-F(g['lower_numerator'],g['lower_denominator']))/2
        self.assertEqual(expected,self.R['joint_upper_q_exact'])
        # Full source budgets are fixed even in a source-only audit.
        self.assertEqual(sum(x['family']=='GST' for x in self.rows),4657)
        with self.assertRaises(ValueError):region(self.rows[:4657])

    def test_density_and_prefix_word_order(self):
        expressions=['{}','GxGy','GyGx','(GxGyGi)^7Gx','Gy(GxGi)^4Gy','(Gy)^4']
        rows=[{'expression':s,'word':parse_word(s)} for s in expressions]
        for model in self.models.values():
            x=model['x'];p=GeneralEvaluator(rows)(x);ind=v.independent_probabilities(rows,x)
            slow=np.array([density_probability(r['word'],x) for r in rows])
            np.testing.assert_allclose(p,ind,atol=1e-12,rtol=0)
            np.testing.assert_allclose(p,slow,atol=1e-12,rtol=0)
        x=self.models['joint']['x']
        # Noncommuting deposited gate conventions must not be silently reversed.
        self.assertGreater(abs(density_probability(parse_word('GxGyGiGx'),x)-density_probability(tuple(reversed(parse_word('GxGyGiGx'))),x)),1e-8)

    def test_unrestricted_dependence_invalidates_coverage(self):
        # X_1=...=X_N=B with B~Bernoulli(1/2): empirical mean is always 0 or 1.
        n=10000;r=F(radius_numerator(n),GRID)
        self.assertLess(r,F(1,2))
        failure_probability=sum(F(1,2) for mean in [F(0),F(1)] if abs(mean-F(1,2))>r)
        self.assertEqual(failure_probability,1)

if __name__=='__main__':unittest.main()
