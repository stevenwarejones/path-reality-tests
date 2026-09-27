"""Independent math, physical propagation, and assumption failures for resource gain."""
import json,sys,unittest
from pathlib import Path
from fractions import Fraction as F
import numpy as np
D=Path(__file__).resolve().parents[1]/'studies/unmeasured-dynamics-gain'
sys.path.insert(0,str(D))
import contrast_certificate as cc
from dynamics_common import BASE,operations,baseline_models,region
from dynamics_verify import audit
from sqd_models import PAULI

class DiscriminationContrast(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rows=json.loads((BASE/'results/observations.json').read_text())
        cls.models={s:json.loads((D/f'results/contrast-{s}.json').read_text()) for s in ['GST','RB']}
        cls.models['joint']=baseline_models()['joint'];cls.R=region(cls.rows)

    def test_exact_offline_certificate_and_two_physical_crossings(self):
        out=cc.run();self.assertEqual(out,json.loads((D/'results/contrast-certificate.json').read_text()))
        for s in ['GST','RB']:
            q=out['models'][s]['separation_below_joint_lower']
            self.assertGreater(F(q['numerator'],q['denominator']),F(1,20))
            self.assertEqual(out['models'][s]['membership']['violations'][s],0)

    def test_independent_kraus_effect_and_extremal_encoding(self):
        for m in self.models.values():
            x=m['x'];_,_,_,channels=operations(x);l,u,t=x[49:52]
            effect=(l+u)*PAULI[0]/2+(l-u)*(np.sin(t)*PAULI[1]+np.cos(t)*PAULI[3])/2
            bare=effect.copy()
            for g in reversed(cc.TARGET):effect=sum(K.conj().T@effect@K for K in channels[g])
            vals,V=np.linalg.eigh(effect);interval=cc.effective_contrast(x)
            self.assertLessEqual(interval['contrast_lower']['decimal'],vals[1]-vals[0])
            self.assertGreaterEqual(interval['contrast_upper']['decimal'],vals[1]-vals[0])
            probs=[]
            for j in [0,1]:
                rho=np.outer(V[:,j],V[:,j].conj())
                for g in cc.TARGET:rho=sum(K@rho@K.conj().T for K in channels[g])
                probs.append(np.trace(bare@rho).real)
            self.assertAlmostEqual(probs[1]-probs[0],vals[1]-vals[0],places=12)
            self.assertAlmostEqual((probs[1]+1-probs[0])/2,(1+vals[1]-vals[0])/2,places=12)

    def test_spectral_bound_in_larger_dimension(self):
        rng=np.random.default_rng(409)
        for d in [2,3,5]:
            for _ in range(20):
                U,_=np.linalg.qr(rng.normal(size=(d,d))+1j*rng.normal(size=(d,d)))
                eigen=rng.uniform(0,1,d);effect=(U*eigen)@U.conj().T;p=[]
                for _ in range(2):
                    A=rng.normal(size=(d,d))+1j*rng.normal(size=(d,d));rho=A@A.conj().T;rho/=np.trace(rho)
                    p.append(np.trace(effect@rho).real)
                self.assertLessEqual(abs(p[0]-p[1]),np.ptp(eigen)+1e-14)

    def test_context_dependent_detector_invalidates_common_effect_step(self):
        # Two constant detectors have diameter zero yet arbitrarily separated means.
        rho=np.eye(2)/2;high=np.eye(2);low=np.zeros((2,2))
        self.assertEqual(np.ptp(np.linalg.eigvalsh(high)),0)
        self.assertEqual(np.ptp(np.linalg.eigvalsh(low)),0)
        self.assertEqual(np.trace((high-low)@rho),1)

    def test_suffix_order_matters(self):
        x=self.models['joint']['x'];right=cc.effective_contrast(x)['contrast_lower']['decimal']
        wrong=cc.effective_contrast(x,tuple(reversed(cc.TARGET)))['contrast_upper']['decimal']
        self.assertGreater(abs(right-wrong),1e-6)

    def test_both_complete_source_deletions_are_essential_for_witnesses(self):
        for source,m in self.models.items():
            if source=='joint':continue
            with self.assertRaisesRegex(ValueError,'retained constraint'):audit(self.rows,m,self.R,['GST','RB'])
            bad={**m,'probabilities':m['probabilities'].copy()}
            i=next(i for i,r in enumerate(self.rows) if r['family']==source and self.R['hi'][i]<.9)
            bad['probabilities'][i]=1
            with self.assertRaisesRegex(ValueError,'retained constraint'):audit(self.rows,bad,self.R,[source])

    def test_row_hash_and_count_corruption_fail(self):
        for field,value in [('word_sha256','0'*64),('k',1)]:
            rows=[dict(r) for r in self.rows];rows[4699][field]=value
            with self.assertRaisesRegex(ValueError,'input mismatch'):cc.validate_inputs(rows)

    def test_interval_rounding_is_outward(self):
        from mpmath import iv
        lo,hi=cc.enclosing_grid(iv.mpf(1)/3)
        self.assertLessEqual(lo,F(1,3));self.assertGreaterEqual(hi,F(1,3))

if __name__=='__main__':unittest.main()
