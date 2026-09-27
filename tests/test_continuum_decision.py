"""Count-to-exclusion regression tests, including physical calibration boundaries."""
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[1]/'studies'/'continuum-finite-models'
sys.path.insert(0,str(ROOT))
try:
    spec=importlib.util.spec_from_file_location('continuum_decision',ROOT/'decision.py')
    d=importlib.util.module_from_spec(spec); spec.loader.exec_module(d)
    import model as m
finally:
    sys.path.pop(0)


def counts(probabilities,n):
    result=np.floor(np.asarray(probabilities)*n).astype(int)
    # Deterministic largest-remainder rounding, preserving all eligible trials.
    order=np.argsort(-(np.asarray(probabilities)*n-result))
    result[order[:n-int(result.sum())]]+=1
    return result.tolist()


class ContinuumDecision(unittest.TestCase):
    def test_tiny_tail_and_endpoints(self):
        for amplitude in (1e-5,1e-10,1e-100):
            c=np.array([np.sqrt(1-amplitude**2),amplitude])
            projected,tail,error=m.project_normalize(c,[True,False])
            self.assertGreater(error,0)
            self.assertAlmostEqual(error/amplitude,1,places=9)
            np.testing.assert_allclose(error,np.linalg.norm(c-projected),rtol=2e-11,atol=0)
        projected,tail,error=m.project_normalize([1,0],[True,False])
        self.assertEqual(error,0); self.assertEqual(tail,0)
        projected,tail,error=m.project_normalize([0,1],[True,False])
        self.assertIsNone(projected); self.assertEqual(tail,1); self.assertEqual(error,2)

    def test_boundary_calibration_centers_are_not_silently_clipped(self):
        for eta,lo,hi in ((0,0,0),(100,100,0),(1,0,1),(99,100,1)):
            cal={'detection':[eta,100],'phase0':[lo,100],'phase_pi':[hi,100]}
            box=d.calibration_box(cal)
            self.assertFalse(box['empty'])
            self.assertGreaterEqual(box['contrast'][0],0)
            self.assertLessEqual(box['contrast'][1],box['efficiency'][1])
        # Statistically incompatible data return inconclusive, not blanket exclusions.
        cal={'detection':[0,10**7],'phase0':[10**7,10**7],'phase_pi':[0,10**7]}
        result=d.exclude_candidates(cal,[[0,0,10000]]*24,[8,16])
        self.assertEqual(result['status'],'inconclusive-calibration')
        self.assertEqual(result['rejected'],[])

    def test_strict_threshold_equality_and_failure_coordinate(self):
        box=(np.array([.2,.3,.1]),np.array([.4,.5,.3])); r=.05
        upper=box[1]+r
        self.assertFalse(d.outside_enlarged_box(upper,box,r).any())
        shifted=upper.copy(); shifted[2]=np.nextafter(upper[2],np.inf)
        np.testing.assert_array_equal(d.outside_enlarged_box(shifted,box,r),[False,False,True])
        lower=box[0]-r
        self.assertFalse(d.outside_enlarged_box(lower,box,r).any())
        lower[0]=np.nextafter(lower[0],-np.inf)
        self.assertTrue(d.outside_enlarged_box(lower,box,r)[0])

    def test_full_loss_null_and_unphysical_count_inputs(self):
        n=100000
        cal={k:[0,n] for k in ('detection','phase0','phase_pi')}
        result=d.exclude_candidates(cal,[[0,0,n]]*24,[8,16,32,256])
        self.assertEqual(result['rejected'],[])
        cal={'detection':[80000,n],'phase0':[76000,n],'phase_pi':[4000,n]}
        result=d.exclude_candidates(cal,[[0,0,n]]*24,[8,16])
        self.assertEqual(result['rejected'],[8,16])
        self.assertTrue(any(w['outcome']=='failure' for w in result['models'][0]['witnesses']))
        for row in ([1,-1,10],[True,0,10],[.5,1,10],[0,0,0],[1,2]):
            with self.assertRaises(ValueError): d.exclude_candidates(cal,[row]*24,[8])

    def test_simultaneous_candidates_and_identical_rule(self):
        n=1000000
        cal={'detection':[800000,n],'phase0':[760000,n],'phase_pi':[40000,n]}
        records=[counts(m.readout(-q,.8,.9),n) for j,t,q in d.DEFAULT_SETTINGS]
        candidates=[8,16,64,256,1024]
        result=d.exclude_candidates(cal,records,candidates,phase_offset=.005,fractional_scale=.001)
        self.assertIn(8,result['rejected']); self.assertIn(1024,result['retained'])
        for candidate in candidates:
            single=d.exclude_candidates(cal,records,[candidate],phase_offset=.005,fractional_scale=.001)
            self.assertEqual(candidate in result['rejected'],bool(single['rejected']))
        with tempfile.TemporaryDirectory() as td:
            source=Path(td)/'counts.json'; target=Path(td)/'decision.json'
            source.write_text(json.dumps(dict(calibration=cal,main_counts=records,candidates=candidates,
                certificates=dict(phase_offset=.005,fractional_scale=.001))))
            subprocess.run([sys.executable,str(ROOT/'decision.py'),str(source),'--output',str(target)],check=True,cwd=td)
            self.assertEqual(json.loads(target.read_text())['rejected'],result['rejected'])

    def test_advertised_budget_guarantees_actual_decision_on_good_events(self):
        result=json.loads((ROOT/'results/sensitivity.json').read_text())
        for sites in (8,32,96):
            row=next(r for r in result['rows'] if r['sites']==sites and r['efficiency']==.8
                     and r['visibility']==.9 and r['phase_offset']==.005
                     and r['fractional_scale']==.001 and r['calibration_radius']==.001)
            n=row['trials_per_setting']; nc=row['calibration_trials_per_stratum']
            rc=m.radius(nc,3,.025)
            cal={'detection':[round((.8+.99*rc)*nc),nc],
                 'phase0':[round((.76-.99*rc)*nc),nc],
                 'phase_pi':[round((.04+.99*rc)*nc),nc]}
            rb=m.radius(n,72,.1)
            records=[]
            for j,t,q in d.DEFAULT_SETTINGS:
                p=m.readout(-q,.8,.9)
                null=m.readout(m.phase_gap(m.Ring(),sites,0,j,t)-q,.8,.9)
                shift=np.sign(null[0]-p[0])*min(rb/2,p[0]/2,p[1]/2)
                p=p+np.array([shift,-shift,0])
                records.append(counts(p,n))
            observed=d.exclude_candidates(cal,records,[sites],phase_offset=.005,fractional_scale=.001)
            self.assertEqual(observed['rejected'],[sites])
            for (j,t,q),record in zip(d.DEFAULT_SETTINGS,records):
                self.assertLess(np.max(np.abs(np.array(record)/n-m.readout(-q,.8,.9))),rb)
            # All calibration perturbations lie in the claimed simultaneous event.
            for key,true in (('detection',.8),('phase0',.76),('phase_pi',.04)):
                self.assertLess(abs(cal[key][0]/nc-true),rc)


if __name__=='__main__': unittest.main()
