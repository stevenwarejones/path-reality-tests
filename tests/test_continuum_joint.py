"""Shared-scale decisions, interval coverage and advertised power certificates."""
import importlib.util
from pathlib import Path
import sys
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[1]/'studies/continuum-finite-models'
sys.path.insert(0,str(ROOT))
try:
    import joint as jnt
    import decision as d
    import model as m
finally:
    sys.path.pop(0)


def counts(p,n):
    expected=np.array(p)*n; row=np.floor(expected).astype(int)
    row[np.argsort(-(expected-row))[:n-int(row.sum())]]+=1
    return row.tolist()


class SharedNuisance(unittest.TestCase):
    def test_shared_scale_recovers_information_lost_in_boxes(self):
        n=10**8; cal={'detection':[8*n//10,n],'phase0':[76*n//100,n],'phase_pi':[4*n//100,n]}
        records=[counts(m.readout(-q,.8,.9),n) for _,_,q in d.DEFAULT_SETTINGS]
        base=d.exclude_candidates(cal,records,[128,192,256],fractional_scale=.001)
        shared=jnt.exclude_shared_candidates(cal,records,[128,192,256],fractional_scale=.001)
        self.assertNotIn(192,base['rejected']); self.assertNotIn(256,base['rejected'])
        self.assertEqual(shared['rejected'],[128,192,256])

    def test_null_coverage_endpoints_loss_and_finite_search(self):
        ring=m.Ring(); n=10**7
        for scale in (-.001,0,.001):
            for phase in (-.005,.005):
                cal={'detection':[8*n//10,n],'phase0':[76*n//100,n],'phase_pi':[4*n//100,n]}
                records=[counts(m.readout(m.phase_gap(ring,128,0,j,t)-q+phase+
                    scale*t*float(ring.lattice(j,128))/ring.hbar,.8,.9),n)
                    for j,t,q in d.DEFAULT_SETTINGS]
                result=jnt.exclude_shared_candidates(cal,records,[128],phase_offset=.005,
                                                     fractional_scale=.001,max_cells=128)
                self.assertEqual(result['rejected'],[])
        cal={k:[0,n] for k in ('detection','phase0','phase_pi')}
        result=jnt.exclude_shared_candidates(cal,[[0,0,n]]*24,[128],fractional_scale=.001,max_cells=1)
        self.assertEqual(result['rejected'],[])
        self.assertFalse(jnt.certify_shared_gap(128,.01,max_cells=1)['certified'])

    def test_single_momentum_scale_overlap_is_exact(self):
        ring=m.Ring(); sites=128; mode=3
        ec=float(ring.continuum(mode)); ea=float(ring.lattice(mode,sites)); delta=ea-ec
        sc=delta/(2*ec); sa=-delta/(2*ea)
        self.assertLess(abs(sc),.001); self.assertLess(abs(sa),.001)
        for t in (.5,2,8,32):
            for q in (0,.3,np.pi/2):
                np.testing.assert_allclose(m.readout(sc*t*ec-q,.8,.9),
                    m.readout(t*delta+sa*t*ea-q,.8,.9),atol=2e-15,rtol=0)
        # The SAME fitted scales cannot match a second momentum's curvature.
        other=float(ring.lattice(1,sites))*(1+sa)-float(ring.continuum(1))*(1+sc)
        self.assertGreater(abs(other),1e-4)

    def test_partition_power_certificate_drives_actual_counts(self):
        ring=m.Ring()
        for sites,target in ((128,.01),(192,.001),(256,.001)):
            cert=jnt.certify_shared_gap(sites,target,max_cells=4096)
            self.assertTrue(cert['certified'])
            n=m.required_trials(target,72); nc,_=m.calibration_trials(.001)
            rc=m.radius(nc,3,.025); rb=m.radius(n,72,.1)
            cal={'detection':[round((.8+.99*rc)*nc),nc],
                 'phase0':[round((.76-.99*rc)*nc),nc],
                 'phase_pi':[round((.04+.99*rc)*nc),nc]}
            for scale in (-.001,0,.001):
                records=[]
                for mode,t,q in d.DEFAULT_SETTINGS:
                    p=m.readout(scale*t*float(ring.continuum(mode))-q,.8,.9)
                    perturb=min(.9*rb,p[0]/2,p[1]/2)
                    records.append(counts(p+np.array([perturb,-perturb,0]),n))
                result=jnt.exclude_shared_candidates(cal,records,[sites],fractional_scale=.001,
                    partitions={sites:cert['decision_partition']},max_cells=1)
                self.assertEqual(result['rejected'],[sites])
                for (mode,t,q),row in zip(d.DEFAULT_SETTINGS,records):
                    p=m.readout(scale*t*float(ring.continuum(mode))-q,.8,.9)
                    self.assertLess(max(abs(np.array(row)/n-p)),rb)
        with self.assertRaises(ValueError):
            jnt.exclude_shared_candidates({k:[0,1000] for k in ('detection','phase0','phase_pi')},
                [[0,0,1000]]*24,[256],fractional_scale=.001,partitions={256:[[0.],[-.0005,.001]]})

    def test_calibration_phase_uncertainty_is_propagated(self):
        n=10**8; eta=.8; visibility=.9; phase=.2
        cal={'detection':[round(eta*n),n],
             'phase0':[round(m.readout(phase,eta,visibility)[0]*n),n],
             'phase_pi':[round(m.readout(np.pi+phase,eta,visibility)[0]*n),n]}
        records=[counts(m.readout(-q,eta,visibility),n) for _,_,q in d.DEFAULT_SETTINGS]
        result=d.exclude_candidates(cal,records,[128],calibration_phase_radius=phase)
        self.assertLessEqual(result['calibration']['contrast'][0],eta*visibility)
        self.assertGreaterEqual(result['calibration']['contrast'][1],eta*visibility)


if __name__=='__main__': unittest.main()
