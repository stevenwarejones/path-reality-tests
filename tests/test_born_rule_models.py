"""Independent physics, calibration and integrity checks for the new study."""
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np
import h5py

STUDY=Path(__file__).resolve().parents[1]/'studies/born-rule-identifiability'
sys.path.insert(0,str(STUDY))
import br_core as b
import br_combinations as q
import br_io
from br_sources import verify
sys.path.remove(str(STUDY))


class BornOpticalTests(unittest.TestCase):
    def setUp(self):
        self.rng=np.random.default_rng(41)
        self.g=np.array([[.18,-.10,.12],[-.10,.36,-.08],[.12,-.08,.16]])

    def test_fixed_leakage_preserves_third_order_identity(self):
        for _ in range(50):
            a=self.rng.normal(size=(3,3))+1j*self.rng.normal(size=(3,3))
            g=a.conj().T@a
            y=b.table(g,background=2.3,leakage=.013*np.exp(.7j))
            self.assertAlmostEqual(y@b.SIGNS,0.,places=11)

    def test_quadratic_detector_creates_ordinary_signal(self):
        g=np.ones((3,3));y=b.table(g,quadratic=.01)
        self.assertAlmostEqual(y@b.SIGNS,.36,places=12)

    def test_psd_completion_sharp_boundary(self):
        u=b.imaginary_limit(self.g)
        for sign in [-1,1]:
            self.assertGreaterEqual(np.linalg.eigvalsh(b.complete_gram(self.g,sign*u)).min(),-1e-12)
            bad=self.g.astype(complex);bad[0,1]+=sign*1.001j*u;bad[1,0]-=sign*1.001j*u
            self.assertLess(np.linalg.eigvalsh(bad).min(),0)
            with self.assertRaises(ValueError):b.complete_gram(self.g,sign*1.001*u)

    def test_phase_formula_independent_quadratic_form(self):
        for _ in range(40):
            u=self.rng.uniform(-1,1)*b.imaginary_limit(self.g);d=self.rng.uniform(-1,1)
            g=b.complete_gram(self.g,u)
            got=b.table(g,d)[7]-b.table(g)[7]
            self.assertAlmostEqual(got,b.phase_shift(self.g,u,d),places=13)

    def test_envelope_contains_domain_and_extrema(self):
        for radius in [0,.01,.3,math.pi/2]:
            lo,hi=b.phase_envelope(self.g,radius);u=b.imaginary_limit(self.g)
            values=[b.phase_shift(self.g,v,d) for v in [-u,u] for d in np.linspace(-radius,radius,20001)]
            self.assertGreaterEqual(min(values),lo-1e-12);self.assertLessEqual(max(values),hi+1e-12)
            self.assertAlmostEqual(min(values),lo,places=8);self.assertAlmostEqual(max(values),hi,places=8)

    def test_phase_cycling_and_error_bound_for_complex_psd(self):
        for _ in range(50):
            a=self.rng.normal(size=(3,3))+1j*self.rng.normal(size=(3,3));g=a.conj().T@a
            d=self.rng.uniform(-3,3)
            self.assertAlmostEqual(b.cancellation_witness(g,d,theta=.1,background=4.),.1,places=11)
            error=self.rng.uniform(-.2,.2)
            self.assertLessEqual(abs(b.cancellation_witness(g,d,error=error)),b.phase_error_bound(g,abs(error))+1e-11)

    def test_analytic_leakage_response_envelope(self):
        g=b.complete_gram(self.g,.5*b.imaginary_limit(self.g));bound=b.leakage_response_bound(g,.014,.002,1.02)
        for _ in range(100):
            scale=self.rng.uniform(.98,1.02);leak=.014*self.rng.uniform(0,1)*np.exp(1j*self.rng.uniform(-math.pi,math.pi))
            w=b.cancellation_witness(g*scale,self.rng.uniform(-math.pi,math.pi),leakage=leak,quadratic=.002)
            self.assertLessEqual(abs(w),bound+1e-12)

    def test_background_and_normalization(self):
        g=self.g;y=b.table(g);s=b.statistics(y)
        shifted=b.statistics(y+3.1)
        for key in ['peres','epsilon_V','denominator_V']:
            self.assertAlmostEqual(float(s[key]),float(shifted[key]),places=12)
        scaled=b.statistics(y*4)
        self.assertAlmostEqual(float(s['epsilon_V']/s['denominator_V']),float(scaled['epsilon_V']/scaled['denominator_V']))
        with self.assertRaises(ValueError):b.statistics(np.ones(8))

    def test_hac_keeps_cycle_gaps(self):
        x=np.array([1.,3.,2.]);ids=np.array([1,3,5])
        self.assertTrue(np.allclose(b.hac_covariance(x,ids,1),b.hac_covariance(x,ids,0)))
        self.assertFalse(np.allclose(b.hac_covariance(x,[1,2,3],1),b.hac_covariance(x,ids,1)))

    def test_budget_tail_certificate_and_impossible_gap(self):
        gap=.01;a=.01;power=.9;n=b.budget(gap,alpha=a,power=power)['per_setting']
        threshold=math.sqrt(3.5*math.log(2/a)/(2*n))
        self.assertLessEqual(math.exp(-2*n*(gap-threshold)**2/3.5),1-power)
        self.assertIsNone(b.budget(0))

    def test_committed_archive_pairs_and_removal(self):
        data=json.loads((STUDY/'results/identifiability.json').read_text())
        for d in data['datasets'].values():
            self.assertLess(d['max_table_residual_V'],1e-10)
            self.assertEqual(d['local_whitened_rank'],1)
            self.assertAlmostEqual(d['added_witness_ordinary_V'],0,places=11)
            self.assertGreater(abs(d['added_witness_alternative_V']),1e-4)
        self.assertTrue(next(r for r in data['joint'] if r['phase_radius_rad']==.03)['contains_zero'])


class BornQubitTests(unittest.TestCase):
    def test_binary_rule_positivity_complements_and_limits(self):
        p=np.linspace(0,1,2001)
        for theta in [-1,-.1,0,.1,1]:
            f=q.probability(p,theta)
            self.assertTrue(np.all((f>=0)&(f<=1)))
            self.assertTrue(np.allclose(f+q.probability(1-p,theta),1))
        self.assertTrue(np.array_equal(q.probability(p,0),p))
        with self.assertRaises(ValueError):q.probability(.5,1.01)

    def test_all_integer_depth_equivalence_and_circuit_matrix(self):
        sx=np.array([[1,-1j],[-1j,1]])/math.sqrt(2)
        zero=np.array([1.,0.]);phi=np.linspace(-math.pi,math.pi,129)
        for theta in [-1,-.1,.01,.9,1]:
            warp=q.warped_angle(phi,theta)
            for n in range(101):
                self.assertTrue(np.allclose(q.probability(q.ideal(phi,n),theta),q.ideal(warp,n),atol=2e-15))
        # Independently evaluate the actual SX Rz(phi+pi) SX^n circuit.
        for f in [.17,.82,2.4]:
            rz=np.diag([np.exp(-1j*(f+math.pi)/2),np.exp(1j*(f+math.pi)/2)])
            for n in [1,2,3,5,9]:
                state=np.linalg.matrix_power(sx,n)@rz@sx@zero
                self.assertAlmostEqual(abs(state[0])**2,float(q.ideal(f,n)),places=12)

    def test_added_orthogonal_preparation_breaks_pair(self):
        theta=.01;x=1/math.sqrt(2)
        radius=2*float(q.h(x,theta))**2
        self.assertAlmostEqual(radius,(1+theta/4)**2,places=13)
        self.assertGreater(radius,1.)
        angle=float(q.warped_angle(math.pi/4,theta))
        self.assertAlmostEqual(math.cos(angle)**2+math.sin(angle)**2,1.)

    def test_mixtures_are_not_silently_density_matrix_equivalent(self):
        # Same z component but different pure-state ensembles can have different
        # polynomial moments; alternate operational state is the ensemble class.
        p1=np.array([.2,.8]);p2=np.array([.5,.5]);weights=np.array([.25,.75])
        mean=float(weights@p1);expected=float(weights@q.probability(p1,.5))
        self.assertNotAlmostEqual(expected,float(q.probability(mean,.5)),places=8)

    def test_count_decision_does_not_exclude_origin_and_requires_calibration(self):
        self.assertFalse(q.quadrature_decision([50,50,50,50],100)[0])
        n=10000000;c=int(round(n*(1+float(q.h(1/math.sqrt(2),.01)))/2))
        self.assertTrue(q.quadrature_decision([c,n-c,c,n-c],n)[0])
        self.assertFalse(q.quadrature_decision([c,n-c,c,n-c],n,radius=1.01)[0])


class BornIntegrityTests(unittest.TestCase):
    def fixture(self,root,bad=None):
        path=root/'fixture.h5'; order=np.array([2,3,7,5,1,6,4,0]);pd=(order+1.).reshape(8,1,1)
        times=np.zeros((8,1,6));times[:,:,:5]=[2020,1,1,0,0];times[:,0,5]=np.arange(8)
        if bad=='duplicate':order[0]=order[1]
        if bad=='time':times[1]=times[0]
        with h5py.File(path,'w') as h:
            h['PD signal']=pd;h['time']=times;h['shuttercombination']=order.reshape(8,1,1)
        data=path.read_bytes();e={'name':path.name,'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest(),'role':'primary','samples_per_setting':1}
        return path,e

    def test_reader_groups_masks_and_rejects_bad_schema(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);path,e=self.fixture(root);y,_=br_io.read(root,e)
            self.assertTrue(np.array_equal(y[0],np.arange(1,9)))
            for bad in ['duplicate','time']:
                path,e=self.fixture(root,bad)
                with self.assertRaises(ValueError):br_io.read(root,e)

    def test_hash_failure_does_not_repair_source(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);path,e=self.fixture(root);data=path.read_bytes();path.write_bytes(data+b'x')
            with self.assertRaises(ValueError):verify(path,e)
            self.assertEqual(path.read_bytes(),data+b'x')

    def test_source_provenance_gates(self):
        d=json.loads((STUDY/'results/portfolio-audit.json').read_text())
        self.assertEqual(d['ibm']['lima1-5s']['evidence_type'],'simulation')
        self.assertFalse(d['ibm']['lagos_bench']['depth_column_present'])
        self.assertGreater(d['photon']['PNR SNSPD']['event_reconciliation'][0]['triggers'],100000)
        self.assertFalse(d['photon']['8-bin TMD']['same_data_calibration_is_independent'])
