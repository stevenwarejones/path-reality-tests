"""Physical, inference and acquisition-failure tests; no network or source data."""
import copy
import importlib.util
import itertools
import json
import pickle
import sys
from pathlib import Path
import tempfile
import unittest

import numpy as np

STUDY = Path(__file__).resolve().parents[1]/'studies/calibrated-quantum-memory'


def load(name):
    spec = importlib.util.spec_from_file_location('memory_'+name, STUDY/(name+'.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


A = load('audit')
C = load('certificate')
S = load('acquisition_audit')
T = load('identifiability')
# Standalone study modules use sibling imports.
sys.modules['identifiability']=T
D = load('delay_candidates')
sys.modules['delay_candidates']=D
J = load('joint_statistics')


class CountAuditTests(unittest.TestCase):
    def rows(self):
        return {k: {'00': 40, '01': 10, '10': 20, '11': 30} for k in A.LABELS}

    def test_regrouping_is_lossless_involution(self):
        rng = np.random.default_rng(1)
        rows = {k: dict(zip(A.OUTCOMES, map(int, rng.multinomial(8000, [.2,.3,.1,.4]))))
                for k in sorted(A.LABELS)}
        archived = A.regroup(rows)
        self.assertEqual(A.regroup(archived), rows)
        self.assertEqual(sum(map(lambda c: sum(c.values()), archived.values())), 8000*324)
        self.assertGreater(len({sum(c.values()) for c in archived.values()}), 1)

    def test_schema_rejects_missing_outcome_and_row(self):
        rows = self.rows()
        del rows[next(iter(rows))]['11']
        with self.assertRaises(ValueError): A.regroup(rows)
        rows = self.rows()
        del rows[next(iter(rows))]
        with self.assertRaises(ValueError): A.regroup(rows)

    def test_schema_rejects_floats_negative_and_boolean_counts(self):
        for bad in (-1, .2, True):
            rows = self.rows()
            rows[next(iter(rows))]['00'] = bad
            with self.assertRaises(ValueError): A.regroup(rows)

    def test_past_diagnostic_does_not_reject_constant_marginal(self):
        result = A.past_bound(self.rows(), .05, 324)
        self.assertTrue(all(r['epsilon_lower'] == 0 for r in result))

    def test_analytic_distance_to_shared_marginal(self):
        rows = self.rows()
        rows['xp,x,xp,x'] = {'00': 900000, '01': 0, '10': 100000, '11': 0}
        rows['xp,x,xm,x'] = {'00': 100000, '01': 0, '10': 900000, '11': 0}
        bound = max(r['epsilon_lower'] for r in A.past_bound(rows,.05,324))
        self.assertGreater(bound, .397)
        self.assertLess(bound, .4)

    def test_hash_mismatch_fails_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory)/'test.json'; p.write_text('{}')
            with self.assertRaises(ValueError):
                A.verify_file(p, {'bytes':2, 'sha256':'0'*64})


class MemoryCertificateTests(unittest.TestCase):
    def test_kraus_complete_positivity_and_normalization(self):
        for lam in (-1/3,0.,.4,.9,1.):
            ks = C.depolarizing_kraus(lam)
            np.testing.assert_allclose(sum(k.conj().T@k for k in ks), C.I, atol=1e-14)
            j = C.choi(ks)
            self.assertGreaterEqual(np.linalg.eigvalsh(j).min(), -1e-14)
            np.testing.assert_allclose(np.trace(j.reshape(2,2,2,2),axis1=1,axis2=3), C.I,atol=1e-14)

    def test_physical_swap_and_bypass_are_equal_on_full_operator_basis(self):
        self.assertLess(C.physical_checks()['maximum_residual'],1e-14)

    def test_reset_probe_distinguishes_exact_pair(self):
        for lam in (.1,.6,1.):
            output=C.channel(C.STATES[0],C.depolarizing_kraus(lam))
            self.assertAlmostEqual(np.trace(C.STATES[0]@output).real,(1+lam)/2)
        self.assertAlmostEqual(C.calibration_probabilities(C.depolarizing_kraus(0))[0,0,0],.5)

    def test_single_outcome_bound_cannot_be_applied_to_hidden_flags(self):
        for rho in C.STATES:
            average = sum(p@rho@p.conj().T/4 for p in [C.I,*C.P])
            corrected = sum(p.conj().T@(p@rho@p.conj().T)@p/4 for p in [C.I,*C.P])
            np.testing.assert_allclose(average,C.I/2,atol=1e-14)
            np.testing.assert_allclose(corrected,rho,atol=1e-14)

    def test_continuum_of_classical_memory_models_obeys_eb_bound(self):
        # Arbitrary measurement direction and different conditional output direction.
        rng=np.random.default_rng(6)
        for _ in range(100):
            v=rng.normal(size=3);v/=np.linalg.norm(v)
            w=rng.normal(size=3);w/=np.linalg.norm(w)
            m=(C.I+sum(v[a]*C.P[a] for a in range(3)))/2
            out=(C.I+sum(w[a]*C.P[a] for a in range(3)))/2
            f=np.mean([np.trace(rho@(np.trace(m@rho)*out+
                       np.trace((C.I-m)@rho)*(C.I-out))).real for rho in C.STATES])
            self.assertAlmostEqual(f,.5+np.dot(v,w)/6)
            self.assertLessEqual(f,2/3+1e-14)

    def test_complete_source_removal_witnesses_and_joint_feasibility(self):
        data=json.loads((STUDY/'results/synthetic-certificate.json').read_text())
        for key,value in data['source_removal_witnesses'].items():
            if type(value) is bool: self.assertTrue(value,key)
        result=C.infer(data['calibration_successes'],data['dynamics_successes'])
        self.assertTrue(result['reject'])
        self.assertGreater(result['margin'],.04)
        self.assertAlmostEqual(result['margin'],data['certificate']['margin'])

    def test_unsupported_tables_fail_closed(self):
        for cal,dyn in [(np.zeros((3,3,1),dtype=int),np.zeros(6,dtype=int)),
                        (np.zeros((3,3,2)),np.zeros(6,dtype=int)),
                        (np.full((3,3,2),8001),np.zeros(6,dtype=int))]:
            with self.assertRaises(ValueError): C.infer(cal,dyn)

    def test_affine_inconsistent_calibration_is_not_memory_rejection(self):
        cal=np.full((3,3,2),4000,dtype=int)
        cal[0,0,:]=8000; cal[0,1,:]=0
        result=C.infer(cal,np.full(6,8000,dtype=int))
        self.assertEqual(result['status'],'calibration_affine_inconsistent')
        self.assertFalse(result['reject'])

    def test_measured_record_never_claims_memory_certificate(self):
        data=json.loads((STUDY/'results/nmn-audit.json').read_text())
        self.assertFalse(data['quantum_memory_certified'])
        self.assertEqual(data['simultaneous_rows'],3240)
        uq=[r for r in data['runs'] if r['source']=='NMN_lab_rslts.json'][0]
        self.assertAlmostEqual(uq['strongest_past_marginal_contrast']['epsilon_lower'],
                               .057844137838283474)




class AcquisitionSemanticsTests(unittest.TestCase):
    def test_safe_primitive_decoder_and_forbidden_object(self):
        values = [0., .25, .75, 1.]*600
        self.assertEqual(S.primitive_float_list(pickle.dumps(values, protocol=4)), values)
        for raw in [pickle.dumps(np.array(values), protocol=4),
                    pickle.dumps({'probabilities': values}, protocol=4),
                    pickle.dumps([[.5,.5]], protocol=4),
                    pickle.dumps([.5,.5], protocol=4)+b'junk',
                    pickle.dumps([float('nan'),.5], protocol=4)]:
            with self.assertRaises(ValueError): S.primitive_float_list(raw)

    def test_mapping_denominators_do_not_resolve_branch_labels(self):
        rows = {k:dict(zip(A.OUTCOMES,[17+i%4,13-i%4,7,3]))
                for i,k in enumerate(sorted(A.LABELS))}
        archived = A.regroup(rows)
        a = S.alternative_regroup(archived,0,'1')
        b = S.alternative_regroup(archived,0,'0')
        self.assertEqual(a, rows)
        self.assertEqual({sum(v.values()) for v in a.values()}, {40})
        self.assertEqual({sum(v.values()) for v in b.values()}, {40})
        for k in a:
            parts=k.split(',');parts[2]=A.opposite(parts[2])
            self.assertEqual(b[k],a[','.join(parts)])

    def test_outcome_axis_cannot_be_silently_swapped(self):
        labels = {'qubit':['q0'],'initial_state':['0','1','+','i+'],'axis':['x','y','z']}
        states = np.zeros((1,5,4,3),dtype=int)
        self.assertEqual(S.qpt_counts(states,labels,5), [[0,0,0]]*4)
        for bad in [dict(labels,axis=['z','y','x']),dict(labels,qubit=['q2'])]:
            with self.assertRaises(ValueError): S.qpt_counts(states,bad,5)
        with self.assertRaises(ValueError): S.qpt_counts(states[:,:4],labels,5)
        states[0,0,0,0]=2
        with self.assertRaises(ValueError): S.qpt_counts(states,labels,5)

    def test_calibration_transfer_fails_on_changed_operation(self):
        cal={'threshold':.001,'length':2000,'integration_weights_angle':.1}
        S.assert_readout_match(cal,dict(cal))
        for key,value in [('threshold',.002),('length',1999),('integration_weights_angle',.2)]:
            with self.assertRaises(ValueError): S.assert_readout_match(cal,dict(cal,**{key:value}))

    def test_measured_snapshot_does_not_promote_readout_to_instrument(self):
        result=json.loads((STUDY/'results/acquisition-audit.json').read_text())
        self.assertFalse(result['quantum_memory_certified'])
        qm=result['quantum_memory_data']
        self.assertFalse(qm['finite_data_instrument_region_built'])
        self.assertEqual(qm['readout_calibration'][0]['counts'],[[4005,995],[744,4256]])
        self.assertEqual(qm['terminal_dynamics'][0]['parents'],[6047])
        self.assertTrue(all(not row['intermediate_flag_axis_present'] for row in qm['terminal_dynamics']))


class IBMIdentifiabilityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.witnesses=json.loads((STUDY/'results/ibm-instrument-witnesses.json').read_text())['witnesses']

    def test_exact_physicality_sharing_and_fixed_instrument_endpoint(self):
        from fractions import Fraction
        for witness in self.witnesses:
            endpoint=T.check_instruments(witness)
            self.assertGreater(endpoint,Fraction(3,1000))
            self.assertLess(endpoint,Fraction(301,100000))
            bad=copy.deepcopy(witness)
            beyond=endpoint+Fraction(1,10**10)
            bad['memory_weight']=[beyond.numerator,beyond.denominator]
            with self.assertRaises(ValueError):T.check_instruments(bad)

    def test_certificate_mutations_fail_closed(self):
        mutations=[lambda w:w['choi'][0].pop(),
                   lambda w:w['choi'].pop(),
                   lambda w:w['choi'][0][0]['real'].__setitem__(1,123),
                   lambda w:w['choi'][0][0]['real'].__setitem__(0,-1),
                   lambda w:w['choi'][0][0]['real'].__setitem__(0,1.5),
                   lambda w:w['settings'].pop()]
        for mutation in mutations:
            bad=copy.deepcopy(self.witnesses[0]);mutation(bad)
            with self.assertRaises(ValueError):T.check_instruments(bad)

    def test_directed_binomial_membership_interior_endpoints_and_rejection(self):
        for k,p in [(50,50),(0,1),(100,99)]:
            self.assertGreaterEqual(T.binomial_membership(100,k,p,100,11664),1)
        for k,p in [(50,1),(0,99),(100,1)]:
            with self.assertRaises(ValueError):T.binomial_membership(100,k,p,100,11664)

    def test_independent_kraus_swap_circuit_on_full_operator_basis(self):
        swap=np.array([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]])
        basis=[np.eye(4)[j].reshape(2,2) for j in range(4)]
        def apply_choi(j,x):
            return np.einsum('ij,iajb->ab',x,j.reshape(2,2,2,2))
        for witness in self.witnesses:
            w=float(T.F(*witness['memory_weight']))
            for pair in witness['choi']:
                for block in pair:
                    re,im=T.arrays(block)
                    j=(np.asarray(re,float)+1j*np.asarray(im,float))/witness['denominator']
                    q=np.trace(j).real/2
                    k=(j-w*q*T.IDENTITY_CHOI)/(1-w)
                    vals,vecs=np.linalg.eigh(k)
                    self.assertGreater(vals.min(),0)
                    kraus=[np.sqrt(vals[i])*vecs[:,i].reshape(2,2).T for i in range(4)]
                    for x in basis:
                        joint=swap@np.kron(x,T.I/2)@swap
                        evolved=sum(np.kron(a,T.I)@joint@np.kron(a.conj().T,T.I) for a in kraus)
                        returned=swap@evolved@swap
                        output=np.trace(returned.reshape(2,2,2,2),axis1=1,axis2=3)
                        direct=sum(a@x@a.conj().T for a in kraus)
                        np.testing.assert_allclose((1-w)*direct+w*output,apply_choi(j,x),atol=3e-14)

    def test_quantum_comb_normalization_causality_and_temporal_entanglement(self):
        w=.003
        phi=np.outer([1,0,0,1],[1,0,0,1])
        direct=np.kron(phi,phi)
        bypass=np.zeros((16,16))
        for i,bits in enumerate(itertools.product(range(2),repeat=4)):
            for j,other in enumerate(itertools.product(range(2),repeat=4)):
                a,b,o,c=bits;aa,bb,oo,cc=other
                bypass[i,j]=phi[2*a+c,2*aa+cc]*(b==bb)*(o==oo)/2
        comb=(1-w)*direct+w*bypass
        self.assertGreaterEqual(np.linalg.eigvalsh(comb).min(),-1e-14)
        self.assertAlmostEqual(np.trace(comb),4)
        marginal=np.trace(comb.reshape([2]*8),axis1=3,axis2=7).reshape(8,8)
        prefix=(1-w)*phi+w*np.eye(4)/2
        np.testing.assert_allclose(marginal,np.kron(prefix,T.I),atol=1e-14)
        np.testing.assert_allclose(np.trace(prefix.reshape(2,2,2,2),axis1=1,axis2=3),T.I)
        pt=(comb/4).reshape([2]*8).transpose(4,5,2,3,0,1,6,7).reshape(16,16)
        v=np.zeros(16);v[1]=1/np.sqrt(2);v[8]=-1/np.sqrt(2)
        self.assertAlmostEqual(v@pt@v,-w/8)

    def test_probe_and_statistical_scope_are_preserved(self):
        result=json.loads((STUDY/'results/ibm-identifiability.json').read_text())
        self.assertFalse(result['empirical_quantum_memory_excluded_or_certified'])
        for witness,report in zip(self.witnesses,result['reports']):
            nums,den=T.probability_numerators(witness,sorted(A.LABELS))
            np.testing.assert_array_equal(nums.sum(axis=1),den)
            probe=T.isolated_probe(witness,sorted(A.LABELS))
            self.assertGreater(probe['total_variation'],.0029)
            self.assertLessEqual(probe['total_variation'],probe['half_diamond_upper'])
            self.assertEqual(probe,report['isolated_instrument_probe'])
            self.assertTrue(report['aggregate_goodness_of_fit_not_certified'])
            self.assertGreater(report['multinomial_deviance_to_saturated'],16000)

    def test_fixed_assignment_ceiling_uses_a_rigorous_tail_upper_bound(self):
        from decimal import Decimal
        counts=np.array([[[7948,0,52,0]]])
        report=T.assignment_ceiling_diagnostic(counts,['zp,z,zp,z'],['example'])
        self.assertTrue(report['rejects_fixed_forward_assignment_point'])
        from scipy.stats import binom
        actual=binom.sf(7947,8000,.976)
        self.assertGreater(float(Decimal(report['binomial_upper_tail_bound'])),actual)
        self.assertLess(float(Decimal(report['binomial_upper_tail_bound'])),1e-32)


class JointStatisticsTests(unittest.TestCase):
    def test_exact_multinomial_pearson_moments_by_enumeration(self):
        from fractions import Fraction as F
        import math
        probs=[F(1,10),F(2,10),F(3,10),F(4,10)]
        for n in (1,2,4):
            mean=second=mass=F(0)
            for a in range(n+1):
                for b in range(n-a+1):
                    for c in range(n-a-b+1):
                        counts=[a,b,c,n-a-b-c]
                        pr=F(math.factorial(n))
                        for k,p in zip(counts,probs):pr*=p**k/math.factorial(k)
                        x=sum((k-n*p)**2/(n*p) for k,p in zip(counts,probs))
                        mass+=pr;mean+=pr*x;second+=pr*x*x
            self.assertEqual(mass,1);self.assertEqual(mean,3)
            self.assertEqual(second-mean*mean,6+(sum(1/p for p in probs)-22)/n)

    def test_directed_statistic_bounds_contain_exact_values(self):
        from decimal import Decimal,localcontext
        from fractions import Fraction as F
        counts=np.array([[2,3,1,4],[4,1,2,3]])
        nums=[1,2,3,4];result=J.pearson_region(counts,nums,10)
        exact=sum((F(int(k))-10*F(p,10))**2/(10*F(p,10)) for row in counts for k,p in zip(row,nums))
        var=2*(6+(sum(F(10,p) for p in nums)-22)/10)
        with localcontext() as ctx:
            ctx.prec=100
            x=Decimal(exact.numerator)/Decimal(exact.denominator)
            vv=Decimal(var.numerator)/Decimal(var.denominator)
            self.assertLessEqual(Decimal(result['pearson_lower']),x)
            self.assertGreaterEqual(Decimal(result['pearson_upper']),x)
            self.assertLessEqual(Decimal(result['variance_lower']),vv)
            self.assertGreaterEqual(Decimal(result['variance_upper']),vv)

    def test_small_exact_experiment_has_nominal_region_coverage(self):
        from fractions import Fraction as F
        import math
        probs=[F(1,100),F(2,100),F(3,100),F(94,100)]
        rejected=F(0);n=4
        for a in range(n+1):
            for b in range(n-a+1):
                for c in range(n-a-b+1):
                    counts=[a,b,c,n-a-b-c];pr=F(math.factorial(n))
                    for k,p in zip(counts,probs):pr*=p**k/math.factorial(k)
                    result=J.pearson_region(np.array([counts]),[1,2,3,94],100)
                    if result['aggregate_region_membership']=='excluded':rejected+=pr
        self.assertGreater(rejected,0)
        self.assertLessEqual(rejected,F(1,40))

    def test_joint_region_rejects_points_without_excluding_class(self):
        result=json.loads((STUDY/'results/ibm-joint-statistics.json').read_text())
        self.assertEqual(len(result['reports']),4)
        for row in result['reports']:
            self.assertEqual(row['aggregate_region_membership'],'excluded')
            self.assertFalse(row['classical_memory_class_excluded'])
            self.assertFalse(row['fitted_degrees_of_freedom_used'])

    def test_delay_probabilities_match_independent_choi_contraction(self):
        catalog=json.loads((STUDY/'results/ibm-delay-candidates.json').read_text())
        labels=sorted(A.LABELS)
        def unitary(vec,den):
            return (den*T.I-1j*sum(k*s for k,s in zip(vec,T.PAULIS)))/np.sqrt(float(den*den+sum(k*k for k in vec)))
        for candidate in catalog['candidates']:
            D.check_candidate(candidate);nums,dens=D.probabilities(candidate,labels)
            probs=np.asarray(nums/dens[...,None],float);S=candidate['spam_denominator']
            witness=candidate['instrument'];matrices=[]
            for pair in witness['choi']:
                matrices.append([(np.asarray(T.arrays(block)[0],float)+1j*np.asarray(T.arrays(block)[1],float))/witness['denominator'] for block in pair])
            for d,pair in enumerate(candidate['rotation_cayley']):
                pre,post=[unitary(v,candidate['rotation_denominator']) for v in pair]
                for row,label in enumerate(labels):
                    a,m,p,z=label.split(',');y=witness['settings'].index(m+','+p)
                    rho=(T.I+sum(k*s/S for k,s in zip(candidate['preparations'][a],T.PAULIS)))/2
                    rho=pre@rho@pre.conj().T
                    bias,*vec=candidate['effects'][z]
                    obs=bias/S*T.I+sum(k*s/S for k,s in zip(vec,T.PAULIS))
                    for b in range(2):
                        output=np.einsum('ij,iajb->ab',rho,matrices[y][b].reshape(2,2,2,2))
                        output=post@output@post.conj().T
                        for c in range(2):
                            expected=np.trace((T.I+(1-2*c)*obs)/2@output).real
                            self.assertAlmostEqual(probs[d,row,2*b+c],expected,places=12)

    def test_new_pair_has_independent_quantum_circuit_verification(self):
        from types import SimpleNamespace
        catalog=json.loads((STUDY/'results/ibm-delay-candidates.json').read_text())
        checker=SimpleNamespace(witnesses=[c['instrument'] for c in catalog['candidates']],assertGreater=self.assertGreater)
        IBMIdentifiabilityTests.test_independent_kraus_swap_circuit_on_full_operator_basis(checker)

    def test_nonphysical_spam_and_invalid_probabilities_fail_closed(self):
        catalog=json.loads((STUDY/'results/ibm-delay-candidates.json').read_text())
        for field,key,value in [('preparations','xp',[10**14,0,0]),('effects','x',[10**12,10**12,0,0])]:
            bad=copy.deepcopy(catalog['candidates'][0]);bad[field][key]=value
            with self.assertRaises(ValueError):D.check_candidate(bad)
        for nums,den in [([1,2,3,3],10),([0,2,3,5],10),([1,2,3,4],10.5)]:
            with self.assertRaises(ValueError):J.pearson_region(np.array([[1,2,3,4]]),nums,den)


if __name__ == '__main__':
    unittest.main()
