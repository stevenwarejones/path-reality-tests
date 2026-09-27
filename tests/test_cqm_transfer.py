"""Independent checks of the forcing model, causal clock and path statistics."""
from itertools import product
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np
from scipy.linalg import expm

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'studies/correlated-qubit-mechanisms'))
from cqm_transfer import (causal_filter, read_forcing, reservoir, transition, forward,
                          duration_predictions, run_duration_counts, summarize)


class TransferChecks(unittest.TestCase):
    def test_exact_markov_transition_and_support_identity(self):
        params=np.array([2.,50.,1700.])
        energy=np.array([0.,.2,1.])
        p01,p11=transition(energy,params)
        for i,e in enumerate(energy):
            u=params[0]+params[1]*e;d=params[2]
            expected=expm(np.array([[-u,u],[d,-d]])*.001)
            np.testing.assert_allclose([p01[i],p11[i]],expected[:,1],atol=1e-14)
        cap=energy.max()
        for a,b in zip(transition(energy,params),transition(np.minimum(energy,cap),params)):
            np.testing.assert_array_equal(a,b)
        self.assertNotEqual(transition(np.array([100.]),params)[0][0],transition(np.array([cap]),params)[0][0])
        with self.assertRaises(ValueError):transition(np.array([-1.]),params)

    def test_forcing_uses_only_available_samples(self):
        t=np.arange(100)*.0001;a=np.arange(100,dtype=float)
        forcing=(t,a,0.,.0001)
        before=reservoir(forcing,np.array([.00305]),.01)
        changed=a.copy();changed[31:]=1e9
        after=reservoir((t,changed,0.,.0001),np.array([.00305]),.01)
        np.testing.assert_array_equal(before,after)
        np.testing.assert_array_equal(causal_filter(a,.0001,.01)[:31],causal_filter(changed,.0001,.01)[:31])
        with self.assertRaises(ValueError):reservoir(forcing,np.array([.02]),.01)
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'forcing.npy'
            voltage=np.full(100,7.);voltage[50:]=20.
            trigger=np.zeros(100);trigger[50:]=2.
            np.save(path,np.array([t,voltage,trigger]))
            clock,centered,edge,step=read_forcing(path,center=True)
            np.testing.assert_array_equal(centered[:50],np.zeros(50))
            np.testing.assert_array_equal(centered[50:],np.full(50,65.))
            self.assertEqual(edge,t[50])

    def test_complete_duration_expectations_by_path_enumeration(self):
        energy=np.linspace(.1,1,8)
        parameters=np.array([[200.,300.,400.],[100.,100.,500.]])
        p,onsets,transitions=forward(np.array([0,0]),energy,parameters)
        p01,p11=transitions[0]
        expected=np.zeros(4);marginals=np.zeros(9)
        for tail in product([0,1],repeat=8):
            states=np.array((0,)+tail,dtype=bool)
            probability=1.
            for j in range(8):
                next_true=p11[j] if states[j] else p01[j]
                probability*=next_true if states[j+1] else 1-next_true
            count,_=run_duration_counts(states)
            expected+=probability*count
            marginals+=probability*states
        np.testing.assert_allclose(duration_predictions(onsets[0],transitions[0]),expected,atol=1e-13)
        np.testing.assert_allclose(p[0],marginals,atol=1e-13)

    def test_boundaries_do_not_create_onsets(self):
        flags=np.zeros((2,2,20),dtype=bool)
        flags[:,0,-1]=True;flags[:,1,0]=True
        r=summarize(flags,np.zeros(19),np.array([[1.,1.,100.],[1.,1.,200.]]),np.arange(20)*.001)
        self.assertEqual(r['total']['exposure_samples'],38)
        self.assertEqual(r['observed_onset_pairs'][10],1)
        self.assertEqual(r['boundary_censored_observed_runs'],[2,2])
        self.assertEqual(np.sum(r['observed_complete_flag_runs']),0)


if __name__=='__main__':unittest.main()
