#!/usr/bin/env python3
"""Prospective six-state recovery certificate with complete reset calibration.

The analytic proof covers arbitrary classical environmental memory under the
explicit causal/anchoring assumptions in theory.md. Demonstrations are synthetic.
"""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.stats import beta

HERE = Path(__file__).resolve().parent
I = np.eye(2, dtype=complex)
P = np.array([[[0, 1], [1, 0]], [[0, -1j], [1j, 0]], [[1, 0], [0, -1]]])
STATES = np.array([(I+s*p)/2 for p in P for s in (1, -1)])
SWAP = np.array([[1, 0, 0, 0], [0, 0, 1, 0],
                 [0, 1, 0, 0], [0, 0, 0, 1]], dtype=complex)


def validate_counts(counts, shots, shape):
    counts = np.asarray(counts)
    if type(shots) is not int or shots <= 0 or counts.shape != shape:
        raise ValueError('invalid shot count or table shape')
    if counts.dtype.kind not in 'iu' or np.any(counts < 0) or np.any(counts > shots):
        raise ValueError('counts must be integers in [0, shots]')
    return counts


def infer(calibration, dynamics, shots=8000, alpha=.05):
    """Calibration has shape (output Pauli, input axis, input sign).

    Dynamics has one matching-output success count per six-state input.
    Bonferroni-Hoeffding events allocate alpha/2 to each complete source.
    Keep the same events when removing a source. No fit or target selection.
    A rejection alone does not establish joint quantum feasibility; callers
    must verify that separately before interpreting it as a memory certificate.
    """
    if not 0 < alpha < 1:
        raise ValueError('alpha must be in (0,1)')
    cal = validate_counts(calibration, shots, (3, 3, 2))
    dyn = validate_counts(dynamics, shots, (6,))
    mean = 2*cal/shots-1
    # +/-1-valued calibration observations.
    radius = math.sqrt(2*math.log(36/(alpha/2))/shots)
    lower, upper = np.maximum(-1., mean-radius), np.minimum(1., mean+radius)
    # All three input-axis pairs constrain the same affine translation t_a.
    tlo = np.max((lower[:, :, 0]+lower[:, :, 1])/2, axis=1)
    thi = np.min((upper[:, :, 0]+upper[:, :, 1])/2, axis=1)
    if np.any(tlo > thi):
        return {'status': 'calibration_affine_inconsistent', 'reject': False}
    tbound = np.maximum(np.abs(tlo), np.abs(thi))
    Tlo = (lower[:, :, 0]-upper[:, :, 1])/2
    Thi = (upper[:, :, 0]-lower[:, :, 1])/2
    Tbound = np.maximum(np.abs(Tlo), np.abs(Thi))
    epsilon = min(1., float((tbound.sum()+Tbound.sum())/2))
    dradius = math.sqrt(math.log(12/(alpha/2))/(2*shots))
    flo = float(np.maximum(0., dyn/shots-dradius).mean())
    fhi = float(np.minimum(1., dyn/shots+dradius).mean())
    null_upper = min(1., 2/3+epsilon)
    margin = flo-null_upper
    return {'status': 'conditional_bound', 'alpha': alpha,
            'calibration_mean_radius': radius, 'dynamics_probability_radius': dradius,
            'reset_epsilon_upper': epsilon, 'recovery_lower': flo,
            'recovery_upper': fhi, 'classical_memory_recovery_upper': null_upper,
            'required_epsilon_lower': max(0., flo-2/3),
            'margin': margin, 'reject': margin > 1e-12,
            'calibration_t_intervals': np.stack([tlo, thi], axis=-1).tolist(),
            'calibration_T_intervals': np.stack([Tlo, Thi], axis=-1).tolist()}


def depolarizing_kraus(lam):
    if not -1/3 <= lam <= 1:
        raise ValueError('unphysical depolarizing parameter')
    return [math.sqrt((1+3*lam)/4)*I] + [math.sqrt((1-lam)/4)*p for p in P]


def channel(rho, kraus):
    return sum(k@rho@k.conj().T for k in kraus)


def choi(kraus):
    # Input x output ordering; unnormalized Choi, trace = 2.
    return sum(np.outer(k.T.reshape(-1), k.T.reshape(-1).conj()) for k in kraus)


def calibration_probabilities(kraus):
    return np.array([[[(1+np.trace(p@channel(rho, kraus)).real)/2
                      for rho in STATES[2*b:2*b+2]] for b in range(3)] for p in P])


def recovery_probabilities(kraus):
    return np.array([np.trace(rho@channel(rho, kraus)).real for rho in STATES])


def swap_recovery(rho, lam):
    """Independent physical construction: SWAP, system reset, SWAP, noise."""
    joint = SWAP@np.kron(rho, I/2)@SWAP.conj().T
    reset = depolarizing_kraus(0)
    joint = sum(np.kron(k, I)@joint@np.kron(k.conj().T, I) for k in reset)
    joint = SWAP@joint@SWAP.conj().T
    system = np.trace(joint.reshape(2, 2, 2, 2), axis1=1, axis2=3)
    return channel(system, depolarizing_kraus(lam))


def cp_interval(k, n, alpha=.05):
    return [0. if k == 0 else float(beta.ppf(alpha/2, k, n-k+1)),
            1. if k == n else float(beta.ppf(1-alpha/2, k+1, n-k))]


def sample_inference(rng, calp, dynp, shots):
    cal = rng.binomial(shots, np.clip(calp, 0, 1))
    dyn = rng.binomial(shots, np.clip(dynp, 0, 1))
    return cal, dyn, infer(cal, dyn, shots)


def physical_checks():
    errors = []
    minimum_eigenvalue = 1.
    for lam in (0., .4, .75, .9, 1.):
        kraus = depolarizing_kraus(lam)
        errors.append(np.max(np.abs(sum(k.conj().T@k for k in kraus)-I)))
        minimum_eigenvalue = min(minimum_eigenvalue, np.linalg.eigvalsh(choi(kraus)).min())
        # Informationally complete inputs, plus complex off-diagonal matrix units:
        # equality on this spanning set proves equality for every terminal input/effect.
        units = [np.eye(4, dtype=complex)[i].reshape(2, 2) for i in range(4)]
        for rho in list(STATES)+units:
            errors.append(np.max(np.abs(swap_recovery(rho, lam)-channel(rho, kraus))))
    return {'maximum_residual': float(max(errors)),
            'minimum_choi_eigenvalue': float(minimum_eigenvalue)}


def demonstration(repetitions=1000):
    shots = 8000
    rng = np.random.default_rng(27092026)
    reset = depolarizing_kraus(0)
    calp = calibration_probabilities(reset)
    cal, dyn, result = sample_inference(rng, calp, np.ones(6), shots)
    # Complete source-only witness tables, not ablations of an optimizer.
    bypass_cal = calibration_probabilities(depolarizing_kraus(1))
    reset_dyn = recovery_probabilities(reset)
    calrad = result['calibration_mean_radius']/2
    drad = result['dynamics_probability_radius']
    witnesses = {
        'quantum_swap_joint_feasible': bool(np.all(np.abs(cal/shots-calp) <= calrad)
                                         and np.all(np.abs(dyn/shots-1) <= drad)),
        'dynamics_only_identity_bypass_feasible': bool(np.all(np.abs(dyn/shots-1) <= drad)),
        'identity_bypass_violates_calibration': bool(np.any(np.abs(cal/shots-bypass_cal) > calrad)),
        'calibration_only_reset_null_feasible': bool(np.all(np.abs(cal/shots-calp) <= calrad)),
        'reset_null_violates_dynamics': bool(np.any(np.abs(dyn/shots-reset_dyn) > drad)),
        'identity_bypass_recovery': 1., 'reset_null_recovery': .5,
    }
    power = []
    for lam in (.4, .75, .9, .97, 1.):
        rejected = sum(sample_inference(rng, calp, np.full(6, (1+lam)/2), shots)[2]['reject']
                       for _ in range(repetitions))
        power.append({'lambda': lam, 'repetitions': repetitions, 'rejections': rejected,
                      'binomial_95_interval': cp_interval(rejected, repetitions)})
    # These simulations are diagnostics, not the proof of composite-null coverage.
    z_dephase = [np.diag([1., 0.]), np.diag([0., 1.])]
    # Four-outcome tetrahedral measure-and-prepare: a classical register is enough.
    tetra = np.array([[1,1,1], [1,-1,-1], [-1,1,-1], [-1,-1,1]])/math.sqrt(3)
    tetra_kraus = [sum(v[j]*P[j] for j in range(3))/math.sqrt(8)+I/math.sqrt(8)
                  for v in tetra]
    rotated = np.array([1., 2., 3.])/math.sqrt(14)
    projector = (I+sum(rotated[j]*P[j] for j in range(3)))/2
    nulls = [
        ('reset_no_memory', calp, recovery_probabilities(reset)),
        ('classical_z_record_survives_reset', calp, recovery_probabilities(z_dephase)),
        ('classical_four_outcome_record', calp, recovery_probabilities(tetra_kraus)),
        ('classical_rotated_record', calp, recovery_probabilities([projector, I-projector])),
        ('failed_reset_identity_bypass', bypass_cal, np.ones(6)),
        ('partial_reset_bypass', calibration_probabilities(depolarizing_kraus(.2)),
         recovery_probabilities(depolarizing_kraus(.2))),
    ]
    controls = []
    for name, cp, dp in nulls:
        rejected = sum(sample_inference(rng, cp, dp, shots)[2]['reject']
                       for _ in range(repetitions))
        controls.append({'name': name, 'recovery_mean': float(dp.mean()),
                         'repetitions': repetitions, 'rejections': rejected,
                         'binomial_95_interval': cp_interval(rejected, repetitions)})
    # This adversary intentionally violates acquisition matching: calibrate reset,
    # use identity in dynamics. A rejection is an assumption failure, not a success.
    mismatch = infer(cal, np.full(6, shots, dtype=int), shots)
    # Averaging a flagged Pauli channel hides the classical correction register.
    flagged_residual = max(float(np.max(np.abs(
        sum(p.conj().T@(p@rho@p.conj().T)@p/4 for p in [I, *P])-rho)))
        for rho in STATES)
    return {'kind': 'synthetic_only', 'shots_per_setting': shots, 'seed': 27092026,
            'calibration_successes': cal.tolist(), 'dynamics_successes': dyn.tolist(),
            'certificate': result, 'source_removal_witnesses': witnesses,
            'physical_checks': physical_checks(), 'power': power, 'null_controls': controls,
            'unmodeled_calibration_mismatch_false_claim': mismatch['reject'],
            'hidden_pauli_flag_recovery_residual': flagged_residual,
            'leakage_and_environment_coupling_covered': False,
            'coverage_basis': 'analytic_Hoeffding_and_full_class_theorem_not_simulations'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = json.dumps(demonstration(), indent=2, sort_keys=True)+'\n'
    target = HERE/'results/synthetic-certificate.json'
    if args.check:
        if target.read_text() != result:
            raise SystemExit('synthetic certificate differs')
        print('Synthetic certificate, physical models, source removal and simulations reproduce')
    else:
        target.parent.mkdir(exist_ok=True)
        target.write_text(result)
        print(target)


if __name__ == '__main__':
    main()
