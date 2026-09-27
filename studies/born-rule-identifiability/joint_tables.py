"""Measured preparation-scan / Viviani combination with explicit ordinary qutrit certificates."""
import argparse
import json
from pathlib import Path
import numpy as np
from br_table_sources import audit, manifest
from br_tables import (directions, h, fit_table, validate_fit, model, probabilities,
                       spectral_boundary, deviance, centered_singular, noise_radius,
                       contextual_response)

HERE = Path(__file__).resolve().parent
KEYS = ['scan', 'viviani']


def plain(value):
    if isinstance(value, dict):
        return {k: plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [plain(v) for v in value]
    if isinstance(value, (np.floating, float)):
        return float(value)
    if isinstance(value, np.integer):
        return int(value)
    return value


def dumps(value, level=0):
    """Readable records with compact numeric arrays (avoid tens of thousands of scalar-only lines)."""
    pad = '  '*level
    if isinstance(value, dict):
        return '{\n'+',\n'.join('  '*(level+1)+json.dumps(k)+': '+dumps(v, level+1)
                                for k, v in value.items())+'\n'+pad+'}'
    if isinstance(value, list):
        compact = json.dumps(value, allow_nan=False)
        if len(compact) <= 160 or all(not isinstance(v, (list, dict)) for v in value):
            return compact
        return '[\n'+',\n'.join('  '*(level+1)+dumps(v, level+1) for v in value)+'\n'+pad+']'
    return json.dumps(value, allow_nan=False)


def protocol():
    return json.loads((HERE/'table-protocol.json').read_text())


def assert_close(actual, expected):
    if isinstance(expected, dict):
        if actual.keys() != expected.keys():
            raise ValueError('Changed result keys')
        for key in expected:
            assert_close(actual[key], expected[key])
    elif isinstance(expected, list):
        if len(actual) != len(expected):
            raise ValueError('Changed result shape')
        for a, b in zip(actual, expected):
            assert_close(a, b)
    elif isinstance(expected, float):
        if not np.isclose(actual, expected, rtol=1e-9, atol=1e-11):
            raise ValueError('Numerical result changed')
    elif actual != expected:
        raise ValueError('Result changed')


def ideal_checks():
    results = {}
    for key in KEYS:
        r, m = directions(key)
        x = r @ m.T
        results[key] = {'ordinary_table': ((1+x)/2).tolist(), 'signed_deformation': {}}
        for theta in [-.01, 0, .01]:
            p = (1+h(x, theta))/2
            value = {'centered_singular_values': centered_singular(p).tolist()}
            if key == 'viviani':
                value['appended_determinant'] = float(np.linalg.det(np.c_[p, np.ones(len(p))]))
                # Independently derived symbolic polynomial, no sympy runtime dependency.
                value['exact_polynomial'] = -3*theta*(3*theta+16)/1024
            results[key]['signed_deformation'][str(theta)] = value
    # Fixed deterministic generic directions: numerical attainability illustration;
    # harmonic argument establishing existence is in joint-tables.md.
    rng = np.random.default_rng(4411)
    r, m = rng.normal(size=(16, 3)), rng.normal(size=(16, 3))
    r /= np.linalg.norm(r, axis=1)[:, None]
    m /= np.linalg.norm(m, axis=1)[:, None]
    results['generic_singular_values'] = {str(t): np.linalg.svd((1+h(r@m.T, t))/2, compute_uv=False).tolist()
                                          for t in [0, -.01, .01]}
    return results


def profile(p, kind, spec):
    return [fit_table(p, kind, t, spec['tangent_coordinate_box'], starts=spec['fit_starts'])
            for t in spec['theta_grid']]


def joint_choice(profiles):
    objectives = [sum(f[i]['frobenius_residual']**2/np.size(f[i]['prediction']) for f in profiles.values())
                  for i in range(len(next(iter(profiles.values()))))]
    return int(np.argmin(objectives)), objectives


def simulations(result, data, spec):
    rng = np.random.default_rng(spec['simulation_seed'])
    out = {}
    reps = spec['simulation_repetitions']
    for kind in KEYS:
        a = np.array(data['records'][kind]['counts_0_1_by_job'])
        shots = a.sum(axis=3)
        n = len(a[0])
        fit = result['records'][kind]['ordinary']
        v = np.array(fit['parameters'])
        stable = model(v, kind)[0]
        deformation = model(v, kind, .01)[0]
        leakage = np.array(result['records'][kind]['leakage'][-1]['prediction'])
        # Two ordinary qubit regimes with correlated preparation/measurement drift.
        # These deliberately violate pooled stationarity; no changed probability law.
        plus, minus = v.copy(), v.copy()
        shift = rng.choice([-1, 1], 2*n+8)*.08
        plus[:2*n+8] += shift
        minus[:2*n+8] -= shift
        drift = np.stack([model(plus if j % 2 else minus, kind)[0] for j in range(len(a))])
        expected = {'stable_qubit': stable, 'deformation_0.01': deformation,
                    'ordinary_qutrit': leakage, 'ordinary_correlated_drift': drift}
        # Preparation label retained by a controller: reset readout with small chance,
        # a legitimate classical-memory competitor, not a two-level independent instrument.
        memory = stable.copy()
        memory += .003*np.linspace(-1, 1, n)[:, None]*np.array([1, -1, 1, -1])[None, :]
        expected['ordinary_preparation_memory'] = np.clip(memory, 0, 1)
        row = {}
        threshold = noise_radius(shots.sum(axis=0), .025)
        derivative = (deformation-stable)/.01
        for label, p in expected.items():
            rejected, estimates = 0, []
            for _ in range(reps):
                k = rng.binomial(shots, p)
                observed = k.sum(axis=0)/shots.sum(axis=0)
                rejected += np.linalg.norm(centered_singular(observed)[3:]) > threshold
                estimates.append(float(np.sum((observed-stable)*derivative)/np.sum(derivative**2)))
            row[label] = {'stable_qubit_class_rejections': int(rejected), 'replications': reps,
                          'theta_projection_mean_known_nuisance': float(np.mean(estimates)),
                          'theta_projection_sd_known_nuisance': float(np.std(estimates, ddof=1))}
        out[kind] = row
    return out


def job_constructions(data, kind):
    """Ordinary drift competitor at the level actually recorded: one physical instrument per job."""
    a = np.array(data['records'][kind]['counts_0_1_by_job'])
    fits = []
    for k, job in enumerate(a):
        fit = fit_table(probabilities(job[None]), kind, leakage_cap=.005, starts=1)
        if fit['max_abs_residual'] > 2e-8:
            fit = fit_table(probabilities(job[None]), kind, leakage_cap=.005, starts=2)
        fits.append(fit)
    return fits


def compute(data):
    spec = protocol()
    result = {'protocol': spec, 'ideal_checks': ideal_checks(), 'records': {}}
    profiles, train_profiles = {}, {}
    for kind in KEYS:
        a = np.array(data['records'][kind]['counts_0_1_by_job'])
        p = probabilities(a)
        cut = len(a)*2//3
        train, holdout = probabilities(a[:cut]), probabilities(a[cut:])
        print('Fit', kind, flush=True)
        fits = profile(p, kind, spec)
        trainfits = profile(train, kind, spec)
        ordinary = fits[spec['theta_grid'].index(0)]
        leaks = [fit_table(p, kind, leakage_cap=cap, starts=spec['fit_starts']) for cap in spec['leakage_caps']]
        trainleak = fit_table(train, kind, leakage_cap=spec['leakage_caps'][-1], starts=spec['fit_starts'])
        best = min(fits, key=lambda x: x['frobenius_residual'])
        sensitivity = {str(box): fit_table(p, kind, best['theta'], box, starts=spec['fit_starts'])
                       for box in spec['sensitivity_boxes']}
        results = {'mean_table': p.tolist(), 'boundary': spectral_boundary(a),
                   'ordinary': ordinary, 'deformation_profile': fits, 'leakage': leaks,
                   'best_separate_grid_theta': best['theta'], 'axis_box_sensitivity': sensitivity,
                   'train_jobs': cut, 'holdout_jobs': len(a)-cut,
                   'train_table': train.tolist(), 'holdout_table': holdout.tolist(),
                   'train_profile': trainfits, 'train_leakage': trainleak,
                   'holdout_train_leakage_max_error': float(abs(holdout-np.array(trainleak['prediction'])).max()),
                   'holdout_train_leakage_deviance': deviance(a[cut:], trainleak['prediction']),
                   'ordinary_context_reset_upper': float(contextual_response(p, np.array(ordinary['prediction'])).max())}
        if kind == 'viviani':
            jobs = a[:, :, :, 0]/a.sum(axis=3)
            results['pooled_determinant'] = float(np.linalg.det(np.c_[p, np.ones(len(p))]))
            results['mean_job_determinant'] = float(np.mean([np.linalg.det(np.c_[q, np.ones(len(q))]) for q in jobs]))
        result['records'][kind] = results
        print('Construct job-specific ordinary models', kind, flush=True)
        results['job_leakage_constructions'] = job_constructions(data, kind)
        profiles[kind], train_profiles[kind] = fits, trainfits
    choice, objectives = joint_choice(profiles)
    train_choice, train_objectives = joint_choice(train_profiles)
    joint = {'grid_theta': spec['theta_grid'], 'equal_acquisition_mse_profile': objectives,
             'selected_grid_theta': spec['theta_grid'][choice],
             'train_selected_grid_theta': spec['theta_grid'][train_choice],
             'train_equal_acquisition_mse_profile': train_objectives,
             'remove_scan_theta': result['records']['viviani']['best_separate_grid_theta'],
             'remove_viviani_theta': result['records']['scan']['best_separate_grid_theta'],
             'add_or_remove_lima_quarter_turn_family': 'no identification gain under its exact common angle warp',
             'shared_calibration_parameters': [], 'holdout': {}}
    for kind in KEYS:
        r = result['records'][kind]
        a = np.array(data['records'][kind]['counts_0_1_by_job'])
        fitted = r['train_profile'][train_choice]
        q = np.array(fitted['prediction'])
        joint['holdout'][kind] = {'prediction_frozen_on_train': q.tolist(),
            'max_abs_error': float(abs(q-np.array(r['holdout_table'])).max()),
            'binomial_deviance_descriptive': deviance(a[r['train_jobs']:], q)}
    joint['ordinary_mean_table_distance_lower'] = max(r['boundary']['mean_table_distance_lower'] for r in result['records'].values())
    joint['ordinary_iid_shot_distance_lower_95pct'] = max(r['boundary']['iid_shot_distance_lower'] for r in result['records'].values())
    joint['ordinary_independent_job_distance_lower_95pct'] = max(r['boundary']['independent_job_distance_lower'] for r in result['records'].values())
    result['joint'] = joint
    # Existing IonQ acquisition tests how pooling changes the apparent table rank.
    a = np.array(data['records']['aria_drift_control']['counts_0_1_by_job'])
    p = probabilities(a)
    jobs = a[:, :, :, 0]/a.sum(axis=3)
    result['aria_drift_control'] = {'boundary': spectral_boundary(a),
        'pooled_determinant': float(np.linalg.det(np.c_[p, np.ones(5)])),
        'job_determinants': [float(np.linalg.det(np.c_[q, np.ones(5)])) for q in jobs],
        'role': 'measured drift stress control; not an independent calibration of Nairobi'}
    result['simulations'] = simulations(result, data, spec)
    return plain(result)


def check(result, data):
    if result['protocol'] != protocol() or data['source'] != manifest():
        raise ValueError('Protocol or source mismatch')
    for kind in KEYS:
        row = result['records'][kind]
        a = np.array(data['records'][kind]['counts_0_1_by_job'])
        p = probabilities(a)
        if not np.allclose(p, row['mean_table'], atol=1e-14, rtol=0):
            raise ValueError('Counts do not reproduce means')
        for name in ['deformation_profile', 'leakage']:
            for fit in row[name]:
                validate_fit(fit, kind, p)
        for fit in row['axis_box_sensitivity'].values():
            validate_fit(fit, kind, p)
        for fit in row['train_profile']+[row['train_leakage']]:
            validate_fit(fit, kind, probabilities(a[:row['train_jobs']]))
        if len(row['job_leakage_constructions']) != len(a):
            raise ValueError('Missing job compatibility constructions')
        for job, fit in zip(a, row['job_leakage_constructions']):
            validate_fit(fit, kind, probabilities(job[None]))
        # Constructive, physical ordinary qutrit explanation must fit every cell.
        if min(f['max_abs_residual'] for f in row['leakage']) > 2e-8:
            raise ValueError('No accurate physical joint-compatibility construction')
        assert_close(row['boundary'], spectral_boundary(a))
    assert_close(result['ideal_checks'], ideal_checks())
    profiles = {kind: result['records'][kind]['deformation_profile'] for kind in KEYS}
    choice, objectives = joint_choice(profiles)
    assert_close(result['joint']['equal_acquisition_mse_profile'], objectives)
    assert_close(result['joint']['selected_grid_theta'], protocol()['theta_grid'][choice])
    train_profiles = {kind: result['records'][kind]['train_profile'] for kind in KEYS}
    choice, objectives = joint_choice(train_profiles)
    assert_close(result['joint']['train_equal_acquisition_mse_profile'], objectives)
    assert_close(result['joint']['train_selected_grid_theta'], protocol()['theta_grid'][choice])
    for kind in KEYS:
        q = np.array(train_profiles[kind][choice]['prediction'])
        a = np.array(data['records'][kind]['counts_0_1_by_job'])
        cut = result['records'][kind]['train_jobs']
        assert_close(result['joint']['holdout'][kind]['prediction_frozen_on_train'], q.tolist())
        assert_close(result['joint']['holdout'][kind]['max_abs_error'], float(abs(q-probabilities(a[cut:])).max()))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache', type=Path)
    parser.add_argument('--refit', action='store_true', help='Repeat nonconvex optimization and simulations')
    parser.add_argument('--check', action='store_true', help='Validate saved physical certificates, not optimizer identity')
    parser.add_argument('--output-dir', type=Path, default=HERE/'results')
    args = parser.parse_args()
    data = audit(args.cache) if args.cache else json.loads((HERE/'results/table-audit.json').read_text())
    if args.refit or not args.check:
        result = compute(data)
    else:
        result = json.loads((HERE/'results/joint-tables.json').read_text())
    check(result, data)
    from br_table_report import render
    report = render(result, data)
    if args.check:
        stored = json.loads((HERE/'results/table-audit.json').read_text())
        if data != stored:
            raise ValueError('Source audit changed')
        if not args.refit and report != (HERE/'results/joint-table-report.md').read_text():
            raise ValueError('Report changed')
        print('Measured joint-table sources, physical certificates and bounds passed')
    else:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        for name, obj in [('table-audit.json', data), ('joint-tables.json', result)]:
            (args.output_dir/name).write_text(dumps(obj)+'\n')
        (args.output_dir/'joint-table-report.md').write_text(report)


if __name__ == '__main__':
    main()
