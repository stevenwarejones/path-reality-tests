#!/usr/bin/env python3
"""Matched retrospective information comparison and finite-window generative checks."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from sklearn.ensemble import HistGradientBoostingClassifier

from cqm_analysis import close
from cqm_mechanical import (HERE, PERIOD, PHASE_ORIGIN, SAMPLE_PERIOD,
                            iter_acquisitions)
from cqm_sources import acquire

FIRST = 250000
STRIDE = 330
WINDOW = 33
SEGMENT = 16666
HISTORY_WIDTHS = (33, 333, 3333, 33333)
FEATURE_NAMES = ['elapsed_s'] + [f'{q}_{name}' for q in ('A', 'B') for name in
    ([f'nonG_last_{w}' for w in HISTORY_WIDTHS] +
     [f'nonG_segment_{i}' for i in range(15)] + ['log_ground_age'])] + ['phase_sin', 'phase_cos']
LEARNER = dict(max_iter=100, max_leaf_nodes=15, min_samples_leaf=100,
               learning_rate=.1, l2_regularization=1., early_stopping=False, random_state=19)


def local_marks(states, starts):
    """0 or six onset/occupation marks, calculated only in the future 33 readouts."""
    future = states[starts[:, None]+np.arange(1, WINDOW+1)] != 0
    triples = future[:, :-2] & future[:, 1:-1] & future[:, 2:]
    exists = triples.any(axis=1)
    onset_bin = np.argmax(triples, axis=1)//11
    long = future.sum(axis=1) >= 12
    return np.where(exists, 1+2*onset_bin+long, 0).astype(np.int8)


def run_features(states, lag, starts=None):
    if starts is None:
        starts = np.arange(FIRST, len(states[0])-WINDOW, STRIDE)
    starts = np.asarray(starts)
    if len(states) != 2 or len(states[0]) != len(states[1]) or np.any(starts < FIRST) or np.any(starts+WINDOW >= len(states[0])):
        raise ValueError('invalid history/target boundary')
    eligible = (states[0][starts] == 0) & (states[1][starts] == 0)
    starts = starts[eligible]
    features = [starts*SAMPLE_PERIOD]
    activity = np.zeros(len(starts))
    marks = []
    for state in states:
        nong = state != 0
        cumulative = np.r_[0, np.cumsum(nong, dtype=np.int64)]
        for width in HISTORY_WIDTHS:
            features.append((cumulative[starts]-cumulative[starts-width])/width)
        activity += (cumulative[starts]-cumulative[starts-3333])/3333/2
        for i in range(15):
            right = starts-i*SEGMENT
            features.append((cumulative[right]-cumulative[right-SEGMENT])/SEGMENT)
        last = np.maximum.accumulate(np.where(nong, np.arange(len(state)), -FIRST))
        features.append(np.log1p(np.minimum(starts-last[starts], FIRST)))
        marks.append(local_marks(state, starts))
    phase = ((starts*SAMPLE_PERIOD + lag + PHASE_ORIGIN) % PERIOD)/PERIOD
    features += [np.sin(2*np.pi*phase), np.cos(2*np.pi*phase)]
    cells = (np.floor(phase*36).astype(int)*4 + np.searchsorted([.01,.05,.2], activity, side='right'))
    return np.column_stack(features).astype(np.float32), np.column_stack(marks), cells.astype(np.int16)


def extract(cache):
    features, marks, cells, runs = [], [], [], []
    for run, states, lag in iter_acquisitions(cache):
        x, m, c = run_features(states, lag)
        features.append(x)
        marks.append(m)
        cells.append(c)
        runs.append(np.full(len(m), run, dtype=np.int16))
    return dict(x=np.concatenate(features), marks=np.concatenate(marks),
                cells=np.concatenate(cells), runs=np.concatenate(runs))


def binary_scores(y, p):
    p = np.clip(p, 1e-12, 1-1e-12)
    return {'windows': len(y), 'observed': int(y.sum()), 'predicted': float(p.sum()),
            'log_loss': float(-np.mean(y*np.log(p)+(1-y)*np.log1p(-p))),
            'brier': float(np.mean((y-p)**2))}


def compare_predictors(data):
    x, runs = data['x'], data['runs']
    y = np.all(data['marks'] > 0, axis=1)
    train = runs < 1550
    selectors = {'history_only': list(range(len(FEATURE_NAMES)-2)),
                 'phase_only': [0, len(FEATURE_NAMES)-2, len(FEATURE_NAMES)-1],
                 'combined': list(range(len(FEATURE_NAMES)))}
    output, scores = {}, {}
    for name in ['constant']+list(selectors):
        print('Fit predictor:', name, flush=True)
        if name == 'constant':
            p_train = np.full(train.sum(), y[train].mean())
            p_test = np.full((~train).sum(), y[train].mean())
        else:
            columns = selectors[name]
            fit = HistGradientBoostingClassifier(**LEARNER).fit(x[train][:, columns], y[train])
            p_train = fit.predict_proba(x[train][:, columns])[:, 1]
            p_test = fit.predict_proba(x[~train][:, columns])[:, 1]
        threshold = float(np.quantile(p_train, .15))
        keep = p_test <= threshold
        block = []
        for b in range(31, 62):
            select = runs[~train]//50 == b
            block.append({'acquisitions': [b*50, b*50+49], **binary_scores(y[~train][select], p_test[select])})
        output[name] = {'validation': binary_scores(y[~train], p_test), 'blocks': block,
                        'policy_training_threshold': threshold,
                        'policy_training_coverage': float(np.mean(p_train <= threshold)),
                        'policy_validation_coverage': float(keep.mean()),
                        'policy_validation': binary_scores(y[~train][keep], p_test[keep]) if keep.any() else None}
        scores[name] = p_test
    differences = {}
    for name in ['history_only', 'phase_only']:
        delta = [a['log_loss']-b['log_loss'] for a,b in zip(output[name]['blocks'], output['combined']['blocks'])]
        differences[name+'_minus_combined'] = {'pooled_log_loss_difference': output[name]['validation']['log_loss']-output['combined']['validation']['log_loss'],
                'block_differences': delta, 'blocks_combined_better': sum(v>0 for v in delta),
                'scope': 'Paired descriptive differences, not independent-block confidence intervals.'}
    return {'predictors': output, 'differences': differences}


def mark_tables(data):
    code = data['marks'][:, 0]*7+data['marks'][:, 1]
    result = []
    for selected in [data['runs'] < 1550, data['runs'] >= 1550]:
        result.append(np.bincount(data['cells'][selected].astype(int)*49+code[selected], minlength=144*49).reshape(144,7,7))
    return result


def smooth_tables(train):
    global_distribution = train.sum(0)+.5
    global_distribution /= global_distribution.sum()
    return train + 100*global_distribution


def environment_mixture(counts):
    """Two shared latent states; local categorical marks independent given state."""
    a0 = counts.sum(2)/counts.sum((1,2))[:, None]
    b0 = counts.sum(1)/counts.sum((1,2))[:, None]
    best = None
    for strength in [.2,.5,.8]:
        # Split components towards/away from the non-event state, with positive mass.
        v = np.array([-1.,1.,1.,1.,1.,1.,1.])
        a = np.stack([a0*(1-strength*v), a0*(1+strength*v)], axis=1)
        b = np.stack([b0*(1-strength*v), b0*(1+strength*v)], axis=1)
        a /= a.sum(2, keepdims=True)
        b /= b.sum(2, keepdims=True)
        weights = np.full((len(counts),2), .5)
        for _ in range(300):
            component = weights[:,:,None,None]*a[:,:,:,None]*b[:,:,None,:]
            responsibility = component/np.maximum(component.sum(1,keepdims=True), 1e-300)
            assigned = counts[:,None,:,:]*responsibility
            mass = assigned.sum((2,3))
            weights = mass/mass.sum(1,keepdims=True)
            a = assigned.sum(3)/mass[:,:,None]
            b = assigned.sum(2)/mass[:,:,None]
        p = (weights[:,:,None,None]*a[:,:,:,None]*b[:,:,None,:]).sum(1)
        objective = (counts*np.log(np.maximum(p,1e-300))).sum((1,2))
        if best is None:
            best = p, objective
        else:
            improve = objective > best[1]
            best[0][improve] = p[improve]
            best[1][improve] = objective[improve]
    return best[0]


def joint_component(counts):
    """Local product plus a common mark on the six nonzero diagonal categories."""
    mass = counts.sum((1,2))
    a = counts.sum(2)/mass[:,None]
    b = counts.sum(1)/mass[:,None]
    weights = np.full(len(counts), .01)
    d = np.full((len(counts),6), 1/6)
    for _ in range(300):
        local = (1-weights[:,None,None])*a[:,:,None]*b[:,None,:]
        common = np.zeros_like(local)
        common[:,np.arange(1,7),np.arange(1,7)] = weights[:,None]*d
        assigned = counts*common/np.maximum(local+common, 1e-300)
        residual = counts-assigned
        common_mass = assigned.sum((1,2))
        local_mass = residual.sum((1,2))
        weights = common_mass/mass
        a = residual.sum(2)/local_mass[:,None]
        b = residual.sum(1)/local_mass[:,None]
        d = assigned[:,np.arange(1,7),np.arange(1,7)]/common_mass[:,None]
    p = (1-weights[:,None,None])*a[:,:,None]*b[:,None,:]
    p[:,np.arange(1,7),np.arange(1,7)] += weights[:,None]*d
    return p


def generative_checks(data):
    train, test = mark_tables(data)
    retained = smooth_tables(train)
    total = retained.sum((1,2))
    local = (retained.sum(2)/total[:,None])[:,:,None]*(retained.sum(1)/total[:,None])[:,None,:]
    models = {'independent_local': local, 'shared_environment': environment_mixture(retained),
              'additional_joint_component': joint_component(retained),
              'joint_table_benchmark': retained/total[:,None,None]}
    a,b = np.indices((7,7))
    masks = {'paired': (a>0)&(b>0),
             'same_onset_bin': (a>0)&(b>0)&((a-1)//2 == (b-1)//2),
             'both_at_least_36us_nonG': (a>0)&(b>0)&(a%2==0)&(b%2==0)}
    result = {}
    for name,p in models.items():
        if not np.isfinite(p).all() or np.any(p<=0) or not np.allclose(p.sum((1,2)),1,atol=1e-10):
            raise ValueError('invalid generative probabilities')
        expected = test.sum((1,2))[:,None,None]*p
        scores = {'log_loss': float(-(test*np.log(p)).sum()/test.sum()), 'structure': {}}
        for key,mask in masks.items():
            scores['structure'][key] = {'observed': int(test[:,mask].sum()), 'predicted': float(expected[:,mask].sum())}
        block_losses=[]
        for block in range(31,62):
            selected = data['runs']//50 == block
            cells = data['cells'][selected]
            marks = data['marks'][selected]
            block_losses.append(float(-np.log(p[cells,marks[:,0],marks[:,1]]).mean()))
        scores['block_log_losses'] = block_losses
        result[name] = scores
    return {'models': result, 'training_windows': int(train.sum()), 'validation_windows': int(test.sum()),
            'scope': 'Finite-window observed mark laws, not full physical or whole-trace compatibility witnesses.'}


def array_digest(value):
    return hashlib.sha256(value.tobytes()).hexdigest()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache', type=Path, required=True)
    parser.add_argument('--features', type=Path, required=True, help='Reconstructible scratch cache, outside repository')
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    if args.features.resolve().is_relative_to(HERE.parents[1]):
        raise ValueError('feature cache must be outside repository')
    acquire(args.cache, verify_only=True, mechanical_prediction=True)
    # Re-extract on every invocation: no stale cached features may stand in for source checks.
    data = extract(args.cache)
    np.savez_compressed(args.features, **data)
    result = {'status': 'Retrospective information comparison using deposited, not verified causal, phase',
              'feature_names': FEATURE_NAMES, 'learner': LEARNER,
              'extraction_digests': {k:array_digest(v) for k,v in data.items()},
              **compare_predictors(data), 'generative_checks': generative_checks(data)}
    output = HERE/'results/complementarity.json'
    if args.check:
        close(result, json.loads(output.read_text()))
        print('Complementarity and generative checks reproduced')
    else:
        output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
        print(output)
