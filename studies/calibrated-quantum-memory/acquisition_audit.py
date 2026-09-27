#!/usr/bin/env python3
"""Reconstruct selected measured schemas without executing source code/pickles.

This is a provenance/observability audit, not a memory certificate or fit.
"""
import argparse
from collections import Counter
from datetime import datetime
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import pickletools
import urllib.request

import numpy as np

HERE = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location('nmn_source_audit', HERE/'audit.py')
audit = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(audit)
QM = 'quantum_memory_data'
WHITE = 'process-tensor-network-tomography'


def primitive_float_list(data):
    """Decode ONLY a protocol-4 list of finite floats via a strict grammar.

    No unpickler, imports, GLOBAL, REDUCE, BUILD, or object construction.
    Frames are serialization boundaries; ignore them only after opcode parsing.
    """
    ops = list(pickletools.genops(data))
    if not ops or ops[-1][2] != len(data)-1:
        raise ValueError('trailing or incomplete pickle')
    tokens = [(op.name, arg) for op, arg, _ in ops if op.name != 'FRAME']
    if tokens[:3] != [('PROTO', 4), ('EMPTY_LIST', None), ('MEMOIZE', None)]:
        raise ValueError('not the declared primitive list grammar')
    values, i = [], 3
    while i < len(tokens)-1:
        if tokens[i] != ('MARK', None):
            raise ValueError('expected flat batch')
        i += 1
        start = len(values)
        while i < len(tokens) and tokens[i][0] == 'BINFLOAT':
            value = tokens[i][1]
            if not math.isfinite(value):
                raise ValueError('nonfinite float')
            values.append(value)
            i += 1
        if len(values) == start or i >= len(tokens) or tokens[i] != ('APPENDS', None):
            raise ValueError('invalid float batch')
        i += 1
    if tokens[i:] != [('STOP', None)]:
        raise ValueError('invalid termination')
    return values


def alternative_regroup(rows, bit, branch):
    audit.validate_run(rows)
    if bit not in (0, 1) or branch not in ('0', '1'):
        raise ValueError('invalid alternative')
    result = {}
    for label, counts in rows.items():
        a, m, p, z = label.split(',')
        partner = rows[','.join((a, m, audit.opposite(p), z))]
        result[label] = {o: (partner if o[bit] == branch else counts)[o]
                         for o in audit.OUTCOMES}
    return result


def mapping_audit(source_dir):
    manifest = json.loads((HERE/'sources.json').read_text())
    records = []
    for entry in manifest['nmn_files']:
        raw = audit.verify_file(source_dir/entry['name'], entry)
        if not entry['name'].endswith('.json'):
            continue
        for delay, rows in sorted(json.loads(raw).items()):
            candidates = {'literal_archive_rows': rows}
            for bit in (0, 1):
                for branch in ('0', '1'):
                    candidates[f'opposite_if_character_{bit}_is_{branch}'] = alternative_regroup(rows, bit, branch)
            diagnostics = {}
            for name, mapped in candidates.items():
                totals = [sum(r.values()) for r in mapped.values()]
                if sum(totals) != sum(sum(r.values()) for r in rows.values()):
                    raise ValueError('alternative lost counts')
                diagnostics[name] = dict(min_total=min(totals), max_total=max(totals),
                                         distinct_totals=len(set(totals)))
            trace_labels = ['xp,x,xp,x', 'xp,x,xm,x']
            records.append(dict(source=entry['name'], delay=delay, alternatives=diagnostics,
                                trace_archived={k:rows[k] for k in trace_labels},
                                trace_inferred_physical={k:audit.regroup(rows)[k] for k in trace_labels}))
    return dict(status='unresolved_historical_transformation',
                note='Equal denominators do not prove semantics; opposite branch choices relabel circuits.',
                runs=records)


def qpt_counts(state, labels, n):
    if labels != {'qubit':['q0'], 'initial_state':['0','1','+','i+'], 'axis':['x','y','z']}:
        raise ValueError('unexpected axis/bit ordering')
    if state.shape != (1, n, 4, 3) or not np.issubdtype(state.dtype, np.integer):
        raise ValueError('unexpected shot dimensions/type')
    if not np.isin(state, [0, 1]).all():
        raise ValueError('nonbinary terminal outcome')
    return state[0].sum(axis=0).astype(int).tolist()


def differences(a, b, prefix=''):
    if isinstance(a, dict) and isinstance(b, dict):
        return [v for key in sorted(a.keys() | b.keys())
                for v in differences(a.get(key), b.get(key), prefix+'/'+key)]
    return [] if a == b else [prefix]


def assert_readout_match(calibration, dynamics):
    if calibration != dynamics:
        raise ValueError('readout operation mismatch; do not transfer calibration')


def quantum_audit(root):
    import h5py
    base = root/QM
    dirs = {prefix:next((base/'data/non_Markovian_152ns/2026-02-21').glob(prefix+'*'))
            for prefix in ('#6047', '#6048', '#6099')}
    caldir = dirs['#6047']
    calnode = json.loads((caldir/'node.json').read_text())
    caldata = json.loads((caldir/'data.json').read_text())
    calstate = json.loads((caldir/'quam_state/state.json').read_text())
    n = calnode['data']['parameters']['model']['num_runs']
    cal = []
    with h5py.File(caldir/'ds.h5', 'r') as h, np.load(caldir/'arrays.npz', allow_pickle=False) as arrays:
        qubits = [x.decode() for x in h['qubit'][:]]
        if qubits != ['q0','q2'] or n != 5000:
            raise ValueError('calibration dimensions changed')
        for key in ['I_g','Q_g','I_e','Q_e']:
            if h[key].shape != (2,n) or not np.isfinite(h[key][:]).all():
                raise ValueError('invalid IQ samples')
        for i, qubit in enumerate(qubits):
            fit = caldata['results'][qubit]
            counts = []
            for prep in ['g','e']:
                rotated = math.cos(fit['angle'])*h['I_'+prep][i] - math.sin(fit['angle'])*h['Q_'+prep][i]
                one = int(np.count_nonzero(rotated > fit['threshold']))
                counts.append([n-one, one])
            conf = arrays[f'results.{qubit}.confusion_matrix']
            if not np.allclose(np.array(counts)/n, conf, atol=1e-14, rtol=0):
                raise ValueError('IQ classification fails to reproduce saved confusion matrix')
            operation = calstate['qubits'][qubit]['resonator']['operations']['readout']
            residual = abs(fit['threshold']*operation['length']/2**12-operation['threshold'])
            if residual > 1e-15:
                raise ValueError('threshold unit conversion mismatch')
            cal.append(dict(qubit=qubit, inputs=['g','e'], outcomes=[0,1], counts=counts,
                            shots_per_input=n, confusion_matrix=conf.tolist(),
                            threshold_voltage=fit['threshold'], threshold_demod=operation['threshold'],
                            threshold_conversion_residual=residual))
    dynamics = []
    for prefix in ['#6048', '#6099']:
        directory = dirs[prefix]
        node = json.loads((directory/'node.json').read_text())
        data = json.loads((directory/'data.json').read_text())
        state = json.loads((directory/'quam_state/state.json').read_text())
        model = node['data']['parameters']['model']
        with h5py.File(directory/'ds.h5', 'r') as h:
            labels = {k:[v.decode() for v in h[k][:]] for k in ['qubit','initial_state','axis']}
            counts = qpt_counts(h['state'][:], labels, model['num_runs'])
            dimensions = list(h['state'].shape)
        bloch = 1-2*np.array(counts)/model['num_runs']
        saved = [data['results']['q0'][s]['raw']['Bloch vector'] for s in labels['initial_state']]
        if not np.allclose(bloch, saved, atol=1e-14, rtol=0):
            raise ValueError('raw counts do not reproduce saved Bloch vectors')
        a, b = calstate['qubits']['q0'], state['qubits']['q0']
        assert_readout_match(a['resonator']['operations']['readout'], b['resonator']['operations']['readout'])
        gap = (datetime.fromisoformat(node['metadata']['run_start']) -
               datetime.fromisoformat(calnode['metadata']['run_end'])).total_seconds()
        dynamics.append(dict(path=str(directory.relative_to(base)), id=node['id'], parents=node['parents'],
                             run_start=node['metadata']['run_start'], run_end=node['metadata']['run_end'],
                             seconds_after_calibration_end=gap, interaction_ns=model['interaction_time_in_ns'],
                             labels=labels, state_shape=dimensions, terminal_one_counts=counts,
                             shots_per_setting=model['num_runs'], raw_bloch_reproduces=True,
                             readout_operation_equal=True, q0_changed_fields=differences(a,b),
                             intermediate_flag_axis_present=False))
    return dict(calibration_path=str(caldir.relative_to(base)), calibration_id=calnode['id'],
                calibration_run_start=calnode['metadata']['run_start'], calibration_run_end=calnode['metadata']['run_end'],
                readout_calibration=cal, terminal_dynamics=dynamics,
                finite_data_instrument_region_built=False,
                classifier_caveat='Saved classifier reproduces training IQ labels; independent validation not verified.',
                match_status='Adjacent readout calibration and terminal QPT; not a matched intermediate instrument.')


def white_audit(root, entries):
    records = []
    for entry in entries:
        if entry['repository'] != WHITE or not entry['path'].endswith('.pickle'):
            continue
        raw = (root/WHITE/entry['path']).read_bytes()
        ops = list(pickletools.genops(raw))
        record = dict(path=entry['path'], role=entry['role'],
                      opcodes=dict(sorted(Counter(op.name for op, _, _ in ops).items())),
                      executed=False)
        if entry['role'] == 'terminal_probability_vector':
            values = primitive_float_list(raw)
            if not all(0 <= v <= 1 for v in values):
                raise ValueError('not a probability vector')
            record.update(length=len(values), min=min(values), max=max(values),
                          all_on_1_over_1024_grid=all(v*1024 == round(v*1024) for v in values),
                          original_counts_or_denominators_stored=False,
                          setting_or_outcome_keys_stored=False)
        else:
            record['decoded'] = False
        records.append(record)
    return dict(records=records,
                warning='Lattice is not a verified shot denominator. Do not reshape 11999 values or impute zeros without keys.',
                matched_intermediate_calibration_verified=False)


def run(source_root, nmn_dir, download=False):
    if source_root.resolve().is_relative_to(HERE.parents[1]):
        raise ValueError('keep source downloads outside repository')
    manifest = json.loads((HERE/'acquisition-sources.json').read_text())
    for entry in manifest['files']:
        path = source_root/entry['repository']/entry['path']
        if download and not path.exists():
            with urllib.request.urlopen(entry['url'], timeout=60) as response:
                raw = response.read()
            if len(raw) != entry['bytes'] or hashlib.sha256(raw).hexdigest() != entry['sha256']:
                raise ValueError('source download integrity failure')
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(raw)
        audit.verify_file(path, entry)
    return dict(kind='measured_acquisition_schema_and_matching_audit', quantum_memory_certified=False,
                nmn=mapping_audit(nmn_dir), quantum_memory_data=quantum_audit(source_root),
                white=white_audit(source_root, manifest['files']))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root', type=Path, required=True)
    parser.add_argument('--nmn-dir', type=Path, required=True)
    parser.add_argument('--download', action='store_true')
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = json.dumps(run(args.source_root, args.nmn_dir, args.download), indent=2, sort_keys=True)+'\n'
    target = HERE/'results/acquisition-audit.json'
    if args.check:
        if target.read_text() != result:
            raise SystemExit('acquisition audit snapshot differs')
        print('Measured IQ labels, QPT counts, source matches, NMN alternatives and safe schemas reproduce')
    else:
        target.write_text(result)
        print(target)


if __name__ == '__main__':
    main()
