#!/usr/bin/env python3
"""Extract a local fixed-background timing cache, outside the git repository."""
import argparse
import importlib.util
import json
from pathlib import Path
import zipfile

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
AUDIT = HERE.parent / 'nist-bell-causal-audit'
SPEC = importlib.util.spec_from_file_location('timing_raw_decoder', AUDIT / 'reconstruct.py')
decoder = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(decoder)


def outside_repository(path):
    path = Path(path).resolve()
    if path.is_relative_to(ROOT):
        raise ValueError('event caches must remain outside the repository')
    return path


def extract(hdf_path, archive_path, output, side='bob', chunk_records=1_000_000):
    output = outside_repository(output)
    manifest = json.loads((AUDIT / 'manifest.json').read_text())
    for key, path in [('hdf5', hdf_path), (side, archive_path)]:
        spec = manifest['files'][key]
        if Path(path).stat().st_size != spec['bytes'] or decoder.sha256(path) != spec['sha256']:
            raise ValueError(f'{key}: source hash mismatch')
    exposures = np.zeros((32, 2), dtype=np.int64)
    row_parts, time_parts, state_parts = [], [], []
    excluded = dict(timestamp_only=0, setting_only=0, both=0)
    base = covered = selected_excluded_events = 0
    previous = None
    jumps = []
    with h5py.File(hdf_path) as h, zipfile.ZipFile(archive_path) as z:
        members = [m for m in z.infolist() if not m.is_dir()]
        if len(members) != 1 or members[0].file_size % 11:
            raise ValueError('unexpected raw archive layout')
        member = members[0].filename
        n = len(h[side + '/settings'])
        offset = decoder.discover_offset(z, member, h[side + '/settings'][:64])
        config = {key: int(h[f'config/{side}/{key}'][()]) for key in ('pk', 'radius', 'bitoffset')}
        expected_jumps = h[side + '/badSyncInfo'][[0, 3], :].T.astype('i8').tolist()
        ranges = decoder.excursion_ranges(expected_jumps, n)
        with z.open(member) as stream:
            for a, next_tag, tail in decoder.closed_blocks(stream, chunk_records):
                if tail:
                    tail_info = {'records': len(a), 'syncs': int((a['ch'] == 6).sum())}
                    continue
                settings, _, legacy, previous, bad, events = decoder.decode_block(
                    a, next_tag, previous, config, include_events=True)
                jumps.extend([[base+i-offset, delta] for i, delta in bad])
                lo, hi = max(base, offset), min(base+len(settings), offset+n)
                if hi > lo:
                    start, stop = lo-offset, hi-offset
                    sl = slice(lo-base, hi-base)
                    old = h[side + '/settings'][start:stop]
                    if not np.array_equal(settings[sl], old):
                        raise ValueError('raw/stored local settings mismatch')
                    if not np.array_equal(legacy[sl], h[side + '/clicks'][start:stop]):
                        raise ValueError('raw/stored archived click-word mismatch')
                    rows = np.arange(start, stop, dtype=np.int64)
                    clock_bad = np.zeros(len(rows), dtype=bool)
                    for left, right in ranges:
                        clock_bad |= (rows >= left) & (rows < right)
                    setting_bad = ~np.isin(old, (1, 2))
                    valid = ~clock_bad & ~setting_bad
                    excluded['timestamp_only'] += int((clock_bad & ~setting_bad).sum())
                    excluded['setting_only'] += int((~clock_bad & setting_bad).sum())
                    excluded['both'] += int((clock_bad & setting_bad).sum())
                    block = rows * 32 // n
                    exposures += np.bincount(2*block[valid] + old[valid]-1,
                                             minlength=64).reshape(32, 2)
                    event_rows = events['row'] + base - offset
                    # Fixed *baseline* pulse membership, with NO phase-radius cut.
                    eligible = ((event_rows >= start) & (event_rows < stop) &
                                np.isin(events['pulse']-config['bitoffset'], (4, 5, 6)))
                    er = event_rows[eligible]
                    good_events = valid[er-start]
                    selected_excluded_events += int((~good_events).sum())
                    row_parts.append(er[good_events])
                    time_parts.append((events['phase'][eligible]-config['pk'])[good_events])
                    state_parts.append((old[er-start][good_events]-1).astype('u1'))
                    covered += hi-lo
                base += len(settings)
    if covered != n or jumps != expected_jumps:
        raise ValueError('incomplete overlap or timestamp excursion mismatch')
    rows = np.concatenate(row_parts)
    times = np.concatenate(time_parts)
    states = np.concatenate(state_parts)
    # Canonical record order. All selected events are retained, including multiclicks.
    order = np.argsort(rows, kind='stable')
    rows, times, states = rows[order], times[order], states[order]
    unique, starts, multiplicity = np.unique(rows, return_index=True, return_counts=True)
    event_exposures = np.bincount(2*(unique*32//n) + states[starts], minlength=64).reshape(32, 2)
    if np.any(event_exposures > exposures) or exposures.sum()+sum(excluded.values()) != n:
        raise ValueError('trial conservation failed')
    summary = dict(schema_version=1, data_kind='real_local_background_no_remote_comparison',
                   side=side, hdf_sha256=manifest['files']['hdf5']['sha256'],
                   archive_sha256=manifest['files'][side]['sha256'],
                   archive_member=member, rows=n, overlap_verified=covered,
                   raw_prefix=offset, raw_suffix=base-offset-n, censored_tail=tail_info,
                   exclusions=excluded, retained_trials=int(exposures.sum()),
                   selected_events=int(len(rows)), event_trials=int(len(unique)),
                   multiclick_trials=int((multiplicity > 1).sum()),
                   max_selected_events_per_trial=int(multiplicity.max(initial=0)),
                   selected_events_in_excluded_rows=selected_excluded_events,
                   pulse_bits=[4, 5, 6], phase_radius_cut=False, phase_peak=config['pk'],
                   no_click_trials=int(exposures.sum()-len(unique)),
                   feature='first event in original record order among baseline-selected pulses',
                   background_by_segment=[])
    for name, mask in [('training', unique*32//n < 8), ('test', unique*32//n >= 8)]:
        for state in (0, 1):
            t = times[starts][mask & (states[starts] == state)]
            segment = slice(0, 8) if name == 'training' else slice(8, 32)
            summary['background_by_segment'].append(dict(segment=name, receiver_setting=state,
                trials=int(exposures[segment, state].sum()), event_trials=int(len(t)),
                first_event_phase_quantiles_bins=np.quantile(t, [.01, .25, .5, .75, .99]).tolist()))
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output, event_row=rows, event_phase=times, event_state=states,
                        exposures=exposures, source_rows=np.array(n),
                        provenance=np.array(json.dumps(summary, sort_keys=True)))
    print(json.dumps(summary, indent=2))
    return summary


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--hdf5', type=Path, required=True)
    p.add_argument('--archive', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--side', choices=['alice', 'bob'], default='bob')
    p.add_argument('--chunk-records', type=int, default=1_000_000)
    a = p.parse_args()
    extract(a.hdf5, a.archive, a.output, a.side, a.chunk_records)


if __name__ == '__main__':
    main()
