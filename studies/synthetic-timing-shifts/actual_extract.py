#!/usr/bin/env python3
"""Reconcile raw records and extract fixed actual-setting categorical counts."""
import argparse
import json
from pathlib import Path
import zipfile

import h5py
import numpy as np
from empirical_extract import AUDIT, decoder


def summarize_rows(local, remote, bad, category, detector, half):
    """Keep uncertain rows as compatible completions, including no-event rows."""
    known = np.zeros((2, 2, 2, 4), dtype=np.int64)
    exposure = np.zeros((2, 2, 2), dtype=np.int64)
    uncertainty = np.zeros((2, 2, 2, 2), dtype=np.int64)
    bad = bad | ~np.isin(local, (1, 2)) | ~np.isin(remote, (1, 2))
    for h in (0, 1):
        for x in (0, 1):
            for y in (0, 1):
                context = (half == h) & (remote == x+1) & (local == y+1)
                trusted = context & ~bad
                exposure[h, x, y] = trusted.sum()
                known[h, x, y, 0] = (trusted & (category > 0)).sum()
                for f in (1, 2, 3):
                    known[h, x, y, f] = (trusted & (category == f)).sum()
                compatible = ((half == h) & bad &
                    ((remote == x+1) | ~np.isin(remote, (1, 2))) &
                    ((local == y+1) | ~np.isin(local, (1, 2))))
                uncertainty[h, x, y] = [compatible.sum(), (compatible & detector).sum()]
    return known, exposure, uncertainty, int(bad.sum())


def extract(hdf_path, archive_path, side, chunk_records=1_000_000):
    manifest = json.loads((AUDIT / 'manifest.json').read_text())
    for key, path in [('hdf5', hdf_path), (side, archive_path)]:
        spec = manifest['files'][key]
        if Path(path).stat().st_size != spec['bytes'] or decoder.sha256(path) != spec['sha256']:
            raise ValueError(f'{key}: source hash mismatch')
    other = 'alice' if side == 'bob' else 'bob'
    counts = np.zeros((2, 2, 2, 4), dtype=np.int64)
    exposures = np.zeros((2, 2, 2), dtype=np.int64)
    unknown = np.zeros((2, 2, 2, 2), dtype=np.int64)
    base = covered = uncertain_rows = selected_events = multiclick_rows = 0
    previous = None
    jumps = []
    tail_info = None
    with h5py.File(hdf_path) as h, zipfile.ZipFile(archive_path) as z:
        members = [m for m in z.infolist() if not m.is_dir()]
        if len(members) != 1 or members[0].file_size % 11:
            raise ValueError('unexpected raw archive layout')
        member = members[0].filename
        n = len(h[side+'/settings'])
        if len(h[other+'/settings']) != n:
            raise ValueError('unequal source lengths')
        offset = decoder.discover_offset(z, member, h[side+'/settings'][:64])
        config = {k: int(h[f'config/{side}/{k}'][()]) for k in ('pk', 'radius', 'bitoffset')}
        expected = h[side+'/badSyncInfo'][[0, 3], :].T.astype('i8').tolist()
        ranges = []
        for station in (side, other):
            ranges += decoder.excursion_ranges(h[station+'/badSyncInfo'][[0, 3], :].T.astype('i8').tolist(), n)
        with z.open(member) as stream:
            for a, next_tag, tail in decoder.closed_blocks(stream, chunk_records):
                if tail:
                    tail_info = dict(records=len(a), syncs=int((a['ch'] == 6).sum()))
                    continue
                setting, _, legacy, previous, bad, events = decoder.decode_block(
                    a, next_tag, previous, config, include_events=True)
                jumps.extend([[base+i-offset, delta] for i, delta in bad])
                lo, hi = max(base, offset), min(base+len(setting), offset+n)
                if hi > lo:
                    start, stop = lo-offset, hi-offset
                    local = h[side+'/settings'][start:stop]
                    remote = h[other+'/settings'][start:stop]
                    if not np.array_equal(setting[lo-base:hi-base], local):
                        raise ValueError('raw/stored settings mismatch')
                    if not np.array_equal(legacy[lo-base:hi-base], h[side+'/clicks'][start:stop]):
                        raise ValueError('raw/stored legacy words mismatch')
                    rows = np.arange(start, stop, dtype=np.int64)
                    clock_bad = np.zeros(len(rows), dtype=bool)
                    for left, right in ranges:
                        clock_bad |= (rows >= left) & (rows < right)
                    er = events['row'] + base-offset
                    overlap = (er >= start) & (er < stop)
                    detector = np.zeros(len(rows), dtype=bool)
                    detector[er[overlap]-start] = True
                    eligible = overlap & np.isin(events['pulse']-config['bitoffset'], (4, 5, 6))
                    selected_events += int(eligible.sum())
                    unique, first, multiplicity = np.unique(er[eligible], return_index=True, return_counts=True)
                    multiclick_rows += int((multiplicity > 1).sum())
                    phase = (events['phase'][eligible]-config['pk'])[first]
                    category = np.zeros(len(rows), dtype=np.uint8)
                    category[unique-start] = np.where(phase < 0, 1, np.where(phase < 2, 2, 3))
                    k, e, u, b = summarize_rows(local, remote, clock_bad, category, detector,
                                               (rows >= n//2).astype('u1'))
                    counts += k
                    exposures += e
                    unknown += u
                    uncertain_rows += b
                    covered += hi-lo
                base += len(setting)
    if covered != n or jumps != expected or exposures.sum()+uncertain_rows != n:
        raise ValueError('source reconciliation or row conservation failed')
    if not np.array_equal(counts[..., 0], counts[..., 1:].sum(axis=-1)):
        raise ValueError('category conservation failed')
    return dict(schema_version=1, side=side, rows=n, half_rows=[n//2, n-n//2],
        hdf_sha256=manifest['files']['hdf5']['sha256'],
        archive_sha256=manifest['files'][side]['sha256'], overlap_verified=covered,
        raw_prefix=offset, raw_suffix=base-offset-n, censored_tail=tail_info,
        source_reconciliation='all local settings, legacy words and timestamp excursions exact',
        pulse_bits=[4, 5, 6], phase_radius_cut=False, phase_peak=config['pk'],
        selected_events_including_uncertain=selected_events, multiclick_rows=multiclick_rows,
        uncertain_rows=uncertain_rows, trusted_exposures=exposures.tolist(),
        known_event_counts=counts.tolist(), compatible_uncertainty=unknown.tolist(),
        axes=dict(known=['half', 'remote_setting', 'local_setting', 'feature'],
                  features=['any', 'early', 'core', 'late'],
                  uncertainty=['half', 'remote_setting', 'local_setting', 'envelope'],
                  envelopes=['unrestricted', 'event_supported']))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--hdf5', type=Path, required=True)
    p.add_argument('--archive', type=Path, required=True)
    p.add_argument('--side', choices=['alice', 'bob'], required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--chunk-records', type=int, default=1_000_000)
    p.add_argument('--check', action='store_true')
    a = p.parse_args()
    result = json.dumps(extract(a.hdf5, a.archive, a.side, a.chunk_records), indent=2, sort_keys=True)+'\n'
    if a.check:
        if a.output.read_text() != result:
            raise SystemExit('actual counts differ')
    else:
        a.output.parent.mkdir(parents=True, exist_ok=True)
        a.output.write_text(result)
    print(a.side+': full-source actual counts '+('verified' if a.check else 'written'))


if __name__ == '__main__':
    main()
