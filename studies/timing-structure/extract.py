#!/usr/bin/env python3
"""Shared streaming extraction and fixed, halo-safe timing-structure aggregates."""
import argparse
import json
from pathlib import Path
import zipfile

import h5py
import numpy as np
from common import AUDIT, ROOT, decoder, LAGS, GROUPS, FEATURES, HISTORY, CLOCK_THRESHOLD

CACHE_DTYPE = np.dtype([('category', 'u1'), ('detector', 'u1'), ('clock', 'u1')])


def build_cache(hdf_path, archive_path, side, directory, chunk_records=1_000_000):
    directory = Path(directory).resolve()
    if directory.is_relative_to(ROOT):
        raise ValueError('record caches must remain outside the repository')
    directory.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((AUDIT/'manifest.json').read_text())
    for key, path in [('hdf5', hdf_path), (side, archive_path)]:
        spec = manifest['files'][key]
        if Path(path).stat().st_size != spec['bytes'] or decoder.sha256(path) != spec['sha256']:
            raise ValueError(f'{key}: source hash mismatch')
    path = directory/f'{side}-structure.npy'
    base = covered = selected = multiclick = 0
    previous = previous_sync = None
    jumps, tail_info = [], None
    with h5py.File(hdf_path) as h, zipfile.ZipFile(archive_path) as z:
        members = [m for m in z.infolist() if not m.is_dir()]
        if len(members) != 1 or members[0].file_size % 11:
            raise ValueError('unexpected archive layout')
        n = len(h[side+'/settings'])
        data = np.lib.format.open_memmap(path, mode='w+', dtype=CACHE_DTYPE, shape=(n,))
        data[:] = 0
        member = members[0].filename
        offset = decoder.discover_offset(z, member, h[side+'/settings'][:64])
        config = {k: int(h[f'config/{side}/{k}'][()]) for k in ('pk', 'radius', 'bitoffset')}
        expected = h[side+'/badSyncInfo'][[0, 3], :].T.astype('i8').tolist()
        with z.open(member) as stream:
            for a, next_tag, tail in decoder.closed_blocks(stream, chunk_records):
                if tail:
                    tail_info = dict(records=len(a), syncs=int((a['ch'] == 6).sum()))
                    continue
                setting, _, legacy, previous, bad, events = decoder.decode_block(
                    a, next_tag, previous, config, include_events=True)
                tags = a['ttag'][a['ch'] == 6]
                # clock[j] records the interval ending at sync j; analysis uses j=i-1.
                durations = np.diff(np.r_[tags[0] if previous_sync is None else previous_sync, tags])
                previous_sync = int(tags[-1])
                jumps.extend([[base+i-offset, d] for i, d in bad])
                lo, hi = max(base, offset), min(base+len(setting), offset+n)
                if hi > lo:
                    start, stop = lo-offset, hi-offset
                    sl = slice(lo-base, hi-base)
                    if not np.array_equal(setting[sl], h[side+'/settings'][start:stop]):
                        raise ValueError('raw/stored setting mismatch')
                    if not np.array_equal(legacy[sl], h[side+'/clicks'][start:stop]):
                        raise ValueError('raw/stored legacy word mismatch')
                    data['clock'][start:stop] = durations[sl] > CLOCK_THRESHOLD
                    er = events['row']+base-offset
                    overlap = (er >= start) & (er < stop)
                    data['detector'][er[overlap]] = 1
                    eligible = overlap & np.isin(events['pulse']-config['bitoffset'], (4, 5, 6))
                    unique, first, multiplicity = np.unique(er[eligible], return_index=True, return_counts=True)
                    phase = (events['phase'][eligible]-config['pk'])[first]
                    data['category'][unique] = np.where(phase < 0, 1, np.where(phase < 2, 2, 3))
                    selected += int(eligible.sum())
                    multiclick += int((multiplicity > 1).sum())
                    covered += hi-lo
                base += len(setting)
    if covered != n or jumps != expected:
        raise ValueError('incomplete reconciliation')
    data.flush()
    del data
    provenance = dict(schema_version=1, side=side, rows=n, overlap_verified=covered,
        hdf_sha256=manifest['files']['hdf5']['sha256'], archive_sha256=manifest['files'][side]['sha256'],
        cache_sha256=decoder.sha256(path), selected_events=selected, multiclick_rows=multiclick,
        raw_prefix=offset, raw_suffix=base-offset-n, censored_tail=tail_info,
        pulse_bits=[4, 5, 6], phase_radius_cut=False, peak=config['pk'])
    (directory/f'{side}-structure.json').write_text(json.dumps(provenance, indent=2, sort_keys=True)+'\n')
    return provenance


def prior_any(values, window=HISTORY):
    """For halo+core input, sum previous window entries, excluding current entry."""
    values = np.asarray(values)
    if len(values) <= window:
        raise ValueError('insufficient history halo')
    cs = np.r_[np.int64(0), np.cumsum(values, dtype=np.int64)]
    return cs[window:len(values)] > cs[:len(values)-window]


def summarize(local, remote, bad, category, detector, half, membership=None, state_unknown=None):
    """A fixed source block; returned uncertainty uses compatibility, not deletion."""
    m = len(local)
    membership = np.ones(m, dtype=bool) if membership is None else membership
    state_unknown = np.zeros(m, dtype=bool) if state_unknown is None else state_unknown
    bad = bad | ~np.isin(local, (1, 2)) | ~np.isin(remote, (1, 2)) | state_unknown
    trusted = ~bad & membership
    counts = np.zeros((2, 2, 2, 4), dtype=np.int64)
    exposures = np.zeros((2, 2, 2), dtype=np.int64)
    unknown = np.zeros((2, 2, 2, 2), dtype=np.int64)
    future_unknown = np.zeros((2, 2, 2), dtype=np.int64)
    if trusted.any():
        index = 16*half[trusted].astype('i8')+8*(remote[trusted]-1)+4*(local[trusted]-1)+category[trusted]
        histogram = np.bincount(index, minlength=32).reshape(2, 2, 2, 4)
        exposures[:] = histogram.sum(axis=-1)
        counts[..., 0] = histogram[..., 1:].sum(axis=-1)
        counts[..., 1:] = histogram[..., 1:]
    for h in (0, 1):
        for y in (0, 1):
            compatible = ((half == h) & bad & (membership | state_unknown) &
                          ((local == y+1) | ~np.isin(local, (1, 2))))
            # Counts each physical row once, even if its remote label is ambiguous.
            future_unknown[h, y] = [compatible.sum(), (compatible & detector).sum()]
            for x in (0, 1):
                context = compatible & ((remote == x+1) | ~np.isin(remote, (1, 2)))
                unknown[h, x, y] = [context.sum(), (context & detector).sum()]
    return counts, exposures, unknown, future_unknown


def aggregate(hdf_path, directory, side, chunk_rows=1_000_000):
    if not isinstance(chunk_rows, int) or chunk_rows < 1:
        raise ValueError('chunk_rows must be positive')
    directory = Path(directory)
    provenance = json.loads((directory/f'{side}-structure.json').read_text())
    path = directory/f'{side}-structure.npy'
    if decoder.sha256(path) != provenance['cache_sha256'] or decoder.sha256(hdf_path) != provenance['hdf_sha256']:
        raise ValueError('cache/HDF hash mismatch')
    data = np.load(path, mmap_mode='r')
    n = provenance['rows']
    if data.dtype != CACHE_DTYPE or data.shape != (n,) or n <= 2*HISTORY:
        raise ValueError('invalid cache shape/type')
    other = 'alice' if side == 'bob' else 'bob'
    groups = [np.zeros((len(GROUPS),)+shape, dtype=np.int64) for shape in
              [(2, 2, 2, 4), (2, 2, 2), (2, 2, 2, 2), (2, 2, 2)]]
    lags = [np.zeros((len(LAGS),)+shape, dtype=np.int64) for shape in
            [(2, 2, 2, 4), (2, 2, 2), (2, 2, 2, 2), (2, 2, 2)]]
    group_membership = np.zeros((len(GROUPS), 2), dtype=np.int64)
    group_unknown = np.zeros((len(GROUPS), 2), dtype=np.int64)
    with h5py.File(hdf_path) as h:
        ranges = []
        for station in (side, other):
            ranges += decoder.excursion_ranges(h[station+'/badSyncInfo'][[0, 3], :].T.astype('i8').tolist(), n)
        for start in range(HISTORY, n-HISTORY, chunk_rows):
            stop = min(n-HISTORY, start+chunk_rows)
            rows = np.arange(start, stop, dtype=np.int64)
            halo_rows = np.arange(start-HISTORY, stop+HISTORY, dtype=np.int64)
            local_halo = h[side+'/settings'][start-HISTORY:stop+HISTORY]
            remote_halo = h[other+'/settings'][start-HISTORY:stop+HISTORY]
            clock_bad = np.zeros(len(halo_rows), dtype=bool)
            for left, right in ranges:
                clock_bad |= (halo_rows >= left) & (halo_rows < right)
            bad_halo = clock_bad | ~np.isin(local_halo, (1, 2)) | ~np.isin(remote_halo, (1, 2))
            m = stop-start
            local = local_halo[HISTORY:HISTORY+m]
            remote = remote_halo[HISTORY:HISTORY+m]
            base_bad = bad_halo[HISTORY:HISTORY+m]
            half = (rows >= n//2).astype('u1')
            category = np.asarray(data['category'][start:stop])
            detector = np.asarray(data['detector'][start:stop], dtype=bool)
            recovery = prior_any(data['detector'][start-HISTORY:stop])
            recovery_unknown = prior_any(bad_halo[:HISTORY+m])
            clock = data['clock'][start-1:stop-1].astype(bool)
            clock_unknown = bad_halo[HISTORY-2:HISTORY+m-2] | bad_halo[HISTORY-1:HISTORY+m-1]
            for g in range(len(GROUPS)):
                member = np.ones(m, bool) if g == 0 else ((clock if g < 3 else recovery) == ((g-1) % 2))
                unknown_state = np.zeros(m, bool) if g == 0 else (clock_unknown if g < 3 else recovery_unknown)
                for target, values in zip(groups, summarize(local, remote, base_bad, category, detector, half, member, unknown_state)):
                    target[g] += values
                for hh in (0, 1):
                    group_membership[g, hh] += ((half == hh) & member & ~unknown_state).sum()
                    group_unknown[g, hh] += ((half == hh) & unknown_state).sum()
            for j, lag in enumerate(LAGS):
                sl = slice(HISTORY+lag, HISTORY+lag+m)
                lag_remote = remote_halo[sl]
                bad = base_bad | bad_halo[sl]
                for target, values in zip(lags, summarize(local, lag_remote, bad, category, detector, half)):
                    target[j] += values
    if not np.array_equal(groups[0][0], lags[0][LAGS.index(0)]):
        raise ValueError('pooled and lag-zero counts disagree')
    def pack(arrays):
        return dict(zip(('counts', 'exposures', 'unknown', 'future_unknown'), [a.tolist() for a in arrays]))
    return dict(schema_version=1, provenance=provenance, interior_rows=n-2*HISTORY,
        boundary_rows=2*HISTORY, half_rows=[n//2-HISTORY, n-n//2-HISTORY],
        groups=pack(groups), lags=pack(lags), group_membership=group_membership.tolist(),
        group_unknown=group_unknown.tolist(),
        axes=dict(groups=list(GROUPS), lags=list(LAGS), features=list(FEATURES),
                  counts=['group_or_lag', 'half', 'remote_setting', 'local_setting', 'feature'],
                  unknown=['group_or_lag', 'half', 'remote_setting', 'local_setting', 'envelope'],
                  future_unknown=['lag', 'half', 'local_setting', 'envelope'],
                  envelopes=['unrestricted', 'event_supported']))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--hdf5', type=Path, required=True)
    p.add_argument('--archive', type=Path)
    p.add_argument('--cache-dir', type=Path, required=True)
    p.add_argument('--side', choices=['alice', 'bob'], required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--chunk-records', type=int, default=1_000_000)
    p.add_argument('--chunk-rows', type=int, default=1_000_000)
    p.add_argument('--from-cache', action='store_true')
    p.add_argument('--check', action='store_true')
    a = p.parse_args()
    if not a.from_cache:
        if a.archive is None:
            p.error('--archive required unless --from-cache')
        build_cache(a.hdf5, a.archive, a.side, a.cache_dir, a.chunk_records)
    result = json.dumps(aggregate(a.hdf5, a.cache_dir, a.side, a.chunk_rows), sort_keys=True, separators=(',', ':'))+'\n'
    if a.check:
        if a.output.read_text() != result:
            raise SystemExit('structure aggregates differ')
    else:
        a.output.parent.mkdir(parents=True, exist_ok=True)
        a.output.write_text(result)
    print(a.side+': timing-structure aggregates '+('verified' if a.check else 'written'))


if __name__ == '__main__':
    main()
