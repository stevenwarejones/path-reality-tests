"""Stream NIST's 11-byte compressed records; independently preserve every click.

No coincidence matching, remote-outcome selection, or deletion of interior syncs.
The phase convention reproduces the archived builder (previous sync interval).
See reconstruction.md for the limited meaning of the recovered row alignment.
"""
import hashlib
import json
from pathlib import Path
import zipfile

import h5py
import numpy as np

DTYPE = np.dtype([('ch', 'u1'), ('ttag', '<i8'), ('xfer', '<u2')])
VARIANTS = ('nominal', 'narrow', 'wide')
LASER_PULSES_PER_SYNC = 800


def excursion_ranges(jumps, n):
    """Opposite timestamp offsets cancel across two nominal sync intervals."""
    if len(jumps) % 2:
        raise ValueError('unpaired timestamp jump')
    ranges = []
    for (start, d0), (end, d1) in zip(jumps[::2], jumps[1::2]):
        if not (start < end and d0 < 0 < d1 and d0+d1 == 2*LASER_PULSES_PER_SYNC):
            raise ValueError('unrecognized paired timestamp excursion')
        lo, hi = max(0, start), min(n, end+2)
        if hi > lo:
            ranges.append([lo, hi])
    return ranges



def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def closed_blocks(stream, chunk_records=1_000_000):
    """Yield complete sync intervals, retaining the final open interval as a tail."""
    if not isinstance(chunk_records, int) or chunk_records < 1:
        raise ValueError('chunk_records must be positive')
    carry = np.empty(0, DTYPE)
    first = True
    while True:
        data = stream.read(DTYPE.itemsize*chunk_records)
        if not data:
            break
        if len(data) % DTYPE.itemsize:
            raise ValueError('truncated 11-byte record')
        a = np.concatenate((carry, np.frombuffer(data, DTYPE)))
        sy = np.flatnonzero(a['ch'] == 6)
        if first and len(sy):
            # This audited run starts with a sync. Refuse to silently drop a prefix.
            if sy[0] != 0:
                raise ValueError('records before initial sync need explicit handling')
            first = False
        if len(sy) < 2:
            carry = a
            continue
        stop = int(sy[-1])
        carry = a[stop:].copy()
        yield a[:stop], int(carry['ttag'][0]), False
    if len(carry):
        yield carry, None, True  # censored tail, never fabricated as a no-click trial


def decode_block(a, next_sync_tag, previous_delta, config):
    sy = np.flatnonzero(a['ch'] == 6)
    if not len(sy) or sy[0] != 0:
        raise ValueError('block must begin at a sync')
    tags = np.r_[a['ttag'][sy], np.int64(next_sync_tag)]
    dt = np.diff(tags)
    rounded = np.floor(dt/(129102/800.) + .5)
    bad = [[int(i), int(rounded[i])] for i in np.flatnonzero(rounded != 800)]
    # Match the source's treatment of enormous paired timestamp jumps.
    # All other irregular intervals require a separately specified decoder.
    for j in np.flatnonzero(rounded != 800):
        if abs(rounded[j]) <= 80_000_000:
            raise ValueError('unhandled short/long sync interval')
        replacement = int(dt[j-1]) if j else previous_delta
        if replacement is None or replacement <= 0:
            raise ValueError('jump without a valid preceding interval')
        dt[j] = replacement
    if np.any(dt <= 0):
        raise ValueError('nonpositive corrected sync duration')
    period = dt/np.floor(dt/(129102/800.)+.5)
    initial = period[0] if previous_delta is None else previous_delta/np.floor(previous_delta/(129102/800.)+.5)
    period = np.r_[initial, period[:-1]]
    settings = np.zeros(len(sy), dtype='u1')
    for ch, bit in ((2, 1), (4, 2)):
        idx = np.searchsorted(sy, np.flatnonzero(a['ch'] == ch), side='right')-1
        np.bitwise_or.at(settings, idx, bit)
    det = np.flatnonzero(a['ch'] == 0)
    idx = np.searchsorted(sy, det, side='right')-1
    delay = a['ttag'][det]-tags[idx]
    phases = delay % period[idx]
    pulses = np.floor(delay/period[idx]).astype('i8')
    bits = pulses-config['bitoffset']
    words = {}
    for variant, shift in (('nominal', 0), ('narrow', -1), ('wide', 1)):
        radius = config['radius']+shift
        if radius <= 0:
            raise ValueError('nonpositive phase radius')
        ok = (abs(phases-config['pk']) < radius) & (bits >= 0) & (bits < 16)
        values = (1 << bits[ok]).astype('u2')
        clicks = np.zeros(len(sy), dtype='u2')
        np.bitwise_or.at(clicks, idx[ok], values)
        words[variant] = clicks
        if variant == 'nominal':
            legacy = np.zeros(len(sy), dtype='u2')
            # Intentionally reproduce archived repeated-index behavior, for diagnosis.
            legacy[idx[ok]] += values
    return settings, words, legacy, int(dt[-1]), bad


def discover_offset(archive, member, pattern):
    """Unique 64-setting fingerprint in the first 500,000 records; then full check."""
    with archive.open(member) as stream:
        a = np.frombuffer(stream.read(11*500_000), DTYPE)
    sy = np.flatnonzero(a['ch'] == 6)
    settings = np.zeros(len(sy), 'u1')
    for ch, bit in ((2, 1), (4, 2)):
        idx = np.searchsorted(sy, np.flatnonzero(a['ch'] == ch), side='right')-1
        if np.any(idx < 0):
            raise ValueError('setting before first sync')
        np.bitwise_or.at(settings, idx, bit)
    hits = np.arange(max(0, len(settings)-len(pattern)+1))
    for j, value in enumerate(pattern):
        hits = hits[settings[hits+j] == value]
    if len(hits) != 1:
        raise ValueError('prefix fingerprint not unique')
    return int(hits[0])


def audit_side(hdf_path, archive_path, side, chunk_records=1_000_000):
    patches = {v: [] for v in VARIANTS}
    fixtures = []
    uncertain_detector_rows = []
    with h5py.File(hdf_path) as h, zipfile.ZipFile(archive_path) as z:
        members = [i for i in z.infolist() if not i.is_dir()]
        if len(members) != 1 or members[0].file_size % 11:
            raise ValueError('unexpected archive layout')
        member = members[0].filename
        n = len(h[side+'/settings'])
        offset = discover_offset(z, member, h[side+'/settings'][:64])
        config = {k: int(h[f'config/{side}/{k}'][()]) for k in ('pk','radius','bitoffset')}
        expected_jumps = {who: h[who+'/badSyncInfo'][[0,3],:].T.astype('i8').tolist() for who in ('alice','bob')}
        ranges = sum((excursion_ranges(expected_jumps[who], n) for who in expected_jumps), [])
        base = 0
        previous = None
        mismatches = {'settings': 0, 'legacy_clicks': 0}
        channels = np.zeros(256, dtype='i8')
        jumps = []
        covered = 0
        tail = None
        with z.open(member) as stream:
            for a, next_tag, is_tail in closed_blocks(stream, chunk_records):
                channels += np.bincount(a['ch'], minlength=256)
                if is_tail:
                    tail = {'records': len(a), 'syncs': int((a['ch']==6).sum())}
                    continue
                settings, words, legacy, next_previous, bad = decode_block(a, next_tag, previous, config)
                jumps.extend([[base+i-offset, delta] for i,delta in bad])
                lo, hi = max(base, offset), min(base+len(settings), offset+n)
                if hi > lo:
                    sl = slice(lo-base, hi-base)
                    old_s = h[side+'/settings'][lo-offset:hi-offset]
                    old_c = h[side+'/clicks'][lo-offset:hi-offset]
                    other = 'bob' if side == 'alice' else 'alice'
                    other_s = h[other+'/settings'][lo-offset:hi-offset]
                    uncertain = ~np.isin(old_s, [1,2]) | ~np.isin(other_s, [1,2])
                    for left, right in ranges:
                        uncertain[max(0,left-(lo-offset)):max(0,min(hi-lo,right-(lo-offset)))] = True
                    sy = np.flatnonzero(a['ch'] == 6)
                    present = np.zeros(len(sy), dtype=bool)
                    present[np.searchsorted(sy,np.flatnonzero(a['ch']==0),side='right')-1] = True
                    for word in words.values():
                        if np.any((word != 0) & ~present):
                            raise ValueError('click without recorded detector event')
                    uncertain_detector_rows.extend((np.flatnonzero(uncertain & present[sl])+lo-offset).tolist())
                    mismatches['settings'] += int(np.count_nonzero(settings[sl] != old_s))
                    mismatches['legacy_clicks'] += int(np.count_nonzero(legacy[sl] != old_c))
                    for variant in VARIANTS:
                        got = words[variant][sl]
                        diff = np.flatnonzero(got != old_c)
                        patches[variant].extend([[int(lo-offset+j), int(old_c[j]), int(got[j])] for j in diff])
                    # Small, provenance-located multiclick fixture, including one preceding
                    # sync so the next interval's phase uses a directly observed duration.
                    diff = np.flatnonzero(words['nominal'][sl] != old_c)
                    if len(diff) and not fixtures:
                        row = lo-base+int(diff[0])
                        sy = np.flatnonzero(a['ch']==6)
                        if row > 0 and row+1 < len(sy):
                            events = a[sy[row-1]:sy[row+1]].copy()
                            fixtures.append(dict(side=side, raw_sync_index=base+row,
                                hdf_row=base+row-offset, config=config,
                                records_hex=events.tobytes().hex(), next_sync_tag=int(a['ttag'][sy[row+1]]),
                                nominal_word=int(words['nominal'][row]),
                                legacy_word=int(legacy[row]), settings=int(settings[row])))
                    covered += hi-lo
                previous = next_previous
                base += len(settings)
        if jumps != expected_jumps[side]:
            raise ValueError('raw timestamp jumps disagree with archived support ranges')
        if covered != n or any(mismatches.values()):
            raise ValueError(f'full reconstruction does not reconcile: {covered}/{n}, {mismatches}')
        return dict(side=side, archive_sha256=sha256(archive_path), member=member,
                    member_bytes=members[0].file_size, hdf_rows=n,
                    closed_raw_intervals=base, prefix_intervals=offset,
                    suffix_intervals=base-offset-n, censored_tail=tail,
                    covered_rows=covered, mismatches=mismatches, timestamp_jumps=jumps,
                    uncertain_detector_rows=uncertain_detector_rows,
                    channel_records={str(i):int(v) for i,v in enumerate(channels) if v},
                    phase_config=config, patches=patches, fixtures=fixtures)


def main():
    import argparse
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--hdf5', type=Path, required=True)
    p.add_argument('--alice', type=Path, required=True)
    p.add_argument('--bob', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args=p.parse_args()
    manifest=json.loads((Path(__file__).parent/'manifest.json').read_text())
    for key, path in [('hdf5',args.hdf5),('alice',args.alice),('bob',args.bob)]:
        if sha256(path) != manifest['files'][key]['sha256']:
            raise SystemExit(f'{key}: input checksum mismatch')
    result={side:audit_side(args.hdf5,getattr(args,side),side) for side in ('alice','bob')}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print('Full raw reconstruction reconciled settings and archived update behavior for both sides')


if __name__ == '__main__':
    main()
