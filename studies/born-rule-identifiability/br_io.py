"""Audited HDF5 schema, mask grouping and whole-cycle source exclusions."""
from datetime import datetime
from pathlib import Path
import h5py
import numpy as np
from br_sources import verify

EXCLUSIONS = {
    23: sorted(set([204,222,226,*range(1,81),45,65,156,226,333,366,121,227,341,46,63,252,335,376])),
    30: sorted(set([186,371,47,67,149,271,335,393,272,444,446,40,424,439,301,410,51]))}


def epoch_rows(a):
    origin = datetime(*a[0, :5].astype(int))
    return np.array([(datetime(*row[:5].astype(int))-origin).total_seconds()+row[5] for row in a])


def read(cache, entry):
    path = verify(Path(cache)/entry['name'], entry)
    with h5py.File(path, 'r') as h:
        inventory = {}
        for key in h:
            d = h[key]
            if not isinstance(d, h5py.Dataset):
                raise ValueError('Unexpected HDF5 group')
            arr = d[()]
            inventory[key] = {'shape': list(arr.shape), 'dtype': str(arr.dtype),
                'attributes': {k: str(v) for k,v in d.attrs.items()},
                'nonfinite': int((~np.isfinite(arr)).sum())}
            if key != 'temperature-housing' and not np.isfinite(arr).all():
                raise ValueError('Nonfinite archive signal/time/settings')
        settings = h['shuttercombination'][:].reshape(-1)
        if not np.equal(settings, settings.astype(int)).all():
            raise ValueError('Noninteger shutter encoding')
        settings = settings.astype(int)
        expected = [0,2,4,6] if entry['role'] == 'contrast' else list(range(8))
        width = len(expected)
        if len(settings) % width:
            raise ValueError('Incomplete cycle')
        ss = settings.reshape(-1, width)
        if not np.all(np.sort(ss, axis=1) == expected):
            raise ValueError('Cycle not a permutation of declared settings')
        pd = h['PD signal'][:]; times = h['time'][:]
        if pd.shape != (len(settings), entry['samples_per_setting'], 1):
            raise ValueError('PD shape changed')
        if times.shape != (*pd.shape[:2], 6):
            raise ValueError('Time shape changed')
        seconds = epoch_rows(times.reshape(-1, 6))
        if np.any(np.diff(seconds) <= 0):
            raise ValueError('Timestamps not strictly ordered')
        y = pd.mean(axis=1).reshape(-1, width)
        y = np.take_along_axis(y, np.argsort(ss, axis=1), axis=1)
        summary = {'datasets': inventory, 'cycles': len(y), 'samples_per_setting': pd.shape[1],
            'shutter_codes': expected, 'elapsed_hours': float((seconds[-1]-seconds[0])/3600),
            'cycle_permutations_valid': True, 'timestamps_strictly_increasing': True,
            'distinct_acquisition_orders': len(set(map(tuple, ss.tolist())))}
    return y, summary


def retained(y, temperature, index_base=1, apply_cuts=True):
    ids = np.arange(len(y))+1
    keep = np.ones(len(y), bool) if not apply_cuts else ~np.isin(np.arange(len(y))+index_base, EXCLUSIONS[temperature])
    return y[keep], ids[keep]
