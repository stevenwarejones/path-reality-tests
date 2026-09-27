"""Frozen representation and family constants; no measured associations used."""
import importlib.util
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
AUDIT = HERE.parent/'nist-bell-causal-audit'
PREVIOUS = HERE.parent/'timing-structure'

def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

decoder = load('fine_decoder', AUDIT/'reconstruct.py')
SIDES = ('alice', 'bob')
PEAKS = {'alice':90, 'bob':125}
OFFSETS = {'alice':28, 'bob':37}
PHASE_CELLS = 162  # legal corrected duration 129022..129182 / 800 < 162
CELLS = 3*PHASE_CELLS
LEVELS = ('coarse', 'pulse', 'width8', 'width2', 'width1')
EPSILONS = (0., .001, .01)
ENVELOPES = ('unrestricted', 'event_supported')
LAMBDAS = np.array([.005,.01,.02,.04,.08,.16,.32,.64,1.,2.])
ALPHA = .005  # each of the two families
BLOCKS = 16
SEED = 20260927
REPETITIONS = 400

def mappings(side):
    phase = np.arange(CELLS) % PHASE_CELLS
    pulse = np.arange(CELLS) // PHASE_CELLS
    cat = np.where(phase < PEAKS[side], 0, np.where(phase < PEAKS[side]+2, 1, 2))
    keys = [cat[:,None], np.c_[cat,pulse], np.c_[cat,pulse,phase//8],
            np.c_[cat,pulse,phase//2], np.c_[cat,pulse,phase]]
    return [np.unique(k, axis=0, return_inverse=True)[1] for k in keys]

MAPS = {s:mappings(s) for s in SIDES}
CDF_LABELS = 2*4*CELLS
PARTITION_LABELS = sum(4*(1+sum(int(m.max())+1 for m in MAPS[s])) for s in SIDES)

def grouped(counts, mapping):
    counts = np.asarray(counts)
    out = np.zeros(counts.shape[:-1]+(int(mapping.max())+1,), dtype=counts.dtype)
    for i in range(out.shape[-1]):
        out[...,i] = counts[...,mapping == i].sum(axis=-1)
    return out

def json_text(value):
    import json
    return json.dumps(value, sort_keys=True, indent=2, allow_nan=False)+'\n'


def artifact_json(value):
    """One independent case per line; full numerical precision and compact arrays."""
    import json
    lines=['{']
    keys=sorted(value)
    for index,key in enumerate(keys):
        item=value[key]; suffix=',' if index+1<len(keys) else ''
        if isinstance(item,list):
            lines.append('  '+json.dumps(key)+': [')
            lines.extend('    '+json.dumps(v,sort_keys=True,separators=(',',':'),allow_nan=False)+(',' if j+1<len(item) else '') for j,v in enumerate(item))
            lines.append('  ]'+suffix)
        else:
            lines.append('  '+json.dumps(key)+': '+json.dumps(item,sort_keys=True,separators=(',',':'),allow_nan=False)+suffix)
    return '\n'.join(lines+['}',''])
