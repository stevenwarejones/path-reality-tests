"""Fixed protocol constants and explicit imports from the audited baseline."""
import importlib.util
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BASE = HERE.parent / 'synthetic-timing-shifts'
AUDIT = HERE.parent / 'nist-bell-causal-audit'


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


decoder = load('structure_decoder', AUDIT / 'reconstruct.py')
LAGS = (-64, -16, -4, -1, 0, 1, 4, 16, 64)
GROUPS = ('pooled', 'clock_0', 'clock_1', 'recovery_0', 'recovery_1')
FEATURES = ('any', 'early', 'core', 'late')
ENVELOPES = ('unrestricted', 'event_supported')
HISTORY = 64
CLOCK_THRESHOLD = 129102


def artifact_json(result):
    """One generated comparison per line, with ordinary readable metadata."""
    import json
    lines = ['{']
    for index, (key, value) in enumerate(sorted(result.items())):
        suffix = ',' if index+1 < len(result) else ''
        if isinstance(value, list):
            lines.append('  '+json.dumps(key)+': [')
            lines.extend('    '+json.dumps(item, sort_keys=True, separators=(',', ':'))+
                         (',' if j+1 < len(value) else '') for j, item in enumerate(value))
            lines.append('  ]'+suffix)
        else:
            encoded = json.dumps(value, indent=2, sort_keys=True).splitlines()
            lines.append('  '+json.dumps(key)+': '+encoded[0])
            lines.extend('  '+line for line in encoded[1:])
            lines[-1] += suffix
    return '\n'.join(lines+['}', ''])
