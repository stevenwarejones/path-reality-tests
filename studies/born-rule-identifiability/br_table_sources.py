"""Read counts by header, retaining jobs and randomized circuit order; never run source code."""
import csv
import hashlib
import io
import json
import re
import zipfile
from pathlib import Path
import numpy as np
from br_sources import verify

HERE = Path(__file__).resolve().parent
DESIGNS = {
    'scan': (7470893, 'wyniki-nairobi.zip', 9, 115, 8, 100000),
    'viviani': (21775462, 'results_nairobi.zip', 5, 9, 15, 100000),
    'aria_drift_control': (21775462, 'results-aria-viv.zip', 5, 20, 5, 1000),
}


def manifest():
    return json.loads((HERE/'table-manifest.json').read_text())


def read_archive(path, nprep, expected_jobs, repetitions, shots):
    with zipfile.ZipFile(path) as z:
        names = [n for n in z.namelist() if not n.startswith('__MACOSX/')]
        jobs = sorted([n for n in names if re.search(r'(testy|tests|results)_\d+\.csv$', n)],
                      key=lambda n: int(re.search(r'_(\d+)\.csv$', n)[1]))
        ids = [int(re.search(r'_(\d+)\.csv$', n)[1]) for n in jobs]
        if len(ids) != expected_jobs or len(set(ids)) != expected_jobs:
            raise ValueError('Unexpected number or duplicate jobs')
        metadata_names = [n for n in names if n.endswith('data.json')]
        metadata = json.loads(z.read(metadata_names[0])) if metadata_names else {}
        jl_names = [n for n in names if 'job_list' in n and n.endswith('.csv')]
        jl = list(csv.DictReader(io.StringIO(z.read(jl_names[0]).decode()))) if jl_names else []
        counts, order_checks, serial = [], [], []
        for name in jobs:
            rows = list(csv.DictReader(io.StringIO(z.read(name).decode('utf-8-sig'))))
            a = np.zeros((nprep, 4, 2), dtype=np.int64)
            reps = np.zeros((nprep, 4), dtype=int)
            order, residuals = [], []
            for row in rows:
                i, j = int(row.get('n', row.get('i'))), int(row['j'])
                if not (0 <= i < nprep and 0 <= j < 4):
                    raise ValueError('Bad setting')
                values = [float(row.get(k, row.get(k+' 0')) or 0) for k in ['0', '1']]
                if any(v < 0 or v != int(v) for v in values) or sum(values) != shots:
                    raise ValueError('Invalid counts or shots')
                a[i, j] += np.array(values, dtype=np.int64)
                reps[i, j] += 1
                order.append([i, j])
                residuals.append(values[0]/shots)
            if not np.all(reps == repetitions):
                raise ValueError('Unbalanced setting repetitions')
            counts.append(a)
            p = a[:, :, 0]/a.sum(2)
            resid = np.array([v-p[i, j] for v, (i, j) in zip(residuals, order)])
            # Circuit-order diagnostic only: the archive has no per-shot timestamps.
            serial.append(float(np.dot(resid[:-1], resid[1:]) / max(np.dot(resid, resid), 1e-30)))
            order_checks.append({'sha256': hashlib.sha256(json.dumps(order).encode()).hexdigest(),
                                 'adjacent_equal_settings': sum(x == y for x, y in zip(order, order[1:])),
                                 'circuit_rows': len(rows)})
        arr = np.array(counts)
        if len(set(x['sha256'] for x in order_checks)) != expected_jobs:
            raise ValueError('Repeated circuit ordering; investigate provenance')
        return {'metadata': metadata, 'actual_jobs': len(jobs), 'job_indices': ids,
                'missing_indices_within_recorded_range': sorted(set(range(max(ids)+1))-set(ids)),
                'job_list_rows': len(jl), 'job_list_complete': len(jl) == len(jobs),
                'job_list_ids': [r['job_id'] for r in jl],
                'setting_order_audit': order_checks, 'lag1_circuit_residual_diagnostic': serial,
                'counts_0_1_by_job': arr.tolist(), 'shots_per_circuit': shots,
                'repetitions_per_setting_job': repetitions,
                'shot_timestamps_available': False,
                'settings_from': 'Each count row, not the incomplete job list',
                'evidence_type': 'hardware_counts'}


def audit(cache):
    for e in manifest()['files']:
        verify(cache/e['cache_name'], e)
    records = {key: read_archive(cache/str(rec)/name, n, j, r, s)
               for key, (rec, name, n, j, r, s) in DESIGNS.items()}
    # The main acquisitions are ten months apart. They share no calibration parameters.
    return {'source': manifest(), 'records': records,
            'selection': 'Nairobi 2022 angle scan + Nairobi 2023 Viviani; Aria is a drift stress control',
            'cross_acquisition_calibration_shared': False}
