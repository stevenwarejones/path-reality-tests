"""Public CAL science-image audit. Pixel diagnostics are NOT a heating limit."""
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import urllib.request
import zipfile
import numpy as np
from PIL import Image
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
MANIFEST = HERE/'cal-manifest.json'
RESULT = HERE/'results/cal-audit.json'


def encode(x):
    return json.dumps(x, indent=2, sort_keys=True, allow_nan=False)+'\n'


def verify(p, entry):
    raw = p.read_bytes()
    if len(raw) != entry['bytes'] or hashlib.sha256(raw).hexdigest() != entry['sha256']:
        raise ValueError(f'Hash/length mismatch: {p.name}')


def fetch(directory):
    directory.mkdir(parents=True, exist_ok=True)
    for e in json.loads(MANIFEST.read_text())['files']:
        p = directory/e['name']
        if not p.exists():
            # NASA's public resolver returns an expiring transport URL. Only the
            # versioned resolver and content hash belong in the pinned manifest.
            target = urllib.request.urlopen(e['url'], timeout=60).read().decode().strip().strip('"')
            if not target.startswith('https://'):
                raise ValueError('Unexpected NASA resolver result')
            data = urllib.request.urlopen(target, timeout=180).read()
            temporary = p.with_suffix(p.suffix+'.part')
            temporary.write_bytes(data)
            verify(temporary, e)
            temporary.replace(p)
        verify(p, e)


def summarize(directory):
    manifest = json.loads(MANIFEST.read_text())
    for e in manifest['files']:
        verify(directory/e['name'], e)
    archive = next(e for e in manifest['files'] if e['name'].endswith('.zip'))
    shots = []
    # Fixed after visual source discovery, before any field fitting. Pixel axes
    # are deliberately not relabelled as physical x/z without a calibration join.
    rows, cols = np.mgrid[760:980, 1530:1640]
    xx = (cols-1585)/55
    yy = (rows-870)/110
    with zipfile.ZipFile(directory/archive['name']) as z:
        table = list(csv.DictReader(io.StringIO(z.read(next(n for n in z.namelist() if n.endswith('.csv'))).decode())))
        for name in sorted(n for n in z.namelist() if n.endswith('image0.png')):
            suffix = name.split('/L1/')[1].rsplit('/', 1)[0]
            path_matches = [r for r in table if r['new_path'].endswith('/L1/'+suffix)]
            matched = [r for r in path_matches if r['image0_time'] != 'not_applicable' and r['image1_time'] != 'not_applicable']
            if len(matched) != 1:
                raise ValueError('Nonunique science-image timestamp join')
            a, b = [np.array(Image.open(io.BytesIO(z.read(n))), dtype=float)
                    for n in (name, name.replace('image0.png', 'image1.png'))]
            if a.shape != (2048, 2048) or b.shape != a.shape:
                raise ValueError('Unexpected camera dimensions')
            aa, bb = a[760:980,1530:1640], b[760:980,1530:1640]
            if min(aa.min(),bb.min()) <= 0:
                raise ValueError('ROI has nonpositive intensity; no silent log floor')
            od = np.log(bb/aa)
            # Descriptive single-Gaussian + planar fringe background. It is not
            # the published two-spin TF/thermal PSF model and has no fit CI.
            def residual(p):
                amp,cx,cy,sx,sy,b0,bx,by=p
                model=amp*np.exp(-.5*((cols-cx)/sx)**2-.5*((rows-cy)/sy)**2)+b0+bx*xx+by*yy
                return (model-od).ravel()
            fit=least_squares(residual,[.2,1583,900,8,20,0,0,0],
                              bounds=([0,1550,800,2,2,-1,-1,-1],[2,1610,950,50,90,1,1,1]),
                              ftol=1e-9,xtol=1e-9,gtol=1e-9,max_nfev=200)
            p=fit.x
            shots.append({'image':name,'nominal_TOF_from_path_ms':float(name.split('/')[-2].split('T')[0].replace('p','.')),
                          'path_rows_without_image_times':len(path_matches)-len(matched),
                          'image0_time_label':matched[0]['image0_time'],'image1_time_label':matched[0]['image1_time'],
                          'roi_intensity_min':float(min(aa.min(),bb.min())),
                          'roi_saturated_pixels':int(np.count_nonzero((aa>=65535)|(bb>=65535))),
                          'log_ratio_min':float(od.min()),'log_ratio_max':float(od.max()),
                          'gaussian_diagnostic':[float(v) for v in p],
                          'fit_success':bool(fit.success),'residual_rms':float(np.sqrt(np.mean(fit.fun**2)))})
    return {'status':'science arrays and timestamp labels joined; physical response calibration NOT established',
            'source':archive,'shots':shots,'images_loaded':2*len(shots),'roi_rows':[760,980],'roi_columns':[1530,1640],
            'diagnostic_parameter_order':['amplitude','column_center_px','row_center_px','column_sigma_px','row_sigma_px','background','column_gradient','row_gradient'],
            'time_scale':'YYYY-day-of-year timestamp labels; clock scale/offset and release-to-image timing unverified',
            'published_52_pK_reproduced':False,'image_calibration_obtained':False,'matched_ISS_trajectory_obtained':False,
            'missing':['object-plane length per pixel and uncertainty for this 2019 sequence',
                       'dark/flat/saturation calibration and PSF with uncertainty',
                       'shot selection, spin components, effective lens waveform and initial state',
                       'release/image time definition and matching ISS position/attitude with errors'],
            'scope':'Gaussian fits are diagnostic pixel summaries only; no physical radii, temperature, field fitting, confidence or modulation claim'}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data-dir',type=Path);p.add_argument('--fetch',action='store_true');p.add_argument('--check',action='store_true')
    a=p.parse_args()
    if a.fetch:
        if a.data_dir is None:p.error('--fetch requires --data-dir')
        fetch(a.data_dir)
    if a.data_dir is None:
        if not a.check:p.error('Writing requires actual source files via --data-dir')
        result=json.loads(RESULT.read_text())
        assert result['source'] in json.loads(MANIFEST.read_text())['files']
        assert len(result['shots'])==9 and result['images_loaded']==18
        assert not any(result[k] for k in ['published_52_pK_reproduced','image_calibration_obtained','matched_ISS_trajectory_obtained'])
        print('CAL offline audit consistency checked; no image re-read')
        return
    result=summarize(a.data_dir)
    if a.check:
        old=json.loads(RESULT.read_text())
        # Fit solvers/BLAS vary slightly across platforms; structural values exact,
        # descriptive fitted pixels within 0.002 px/parameter and residual 1e-7.
        for s,t in zip(result['shots'],old['shots']):
            if not np.allclose(s.pop('gaussian_diagnostic'),t.pop('gaussian_diagnostic'),atol=.002,rtol=1e-5):
                raise AssertionError('Diagnostic fit changed')
            for key in ['log_ratio_min','log_ratio_max']:
                if not np.isclose(s.pop(key),t.pop(key),atol=1e-12,rtol=1e-12):
                    raise AssertionError('Pixel log ratio changed')
            if not np.isclose(s.pop('residual_rms'),t.pop('residual_rms'),atol=1e-7,rtol=1e-5):
                raise AssertionError('Diagnostic residual changed')
        if result != old:raise AssertionError('CAL source summary changed')
        print('CAL full-source image, timestamp and diagnostic checks passed')
    else:
        RESULT.write_text(encode(result))
        print('Wrote diagnostic CAL source audit; no calibrated heating bound')

if __name__=='__main__':main()
