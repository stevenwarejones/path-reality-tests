"""Auditable conditional reconstruction of the ONE deposited Fig S1 example.

Image grayscales are not assumed to be photon counts. The primary outputs are
unnormalized pointer differences and division by the measured reference modulus.
A fitted complex scale is a descriptive projection, never a calibration estimate.
"""
import hashlib
import json
from pathlib import Path
import numpy as np
from openpyxl import load_workbook

ROOT=Path(__file__).resolve().parent
WEN=ROOT.parent/'wen-2026-propagator'
CHANNELS=('Plus','Minus','R','L','PSI')


def verified_sheet(name):
    manifest=json.loads((WEN/'manifest.json').read_text())
    entry=next(e for e in manifest['files'] if e['name']==name)
    path=WEN/'raw/dataset'/name
    actual=hashlib.sha256(path.read_bytes()).hexdigest()
    if actual!=entry['sha256']:
        raise ValueError(f'source hash mismatch: {name}')
    book=load_workbook(path,read_only=True,data_only=True)
    try:
        rows=list(book.active.values)
    finally:
        book.close()
    return rows,actual


def load_example():
    frames={}; hashes={}; x=y=None
    for ch in CHANNELS:
        rows,digest=verified_sheet(f'FigS1B-{ch}.xlsx')
        xx=np.array([r[0] for r in rows[1:322]],float)
        yy=np.array(rows[0][1:102],float)
        frame=np.array([r[1:102] for r in rows[1:322]],float)
        if frame.shape!=(321,101) or not np.isfinite(frame).all() or (frame<0).any():
            raise ValueError('invalid image')
        if x is not None and (not np.array_equal(x,xx) or not np.array_equal(y,yy)):
            raise ValueError('unaligned image coordinates')
        x,y=xx,yy; frames[ch]=frame; hashes[ch]=digest
    rows,digest=verified_sheet('FigS1C.xlsx'); hashes['slice']=digest
    target_x=np.array([r[3] for r in rows[1:18]],float)
    target=np.array([r[6]+1j*r[4] for r in rows[1:18]])
    se=np.array([[r[7],r[5]] for r in rows[1:18]],float)
    indices=[]
    for value in target_x:
        match=np.flatnonzero(x==value)
        if len(match)!=1:
            raise ValueError('slice coordinate has no unique image row')
        indices.append(int(match[0]))
    return x,y,frames,target_x,np.array(indices),target,se,hashes


def reduce_frames(frames,y,indices,roi):
    masks={'all_y':np.ones(len(y),bool),'central_y':y==0,'strip_y_1':abs(y)<=1}
    if roi not in masks or not masks[roi].any():
        raise ValueError('unknown/empty ROI')
    # Equal-weight means at the exact registered x row; no interpolation or x fit.
    return {ch:frame[indices][:,masks[roi]].mean(axis=1) for ch,frame in frames.items()}


def reconstruct(signals,gains=None,backgrounds=None):
    gains={ch:1. for ch in CHANNELS} if gains is None else gains
    backgrounds={ch:0. for ch in CHANNELS} if backgrounds is None else backgrounds
    if any(gains[ch]<=0 or backgrounds[ch]<0 for ch in CHANNELS):
        raise ValueError('invalid calibration')
    corrected={ch:(np.asarray(signals[ch])-backgrounds[ch])/gains[ch] for ch in CHANNELS}
    if any((v<0).any() for v in corrected.values()) or (corrected['PSI']<=0).any():
        raise ValueError('calibration produces negative signal or zero reference')
    sx=corrected['Plus']-corrected['Minus']; sy=corrected['R']-corrected['L']
    return (sx+1j*sy)/np.sqrt(corrected['PSI'])


def conditional_box(signals,gain_radius=0.,background_max=0.):
    """Guaranteed coordinate box before unknown common scale/reference phase.

The same channel gain/background may be shared across coordinates: these outer
bounds allow more combinations and therefore do not assert joint attainability.
All endpoints are deterministic sensitivity inputs, not confidence limits.
"""
    if not 0<=gain_radius<1 or background_max<0:
        raise ValueError('invalid sensitivity budget')
    lows={ch:(v-background_max)/(1+gain_radius) for ch,v in signals.items()}
    highs={ch:v/(1-gain_radius) for ch,v in signals.items()}
    if any((v<0).any() for v in lows.values()) or (lows['PSI']<=0).any():
        raise ValueError('box includes negative signals or vanishing reference')
    denom=(np.sqrt(lows['PSI']),np.sqrt(highs['PSI']))
    out=[]
    for plus,minus in [('Plus','Minus'),('R','L')]:
        lower=lows[plus]-highs[minus]; upper=highs[plus]-lows[minus]
        corners=np.array([lower/denom[0],lower/denom[1],upper/denom[0],upper/denom[1]])
        out.extend([corners.min(axis=0),corners.max(axis=0)])
    return np.array(out).T


def projective_comparison(candidate,target):
    norm=float(np.vdot(candidate,candidate).real)
    if norm<=0 or np.linalg.norm(target)==0:
        raise ValueError('zero vector')
    scale=np.vdot(candidate,target)/norm
    residual=float(np.linalg.norm(scale*candidate-target)/np.linalg.norm(target))
    return {'fitted_complex_scale':[float(scale.real),float(scale.imag)],
            'relative_l2_residual':residual,
            'interpretation':'descriptive best common complex scale; not an independent calibration or significance test'}


def analyze():
    x,y,frames,xp,indices,target,se,hashes=load_example()
    rois={}
    for roi in ('all_y','central_y','strip_y_1'):
        signals=reduce_frames(frames,y,indices,roi)
        z=reconstruct(signals)
        denom=signals['Plus']+signals['Minus']
        mismatch=(denom-signals['R']-signals['L'])/denom
        cases={}
        for gr,bg in [(0.,0.),(.01,0.),(.01,1.),(.05,1.)]:
            cases[f'gain_{gr:g}_background_{bg:g}']=conditional_box(signals,gr,bg).tolist()
        rois[roi]={'signals':{ch:v.tolist() for ch,v in signals.items()},
                   'reconstructed_before_scale':np.column_stack([z.real,z.imag]).tolist(),
                   'comparison':projective_comparison(z,target),
                   'pair_total_relative_difference':mismatch.tolist(),
                   'conditional_coordinate_boxes':cases}
    # A nonconstant output phase changes K and psi together but leaves all
    # existing local pointer bilinears and reference intensities invariant.
    z=np.array(rois['all_y']['reconstructed_before_scale']); K=z[:,0]+1j*z[:,1]
    psi=np.sqrt(np.array(rois['all_y']['signals']['PSI']))
    theta=np.pi*xp/16
    C=np.conj(psi)*K
    completions={}
    for sign in (-1,1):
        phase=np.exp(1j*sign*theta); kp=phase*K; psip=phase*psi
        np.testing.assert_allclose(np.conj(psip)*kp,C,rtol=1e-13,atol=1e-12)
        np.testing.assert_allclose(abs(psip)**2,abs(psi)**2,rtol=1e-13,atol=1e-12)
        completions[str(sign)]={'reference_phase_radians':(sign*theta).tolist(),
                               'K_before_common_scale':np.column_stack([kp.real,kp.imag]).tolist()}
    return {'kind':'conditional real-data reconstruction; one image set compared with repeat-mean slice',
            'source_hashes':hashes,'image_shape':[321,101],'x_coordinates':xp.tolist(),
            'deposited_complex_slice':np.column_stack([target.real,target.imag]).tolist(),
            'deposited_se_re_im_not_used_as_confidence_region':se.tolist(),
            'roi_results':rois,'phase_gauge_completions':completions,
            'unbounded_without_calibration':['common amplitude scale without lower/upper reference normalization',
              'spatial reference phase without a phase-sensitive reference acquisition'],
            'scope':'ROI, gains, backgrounds, normalization, shared reference and repeated acquisitions are not recovered by fitting the deposited slice.'}
