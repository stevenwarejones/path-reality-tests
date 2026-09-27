"""Execute the feasibility boundary; --check never updates committed artifacts."""
import argparse
import json
from pathlib import Path
import sys
import numpy as np
from scipy.optimize import brentq
from noise_kernels import (sphere_spatial_difference,sphere_difference_lower,band_held_lower,
    quantum_tv_bound,harmonic_visibility_bound,separated_band_gain_upper,line_periodogram_gain,
    mechanical_gain,held_kernel)
from noise_sources import HERE,BASE,summarize,encode,digest
sys.path.insert(0,str(BASE))
from moving_response import assembly_psd_upper_per_a


def calculate(source,protocol):
    b=protocol['benchmark'];T=b['duration_s'];rc=protocol['rc_m']
    G=sphere_spatial_difference(b['radius_m'],b['density_kg_m3'],b['separation_m'],rc)
    lower=sphere_difference_lower(b['radius_m'],b['density_kg_m3'],b['separation_m'],rc)
    time_lower=band_held_lower(1.,T)
    variance=b['target_exponent']/(lower*T*T*time_lower)
    mass=protocol['quantum_domain']['maximum_mass_u'];exposure=protocol['quantum_domain']['maximum_exposure_s']
    def sodium_ratio(sigma2,tail=0.,gain=1.):
        eps=quantum_tv_bound(sigma2,mass,exposure,tail)
        return max(harmonic_visibility_bound(eps,s['optical_mean_probability']*gain,
                        s['optical_visibility'],s['ols_row_l1'])/s['visibility_se']
                   for s in source['sodium']['scans'])
    max_sodium=brentq(lambda factor:sodium_ratio(variance*factor)-1.,1e-4,1e5,xtol=1e-10)*variance
    tail=brentq(lambda eta:sodium_ratio(variance,eta)-1.,0.,.01,xtol=1e-14)
    assembly=assembly_psd_upper_per_a(rc,1.)
    q0=assembly['q0_upper'] # covariance per sigma^2, independent of artificial speed 1
    sensor=source['sensor'];n=protocol['sensor']['samples'];dt=protocol['sensor']['sample_interval_s']
    fmin=min(x['minimum_hz'] for x in sensor['spectra']);fmax=max(x['maximum_hz'] for x in sensor['spectra'])
    leakage=separated_band_gain_upper(fmin,fmax,1.,n,dt)
    f0=sensor['workbook_cells']['N20']; k=sensor['workbook_cells']['L20']
    # For 0 <= nu <= 1, |chi(nu)/chi(0)|² <= [1-(1/f0)^2]^-2, all Q.
    mechanical_upper=(1-(1/f0)**2)**-2
    calibration=sensor['dc_flux_force_calibration_phi0_squared_per_N2']
    minimum_psd=min(x['minimum_psd_phi0_squared_per_hz'] for x in sensor['spectra'])
    gain=protocol['sensor']['calibration_multiplier_stress']
    projection_gain=2*n*dt  # any fixed orthogonal sample-space detrending
    sensor_per_variance=gain*q0*calibration*mechanical_upper*projection_gain
    max_sensor=minimum_psd/sensor_per_variance
    # Off-grid lines: evaluated exactly, never used to certify the continuum.
    atom_results=[]
    for nu in [.1,.137123456789,.5,.999999999,1.]:
        line_gain=max(float(line_periodogram_gain(fmin,nu,n,dt,s)) for s in [False,True])
        atom_results.append({'frequency_hz':nu,'variance_s_minus2':variance,
                             'benchmark_exponent_numerical':G*variance*T*T*float(held_kernel(nu,T)),
                             'force_variance_upper_N2':variance*q0,
                             'raw_sensor_psd_upper_at_min_frequency':gain*variance*q0*calibration*mechanical_upper*line_gain})
    # A continuum witness is the uniform variance density sigma²/.9 on [.1,1].
    # Use deterministic quadrature only for its descriptive exponent.
    nodes,weights=np.polynomial.legendre.leggauss(64)
    uniform_time=float(np.dot(weights,held_kernel(.55+.45*nodes,T))/2)
    margins=[]
    for mscale,tscale,transmission,cal in [(1,1,1,1),(1,1,.5,10),(2,1,1,10),(1,2,1,10)]:
        ep=quantum_tv_bound(variance,mass*mscale,exposure*tscale)
        ratio=max(harmonic_visibility_bound(ep,s['optical_mean_probability']*transmission,
                          s['optical_visibility'],s['ols_row_l1'])/s['visibility_se'] for s in source['sodium']['scans'])
        margins.append({'mass_multiplier':mscale,'exposure_multiplier':tscale,'transmission_multiplier':transmission,
                        'additional_sensor_calibration_multiplier':cal,'visibility_SE_ratio_upper':ratio,
                        'sensor_minimum_observed_psd_fraction_upper':variance/max_sensor*cal})
    return {'status':'Gate A blocked for complementary empirical identification; conditional bounds are implemented',
            'claim_flags':{'empirical_joint_survivor':False,'exclusion':False,'genuine_identification_gain':False,
                           'conditional_ensemble_small_perturbation_family':True,'continuum_temporal_band_covered':True,
                           'arbitrary_frames_covered':False,'outcome_selection':False},
            'benchmark':{'spatial_difference_numerical':G,'spatial_difference_analytic_lower':lower,
                         'temporal_band_analytic_lower':time_lower,'variance_s_minus2':variance,
                         'certified_exponent_lower':lower*variance*T*T*time_lower,
                         'uniform_band_exponent_numerical':G*variance*T*T*uniform_time,
                         'analytic_spatial_conservatism':G/lower},
            'sodium':{'unconditional_probability_change_upper':quantum_tv_bound(variance,mass,exposure),
                      'maximum_visibility_change_over_SE_upper':sodium_ratio(variance),
                      'maximum_unbounded_input_tail_for_one_SE_certificate':tail,
                      'conditional_one_SE_variance_ceiling':max_sodium,
                      'ensemble_bound_for_prior_SSB_variance':quantum_tv_bound(35.6535/1e6/2,mass,exposure),
                      'single_realization_95pct_trace_distance_bound':min(1.,np.sqrt(variance*mass**2*exposure**2/.05)),
                      'ensemble_is_not_single_realization_or_joint_likelihood':True},
            'sensor':{'complete_assembly_q0_upper':q0,'discrete_window_uniform_gain_upper_s':leakage,
                      'dc_flux_force_calibration':calibration,'mechanical_low_band_gain_upper':mechanical_upper,
                      'fixed_orthogonal_detrend_gain_upper_s':projection_gain,
                      'undetrended_leakage_fraction_upper':variance*gain*q0*calibration*mechanical_upper*leakage/minimum_psd,
                      'minimum_measured_psd':minimum_psd,'max_predicted_psd_upper':variance*sensor_per_variance,
                      'fraction_of_minimum_measured_psd_upper':variance/max_sensor,
                      'conditional_measured_total_PSD_variance_ceiling':max_sensor,
                      'field_induced_displacement_rms_upper_m':np.sqrt(variance*q0)/k,
                      'same_field_backaction_ratio_to_rc':np.sqrt(variance*q0)/k/rc,
                      'calibration_stress_multiplier':gain},
            'source_removal':{'sodium_only_certificate_variance':max_sodium,
                              'sensor_only_certificate_variance':max_sensor,
                              'both_certificate_variance':min(max_sodium,max_sensor),
                              'sensor_addition_gain_in_selected_band':max_sodium/min(max_sodium,max_sensor),
                              'actual_admitted_by_one_rejected_by_joint_spectrum':None,
                              'reason':'These are sufficient certificate domains, not empirical acceptance sets. Sensor is redundant in selected band; minimum complementarity bar is NOT met.'},
            'off_grid_atomic_witnesses':atom_results,'uniform_witness_density_s_minus2_per_hz':variance/.9,
            'sensitivity':margins,
            'required_next_input':{'sodium':'validated exposure/mass tail and preparation independence; calibrated joint count/flux likelihood or exact controlled moving response',
                                  'sensor':'matched calibrated low-frequency force or another measured geometry, including control response and mode',
                                  'LPF':'matched calibrated angular science/control data and reconstructed trajectory error envelope',
                                  'CAL':'public PSI-32 images and timestamp labels obtained; object-plane scale/PSF, effective lens/spin fit, clock and ISS trajectory error join unestablished; 52 pK not reproduced'}}


def report(r):
    b,s,m=r['benchmark'],r['sodium'],r['sensor'];rem=r['source_removal']
    return f'''# Moving-noise feasibility boundary

**Gate A is blocked for the requested complementary empirical result.** This is
an implemented response-bound diagnostic, not completion of that research mission.
The selected slow terrestrial field band gives no measured-source identification
gain when the sensor is added. No large scan or joint confidence claim follows.

| Executed quantity | Result | Status |
|---|---:|---|
| Shared field variance | {b['variance_s_minus2']:.8g} s^-2 | Constructed |
| Held-sphere dephasing lower bound | {b['certified_exponent_lower']:.8g} | Analytic continuum certificate |
| Uniform-band dephasing | {b['uniform_band_exponent_numerical']:.8g} | Numerical evaluation |
| Sodium probability perturbation | ≤ {s['unconditional_probability_change_upper']:.8g} | Conditional exact quantum bound |
| Largest visibility perturbation / measured SE | ≤ {s['maximum_visibility_change_over_SE_upper']:.8g} | Descriptive sensitivity, not coverage |
| Sensor response / minimum measured PSD | ≤ {m['fraction_of_minimum_measured_psd_upper']:.8g} | Conditional assembly/calibration bound, 10× gain stress |
| Field-induced displacement / correlation length | ≤ {m['same_field_backaction_ratio_to_rc']:.8g} | Linear-mode self-consistency diagnostic |
| Unbounded input tail allowed before 1-SE certificate fails | {s['maximum_unbounded_input_tail_for_one_SE_certificate']:.8g} | Needed bound; not measured |
| Prior SSB witness under new quantum bound | {s['ensemble_bound_for_prior_SSB_variance']:.1f} | Vacuous; prior gap remains |

The field law is stationary in one Earth-following rigid coordinate system, with
Gaussian spatial covariance at 100 nm and any positive temporal variance measure
on 0.1–1 Hz. All atoms and between-grid spectra in this declared band are covered;
other frequencies and inertial frame laws are not. The benchmark is recomputed
in that same law. This does not reinterpret PR #12's SSB field as terrestrial.

## What the measured sources do and do not add

All 95 sodium scans (3,895 count bins) and all 14 force-sensor spectra are loaded
from hash-verified originals in the full-source route. The deposited fringe fits
and thermal-regression calibration are reproduced independently. The sensor bound
uses the actual discrete Blackman sample count, both sidebands and both window
conventions. The primary sensor bound additionally covers any fixed orthogonal
sample-space detrending, at a quantified sensitivity cost; no frequency grid
supplies the uniform certificate.

| Source subset | Variance ceiling of sufficient certificate (s^-2) |
|---|---:|
| Sodium | {rem['sodium_only_certificate_variance']:.8g} |
| Sensor | {rem['sensor_only_certificate_variance']:.8g} |
| Both | {rem['both_certificate_variance']:.8g} |

The gain is exactly {rem['sensor_addition_gain_in_selected_band']:.1f}. Leaving this
band or exceeding either ceiling is **not exclusion**. An upper-response
certificate cannot establish a rejected spectrum. Thus these results do not meet
the mission's minimum complementarity requirement. They locate a useful boundary:
a full quantum norm estimate can avoid the sodium Markov approximation, but in
the class where it is informative the existing sensor band adds no identification.

## Why this is not an empirical survivor

The mass/exposure support, independent preparation and apparatus response are
conditional inputs. The archived velocity model does not measure its tails;
preparation history and correlated repeated events are not a calibrated joint
likelihood. A 95% single-realization trace-distance estimate is only
{s['single_realization_95pct_trace_distance_bound']:.5g}, much weaker than the
ensemble estimate. Existing optical residuals remain. No simulated recovery is
misrepresented as recovery through an unavailable calibrated experimental chain.

See [theory](../theory.md), [source/candidate audit](../sources.md),
[prior art](../prior-art.md), and [claim status](../README.md).
Numerical records, nuisance stresses, off-grid spectra, and required next inputs
are in [boundary.json](boundary.json). Original PR #12 files remain unchanged.
'''


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data-dir',type=Path);p.add_argument('--check',action='store_true')
    args=p.parse_args();out=HERE/'results';out.mkdir(exist_ok=True)
    path=out/'sources.json'
    source=summarize(args.data_dir) if args.data_dir else json.loads(path.read_text())
    # Canonicalize before analysis: offline/full-source use identical inputs.
    source=json.loads(encode(source))
    protocol=json.loads((HERE/'protocol.json').read_text())
    result=json.loads(encode(calculate(source,protocol)))
    provenance={str(f.relative_to(HERE.parent.parent)):digest(f) for f in
                [HERE/'protocol.json',HERE/'noise_kernels.py',HERE/'noise_sources.py',HERE/'analyze.py',
                 BASE/'manifest.json',BASE/'results/baselines.json',BASE/'moving_response.py',BASE/'response.py',
                 BASE/'io_data.py',BASE/'baselines.py',HERE/'plot.py',HERE/'cal_audit.py',HERE/'cal_response.py',
                 HERE/'cal-manifest.json',HERE/'results/cal-audit.json']}
    provenance['generated_sources_sha256']=__import__('hashlib').sha256(encode(source).encode()).hexdigest()
    artifacts={'sources.json':encode(source),'boundary.json':encode(result),'report.md':report(result),'provenance.json':encode(provenance)}
    for name,content in artifacts.items():
        target=out/name
        if args.check:
            if target.read_text()!=content: raise SystemExit('Reproduction failed: '+name)
        else: target.write_text(content)
    print('Feasibility boundary reproduces ('+('full-source' if args.data_dir else 'offline')+'); empirical gate remains blocked')

if __name__=='__main__':main()
