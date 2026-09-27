"""Joint dephasing feasibility with trajectories, finite windows and assembly bounds."""
import argparse
import hashlib
import json
from pathlib import Path
import tempfile
import numpy as np
from scipy import constants as co
from scipy.integrate import quad
from response import M0,L,velocity_nodes,optical_coefficients,triangular_white_coefficient,layered_force_coefficient
from combinations import cube_response
from io_data import HERE,verify_sources
from moving_response import (moving_sphere_exponent_per_a,assembly_psd_upper_per_a,
    moving_psd_upper_per_a,lateral_ou_factor,lateral_frozen_factor,blackman_dc_gain,
    periodogram_weights,frozen_periodogram_per_a,relative_pair_factor,triangle_path_exponent_per_a)
from motion_sources import source_summary
from moving_plot import generate as plot


def sodium_forward(source,baseline,a,wind_min,mean=158.,sigma=9.):
    """Nonlinear signed-harmonic interval under the local advected Markov map.

    Unknown wind direction: exponent lies between 0 and the unprojected
    separation response at the slowest admitted relative speed. This is a box
    outer bound, not independent frame choices or fitted resolved fringes.
    """
    m=np.array(source['mass_quadrature_u']);mw=np.array(source['normalized_mass_weights'])
    v,vw=velocity_nodes(mean=mean,sigma=sigma);w=mw[:,None]*vw
    speed=wind_min-v
    if np.any(speed<=0):raise ValueError('Relative speed interval contains zero; local advection bound unavailable')
    dmax=a*np.sqrt(np.pi)*1e-7/speed[None,:]*triangular_white_coefficient(m[:,None],v[None,:],1e-7)
    attenuation=np.exp(-dmax)
    rows=[]
    observed={x['scan']:x for x in baseline['sodium']['scans']}
    for name,power in sorted(source['powers_by_scan'].items(),key=lambda item:int(item[0].rsplit('_',1)[1])):
        s0,s1=optical_coefficients(m[:,None],v[None,:],power)
        c=w*s1;den=np.sum(w*s0)
        low=np.sum(np.where(c>=0,c*attenuation,c));high=np.sum(np.where(c>=0,c,c*attenuation))
        vlo=0. if low<=0<=high else 2*min(abs(low),abs(high))/den
        vhi=2*max(abs(low),abs(high))/den
        original=2*abs(np.sum(c))/den
        delta=max(abs(vlo-original),abs(vhi-original))
        rows.append({'scan':name,'baseline_visibility':float(original),'moving_visibility_interval':[float(vlo),float(vhi)],
                     'max_change_over_measured_visibility_se':float(delta/observed[name]['visibility_se_gaussian']),
                     'absolute_harmonic_condition_number':float(np.sum(abs(c))/max(abs(np.sum(c)),1e-300))})
    # A simple optical timescale is the narrowest grating waist / largest beam speed.
    tc=1e-7/(wind_min-max(v)); massmin=min(m)*M0
    return {'rows':rows,'max_exponent':float(dmax.max()),
            'max_visibility_change_over_se':max(r['max_change_over_measured_visibility_se'] for r in rows),
            'validity':{'spatial_memory_s':tc,'memory_over_optical_transit':tc/(540e-6/max(v)),
                        'recoil_parameter':co.hbar*tc/(massmin*(1e-7)**2),
                        'max_cluster_radius_over_rc':(3*max(m)*M0/(4*np.pi*971))**(1/3)/1e-7,
                        'scope':'local advected Markov, paraxial, point cluster; no exact long-memory Talbot-Lau map claimed'}}


def window_gain(f,f0,q,segment_s):
    """Continuous Blackman approximation to the very finely sampled periodogram.

    Computes white-force mechanical convolution / point-frequency response.
    N=2^22 at 100 kHz makes sampling corrections negligible in the narrow band;
    symmetric-versus-periodic DC leakage is evaluated separately exactly.
    """
    coeff={0:.42,-1:-.25,1:-.25,-2:.04,2:.04}
    norm=sum(c*c for c in coeff.values())
    def kernel(delta):
        x=delta*segment_s
        value=sum(c*np.sinc(x-j)*np.exp(-1j*np.pi*(x-j)) for j,c in coeff.items())
        return segment_s*abs(value)**2/norm
    def mech(nu): return f0**4/((f0*f0-nu*nu)**2+(nu*f0/q)**2)
    # Resolve the narrow mechanical pole explicitly and its Blackman main lobe.
    cuts=sorted(set([f0-1,f0-.1,f0-.01,f0-f0/q,f0,f0+f0/q,f0+.01,f0+.1,f0+1,f]))
    value=sum(quad(lambda nu:(kernel(f-nu)+kernel(f+nu))*mech(nu)/mech(f),lo,hi,epsabs=1e-8,epsrel=1e-7,limit=150)[0]
              for lo,hi in zip(cuts[:-1],cuts[1:]))
    # Bound omitted |nu-f| >= delta using Blackman's exact rational envelope.
    delta=1-abs(f-f0);zmin=delta*segment_s
    c=(.18+1.68/zmin**2)/(np.pi*(1-1/zmin**2)*(1-4/zmin**2))
    tail=2*c*c/(5*norm*segment_s**5*delta**5)*(q*q/(1-1/(4*q*q)))/mech(f)
    return value+tail


def calculate(inputs,baseline):
    rc=1e-7;T=.01;C=10.;tau=1e6
    force_upper=baseline['force']['literal_workbook']['cached_upper95_N2_per_hz'];torque_upper=5.7e-34
    traj=inputs['trajectory']; lp=np.array(traj['lpf']['velocity_m_s']);el=np.array(traj['earth_lpf']['velocity_m_s'])
    es=np.array(traj['earth_sodium']['velocity_m_s'])
    vbench=float(np.linalg.norm(es,axis=1).max()+500.)
    lp_axis=lp.mean(axis=0);lp_axis/=np.linalg.norm(lp_axis)
    lp_min=float(np.min(lp@lp_axis))
    # Unknown sensor epochs/orientation are not imputed from the sodium epoch.
    sensor_min=28000.  # declared terrestrial barycentric speed envelope, including rotation
    a=C/moving_sphere_exponent_per_a(T,tau,vbench,1e-7,2200.,1e-7,rc)
    a_static=C/moving_sphere_exponent_per_a(T,tau,0.,1e-7,2200.,1e-7,rc)
    assembly=assembly_psd_upper_per_a(rc,sensor_min)
    qtorque=cube_response(rc)[1]
    torque_per_a=moving_psd_upper_per_a(qtorque,1.928,np.sqrt(3)*.046,rc,lp_min-2000.)
    sodium=sodium_forward(inputs['sodium'],baseline,a,float(np.linalg.norm(es,axis=1).min()-500))
    # Largest admitted a at which at least one scan can change by one measured SE.
    def first_resolvable(wind):
        low,high=a,a*1e10
        for _ in range(36):
            mid=np.sqrt(low*high)
            if sodium_forward(inputs['sodium'],baseline,mid,wind)['max_visibility_change_over_se']<1:low=mid
            else:high=mid
        return high
    sodium_a_scale=first_resolvable(float(np.linalg.norm(es,axis=1).min()-500))
    fgrid=np.array(inputs['sensor']['frequency_grid_hz']);n=2**22;dt=1e-5;seg=n*dt
    fp=baseline['force']['spectra'][0]['f0_hz']
    # Explicitly choose retained-bin offsets, not removed peak bins.
    k=int(round(fp*seg));frequencies=(k+np.array([-12,-6,-4,4,6,12]))/seg
    windows=[]
    for q in [1e5,1e6,2.83e6]:
        for f in frequencies:
            windows.append({'qprime':q,'frequency_hz':float(f),'white_force_window_gain':window_gain(f,fp,q,seg),
                            'blackman_constant_gain_periodic':blackman_dc_gain(n,dt,f),
                            'blackman_constant_gain_symmetric':blackman_dc_gain(n,dt,f,True)})
    maxgain=max(x['white_force_window_gain'] for x in windows)
    comparisons=[]
    for error in [0.,2000.,10000.]:
        tupper=moving_psd_upper_per_a(qtorque,1.928,np.sqrt(3)*.046,rc,lp_min-error)
        fu=force_upper/(assembly['all_cross_terms_upper']*max(1.,maxgain))
        tu=torque_upper/tupper
        comparisons.append({'lpf_velocity_error_radius_m_s':error,'sensor_a_upper_s_minus2':fu,'space_a_upper_s_minus2':tu,
                            'joint_a_upper_s_minus2':min(fu,tu),'benchmark_a':a,
                            'force_fraction_upper':a/fu,'torque_fraction_upper':a/tu,
                            'survives_half_envelopes_and_10pct_benchmark_stress':bool(a*1.1<.5*min(fu,tu)),
                            'sodium_one_se_sensitivity_a':sodium_a_scale,
                            'interpretation':'sufficient envelope survivor using upper response bounds; not calibrated joint coverage'})
    # Two frame declarations; the terrestrial one is an Earth-following coordinate
    # model (non-inertial), not an independently optimized frame for each channel.
    zero_b=periodogram_weights(4096,1.,.003,'blackmanharris',False)
    zero_bd=periodogram_weights(4096,1.,.003,'blackmanharris',True)
    # Frozen static sensor force passes through DC susceptibility 1/k; comparison
    # to the resonance is conservative using gain <= max DC gain in audit bins.
    zero_sensor_a=assembly_psd_upper_per_a(rc,1.)['q0_upper']
    leak=max(x['blackman_constant_gain_symmetric'] for x in windows)
    terrestrial_force=a_static*zero_sensor_a/2*leak
    geo=lp-el;geo_axis=geo.mean(axis=0);geo_axis/=np.linalg.norm(geo_axis)
    geospeed_min=float(np.min(geo@geo_axis))
    geo_position=np.array(traj['lpf']['position_m'])-np.array(traj['earth_lpf']['position_m'])
    terrestrial_velocity=geo-np.cross(np.array([0.,0.,7.292115e-5]),geo_position)
    terrestrial_speeds=np.linalg.norm(terrestrial_velocity,axis=1)
    # Profile only explicitly sampled common inertial velocities. The domain is
    # a sensitivity neighbourhood of SSB, not a physically mandatory cutoff.
    center=es.mean(axis=0); arm=np.cross(center,[0.,0.,1.]);arm/=np.linalg.norm(arm)
    profile=[]
    for speed in [0.,5000.,10000.]:
        for axis in ([np.array([1.,0.,0.])] if speed==0 else np.vstack([np.eye(3),-np.eye(3)])):
            u=speed*axis;vb=center-u;vbmag=np.linalg.norm(vb);ca=float(-u@arm/vbmag)
            aa=C/moving_sphere_exponent_per_a(T,tau,vbmag,1e-7,2200.,1e-7,rc,ca)
            ss=sensor_min-speed;lp_frame=lp-u;axis_frame=lp_frame.mean(axis=0);axis_frame/=np.linalg.norm(axis_frame)
            lv=np.min(lp_frame@axis_frame)-2000.
            ab=assembly_psd_upper_per_a(rc,ss)['all_cross_terms_upper']
            tb=moving_psd_upper_per_a(qtorque,1.928,np.sqrt(3)*.046,rc,lv)
            profile.append({'common_u_icrf_m_s':u.tolist(),'a_required':aa,'force_fraction_upper':aa*ab*max(1.,maxgain)/force_upper,
                            'torque_fraction_upper':aa*tb/torque_upper,'benchmark_cos_angle':ca})
    limits=[]
    for length,omega,speed in [(113e-6,2*np.pi*3532.7,30000.),(.046,2*np.pi*.003,lp_min)]:
        frozen=lateral_frozen_factor(omega,speed,length,rc)
        for memory in [length/speed*.1,length/speed,length/speed*100,1.]:
            finite=memory*lateral_ou_factor(omega,memory,speed,length,rc)
            limits.append({'length_m':length,'omega_rad_s':omega,'speed_m_s':speed,'tau_s':memory,
                           'finite_tau_factor_times_tau_s':finite,'distributional_limit_s':frozen,'ratio':finite/frozen})
    return {'scope':'Gaussian random-unitary ensemble dephasing, with conditional local-Markov interferometry; not outcome selection',
            'benchmark':{'radius_m':1e-7,'density_kg_m3':2200.,'separation_m':1e-7,'time_s':T,'target_exponent':C,
                         'tau_s':tau,'lambda_s_inverse':a*tau,'rc_m':rc,'a_lambda_over_tau_s_minus2':a,
                         'speed_m_s':vbench,'separation_perpendicular_to_velocity':True,
                         'stationary_a':a_static,'motion_renormalization':a/a_static,
                         'frozen_atom_a':C/moving_sphere_exponent_per_a(T,np.inf,vbench,1e-7,2200.,1e-7,rc)},
            'velocity_audit':{'lpf_speed_minmax_m_s':[float(np.linalg.norm(lp,axis=1).min()),float(np.linalg.norm(lp,axis=1).max())],
                              'lpf_projection_speed_lower_m_s':lp_min,
                              'lpf_geocentric_speed_minmax_m_s':[geospeed_min,float(np.linalg.norm(lp-el,axis=1).max())],
                              'lpf_predicted_not_reconstructed':traj['lpf']['prediction_after_2016_04_14'],
                              'sensor_assumed_speed_minmax_m_s':[28000.,32000.],
                              'sodium_earth_speed_minmax_m_s':[float(np.linalg.norm(es,axis=1).min()),float(np.linalg.norm(es,axis=1).max())]},
            'sodium':sodium,'sodium_auxiliary_removal':{
                'no_advection_prescribed_path_exponent':float(a*triangle_path_exponent_per_a(171283.8,158,0.,rc,tau)),
                'mean_velocity_only':sodium_forward(inputs['sodium'],baseline,a,30000.,sigma=.00001)['max_visibility_change_over_se'],
                'sigma_5_to_7_percent':[sodium_forward(inputs['sodium'],baseline,a,28000.,sigma=158*s)['max_visibility_change_over_se'] for s in [.05,.07]]},
            'assembly_upper':assembly,'comparisons':comparisons,'windows':windows,
            'window_scope':'actual sensor Blackman sampling; LPF normalized window-independent PSD upper, not a raw likelihood',
            'terrestrial_frame':{'coordinate_law':'Earth-corotating field: Earth translation plus rotation about ICRF z at 7.292115e-5 rad/s; distinct non-inertial model',
                                  'a_required':a_static,'sensor_frozen_DC_leakage_upper_N2_per_hz':terrestrial_force,
                                  'lpf_instantaneous_speed_minmax_m_s':[float(terrestrial_speeds.min()),float(terrestrial_speeds.max())],
                                  'space_fraction_bound':None,
                                  'space_missing_bound_reason':'corotating trajectory lacks a validated whole-window chord bound; instantaneous speed is not substituted for one',
                                  'sodium_frozen_repeated_path_warning':'no event ergodicity asserted in a static terrestrial field; ensemble map alone insufficient'},
            'finite_window_zero_speed':{'lpf_example_duration_s':4096.,'frequency_hz':.003,
                                        'no_detrend_per_a':frozen_periodogram_per_a(qtorque,zero_b),
                                        'linear_detrend_per_a':0.,
                                        'projection_zero_verified':bool(frozen_periodogram_per_a(qtorque,zero_bd)<1e-20*frozen_periodogram_per_a(qtorque,zero_b)),
                                        'evidence':'response diagnostic, not LPF actual segment length'},
            'common_frame_profile':profile,'profile_scope':'13 sampled velocities on coordinate axes, speeds 0/5/10 km/s about SSB; no continuum or all-frame exclusion',
            'finite_tau_limit_checks':limits,
            'independent_triangle_check':{'exact_prescribed_per_a':triangle_path_exponent_per_a(171283.8,158,30000.,rc,tau),
                'local_markov_per_a':np.sqrt(np.pi)*rc/30000.*triangular_white_coefficient(171283.8,158,rc)},
            'nuisance_margins':{'sensor_additional_PSD_gain_allowed':force_upper/(a*assembly['all_cross_terms_upper']*max(1.,maxgain)),
                'sensor_tenfold_response_and_half_envelope':bool(10*a*assembly['all_cross_terms_upper']*max(1.,maxgain)<.5*force_upper),
                'space_additional_PSD_gain_allowed_2km_s':torque_upper/(a*torque_per_a)},
            'trajectory_removal':{'counterfactual_speed_m_s':.001,
                'axis_aligned_half_differential_torque_fraction':a*qtorque*lateral_frozen_factor(2*np.pi*.003,.001,.046,rc)/torque_upper,
                'nominal_trajectory_uniform_fraction':a*torque_per_a/torque_upper,
                'interpretation':'same field strength; without motion metadata millimetre/s drift remains possible and can produce a large response; no empirical time modulation or exclusion claimed'},
            'relative_pair':{'baseline_m':.376,'along_velocity_frozen_factor_at_3mHz':float(relative_pair_factor(2*np.pi*.003*.376/lp_min)),
                             'independent_mass_factor':.5,'scope':'geometry diagnostic only; physical bound allows all orientations and delayed correlations'},
            'general_spectrum':{'constructive_atom':'C_t=a/2, temporal measure pi*a*delta_0 under C=int S exp(-i omega t)domega/(2pi)',
                                'physical_kernel_LP_enabled':False,'reason':'constructive spectrum suffices here; no continuum interpolation/tail certificate asserted'},
            'claim':'sufficient joint envelope survivor conditional on declared source/mode/response and trajectory-error bounds; no new empirical exclusion'}


def render(result):
    b=result['benchmark'];c=result['comparisons'][1];s=result['sodium'];v=result['velocity_audit']
    lines=['# Moving-field joint dephasing compatibility','',
           '**Result:** motion removes the vanishing point-frequency response of the old stationary construction. Recomputing the benchmark in the same frame still gives a sufficient joint envelope survivor under the explicit assumptions below. This is ensemble dephasing, not outcome selection.','',
           f"For a 100 nm silica sphere, 100 nm transverse separation and 10 ms duration, the SSB-frame witness has tau={b['tau_s']:.0e} s, lambda={b['lambda_s_inverse']:.6g} /s, and lambda/tau={b['a_lambda_over_tau_s_minus2']:.6g} /s². Its exponent is 10 by construction. Motion increases the required lambda/tau by {b['motion_renormalization']:.3g} relative to the stationary benchmark.",'',
           '| Constraint / comparison | Result |','|---|---|',
           f"| Complete sensor assembly, all cross terms bounded, Blackman window stress | <= {c['force_fraction_upper']:.5g} of the quoted force envelope |",
           f"| LPF half-differential torque, any orientation/cross correlation, nominal velocity minus 2 km/s | <= {c['torque_fraction_upper']:.5g} of the published envelope |",
           f"| Sodium signed mass/velocity mixture, orientation outer bound | visibility changes <= {s['max_visibility_change_over_se']:.5g} measured standard errors |",
           f"| Force alone: sufficient allowed a | <= {c['sensor_a_upper_s_minus2']:.5g} /s² |",
           f"| Space alone: sufficient allowed a | <= {c['space_a_upper_s_minus2']:.5g} /s² |",
           f"| Joint: sufficient allowed a | <= {c['joint_a_upper_s_minus2']:.5g} /s² |",
           f"| Sodium one-SE response scale (not a confidence limit) | {c['sodium_one_se_sensitivity_a']:.5g} /s² |",'',
           'The response upper bounds certify a sufficient allowed region; their failure does not certify exclusion. Dataset removal is given by the single-channel columns. No joint 95% interpretation is assigned to the published envelopes or to the sodium distinguishability diagnostic.','',
           '## What the auxiliary records change','',
           f"JPL vectors give nominal LPF speeds {v['lpf_speed_minmax_m_s'][0]:.2f}–{v['lpf_speed_minmax_m_s'][1]:.2f} m/s in the SSB frame and {v['lpf_geocentric_speed_minmax_m_s'][0]:.2f}–{v['lpf_geocentric_speed_minmax_m_s'][1]:.2f} m/s relative to Earth over the padded February 2017 run. The LPF source explicitly labels this a prediction after April 2016; 0, 2 and 10 km/s error radii are sensitivity assumptions, not tracking uncertainties.",'',
           'The sodium mass spectrum, optical powers and author velocity distribution control the nonlinear harmonic calculation. Suppressing advection in the finite-time prescribed-path diagnostic gives large dephasing; the moving calculation leaves it essentially unchanged. The archive supplies a velocity model, not event-resolved velocities. Removing its width and varying it through the paper\'s 5–7% range is saved in the JSON.','',
           'A frozen field can correlate two spatially separated masses at a time lag equal to their separation divided by speed. Equal-time separation is therefore insufficient to assume independent colored force spectra. The implemented torque upper bound includes the worst allowed pair correlation; the JSON also records a constructive delayed-correlation zero.','',
           '## Finite observation and geometry','',
           'The sensor paper specifies 100 kHz sampling, 2^22-sample Blackman blocks and 40–80 averages. The code evaluates mechanical-pole/window convolution at retained bins, exact constant leakage for both Blackman conventions, and Qprime sensitivity. For moving channels, a PSD bound uniform in frequency survives every nonnegative normalized spectral window. Transfer calibration remains a premise: this does not reconstruct missing filtered time streams.','',
           'The sensor bound adds the layered load, Si flexure, magnet and an epoxy mass bound at the spectral-amplitude level, retaining every cross term. It assumes a vertical mode between zero and one when normalized at the load; a factor-ten response stress and the maximum gain allowed by the data are documented in the derivation. The LPF bound is uniform over cube orientation and relative-mass phase, so missing attitude cannot invalidate it within the rigid-source model.','',
           '## Quantifiers and remaining empirical limit','',
           'One common SSB field supplies the benchmark and every response. A separate Earth-corotating model is calculated and explicitly marked non-inertial; no instrument is assigned its own optimizing frame. The 13 inertial frame samples are a sensitivity profile, not an exclusion over all velocities. No high-speed or large-tau cutoff is used to manufacture exclusion.','',
           'The constructive frozen spectral atom and finite-OU witness are response-level survivors. The sodium result is conditional on the local advected Markov/paraxial map and measures change from the audited optical baseline, whose existing residuals remain. It is not a newly calibrated goodness-of-fit acceptance. The LPF envelope remains a published summary, with no reconstructed attitude/control-subtraction likelihood. Thus the implemented creative combination is substantive but does not yet meet the stronger standard of a fully calibrated empirical joint survivor.','',
           '![Moving-field comparisons](moving-diagnostics.svg)','',
           'See [physical derivation and acquisition audit](../moving-frame.md) and [machine-readable comparisons](motion.json).','']
    return '\n'.join(lines)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--data-dir',type=Path);p.add_argument('--check',action='store_true');args=p.parse_args()
    path=HERE/'results/motion.json'
    expected=json.loads(path.read_text()) if path.exists() else None
    if args.data_dir:verify_sources(args.data_dir);inputs=source_summary(args.data_dir)
    elif expected:inputs=expected['inputs']
    else:p.error('First generation requires --data-dir')
    baseline=json.loads((HERE/'results/baselines.json').read_text())['data']
    from analyze import canonical,compare,serial,metadata
    meta=metadata()
    meta['baseline_data_sha256']=hashlib.sha256(serial(baseline).encode()).hexdigest()
    meta['auxiliary_summary_sha256']=hashlib.sha256(serial(inputs).encode()).hexdigest()
    value={'metadata':meta,'inputs':canonical(inputs),'data':canonical(calculate(inputs,baseline))}
    if args.check:
        if expected['metadata']!=value['metadata']:raise AssertionError('Stale motion provenance')
        compare(expected['inputs'],value['inputs']);compare(expected['data'],value['data'])
        if (HERE/'results/moving-report.md').read_text()!=render(expected['data']):raise AssertionError('Stale moving report')
        with tempfile.TemporaryDirectory() as tmp:
            dest=Path(tmp)/'moving.svg';plot(expected['data'],dest)
            if dest.read_bytes()!=(HERE/'results/moving-diagnostics.svg').read_bytes():raise AssertionError('Stale moving plot')
        print('Moving responses, joint witness, auxiliary removals and windows reproduce'+(' (full sources)' if args.data_dir else ' (offline summaries)'))
    else:
        path.write_text(serial(value));saved=json.loads(path.read_text());(HERE/'results/moving-report.md').write_text(render(saved['data']))
        plot(saved['data'],HERE/'results/moving-diagnostics.svg')
        print('Generated moving-field analysis and joint dephasing report')


if __name__=='__main__':main()
