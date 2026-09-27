"""Compare three data combinations and construct long-memory survival sequences.

The selected calculation is a deterministic published-envelope feasibility study,
not a new joint confidence region. Radiation uses a declared response approximation.
"""
import io
from pathlib import Path
import zipfile
import numpy as np
from scipy import constants as co, special as sp
from scipy.integrate import quad
from scipy.stats import chi2

from response import (M0,L,velocity_nodes,optical_coefficients,triangular_white_coefficient,
                      layered_force_coefficient,sphere_static_coefficient,ou_time,ou_spectrum,
                      moving_local_ou_factor)
from io_data import sodium,table


def cube_response(rc,side=.046,mass=1.928):
    """Single-cube two-sided force/torque coefficients divided by lambda.

    Equal to ONE-sided PSD coefficients of (F1-F2)/2, (tau1-tau2)/2
    for independent identical masses. Derived from separable Fourier integrals.
    Valid numerical implementation used only at rc/side <= 1e-3.
    """
    if not 0<rc/side<=1e-3: raise ValueError("Cube quadrature regime not implemented")
    x = side/(2*rc); e = np.exp(-x*x)
    i0 = 2*np.pi*(side*sp.erf(x)+2*rc/np.sqrt(np.pi)*np.expm1(-x*x))
    i2 = 2*np.sqrt(np.pi)/rc*(1-e)
    cross = -.5*i0+rc*rc*i2
    m0 = rc*np.sqrt(np.pi)*sp.erf(x)
    m1 = 2*rc*rc*(1-e)
    m3 = 8*rc**4*(1-(1+x*x)*e)
    derivative = 2*np.sqrt(np.pi)/rc*(side**3*m0/12-side*side*m1/4+m3/6)
    factor = co.hbar**2*rc**3/(np.pi**1.5*M0**2)*(mass/side**3)**2
    return factor*i2*i0*i0, factor*i0*(2*i2*derivative-2*cross*cross)


def radiation_charge_factor(energy_kev,alpha=1.,beta=1.04):
    """XENONnT 2026 supplement Eq S1 and Table S1; rc >= 1.15e-8 m."""
    n = np.array([2,2,6,2,6,10,2,6,10,2,6])
    radius = np.array([.015,.064,.055,.170,.165,.149,.395,.412,.461,1.024,1.231])*1e-10
    energy = np.asarray(energy_kev)*1000*co.e
    scale = energy/(co.hbar*co.c)
    sinc = lambda x: np.sinc(x/np.pi)
    charge = np.full_like(energy,54**2+54,dtype=float)
    for j in range(len(n)):
        charge -= 108*n[j]*sinc(radius[j]*scale)
        charge += n[j]*(n[j]-1)*sinc(alpha*radius[j]*scale)
        for k in range(j):
            charge += 2*n[j]*n[k]*sinc(beta*abs(radius[j]-radius[k])*scale)
    if np.any(charge < -1e-8): raise ArithmeticError("Negative radiation spectrum")
    return np.maximum(0,charge)


def radiation_response(energies,efficiency,rc,tau,alpha=1.):
    """Expected counts/lambda: 1.16 tonne year, identity reconstructed-energy map.

    No resolution/migration model: this is a feasibility approximation, not an
    experimental confidence statement about lambda.
    """
    if rc<1.15e-8: raise ValueError("Outside charge-factor approximation")
    e = np.asarray(energies)
    prefactor = co.hbar*co.e**2/(4*np.pi**2*co.epsilon_0*co.c**3*M0**2*rc**2)
    atoms_seconds = 1160/(131.293*M0)*(365.25*86400)
    omega = e*1000*co.e/co.hbar
    return atoms_seconds*np.trapezoid(prefactor*radiation_charge_factor(e,alpha)/e*
                                     efficiency*ou_spectrum(omega,tau),e)


def xenon_audit(directory):
    with zipfile.ZipFile(Path(directory)/"XENONnT_ER_v2.zip") as z:
        events = np.loadtxt(io.BytesIO(z.read("zenodo/data_unbinned_1to30kev.txt")))
        efficiency = table(z.read("zenodo/efficiency.txt"),2)
        background = table(z.read("zenodo/bkg_model.txt"),2)
        if events.ndim!=1 or not np.isfinite(events).all() or np.any((events<1)|(events>30)):
            raise ValueError("Invalid ER energies")
        if np.any((efficiency[:,1]<0)|(efficiency[:,1]>1)): raise ValueError("Invalid efficiency")
    energy = np.linspace(1,30,2001)
    eff = np.interp(energy,efficiency[:,0],efficiency[:,1])
    upper = float(chi2.ppf(.95,2*(len(events)+1))/2)
    return {"evidence":"observed","events":len(events),"energy_min_kev":float(events.min()),"energy_max_kev":float(events.max()),
            "efficiency_rows":len(efficiency),"background_rows":len(background),
            "signal_counts_upper95_no_background_subtraction":upper,
            "interpretation":"Poisson count bound with nonnegative unknown background; parameter conversion requires response assumptions"},energy,eff


def white_scan_diagnostics(directory,protocol):
    scans,masses,mw = sodium(directory)
    v,vw = velocity_nodes(); weight = mw[:,None]*vw[None,:]
    rates = np.array(protocol["white_rate_grid_s_inverse"])
    output = []
    for rc in protocol["rc_grid_m"]:
        damping = np.exp(-rates[:,None,None]*triangular_white_coefficient(masses[:,None],v[None,:],rc)[None,:,:])
        best = {"increase":0.,"scan":None}; curves = []
        for scan in scans:
            s0,s1 = optical_coefficients(masses[:,None],v[None,:],scan.powers_w)
            vis = 2*abs(np.sum(damping*(weight*s1)[None,:,:],axis=(1,2)))/np.sum(weight*s0)
            j = int(np.argmax(vis-vis[0]))
            if vis[j]-vis[0]>best["increase"]:
                best = {"increase":float(vis[j]-vis[0]),"scan":scan.name,"lambda":float(rates[j]),
                        "zero_visibility":float(vis[0]),"nonzero_visibility":float(vis[j])}
            if scan.name in protocol["display_scans"]: curves.append({"scan":scan.name,"visibility":vis.tolist()})
        output.append({"rc_m":rc,"largest_visibility_increase":best,"display_curves":curves})
    return {"evidence":"observed","interpretation":"fixed-calibration predictions using measured powers/mass weights; no fitted collapse rate",
            "rate_grid":rates.tolist(),"results":output,"joint_confidence_region_enabled":False}


def combinations(directory,protocol, radiation_input=None, original_diagnostic=None):
    if radiation_input is None:
        audit,energy,eff = xenon_audit(directory)
    else:
        audit = radiation_input["audit"]
    taus = [0.,1e-12,1e-8,1e-4,1.,1e4,1e8,1e12]
    selected,runner,escapes = [],[],[]
    for benchmark in protocol["benchmarks"]:
        for rc in protocol["rc_grid_m"]:
            g = sphere_static_coefficient(benchmark["radius_m"],benchmark["density_kg_m3"],benchmark["separation_m"],rc)
            load = 2*layered_force_coefficient(rc)
            _,torque = cube_response(rc)
            for tau in taus:
                required = benchmark["required_log_suppression"]/(g*ou_time(benchmark["time_s"],tau))
                a = protocol["published_force_upper_N2_per_hz"]/(load*ou_spectrum(2*np.pi*3532.7,tau))
                b = 5.7e-34/(torque*ou_spectrum(2*np.pi*.003,tau))
                joint = min(a,b)
                selected.append({"benchmark":benchmark["name"],"rc_m":rc,"tau_s":tau,
                                 "required_lambda":required,"force_only_lambda_upper":float(a),
                                 "space_only_lambda_upper":float(b),"joint_lambda_upper":float(joint),
                                 "benchmark_feasible":bool(required<=joint),
                                 "nuisance_half_to_double_envelope_feasible":[bool(required<=.5*joint),bool(required<=2*joint)]})
                # Physical covariance variance stays finite as tau -> infinity, despite lambda -> infinity.
                if tau>=1:
                    escapes.append({"benchmark":benchmark["name"],"rc_m":rc,"tau_s":tau,"lambda_s_inverse":required,
                                    "force_fraction_of_upper":float(required/a),"space_fraction_of_upper":float(required/b),
                                    "macro_log_suppression":float(required*g*ou_time(benchmark["time_s"],tau)),
                                    "lambda_over_2tau":required/(2*tau)})
    for rc in protocol["rc_grid_m"]:
        for tau in taus:
            a = protocol["published_force_upper_N2_per_hz"]/(2*layered_force_coefficient(rc)*ou_spectrum(2*np.pi*3532.7,tau))
            upper = []
            for alpha in [1.,1.5]:
                key = f"{rc}:{tau}:{alpha}"
                k = radiation_response(energy,eff,rc,tau,alpha) if radiation_input is None else radiation_input["counts_per_lambda"][key]
                upper.append(audit["signal_counts_upper95_no_background_subtraction"]/k)
            runner.append({"rc_m":rc,"tau_s":tau,"force_only_upper":float(a),"radiation_only_upper_alpha_range":upper,
                           "joint_upper_alpha_range":[min(float(a),u) for u in upper],
                           "radiation_efficiency_20pct_reduction_upper":[u/.8 for u in upper],
                           "evidence":"observed","interpretation":"identity-energy response feasibility only; not official 2026 limit"})
    motion = [{"tau_s":tau,"static_T10ms_ratio":ou_time(.01,tau)/.01,
               "force_spectrum_ratio":float(ou_spectrum(2*np.pi*3532.7,tau)),
               "moving_local_approximation":moving_local_ou_factor(158.,1e-7,tau)} for tau in protocol["ou_tau_grid_s"]]
    return {"selected":{"combination":"layered force sensor + LISA Pathfinder published rotational envelope",
                        "evidence":"published_constraint","statistical_interpretation":"deterministic compatibility with quoted envelopes; NOT joint 95% coverage",
                        "rows":selected,"escape_witnesses":escapes,
                        "physical_premise":"stationary-body factorization; velocity through a common preferred noise frame is not calibrated; colored results are conditional envelope calculations only"},
            "runner_up":{"combination":"layered force sensor + XENONnT 1-30 keV public ER subset",
                         "audit":audit,"rows":runner,"official_2026_limit_reproduced":False},
            "original_pairing":{"diagnostic":original_diagnostic if original_diagnostic is not None else white_scan_diagnostics(directory,protocol),"colored_motion_comparison":motion},
            "claim_level":"conditional feasibility and constructive limitation of finite-positive-frequency envelopes",
            "unrestricted_experimental_exclusion_enabled":False}
