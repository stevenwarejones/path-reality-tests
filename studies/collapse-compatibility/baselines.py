"""Source reproductions and explicit discrepancies, separate from new inference."""
import numpy as np
from scipy.optimize import curve_fit,least_squares
from scipy.stats import binned_statistic
from scipy import constants as co,special as sp
from io_data import sodium,force
from response import D,L,M0,velocity_nodes,optical_coefficients,author_sphere_tau_exponent


def harmonic_fit(position,counts):
    x,y = np.asarray(position),np.asarray(counts)
    design = np.array([np.ones(len(x)),np.cos(2*np.pi*x/D),np.sin(2*np.pi*x/D)]).T
    beta = np.linalg.lstsq(design,y,rcond=None)[0]
    residual = y-design@beta
    cov = np.linalg.inv(design.T@design)*(residual@residual)/(len(y)-3)
    amplitude = np.hypot(*beta[1:]); vis = amplitude/beta[0]
    grad = np.array([-vis/beta[0],beta[1]/(beta[0]*amplitude),beta[2]/(beta[0]*amplitude)])
    return {"mean_counts":float(beta[0]),"visibility":float(vis),"visibility_se_gaussian":float(np.sqrt(grad@cov@grad)),
            "phase_rad":float(np.arctan2(-beta[2],beta[1])),"residual_rms_counts":float(np.sqrt(np.mean(residual**2)))}


def author_sinusoid(position,counts):
    means = binned_statistic(position%D,counts,statistic="mean",bins=round(D/15e-9),range=(0,D))[0]
    coeff = np.fft.fft(means)/len(means)
    p0 = [abs(coeff[0]),2*abs(coeff[1]/coeff[0]),np.angle(coeff[1]*np.exp(-2j*np.pi*15e-9/2/D))]
    def model(x,mean,vis,phase): return mean*(1+vis*np.cos(2*np.pi*x/D+phase))
    p,cov = curve_fit(model,position,counts,p0=p0,bounds=([min(counts),0,-np.pi],[max(counts),1,np.pi]))
    return float(p[1]),float(np.sqrt(cov[1,1]))


def sodium_baseline(directory):
    scans,masses,mw = sodium(directory)
    v,vw = velocity_nodes(); weight = mw[:,None]*vw[None,:]
    result = []
    for scan in scans:
        fit = harmonic_fit(scan.position_m,scan.counts)
        av,ase = author_sinusoid(scan.position_m,scan.counts)
        s0,s1 = optical_coefficients(masses[:,None],v[None,:],scan.powers_w)
        _,classical = optical_coefficients(masses[:,None],v[None,:],scan.powers_w,classical=True)
        norm = np.sum(weight*s0)
        result.append({"scan":scan.name,"g2_power_mw":float(scan.powers_w[1]*1000),**fit,
                       "author_visibility":av,"author_visibility_se":ase,
                       "quantum_visibility":float(2*abs(np.sum(weight*s1))/norm),
                       "classical_visibility":float(2*abs(np.sum(weight*classical))/norm),
                       "signed_amplitude_fraction_negative":float(np.sum(weight*(s1 < -1e-15))),
                       "transmission":float(norm),"observed_total_counts":int(scan.counts.sum())})
    return {"evidence":"observed","interpretation":"descriptive fits and fixed-calibration model predictions",
            "scan_count":len(scans),"count_bins":sum(len(s.counts) for s in scans),
            "mean_mass_u":float(masses@mw),"mass_nodes":len(masses),"velocity_nodes":len(v),
            "max_visibility_fit_discrepancy":max(abs(r["visibility"]-r["author_visibility"]) for r in result),
            "rms_quantum_visibility_residual":float(np.sqrt(np.mean([(r["visibility"]-r["quantum_visibility"])**2 for r in result]))),
            "scans":result}


def weighted_line(x,y,error):
    a = np.array([np.ones(len(x)),x]).T
    cov = np.linalg.inv((a/error[:,None]).T@(a/error[:,None]))
    beta = cov@((a/error[:,None]).T@(y/error))
    return beta,cov,float(np.sum(((y-a@beta)/error)**2))


def psd_model(freq,pars,qprime):
    a,b,c,f0,offset = pars
    den = (freq*freq-f0*f0)**2+(freq*f0/qprime)**2
    return a*1e-12+(b*1e-18*f0**4+c*1e-12*(freq*freq-(f0+offset)**2)**2)/den


def fit_psd(spectrum,qprime=1e6,omit_shift=0):
    f,y,error = spectrum.T
    nav = int(round(float(np.median((y/error)**2))))
    if not np.allclose((y/error)**2,nav,rtol=2e-6): raise ValueError("Inconsistent PSD averaging errors")
    center = int(np.argmax(y)); start = center-3+omit_shift
    omitted = np.arange(start,start+6)
    if start<0 or start+6>len(f): raise ValueError("Invalid leakage mask")
    keep = np.ones(len(f),bool); keep[omitted] = False
    x,data = f[keep],y[keep]
    p = np.array([1.,.5,.1,f[center],1.]); sigma = error[keep]; old = np.inf
    for iteration in range(40):
        fit = least_squares(lambda par:(psd_model(x,par,qprime)-data)/sigma,p,
                            bounds=([0,0,0,f[center]-.07,-10],[100,100,100,f[center]+.07,10]),
                            x_scale="jac",ftol=1e-11,xtol=1e-11,gtol=1e-9,max_nfev=2000)
        if not fit.success: raise ArithmeticError("PSD optimizer failed")
        p = fit.x; prediction = psd_model(x,p,qprime); sigma = prediction/np.sqrt(nav)
        chi = float(np.sum(((prediction-data)/sigma)**2))
        if abs(chi-old)/chi<1e-4: break
        old = chi
    else: raise ArithmeticError("Recursive fit did not stabilize")
    return {"B_phi0_squared_per_hz":float(p[1]*1e-18),"f0_hz":float(p[3]),"parameters_scaled":p.tolist(),
            "qprime_assumed":qprime,"nav_inferred_from_errors":nav,"omitted_zero_based_indices":omitted.tolist(),
            "omitted_hz":f[omitted].tolist(),"rows":len(f),"iterations":iteration+1,"chi_squared":chi,
            "degrees_of_freedom":len(data)-5,
            "interpretation":"conditional recursive fit; measured Qprime and exact leakage mask missing"}


def force_baseline(directory):
    rows,spectra,cells,formulas = force(directory)
    high = rows[:,0]>=100
    beta,cov,chi = weighted_line(rows[high,2]*1e7,rows[high,3]*1e19,rows[high,4]*1e19)
    b0,b1 = beta*np.array([1e-19,1e-12]); se = np.sqrt(np.diag(cov))*np.array([1e-19,1e-12])
    kB,k,ek,f0,ef0 = [cells[c] for c in ["K20","L20","M20","N20","O20"]]
    b0p,eb0p,b1p,eb1p = [cells[c] for c in ["D20","E20","D23","E23"]]
    sf = 4*kB*k*b0p/(2*np.pi*f0*b1p)
    terms = [abs(sf)*ek/k,abs(sf)*ef0/f0,4*kB*k*eb0p/(2*np.pi*f0*b1p),
             4*kB*abs(b0p)*ek*eb1p/(2*np.pi*f0*b1p*b1p)]
    corrected = terms[:3]+[abs(sf)*eb1p/b1p]
    psds = []
    for temp,spectrum in spectra.items():
        fit = fit_psd(spectrum); pub = rows[rows[:,0]==temp,3][0]
        fit.update(temperature_mk=temp,published_B=float(pub),relative_B_discrepancy=float(fit["B_phi0_squared_per_hz"]/pub-1))
        psds.append(fit)
    sensitivity = []
    for q in [1e5,1e6,1e7]:
        for shift in [-1,0,1]:
            f = fit_psd(spectra[580],q,shift)
            sensitivity.append({"qprime":q,"mask_shift":shift,"B":f["B_phi0_squared_per_hz"],"chi_squared":f["chi_squared"]})
    return {"evidence":"observed","interpretation":"source reproduction, not independent repeated evidence",
            "weighted_y_only_fit":{"B0":float(b0),"B1":float(b1),"se_B0":float(se[0]),"se_B1":float(se[1]),"chi_squared":chi,"ndf":int(sum(high)-2)},
            "published_orthogonal_fit":{"B0":b0p,"se_B0":eb0p,"B1":b1p,"se_B1":eb1p},
            "B0_discrepancy_in_published_se":float((b0-b0p)/eb0p),
            "literal_workbook":{"force_estimate_N2_per_hz":sf,"linear_error_sum":sum(terms),
                                "upper95_N2_per_hz":1.17*sum(terms),"cached_upper95_N2_per_hz":cells["N29"],"formulas":formulas},
            "sensitivity_only_corrected_derivatives":{"linear_error_sum":sum(corrected),"independent_quadrature_error":float(np.linalg.norm(corrected)),"new_confidence_limit_enabled":False},
            "spectra":psds,"psd_580mk_sensitivity":sensitivity,
            "temperature_fit_inputs":[{"temperature_mk":int(row[0]),"T_over_Q_K":float(row[2]),"B":float(row[3]),"se_B":float(row[4])} for row in rows]}


def macroscopicity_reproduction(directory):
    scans,masses,mw = sodium(directory)
    v,vw = velocity_nodes(39); tau = np.logspace(14,16,400); weight = mw[:,None]*vw[None,:]
    damping = np.exp(-author_sphere_tau_exponent(masses[:,None],v[None,:])[None,:,:]/tau[:,None,None])
    t = L/158; m = masses.mean()*M0; y = co.h*t/(D*m*1e-8)
    t_tilde = t*(m/co.m_e)**2*(1-np.sqrt(np.pi/2)/y*sp.erf(y/np.sqrt(2)))
    prior = t_tilde/tau**2*np.sqrt(1/np.sqrt(1-.5**2*np.exp(-4*t_tilde/tau))-1)
    logpost = np.log(prior/prior.max())
    for scan in scans:
        s0,s1 = optical_coefficients(masses[:,None],v[None,:],scan.powers_w)
        transmission = np.sum(weight*s0)
        vis = 2*abs(np.sum(damping*(weight*s1)[None,:,:],axis=(1,2)))/transmission
        probability = transmission*(1+vis[:,None]*np.cos(2*np.pi*scan.position_m[None,:]/D-2.48))
        if np.any((probability<=0)|(probability>=1)): raise ArithmeticError("Invalid author likelihood")
        logpost += np.sum(scan.counts[None,:]*np.log(probability/(1-probability))+scan.incident_rate_hz*np.log1p(-probability),axis=1)
        logpost -= logpost.max()
    density = np.exp(logpost); widths = np.r_[np.diff(tau),np.diff(tau)[-1]]
    cdf = np.cumsum(density*widths); cdf /= cdf[-1]
    quantile = tau[np.searchsorted(cdf,.05)]
    return {"evidence":"observed","interpretation":"author fixed-calibration Bayesian convention, not a joint confidence region",
            "tau5_s":float(quantile),"log10_tau5":float(np.log10(quantile)),"length_hbar_over_sigma_q_m":1e-8,
            "tau_grid_s":[1e14,1e16,400],"phase_rad":-2.48,
            "density_relative_to_max_at_endpoints":[float(density[0]),float(density[-1])],
            "uncertainty_warning":"incident rate treated as binomial trial count; phase estimated from these scans",
            "CSL_mapping":"rc=hbar_over_sigma_q/sqrt(2); lambda=(m0/me)^2/tau_e"}
