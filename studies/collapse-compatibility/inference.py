"""Conditional statistical calibration; never presented as the experimental likelihood."""
import numpy as np
from scipy.stats import norm,beta


def gaussian_simultaneous_upper(measured,known_sigma,alpha):
    x,s = np.asarray(measured),np.asarray(known_sigma)
    if x.shape!=s.shape or x.ndim!=1 or not len(x) or np.any(s<=0) or not 0<alpha<1: raise ValueError("Invalid Gaussian inputs")
    if not np.isfinite(x).all() or not np.isfinite(s).all(): raise ValueError("Nonfinite inputs")
    return np.maximum(0,x+norm.isf(alpha/len(x))*s)


def calibration_demo(seed,repeats,alpha):
    rng = np.random.default_rng(seed)
    matrix = np.array([[1.,.2],[.1,1.]])
    sigma = np.array([.1,.15])
    result = []
    for label,spectrum in [("zero_coupling",[0.,0.]),("injected",[.7,.4])]:
        truth = matrix@np.array(spectrum)
        samples = truth+rng.normal(size=(repeats,2))*sigma
        upper = np.maximum(0,samples+norm.isf(alpha/2)*sigma)
        failures = int(np.any(upper<truth,axis=1).sum())
        ci = [0. if failures==0 else beta.ppf(.025,failures,repeats-failures+1),
              1. if failures==repeats else beta.ppf(.975,failures+1,repeats-failures)]
        result.append({"scenario":label,"trials":repeats,"true_spectrum":spectrum,
                       "false_exclusions_of_true_response":failures,"binomial_mc_interval95":ci,
                       "recovered_spectrum_from_mean":np.linalg.solve(matrix,samples.mean(axis=0)).tolist()})
    return {"evidence":"simulated","alpha":alpha,"seed":seed,"cases":result,
            "procedure":"two known-variance Gaussian responses, Bonferroni upper bounds",
            "scope":"calibration demonstration only; NOT sodium counts or Blackman periodograms"}

