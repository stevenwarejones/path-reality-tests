"""Optional independent check against the inspected, hash-pinned author module.

Only talbotlau.py is executed, after the complete source archive passes verification.
Its top level contains imports and function definitions, no acquisition or plotting.
"""
import argparse
import types
import zipfile
from pathlib import Path
import numpy as np
from scipy import constants as co
from io_data import verify_sources,digest
from response import optical_coefficients,author_sphere_tau_exponent,L,D,M0


def oracle(directory):
    verify_sources(directory)
    with zipfile.ZipFile(Path(directory)/"Na-Cluster-Interference.zip") as z:
        source = z.read("Na-Cluster-Interference/talbotlau.py")
    module = types.ModuleType("pinned_author_talbotlau")
    exec(compile(source,"pinned_author_talbotlau.py","exec"),module.__dict__)
    errors = []
    for mass in [150000.,170000.,190000.]:
        for velocity in [140.,158.,176.]:
            particle = {"mass":mass*M0}
            particle["alpha_m3_of_m"] = lambda: -4.45e-30*particle["mass"]/co.m_u/22.9897
            particle["sigma_abs_of_m"] = lambda: (.5837*(particle["mass"]/co.m_u/1000)**(2/3)-1.5)*4.24e-20
            particle["radius_of_m"] = lambda: (3*particle["mass"]/(4*np.pi*971))**(1/3)
            powers = [.06162,.019844,.06786]
            gratings = []
            for power,waist in zip(powers,[560e-6,540e-6,575e-6]):
                grating = dict(type="laser",period=D,power=power,w_y=waist,reflectivity=1,**{"ionization probability":1,"classical":False})
                module.updategrating(grating,particle,velocity)
                gratings.append(grating)
            setup = dict(particle=particle,T_1=L/velocity,T_2=L/velocity,tau_e=np.nan,hbar_over_sigma_q=1e-8)
            setup.update({f"grating_{j+1}":g for j,g in enumerate(gratings)})
            s0,s1 = optical_coefficients(np.array([mass]),np.array([velocity]),powers)
            errors.extend([float(abs(s0[0]-module.S(0,setup))),float(abs(s1[0]-module.S(1,setup)))])
            setup["tau_e"] = np.array([1e14,1e15,1e16])
            ours = np.exp(-author_sphere_tau_exponent(mass,velocity)/setup["tau_e"])
            errors.extend(abs(ours-module.R(1,setup)).tolist())
    maximum = float(max(errors))
    if maximum>1e-10: raise ArithmeticError("Author oracle disagreement")
    return {"module_sha256":digest(source),"evaluations":len(errors),"max_absolute_discrepancy":maximum,
            "scope":"independent optics/spherical-reduction implementation versus inspected author functions"}


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir",type=Path,required=True)
    print(oracle(parser.parse_args().data_dir))
