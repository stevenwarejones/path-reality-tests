#!/usr/bin/env python3
"""Construct an explicit physical inner example, using external singleton records."""
from pathlib import Path
import argparse,json,math
from fractions import Fraction as F
import numpy as np
import xarray as xr
from scipy.optimize import least_squares
from sources import reconstruct

HERE=Path(__file__).resolve().parent

def main():
    p=argparse.ArgumentParser();p.add_argument('--cache',type=Path,required=True);args=p.parse_args()
    reconstruct(args.cache) # hash/schema validation before use
    s=xr.open_dataarray(args.cache/'selected/1D_singles.nc')
    v=s.sel(key=30).isel(image_id=1).values
    v=np.nan_to_num(v[~np.all(np.isnan(v),axis=(1,2))]);cts=v.sum(axis=0).astype(int)
    q=cts/len(v);probs=np.stack([q,np.roll(q,1,axis=0),np.roll(q,-1,axis=0)],axis=-1)
    amp=np.sqrt(probs)
    def residual(ph):
        theta=np.c_[np.zeros(28),ph.reshape(28,2)]
        A=(amp*np.exp(1j*theta[:,None,:])).reshape(-1,3)
        G=A.conj().T@A
        z=np.array([G[0,1],G[0,2],G[1,2]])
        return np.r_[z.real,z.imag]
    res=least_squares(residual,np.random.default_rng(728).uniform(-np.pi,np.pi,56),
        xtol=1e-13,ftol=1e-13,gtol=1e-13,max_nfev=1000)
    assert res.success
    ts=[str(F(round(math.tan(t/2)*10**8),10**8)) for t in res.x]
    out={'description':'Exact probabilities and rational unit phases; sqrt probabilities define a contraction completed by loss modes. All identical internal states for the joint-summary example; orthogonal labels for the full-singleton-source example.',
         'single_shots':len(v),'nonzero_cell_counts':[[int(y),int(x),int(cts[y,x])] for y,x in np.argwhere(cts>0)],
         'phase_tangents':ts,'phase_rule':'(1-t^2+2it)/(1+t^2)','many_rows':list(range(8,20))}
    (HERE/'results/physical-witness.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')

if __name__=='__main__':main()
