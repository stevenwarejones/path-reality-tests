#!/usr/bin/env python3
"""Acquire pinned originals externally and reconstruct aggregate statistics."""
from pathlib import Path
import argparse
import hashlib
import json
import urllib.request
import zipfile

HERE=Path(__file__).resolve().parent

def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(2**20),b''):h.update(block)
    return h.hexdigest()

def reconstruct(cache):
    import numpy as np
    import xarray as xr
    manifest=json.loads((HERE/'sources.json').read_text())
    archive=cache/'boson-1.0.4.zip'
    assert digest(archive)==manifest['atomic_archive']['sha256']
    base=cache/'selected'
    base.mkdir(exist_ok=True)
    with zipfile.ZipFile(archive) as z:
        for name,sha in manifest['atomic_members'].items():
            target=base/name
            if not target.exists():
                target.write_bytes(z.read('boson-1.0.4/Boson sampling data/'+name))
            assert digest(target)==sha, name
    a=xr.open_dataarray(base/'1D.nc');s=xr.open_dataarray(base/'1D_singles.nc')
    def shots(d,key):
        v=d.sel(key=key).values
        valid=~np.all(np.isnan(v),axis=(0,2,3))
        v=v[:,valid]
        mask=~np.isnan(v[1,0])
        assert np.array_equal(~np.isnan(v),np.broadcast_to(mask,v.shape))
        v=np.nan_to_num(v)
        assert np.all((v==0)|(v==1))
        assert np.all(v[0]==v[0,0])
        return v,mask
    av,am=shots(a,3);sv,sm=shots(s,30)
    assert len(sv[0])==931 and len(av[0])==2999
    for key,shift in [(31,1),(32,-1)]:
        v,m=shots(s,key)
        assert np.array_equal(np.roll(sv,shift,axis=2),v)
        assert np.array_equal(np.roll(sm,shift,axis=0),m)
    assert np.array_equal(av[0,0],sum(np.roll(sv[0,0],shift,axis=0) for shift in [0,1,-1]))
    rows=np.flatnonzero(sm.any(axis=1));mrows=np.flatnonzero(am.any(axis=1))
    assert rows.tolist()==list(range(9,21)) and mrows.tolist()==list(range(8,20))
    assert np.all(sv[0].sum(axis=(1,2))==1)
    assert np.all(sv[1].sum(axis=(1,2))<=1)
    assert np.all(av[0].sum(axis=(1,2))==3)
    assert np.all(av[1].sum(axis=(1,2))<=3)
    witness_path=HERE/'results/physical-witness.json'
    if witness_path.exists():
        witness=json.loads(witness_path.read_text())
        cell=sv[1].sum(axis=0).astype(int)
        expected=[[int(y),int(x),int(cell[y,x])] for y,x in np.argwhere(cell>0)]
        assert witness['nonzero_cell_counts']==expected, 'Witness does not reproduce singleton cells'
    scan=[]
    for n in [2,3,4,5]:
        v,_=shots(a,n)
        assert np.all(v[0].sum(axis=(1,2))==n)
        assert np.all(v[1].sum(axis=(1,2))<=n)
        scan.append({'n':n,'prepared_shots':len(v[1]),
          'full_survival':int(np.sum(v[1].sum(axis=(1,2))==n)),
          'bunch_shots':int(np.sum(np.any(v[1].sum(axis=2)==n,axis=1)))})
    return {'selected_n':3,'single_shots':len(sv[1]),'many_shots':len(av[1]),
      'bunch_shots':scan[1]['bunch_shots'],'full_survival':scan[1]['full_survival'],
      'reference_row_indices':rows.tolist(),'many_row_indices':mrows.tolist(),
      'reference_row_counts':sv[1].sum(axis=(0,2))[rows].astype(int).tolist(),
      'copies':{'31':1,'32':-1},'scan_family':scan,
      'provenance':'Derived aggregate statistics, not raw events; Young et al., Zenodo 10453016, boson-1.0.4.zip.'}

def main():
    p=argparse.ArgumentParser();p.add_argument('--cache',type=Path,required=True)
    p.add_argument('--download',action='store_true');p.add_argument('--check',action='store_true');args=p.parse_args()
    args.cache.mkdir(parents=True,exist_ok=True)
    src=json.loads((HERE/'sources.json').read_text())['atomic_archive']
    target=args.cache/'boson-1.0.4.zip'
    if args.download and not target.exists():
        temp=target.with_suffix('.part')
        with urllib.request.urlopen(src['url'],timeout=120) as r,temp.open('wb') as f:
            while block:=r.read(2**20):f.write(block)
        assert digest(temp)==src['sha256']
        temp.rename(target)
    counts=reconstruct(args.cache)
    text=json.dumps(counts,indent=2,sort_keys=True)+'\n';path=HERE/'results/counts.json'
    if args.check:assert path.read_text()==text,'Source reconstruction differs'
    else:path.write_text(text)
    print(text)

if __name__=='__main__':main()
