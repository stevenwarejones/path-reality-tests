#!/usr/bin/env python3
"""Read-only audit of neighbouring controls and full-image support of the old witness."""
from pathlib import Path
from collections import Counter
import argparse,ast,hashlib,json,re,urllib.request,zipfile
import numpy as np
import xarray as xr
from scipy.stats import beta
from sources import digest,reconstruct
HERE=Path(__file__).resolve().parent

def audit(cache):
    reconstruct(cache)
    manifest=json.loads((HERE/'sources.json').read_text());base=cache/'controls';base.mkdir(exist_ok=True)
    notebook=cache/'atomic-notebook.ipynb'
    if not notebook.exists():
        with urllib.request.urlopen(manifest['notebook']['url'],timeout=120) as r:notebook.write_bytes(r.read())
    assert digest(notebook)==manifest['notebook']['sha256']
    text='\n'.join(''.join(c['source']) for c in json.loads(notebook.read_text())['cells'])
    many_map=ast.literal_eval(re.search(r'dataDict3NN\s*=\s*(\{.*?\})',text,re.S).group(1))
    run_maps=[ast.literal_eval(re.search(r'dataDict3NNSingle_'+str(j)+r'\s*=\s*(\{.*?\})',text,re.S).group(1)) for j in range(3)]
    with zipfile.ZipFile(cache/'boson-1.0.4.zip') as z:
        for name,sha in manifest['control_members'].items():
            f=base/name
            if not f.exists():f.write_bytes(z.read('boson-1.0.4/Boson sampling data/'+name))
            assert digest(f)==sha
    many=xr.open_dataarray(base/'3NN.nc')
    singles=[xr.open_dataarray(base/f'3NN_{j}.nc') for j in range(3)]
    def read(d,key):
        v=d.sel(key=key).values;ok=~np.all(np.isnan(v),axis=(0,2,3));v=np.nan_to_num(v[:,ok])
        assert np.all((v==0)|(v==1))
        assert np.all(v[0]==v[0,0])
        return v
    records=[]
    for key in many.key.values:
        v=read(many,float(key));assert np.all(v[0].sum(axis=(1,2))==3) and np.all(v[1].sum(axis=(1,2))<=3)
        hists=[];ns=[];b=[];coordinates=[]
        for d in singles:
            z=read(d,float(key));initial=np.argwhere(z[0,0]);assert len(initial)==1
            assert np.all(z[1].sum(axis=(1,2))<=1)
            iy,ix=map(int,initial[0]);N=z.shape[1];ns.append(N)
            coordinates.append([int(d.y.values[iy]),int(d.x.values[ix])])
            ct=z[1].sum(axis=(0,2));hist={int(i-iy):int(n)/N for i,n in enumerate(ct) if n}
            hist['lost']=float(1-ct.sum()/N);hists.append(hist)
            # One specified outcome, displacement y=0; ordinary 95% marginal interval, exploratory only.
            k=int(ct[iy]);l=0 if k==0 else beta.ppf(.025,k,N-k+1)
            u=1 if k==N else beta.ppf(.975,k+1,N-k)
            b.append([round(float(l),8),round(float(u),8)])
        tv=[]
        for i,j in [(0,1),(0,2),(1,2)]:
            keys=set(hists[i])|set(hists[j]);tv.append(round(sum(abs(hists[i].get(k,0)-hists[j].get(k,0)) for k in keys)/2,8))
        is_selected=bool(abs(key-2.45)<1e-9);is_late=bool(abs(key-4.65)<1e-9)
        run_key='2.45' if is_selected else ('4.65' if is_late else '(.5, 3.5, 15)')
        origins=[r[run_key][0] for r in run_maps]
        records.append({'nominal_time_ms':round(float(key),12),'many_shots':v.shape[1],
            'bunch_count':int(np.any(v[1].sum(axis=2)==3,axis=1).sum()),
            'many_run_ids':many_map[run_key],'single_shots':ns,'single_run_ids':origins,'single_input_sites':coordinates,
            'empirical_displacement_row_tv':tv,'center_row_probability_marginal_95_intervals':b,
            'same_nominal_selected_setting':is_selected,
            'independent_input_runs':len(set(origins))==3})
    # Full higher-order image support audit, including all observed output pixels.
    counts=json.loads((HERE/'results/counts.json').read_text())
    w=json.loads((HERE/'results/physical-witness.json').read_text())
    supported=set((y+s,x) for y,x,k in w['nonzero_cell_counts'] for s in [0,1,-1])
    a=xr.open_dataarray(cache/'selected/1D.nc');v=read(a,3);patterns=Counter()
    # Align after materializing: avoid backend-dependent lazy advanced-index ordering.
    iy=np.searchsorted(a.y.values,many.y.values);ix=np.searchsorted(a.x.values,many.x.values)
    assert np.array_equal(a.y.values[iy],many.y.values) and np.array_equal(a.x.values[ix],many.x.values)
    assert np.array_equal(np.take(np.take(v,iy,axis=2),ix,axis=3),read(many,2.45)), 'Different packaging is not an independent record'
    for img in v[1]:patterns[tuple(tuple(map(int,p)) for p in np.argwhere(img))]+=1
    bad={pat:k for pat,k in patterns.items() if any(p not in supported for p in pat)}
    # Every such pixel has identically zero amplitude from all three inputs.
    max_bad=max(bad.items(),key=lambda kv:kv[1])
    return {'records':records,'tv_intervals_status':'Exploratory point-TV and one-outcome marginal intervals; not simultaneous transfer guarantees.',
      'selected_date':'2022-07-02','source_notebook_sha256':manifest['notebook']['sha256'],
      'selected_full_image_support_audit':{'prepared_shots':sum(patterns.values()),
        'distinct_observed_parity_patterns':len(patterns),'witness_zero_probability_patterns':len(bad),
        'shots_in_zero_probability_patterns':sum(bad.values()),
        'example_pattern':[list(p) for p in max_bad[0]],'example_pattern_count':max_bad[1],
        'quantum_summary_witness_passes_full_image_support':False,
        'scope':'Rules out this committed inner example as a full-image fit, not all quantum or cluster models.'}}

def main():
    p=argparse.ArgumentParser();p.add_argument('--cache',type=Path,required=True);p.add_argument('--check',action='store_true');a=p.parse_args()
    result=audit(a.cache);text=json.dumps(result,indent=2,sort_keys=True)+'\n';out=HERE/'results/control-audit.json'
    if a.check:assert out.read_text()==text
    else:out.write_text(text)
    print(json.dumps(result['selected_full_image_support_audit'],indent=2))
    for r in result['records']:
        if r['nominal_time_ms'] in [2.45,4.65,2.428571428571]:print(r)

if __name__=='__main__':main()
