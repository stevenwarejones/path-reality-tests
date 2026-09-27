"""Matched stability records: exact overlap audit and a two-context control."""
import argparse,collections,hashlib,json,sys,urllib.request
from pathlib import Path
from fractions import Fraction as F
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'unmeasured-dynamics-gain'))
from dynamics_common import BASE,region,GRID,save,sources
from dynamics_verify import record

def fetch(cache,download=False):
    cache=Path(cache).resolve()
    if HERE.parents[1] in cache.parents or cache==HERE.parents[1]:raise ValueError('external cache required')
    manifest=json.loads((HERE/'stability-manifest.json').read_text())
    for f in manifest['files']:
        p=cache/f['cache_name']
        if download and not p.exists():
            cache.mkdir(parents=True,exist_ok=True);data=urllib.request.urlopen(f['url'],timeout=60).read()
            if hashlib.sha256(data).hexdigest()!=f['sha256']:raise ValueError('download hash mismatch')
            p.write_bytes(data)
        data=p.read_bytes()
        if len(data)!=f['bytes'] or hashlib.sha256(data).hexdigest()!=f['sha256']:raise ValueError('audit source hash mismatch')
    tree=json.loads((cache/'tree.json').read_text())
    if tree['sha']!=manifest['commit'] or tree['truncated']:raise ValueError('incomplete deposit inventory')
    files=[r['path'] for r in tree['tree'] if r['type']=='blob' and r['path'].startswith('ExperimentalData/')]
    notebook=json.loads((cache/'Figure5.ipynb').read_text());cells=notebook.get('cells',[]) or notebook['worksheets'][0]['cells']
    text='\n'.join(''.join(c.get('source',c.get('input',[]))) for c in cells)
    markers=['generate_fake_data','NonMarkovRBData','theta = 1.25e-2','envSwich','2015_03_30-RB_0320_condensed_cliffs.txt']
    if not all(m in text for m in markers):raise ValueError('analysis semantics changed')
    return {'experimental_data_paths':files,'analysis_markers_confirmed':markers,
            'interpretation':'Figure5 non-Markovian outcomes are simulated, not matched calibration observations; notebook and pickles never executed'}

def run(source_cache=None,audit_cache=None,download=False):
    rows=json.loads((BASE/'results/observations.json').read_text())
    if source_cache:
        raw=sources(source_cache)
        if len(raw)!=len(rows) or any(any(a[k]!=b[k] for k in ['family','row','word_sha256','k','n','length']) for a,b in zip(raw,rows)):raise ValueError('count-source binding mismatch')
    byhash=collections.defaultdict(list)
    for i,r in enumerate(rows):byhash[r['word_sha256']].append(i)
    overlap=[(h,ids) for h,ids in byhash.items() if {rows[i]['family'] for i in ids}=={'GST','RB'}]
    R=region(rows);q=F(1,100);slacks=[];records=[]
    for h,ids in overlap:
        for i in ids:slacks.extend([q-F(R['lo_num'][i],GRID),F(R['hi_num'][i],GRID)-q])
        records.append({'word_sha256':h,'indices':ids,'counts':[[rows[i]['k'],rows[i]['n']] for i in ids]})
    if min(slacks)<=0:raise ValueError('overlap-only context countermodel no longer compatible')
    out={'matched_words':len(overlap),'matched_rows':sum(len(i) for h,i in overlap),'records':records,
         'overlap_only_countermodel':{'rho':'|0><0|','all_native_gates':'identity','GST_effect':'diag(1/100,0)','RB_effect':'diag(1/100,1)',
           'operator_norm_effect_difference':1,'all_word_probability':record(q),'minimum_retained_cell_slack':record(min(slacks)),
           'scope':'Fits every matched-word count cell only; does NOT fit all GST records or retained source-wide groups'},
         'paper_assertions':['Interleaved GST/RB in final March 30 acquisition','Active pi-time and frequency feedback using interleaved probes','Final acquisition still reported statistically significant non-Markovianity'],
         'confirmed_missing_in_count_schema':['timestamps','per-shot order','feedback records','matched independently calibrated detector or state probes'],
         'not_proved':'No global stability upper bound follows from duplicate-word agreement alone; this is not an impossibility theorem using all GST/RB dynamics constraints'}
    if audit_cache:
        inventory=fetch(audit_cache,download)
        old=HERE/'results/stability-inventory.json'
        if old.exists() and json.loads(old.read_text())!=inventory:raise ValueError('inventory mismatch')
        if not old.exists():save(old,inventory)
    return out
if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--source-cache');a.add_argument('--audit-cache');a.add_argument('--download',action='store_true');a.add_argument('--check',action='store_true');v=a.parse_args()
    if v.download and not v.audit_cache:a.error('--download needs --audit-cache')
    out=run(v.source_cache,v.audit_cache,v.download);p=HERE/'results/stability.json'
    if v.check:
        if json.loads(p.read_text())!=out:raise ValueError('stability artifact mismatch')
    else:save(p,out)
    print(out['matched_words'],out['matched_rows'],out['overlap_only_countermodel']['minimum_retained_cell_slack']['decimal'])
