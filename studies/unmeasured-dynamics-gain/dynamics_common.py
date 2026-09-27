"""Reuse the unchanged PR17 confidence construction and physical gate class."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
HERE=Path(__file__).resolve().parent
GAIN=HERE.parent/'gst-rb-combination-gain'
sys.path.insert(0,str(GAIN))
from gain_core import BASE,GRID,region,residual_constraints,operations,GeneralEvaluator,save,sources
from sqd_sources import parse_word
from sqd_general import bounds

def baseline_models():return json.loads((GAIN/'results/witnesses.json').read_text())
def target_hash(w):return hashlib.sha256(bytes(w)).hexdigest()
def rows_from_cache(cache):return sources(cache)

def candidates(rows,R):
    known={r['word'] for r in rows};i0=next(i for i,r in enumerate(rows) if r['word']==())
    lookup={r['word']:i for i,r in enumerate(rows) if r['family']=='GST'}
    out=[]
    for i,r in enumerate(rows):
        if r['family']!='RB':continue
        w=r['word'];a=R['hi'][i0];b=R['hi'][i]
        for name,target,f in [('repeat',w*2,None),('append_xx',w+(1,1),(1,1)),('prepend_xx',(1,1)+w,(1,1)),('append_yy',w+(2,2),(2,2)),('prepend_yy',(2,2)+w,(2,2))]:
            if target in known:continue
            record={'kind':name,'source_index':i,'word':target,'word_sha256':target_hash(target),'length':len(target)}
            if name=='repeat':record['upper']=min(1,(2*np.sqrt(b)+np.sqrt(a))**2)
            elif name.startswith('append'):
                fi=lookup[f];record['flip_index']=fi
                record['lower']=max(0,R['lo'][fi]-(np.sqrt(a*(1-b))+np.sqrt(b*(1-a))) if a+b<=1 else 0)
            out.append(record)
    return out
