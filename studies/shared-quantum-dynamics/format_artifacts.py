"""Keep derived numeric artifacts reviewable without dropping stored precision."""
import argparse,json
from sqd_sources import HERE


def formatted(obj,depth=0):
    indent='  '*depth
    if isinstance(obj,dict):
        if not obj:return '{}'
        return '{\n'+',\n'.join('  '+indent+json.dumps(k)+': '+formatted(v,depth+1) for k,v in sorted(obj.items()))+'\n'+indent+'}'
    if isinstance(obj,list):
        if all(not isinstance(v,(dict,list)) for v in obj):return json.dumps(obj,allow_nan=False)
        if all(isinstance(v,dict) and not any(isinstance(a,(dict,list)) for a in v.values()) for v in obj):
            return '[\n'+',\n'.join('  '+indent+json.dumps(v,sort_keys=True,allow_nan=False) for v in obj)+'\n'+indent+']'
        return '[\n'+',\n'.join('  '+indent+formatted(v,depth+1) for v in obj)+'\n'+indent+']'
    return json.dumps(obj,allow_nan=False)


if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--check',action='store_true');args=a.parse_args()
    for p in (HERE/'results').glob('*.json'):
        text=formatted(json.loads(p.read_text()))+'\n'
        if args.check:
            if text!=p.read_text():raise SystemExit('Artifact formatting differs: '+p.name)
        else:p.write_text(text)
