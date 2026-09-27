"""Read pinned alternate RB spelling as data; never execute deposited scripts."""
import ast,re,urllib.request
from dynamics_common import *
from sqd_sources import external_cache,digest

def audit_extra(cache,rows,download=False):
    cache=external_cache(cache)
    e=json.loads((HERE/'extra-manifest.json').read_text())['files'][0];p=cache/e['cache_name']
    if download and not p.exists():
        b=urllib.request.urlopen(e['url'],timeout=90).read()
        if digest(b)!=e['sha256']:raise ValueError('alternate download checksum')
        p.write_bytes(b)
    b=p.read_bytes()
    if len(b)!=e['bytes'] or digest(b)!=e['sha256']:raise ValueError('alternate source checksum')
    mappings={}
    for line in (cache/'Cliffs_to_prims.py').read_text().splitlines():
        m=re.fullmatch(r"cliff_to_prims\['(Gc[0-9]+)'\] = (\[.*\])",line)
        if m:mappings[m[1]]=parse_word(''.join(ast.literal_eval(m[2])))
    if len(mappings)!=48:raise ValueError('unexpected Clifford alphabet')
    lines=b.decode().splitlines();rb=[r for r in rows if r['family']=='RB']
    if lines[0]!='## Columns = plus count, count total' or len(lines)-1!=len(rb):raise ValueError('alternate count schema')
    for line,r in zip(lines[1:],rb):
        f=line.split()
        if len(f)!=5 or f[3:]!=['0','0'] or not re.fullmatch(r'(?:Gc[0-9]+)+',f[0]):raise ValueError('alternate row schema')
        w=tuple(g for c in re.findall(r'Gc[0-9]+',f[0]) for g in mappings[c])
        if w!=r['word'] or (float(f[1]),float(f[2]))!=(r['k'],r['n']):raise ValueError('Clifford/primitive mismatch')
    return {'rows':len(rb),'independent_new_counts':0}

if __name__=='__main__':
    import argparse
    a=argparse.ArgumentParser();a.add_argument('--cache',required=True);a.add_argument('--download',action='store_true');a.add_argument('--literature',action='store_true');args=a.parse_args()
    rows=sources(args.cache,args.download);print(audit_extra(args.cache,rows,args.download))
    if args.literature:
        cache=external_cache(args.cache)
        for e in json.loads((HERE/'literature-manifest.json').read_text())['files']:
            p=cache/e['cache_name']
            if args.download and not p.exists():
                b=urllib.request.urlopen(e['url'],timeout=90).read()
                if digest(b)!=e['sha256']:raise ValueError('literature download checksum')
                p.write_bytes(b)
            if len(p.read_bytes())!=e['bytes'] or digest(p.read_bytes())!=e['sha256']:raise ValueError('literature checksum')
        print('Pinned public literature supplement verified; excluded from count likelihood.')
