"""Audit runner-up sources as data only. No deposited Python/notebook is executed."""
import argparse,csv,hashlib,io,json,re,urllib.request,zipfile
from pathlib import Path
from sqd_sources import HERE,external_cache,digest


def audit(cache,download=False):
    cache=external_cache(cache);manifest=json.loads((HERE/'candidate-manifest.json').read_text());out={}
    for e in manifest['files']:
        p=cache/e['cache_name']
        if download and not p.exists():
            b=urllib.request.urlopen(e['url'],timeout=90).read()
            if digest(b)!=e['sha256']:raise ValueError('candidate checksum mismatch')
            cache.mkdir(parents=True,exist_ok=True);p.write_bytes(b)
        if digest(p.read_bytes())!=e['sha256']:raise ValueError('candidate checksum mismatch')
    with zipfile.ZipFile(cache/'belem.zip') as z:
        names=sorted(n for n in z.namelist() if re.fullmatch(r'wyniki/wyniki_testy_\d+\.csv',n))
        denominators=set();total=0;rows=0
        for name in names:
            for r in csv.DictReader(io.StringIO(z.read(name).decode())):
                k,n=int(r['0']),int(r['0'])+int(r['1'])
                if not 0<=k<=n:raise ValueError('invalid Belem counts')
                denominators.add(n);rows+=1;total+=n
        out['repeated_operation']={'jobs':len(names),'rows':rows,'shots':total,'denominators':sorted(denominators),'metadata':json.loads(z.read('wyniki/data.json')),
          'match':'Belem May 2023 is not Nairobi December 2022 or October 2023; no shared-gate acquisition bridge'}
    with zipfile.ZipFile(cache/'leakage.zip') as z:
        names=z.namelist();out['leakage_control']={'files':len([n for n in names if not n.endswith('/')]),'csv_files':sum(n.endswith('.csv') for n in names),
          'rb_readme':z.read('v2/Fig3/Fig3d/readme.txt').decode().strip(),'leakage_readme':z.read('v2/Fig3/Fig3e/readme.txt').decode().strip(),
          'status':'measured probability/fit arrays and explicitly named simulation arrays; not transferred to target ion acquisition'}
    for file in ('h1-spam.json','h1-rb.json'):
        x=json.loads((cache/file).read_text());N=x['shots'];raw=x['raw_data'];bad=[];examples=[]
        for key,data in raw.items():
            if any(len(v)!=N for v in data.values()):raise ValueError('wrong shot denominator')
            nums=tuple(map(int,re.search(r'\((\d+), (\d+)\)',key).groups()))
            family='SPAM' if file=='h1-spam.json' else 'SQ_RB'
            exp=x['expected_output'][f'{family}: ({nums[0]}, {nums[1]})']
            for q in range(10):
                # Deposited c/l arrays store bit strings in c[9]..c[0] order.
                observed=sum(s[9-q]==exp[str(q)] for s in data['c'])
                saved=x['survival'][str(q)][str(nums[1])] if family=='SPAM' else x['survival'][str(q)][str(nums[0])][str(nums[1])]
                if observed!=saved:bad.append([key,q,'survival',observed,saved])
                if family=='SQ_RB':
                    unflagged=sum(s[9-q]=='0' for s in data['l']);saved_l=x['leakage_postselect'][str(q)][str(nums[0])][str(nums[1])]
                    if unflagged!=saved_l:bad.append([key,q,'leakage_postselect',unflagged,saved_l])
                if q==0 and len(examples)<2:examples.append({'key':key,'q':q,'correct':observed,'shots':N})
        if bad:raise ValueError('Quantinuum schema/count mismatch '+str(bad[:3]))
        out[file]={'keys':list(x),'circuits':len(raw),'shots_per_circuit':N,'logical_qubits':10,'examples':examples,
          'schema_note':'leakage_postselect stores unflagged counts in this file, NOT leakage probabilities; QASM contains an ancilla gadget',
          'missing':'per-circuit timestamps, persistent physical ion identifiers, matched calibration/gadget response records'}
    return out


if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--cache',required=True);a.add_argument('--download',action='store_true');a.add_argument('--check',action='store_true');args=a.parse_args();out=audit(args.cache,args.download);p=HERE/'results/candidates.json'
    if args.check:
        if out!=json.loads(p.read_text()):raise ValueError('candidate audit mismatch')
        print('Candidate audits reproduce')
    else:p.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
