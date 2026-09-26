"""Reconstruct count tables and report assumption-indexed bounds for one NIST run."""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np

from artifacts import write_json
from inference import score_interval
from reconstruct import VARIANTS, sha256

HERE=Path(__file__).resolve().parent


def counts(a,b,ac,bc):
    arrays=[np.asarray(x) for x in (a,b,ac,bc)]
    if any(x.ndim != 1 or len(x)!=len(arrays[0]) for x in arrays):
        raise ValueError('unaligned record arrays')
    if any(np.any((x<0)|(x>3)) for x in arrays[:2]) or any(np.any((x!=0)&(x!=1)) for x in arrays[2:]):
        raise ValueError('unexpected setting or outcome code')
    a,b,ac,bc=arrays
    codes=(((a.astype('i8')*4+b)*2+ac)*2+bc)
    return np.bincount(codes,minlength=64).reshape(4,4,2,2)


def jump_ranges(reconstruction,n):
    ranges=[]
    for side in ('alice','bob'):
        jumps=reconstruction[side]['timestamp_jumps']
        if len(jumps)%2:
            raise ValueError('unpaired timestamp jump')
        for (start,d0),(end,d1) in zip(jumps[::2],jumps[1::2]):
            if not (start<end and d0<0<d1 and d0+d1==1600):
                raise ValueError('unrecognized paired timestamp excursion')
            lo,hi=max(0,start),min(n,end+2)
            if hi>lo:ranges.append([lo,hi])
    return ranges


def apply_patches(old,start,patches):
    result=old.copy()
    if not len(patches):return result
    indices=patches[:,0]
    lo,hi=np.searchsorted(indices,[start,start+len(old)])
    selected=patches[lo:hi]
    local=selected[:,0]-start
    if np.any(old[local]!=selected[:,1]):
        raise ValueError('patch base does not match archived words')
    result[local]=selected[:,2]
    return result


def reconstruct_counts(hdf_path, reconstruction, protocol):
    slots=[protocol['primary_slots_zero_based']]+protocol['sensitivity_slots_zero_based']
    masks=[sum(1<<i for i in group) for group in slots]
    nb=protocol['chronological_blocks']
    with h5py.File(hdf_path) as f:
        n=len(f['alice/settings'])
        for side in ('alice','bob'):
            if len(f[side+'/settings'])!=n or len(f[side+'/clicks'])!=n:
                raise ValueError('different stream lengths')
            r=reconstruction[side]
            if r['covered_rows']!=n or any(r['mismatches'].values()):
                raise ValueError('unreconciled reconstruction')
        ranges=jump_ranges(reconstruction,n)
        cuts=np.linspace(0,n,nb+1,dtype='i8')
        tables=np.zeros((3,len(masks),nb,4,4,2,2),dtype='i8')
        trusted=np.zeros_like(tables)
        stored=np.zeros((4,4,2,2),dtype='i8')
        patches={}
        for side in ('alice','bob'):
            for v in VARIANTS:
                p=np.array(reconstruction[side]['patches'][v],dtype='i8').reshape(-1,3)
                if len(p) and (np.any(np.diff(p[:,0])<=0) or np.any(p[:,0]<0) or np.any(p[:,0]>=n)
                               or np.any(p[:,1:]<0) or np.any(p[:,1:]>65535)):
                    raise ValueError('invalid patch list')
                patches[side,v]=p
        for block,(begin,end) in enumerate(zip(cuts[:-1],cuts[1:])):
            for start in range(int(begin),int(end),1_000_000):
                stop=min(start+1_000_000,int(end))
                a=f['alice/settings'][start:stop];b=f['bob/settings'][start:stop]
                olda=f['alice/clicks'][start:stop];oldb=f['bob/clicks'][start:stop]
                stored+=counts(a,b,(olda&masks[0])>0,(oldb&masks[0])>0)
                secure=np.ones(len(a),dtype=bool)
                for lo,hi in ranges:
                    secure[max(0,lo-start):max(0,min(len(a),hi-start))]=False
                for vi,v in enumerate(VARIANTS):
                    ca=apply_patches(olda,start,patches['alice',v]);cb=apply_patches(oldb,start,patches['bob',v])
                    for wi,mask in enumerate(masks):
                        oa=(ca&mask)>0;ob=(cb&mask)>0
                        tables[vi,wi,block]+=counts(a,b,oa,ob)
                        trusted[vi,wi,block]+=counts(a[secure],b[secure],oa[secure],ob[secure])
    return dict(schema_version=1,source_sha256=sha256(hdf_path),
                protocol=protocol,n=n,block_edges=cuts.tolist(),slot_groups=slots,
                phase_variants=list(VARIANTS),timestamp_uncertain_ranges=ranges,
                stored_primary=stored.tolist(),tables=tables.tolist(),trusted_tables=trusted.tolist())


def contrast(table, trusted, n, direction, receiver_setting):
    """Descriptive complete-setting arm rates and full-row bounded score."""
    table=np.asarray(table,dtype='i8');trusted=np.asarray(trusted,dtype='i8')
    if direction=='B_to_A':
        table=table.transpose(1,0,3,2);trusted=trusted.transpose(1,0,3,2)
    elif direction!='A_to_B':
        raise ValueError('unknown direction')
    r=receiver_setting+1
    arms=[]
    for s in (1,2):
        arm=table[s,r]
        den=int(arm.sum());clicks=int(arm[:,1].sum())
        arms.append(dict(n=den,clicks=clicks,rate=None if not den else clicks/den))
    gap=None if any(x['rate'] is None for x in arms) else arms[1]['rate']-arms[0]['rate']
    total=0
    for s in (1,2):
        sign=2*(s-1)-1
        total+=2*sign*int(trusted[s,r,:,1].sum()-trusted[s,r,:,0].sum())
    # Both settings must be unambiguous. Unknown rows remain in the denominator.
    known=int(trusted[1:3,1:3].sum())
    return dict(direction=direction,receiver_setting=receiver_setting,arms=arms,
                descriptive_gap=gap,n=int(n),known_rows=known,unknown_rows=int(n-known),score_sum=total)


def selected_contrast(table,direction,receiver):
    table=np.asarray(table,dtype='i8')
    if direction=='B_to_A':table=table.transpose(1,0,3,2)
    arms=[]
    for sender in (1,2):
        row=table[sender,receiver+1,1,:]
        n=int(row.sum());k=int(row[1])
        arms.append(dict(n=n,clicks=k,rate=None if n==0 else k/n))
    gap=None if any(a['rate'] is None for a in arms) else arms[1]['rate']-arms[0]['rate']
    return dict(direction=direction,receiver_setting=receiver,selection='sender clicked',
                arms=arms,descriptive_selected_gap=gap,causal_interpretation=False)


def claim_status(ledger):
    requirements=ledger['requirements']
    causal=[k for k in ledger['causal_required'] if requirements[k]['status']!='established']
    spacelike=causal+[k for k in ledger['spacelike_additional_required'] if requirements[k]['status']!='established']
    return dict(causal_claim_enabled=not causal,spacelike_claim_enabled=not spacelike,
                causal_blockers=causal,spacelike_blockers=spacelike)


def summarize(bundle):
    p=bundle['protocol'];tables=np.array(bundle['tables'],dtype='i8');secure=np.array(bundle['trusted_tables'],dtype='i8')
    edges=bundle['block_edges'];nb=len(edges)-1
    if tables.shape!=(3,3,nb,4,4,2,2) or secure.shape!=tables.shape or np.any(secure>tables) or np.any(secure<0):
        raise ValueError('unexpected count schema')
    sizes=np.diff(edges)
    if np.any(tables.sum(axis=(3,4,5,6))!=sizes[None,None,:]):
        raise ValueError('table totals do not reconstruct all rows')
    comparisons=3*3*(nb+1)*4
    rows=[]
    for vi,v in enumerate(VARIANTS):
        for wi,slots in enumerate(bundle['slot_groups']):
            for block in range(-1,nb):
                t=tables[vi,wi].sum(axis=0) if block==-1 else tables[vi,wi,block]
                st=secure[vi,wi].sum(axis=0) if block==-1 else secure[vi,wi,block]
                n=bundle['n'] if block==-1 else int(sizes[block])
                for direction in ('A_to_B','B_to_A'):
                    for receiver in (0,1):
                        row=contrast(t,st,n,direction,receiver)
                        row.update(phase=v,slots=slots,block='all' if block==-1 else block,
                                   start=0 if block==-1 else int(edges[block]))
                        row['ideal_assignment_interval']=score_interval(row['score_sum'],row['unknown_rows'],n,
                            alpha=p['family_alpha'],comparisons=comparisons,setting_tv=0,leakage=0,start=row['start'])
                        rows.append(row)
    primary=[r for r in rows if r['phase']=='nominal' and r['slots']==p['primary_slots_zero_based'] and r['block']=='all']
    sensitivity=[]
    for row in primary:
        for eps in p['joint_setting_tv_scenarios']:
            for leakage in p['ordinary_leakage_gap_scenarios']:
                sensitivity.append(dict(direction=row['direction'],receiver_setting=row['receiver_setting'],
                    assumed_joint_setting_tv=eps,assumed_leakage_gap=leakage,
                    interval=score_interval(row['score_sum'],row['unknown_rows'],row['n'],
                      alpha=p['family_alpha'],comparisons=comparisons,setting_tv=eps,leakage=leakage)))
    return dict(schema_version=1,n=bundle['n'],comparisons=comparisons,family_alpha=p['family_alpha'],
                **claim_status(json.loads((HERE/'sufficiency.json').read_text())),
                interpretation='Retrospective recorded-row contrasts and assumption-indexed signed-effect intervals; no certified spacelike causal bound.',
                rows=rows,primary_sensitivity=sensitivity,
                selection_diagnostics=[selected_contrast(tables[0,0].sum(axis=0),direction,r)
                                       for direction in ('A_to_B','B_to_A') for r in (0,1)])


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--hdf5',type=Path)
    p.add_argument('--reconstruction',type=Path)
    p.add_argument('--output-dir',type=Path)
    p.add_argument('--check',action='store_true')
    a=p.parse_args()
    if a.check:
        bundle=json.loads((HERE/'results/counts.json').read_text())
        if bundle['protocol'] != json.loads((HERE/'protocol.json').read_text()):
            raise SystemExit('protocol changed without source regeneration')
        result=summarize(bundle)
        if result!=json.loads((HERE/'results/analysis.json').read_text()):
            raise SystemExit('analysis snapshot differs')
        print('396 simultaneous assumption-indexed intervals reproduced from committed counts; raw sources not rechecked')
        return
    if not all((a.hdf5,a.reconstruction,a.output_dir)):
        p.error('provide --hdf5, --reconstruction, and --output-dir')
    manifest=json.loads((HERE/'manifest.json').read_text())
    if sha256(a.hdf5)!=manifest['files']['hdf5']['sha256']:
        raise SystemExit('HDF5 checksum mismatch')
    reconstruction=json.loads(a.reconstruction.read_text())
    for side in ('alice','bob'):
        if reconstruction[side]['archive_sha256']!=manifest['files'][side]['sha256']:
            raise SystemExit('raw archive provenance mismatch')
    protocol=json.loads((HERE/'protocol.json').read_text())
    bundle=reconstruct_counts(a.hdf5,reconstruction,protocol)
    a.output_dir.mkdir(parents=True,exist_ok=True)
    for name,obj in [('counts',bundle),('analysis',summarize(bundle))]:
        write_json(a.output_dir/(name+'.json'),obj)
    print(f'{bundle["n"]:,} paired rows counted, including no-click and ambiguous-setting rows')


if __name__=='__main__':
    main()
