#!/usr/bin/env python3
"""Stream the audited decoder into fine aggregates and an external lossless cache."""
import argparse
import json
from pathlib import Path
import zipfile
import h5py
import numpy as np
from common import (HERE, ROOT, AUDIT, PREVIOUS, decoder, PEAKS, OFFSETS,
                    PHASE_CELLS, CELLS, BLOCKS, MAPS, grouped, json_text)
CACHE_DTYPE = np.dtype([('row','<i8'),('pulse_bit','u1'),('phase','<f8')])


def select(events, config, lo, hi):
    """First eligible event in raw record order, including tied timestamps."""
    er = events['row']
    eligible = ((er >= lo) & (er < hi) &
                np.isin(events['pulse']-config['bitoffset'], (4,5,6)))
    indices = np.flatnonzero(eligible)
    rows, first, multiplicity = np.unique(er[indices], return_index=True, return_counts=True)
    chosen = indices[first]
    phase = events['phase'][chosen]
    if np.any(~np.isfinite(phase) | (phase < 0) | (phase >= PHASE_CELLS)):
        raise ValueError('phase outside frozen enclosing domain')
    bits = events['pulse'][chosen]-config['bitoffset']
    cells = (bits-4)*PHASE_CELLS+np.floor(phase).astype('i8')
    coarse = np.where(phase < config['pk'],0,np.where(phase < config['pk']+2,1,2))
    from_cells = np.where(cells % PHASE_CELLS < config['pk'],0,
                         np.where(cells % PHASE_CELLS < config['pk']+2,1,2))
    if not np.array_equal(coarse,from_cells):
        raise ValueError('coarse map not deterministic')
    return rows, cells, bits, phase, int(eligible.sum()), int((multiplicity>1).sum())


def summarize(rows, local, remote, bad, cell, detector, n):
    """Independent of streaming chunk boundaries; no-event is cell=-1."""
    counts = np.zeros((2,2,2,CELLS), dtype='i8')
    exposure = np.zeros((2,2,2), dtype='i8')
    unknown = np.zeros((2,2,2,2), dtype='i8')
    blocks = np.zeros((BLOCKS,2,CELLS), dtype='i8')
    block_exposure = np.zeros((BLOCKS,2),dtype='i8')
    interior = (rows >= 64) & (rows < n-64)
    bad = bad | ~np.isin(local,(1,2)) | ~np.isin(remote,(1,2))
    trusted = interior & ~bad
    half = (rows >= n//2).astype('i8')
    block = np.clip((rows-64)*BLOCKS//(n-128),0,BLOCKS-1)
    context = 4*half+2*(remote.astype('i8')-1)+local.astype('i8')-1
    exposure[:] = np.bincount(context[trusted],minlength=8).reshape(2,2,2)
    event = trusted & (cell >= 0)
    counts[:] = np.bincount(context[event]*CELLS+cell[event],minlength=8*CELLS).reshape(2,2,2,CELLS)
    block_exposure[:] = np.bincount(2*block[trusted]+local[trusted]-1,minlength=BLOCKS*2).reshape(BLOCKS,2)
    blocks[:] = np.bincount((2*block[event]+local[event]-1)*CELLS+cell[event],minlength=BLOCKS*2*CELLS).reshape(BLOCKS,2,CELLS)
    for hh in (0,1):
        for x in (0,1):
            for y in (0,1):
                compatible = (interior & bad & (half==hh) &
                    ((remote==x+1) | ~np.isin(remote,(1,2))) &
                    ((local==y+1) | ~np.isin(local,(1,2))))
                unknown[hh,x,y] = [compatible.sum(),(compatible & detector).sum()]
    return counts,exposure,unknown,blocks,block_exposure


def compare_coarse(result, baseline):
    fine = np.array(result['counts'])
    categories = grouped(fine, MAPS[result['side']][0])
    features = np.concatenate([fine.sum(axis=-1,keepdims=True),categories],axis=-1)
    for got, want in [(features, baseline['groups']['counts'][0]),
                      (result['exposures'],baseline['groups']['exposures'][0]),
                      (result['unknown'],baseline['groups']['unknown'][0])]:
        if not np.array_equal(got,want):
            raise ValueError('fine-to-coarse baseline reproduction failed')


def extract(hdf_path, archive_path, side, cache_dir, chunk_records=1_000_000):
    cache_dir=Path(cache_dir).resolve()
    if cache_dir.is_relative_to(ROOT):
        raise ValueError('cache must be outside repository')
    cache_dir.mkdir(parents=True,exist_ok=True)
    manifest=json.loads((AUDIT/'manifest.json').read_text())
    for key,path in [('hdf5',hdf_path),(side,archive_path)]:
        spec=manifest['files'][key]
        if Path(path).stat().st_size!=spec['bytes'] or decoder.sha256(path)!=spec['sha256']:
            raise ValueError(key+': source hash mismatch')
    totals=[np.zeros(s,dtype='i8') for s in [(2,2,2,CELLS),(2,2,2),(2,2,2,2),(BLOCKS,2,CELLS),(BLOCKS,2)]]
    base=covered=selected=multi=0
    previous=None
    jumps=[]; caches=[]; tail_info=None
    other='bob' if side=='alice' else 'alice'
    with h5py.File(hdf_path) as h, zipfile.ZipFile(archive_path) as z:
        members=[m for m in z.infolist() if not m.is_dir()]
        if len(members)!=1 or members[0].file_size%11:
            raise ValueError('unexpected raw layout')
        member=members[0].filename
        n=len(h[side+'/settings'])
        if n<=128 or len(h[other+'/settings'])!=n:
            raise ValueError('invalid populations')
        config={k:int(h[f'config/{side}/{k}'][()]) for k in ('pk','radius','bitoffset')}
        if (config['pk'],config['bitoffset'])!=(PEAKS[side],OFFSETS[side]):
            raise ValueError('changed frozen configuration')
        offset=decoder.discover_offset(z,member,h[side+'/settings'][:64])
        expected=h[side+'/badSyncInfo'][[0,3],:].T.astype('i8').tolist()
        ranges=sum([decoder.excursion_ranges(h[s+'/badSyncInfo'][[0,3],:].T.astype('i8').tolist(),n) for s in (side,other)],[])
        with z.open(member) as stream:
            for a,next_tag,tail in decoder.closed_blocks(stream,chunk_records):
                if tail:
                    tail_info=dict(records=len(a),syncs=int((a['ch']==6).sum()))
                    continue
                setting,_,legacy,previous,bad,events=decoder.decode_block(a,next_tag,previous,config,include_events=True)
                jumps.extend([[base+i-offset,d] for i,d in bad])
                lo,hi=max(base,offset),min(base+len(setting),offset+n)
                if hi>lo:
                    start,stop=lo-offset,hi-offset
                    local=h[side+'/settings'][start:stop]; remote=h[other+'/settings'][start:stop]
                    if not np.array_equal(setting[lo-base:hi-base],local) or not np.array_equal(legacy[lo-base:hi-base],h[side+'/clicks'][start:stop]):
                        raise ValueError('raw/stored mismatch')
                    rows=np.arange(start,stop,dtype='i8')
                    clock_bad=np.zeros(len(rows),bool)
                    for left,right in ranges:
                        clock_bad |= (rows>=left)&(rows<right)
                    detector=np.zeros(len(rows),bool)
                    er=events['row']+base-offset
                    inside=(er>=start)&(er<stop)
                    detector[er[inside]-start]=True
                    selected_rows,cells,bits,phase,k,m=select(events,config,lo-base,hi-base)
                    event_rows=selected_rows+base-offset
                    cell=np.full(len(rows),-1,dtype='i8'); cell[event_rows-start]=cells
                    record=np.empty(len(event_rows),CACHE_DTYPE)
                    record['row']=event_rows; record['pulse_bit']=bits; record['phase']=phase
                    caches.append(record)
                    selected+=k; multi+=m; covered+=len(rows)
                    for target,values in zip(totals,summarize(rows,local,remote,clock_bad,cell,detector,n)):
                        target+=values
                base+=len(setting)
    if covered!=n or jumps!=expected:
        raise ValueError('incomplete reconciliation')
    cache=np.concatenate(caches); path=cache_dir/f'{side}-fine-events.npy'; np.save(path,cache)
    result=dict(schema_version=1,side=side,rows=n,interior_rows=n-128,
        counts=totals[0].tolist(),exposures=totals[1].tolist(),unknown=totals[2].tolist(),
        pooled_block_counts=totals[3].tolist(),pooled_block_exposures=totals[4].tolist(),
        provenance=dict(hdf_sha256=manifest['files']['hdf5']['sha256'],archive_sha256=manifest['files'][side]['sha256'],
            cache_sha256=decoder.sha256(path),overlap_verified=covered,raw_prefix=offset,raw_suffix=base-offset-n,
            selected_events=selected,first_event_rows=len(cache),multiclick_rows=multi,censored_tail=tail_info,config=config,
            phase_min=float(cache['phase'].min()),phase_max=float(cache['phase'].max()),
            retained_pair='pulse bit and original float64 decoder phase',phase_radius_cut=False,
            protocol_sha256=decoder.sha256(HERE/'protocol.md')),
        axes=dict(counts=['half','remote','local','cell'],unknown=['half','remote','local','envelope'],
                  pooled_blocks=['chronological block','local','cell'],envelopes=['unrestricted','event_supported']))
    compare_coarse(result,json.loads((PREVIOUS/f'{side}-counts.json').read_text()))
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('hdf5','archive','cache-dir','output'): p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--side',choices=['alice','bob'],required=True)
    p.add_argument('--chunk-records',type=int,default=1_000_000)
    p.add_argument('--check',action='store_true'); a=p.parse_args()
    r=extract(a.hdf5,a.archive,a.side,a.cache_dir,a.chunk_records)
    # One aggregate array per line keeps the committed representation manageable.
    content=json.dumps(r,sort_keys=True,separators=(',',':'))+'\n'
    if a.check:
        if a.output.read_text()!=content: raise SystemExit('fine aggregates differ')
    else: a.output.write_text(content)
    print(a.side+': fine aggregates and exact coarse reproduction '+('verified' if a.check else 'written'))
if __name__=='__main__': main()
