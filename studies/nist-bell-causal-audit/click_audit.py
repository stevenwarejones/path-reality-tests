"""Retrospective click-only analysis, retaining the centered-score reference."""
import argparse
import json
from pathlib import Path
import numpy as np
from artifacts import write_json
from click_inference import effect_interval, LAMBDAS, EVENT_LABELS
from audit import claim_status, contrast

HERE=Path(__file__).resolve().parent


def possible_context(hist, sender, receiver):
    # Codes 0/3 do not specify a unique binary setting. Both latent settings remain possible.
    return int(np.asarray(hist)[np.ix_([0,sender+1,3],[0,receiver+1,3])].sum())


def summarize(bundle):
    tables=np.asarray(bundle['tables']); trusted=np.asarray(bundle['trusted_tables'])
    unknown=np.asarray(bundle['unknown_contexts']); events=np.asarray(bundle['uncertain_event_contexts'])
    nb=len(bundle['block_edges'])-1
    if unknown.shape!=(nb,4,4) or events.shape!=(nb,2,4,4) or np.any(events<0) or np.any(events>unknown[:,None]):
        raise ValueError('invalid uncertain-event counts')
    rows=[]
    for vi,phase in enumerate(bundle['phase_variants']):
        for wi,slots in enumerate(bundle['slot_groups']):
            for block in range(-1,nb):
                def part(array): return array.sum(axis=0) if block==-1 else array[block]
                t=part(tables[vi,wi]);st=part(trusted[vi,wi]);u=part(unknown);e=part(events)
                start=0 if block==-1 else bundle['block_edges'][block]
                n=bundle['n'] if block==-1 else bundle['block_edges'][block+1]-start
                for direction,si in (('A_to_B',1),('B_to_A',0)):
                    oriented=st if si==1 else st.transpose(1,0,3,2)
                    uh=u if si==1 else u.T;eh=e[si] if si==1 else e[si].T
                    for r in (0,1):
                        row=contrast(t,st,n,direction,r)
                        known=[int(oriented[x+1,r+1,:,1].sum()) for x in (0,1)]
                        possible=[possible_context(eh,x,r) for x in (0,1)]
                        arbitrary=[possible_context(uh,x,r) for x in (0,1)]
                        supported=list(zip(known,[k+m for k,m in zip(known,possible)]))
                        unrestricted=list(zip(known,[k+m for k,m in zip(known,arbitrary)]))
                        row.update(phase=phase,slots=slots,block='all' if block==-1 else block,start=start,
                            known_click_counts=known,possible_unknown_clicks=possible,
                            unrestricted_unknown_clicks=arbitrary,click_only_known_mean=4*(known[1]-known[0])/n,
                            event_supported_interval=effect_interval(supported,n,start=start,alpha=bundle['protocol']['family_alpha']),
                            unrestricted_interval=effect_interval(unrestricted,n,start=start,alpha=bundle['protocol']['family_alpha']))
                        rows.append(row)
    primary=[r for r in rows if r['phase']=='nominal' and r['slots']==bundle['protocol']['primary_slots_zero_based'] and r['block']=='all']
    sensitivity=[]
    for row in primary:
        arms=list(zip(row['known_click_counts'],[k+m for k,m in zip(row['known_click_counts'],row['possible_unknown_clicks'])]))
        for eps in bundle['protocol']['joint_setting_tv_scenarios']:
            for leak in bundle['protocol']['click_only_leakage_gap_scenarios']:
                sensitivity.append(dict(direction=row['direction'],receiver_setting=row['receiver_setting'],
                    interval=effect_interval(arms,row['n'],alpha=bundle['protocol']['family_alpha'],per_history_tv=eps,leakage=leak)))
    return dict(schema_version=1,analysis='retrospective click-only v2',n=bundle['n'],family_alpha=bundle['protocol']['family_alpha'],
        lambda_grid=list(LAMBDAS),event_labels=EVENT_LABELS,
        **claim_status(json.loads((HERE/'sufficiency.json').read_text())),rows=rows,primary_sensitivity=sensitivity)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check',action='store_true');parser.add_argument('--directory',type=Path,default=HERE/'results')
    args=parser.parse_args();result=summarize(json.loads((args.directory/'counts.json').read_text()))
    if args.check:
        if result!=json.loads((args.directory/'click-analysis.json').read_text()):raise SystemExit('click analysis snapshot differs')
        print('396 click-only intervals reproduced; physical premises remain conditional')
    else:write_json(args.directory/'click-analysis.json',result)


if __name__=='__main__':main()
