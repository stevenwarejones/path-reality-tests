#!/usr/bin/env python3
"""Execute only the frozen fine-resolution family; all assumptions are retained."""
import argparse
import json
import numpy as np
from common import HERE, PREVIOUS, SIDES, LEVELS, EPSILONS, ENVELOPES, MAPS, CDF_LABELS, PARTITION_LABELS, json_text, artifact_json
from extract import compare_coarse
from inference import cdf,partitions


def native(x):
    if isinstance(x,np.ndarray): return x.tolist()
    if isinstance(x,dict): return {k:native(v) for k,v in x.items()}
    return x


def run():
    cases=[]; baselines=[]; sources={}
    old=json.loads((PREVIOUS/'results.json').read_text())
    for side in SIDES:
        source=json.loads((HERE/f'{side}-counts.json').read_text())
        compare_coarse(source,json.loads((PREVIOUS/f'{side}-counts.json').read_text()))
        sources[side]=source['provenance']
        counts=np.array(source['counts']).sum(axis=0)
        unknown=np.array(source['unknown']).sum(axis=0)
        for y in (0,1):
            for epsilon in EPSILONS:
                for e,envelope in enumerate(ENVELOPES):
                    entry=dict(receiver=side,local_setting=y,epsilon=epsilon,envelope=envelope)
                    try:
                        a=cdf(counts[:,y],unknown[:,y,e],source['interior_rows'],epsilon)
                        b=partitions(counts[:,y],unknown[:,y,e],source['interior_rows'],side,epsilon)
                        entry.update(status='compatible',cdf=native(a),partitions=native(b))
                    except ValueError as error:
                        entry.update(status='incompatible',reason=str(error))
                    cases.append(entry)
        # Historical bound is reproduced from unchanged #13 artifacts, not charged to this new family.
        for v in old['strata']:
            if v['receiver']==side and v['group']=='pooled' and v['period']=='full':
                baselines.append(v)
    return dict(schema_version=1,cdf_labels=CDF_LABELS,partition_labels=PARTITION_LABELS,
                error_budget=.01,levels=list(LEVELS),sources=sources,cases=cases,historical_baselines=baselines)


def fmt(v): return '['+', '.join(f'{z*1e6:.1f}' for z in v)+']'


def report(r):
    chosen=[v for v in r['cases'] if v['epsilon']==0 and v['envelope']=='event_supported']
    okay=[v for v in r['cases'] if v['status']=='compatible']
    detections=sum(v['cdf']['k_grid'][0]>0 or any(t[0]>0 for t in v['partitions']['tv']) for v in okay)
    gains=sum(any(t[0]>0 for t in v['partitions']['gain'][1:]) for v in okay)
    lines=['# What finer archived timing buys', '',
        f'**{detections}/{len(okay)} assumption-specific cases detect a difference; {gains} resolve positive information loss J.**',
        f'{len(r["cases"])-len(okay)} confidence/model intersections are incompatible (flagged, not counted as discoveries).',
        'These are overlapping sensitivity cases, not independent replications. A null result is not equality.', '',
        'Both receivers use the same 107,109,468 interior rows, including no-events.',
        'All source rows reconcile and the fine extraction exactly reproduces the historical',
        'coarse counts, exposures and uncertainty totals. No new experiment or physical calibration.', '',
        '## Measured comparison', '',
        'All bounds below are simultaneous full-trial ppm, epsilon=0, event-supported',
        'completion. Other epsilon values (.001,.01) and unrestricted completions are',
        'in results.json. CDF and partition families jointly spend .01 for this extension.', '',
        '| Receiver/local setting | Representation | Outcome resolution | Effect bound (ppm) | Added information versus coarse |',
        '|---|---|---|---|---|']
    for v in chosen:
        label=f'{v["receiver"]}/{v["local_setting"]}'
        if v['status']!='compatible':
            lines.append(f'| {label} | all | — | incompatible | — |'); continue
        historical=next(b for b in r['historical_baselines'] if b['receiver']==v['receiver'] and b['local_setting']==v['local_setting'] and b['epsilon']==0 and b['envelope']=='event_supported')
        lines.append(f'| {label} | Original #13 | 4 outcomes | TV {fmt(historical["weighted_tv"])} | Historical reference; original error family |')
        lines.append(f'| {label} | Cumulative | 486 fixed thresholds | K_grid {fmt(v["cdf"]["k_grid"])} | Pulse and within-category cumulative contrasts |')
        lines.append(f'| {label} | Cumulative interpolation | All retained phases | K_full {fmt(v["cdf"]["k_full"])} | Conservative monotone bridge between thresholds |')
        for j,level in enumerate(LEVELS):
            bins=int(MAPS[v['receiver']][j].max())+2
            lines.append(f'| {label} | {level} | {bins} outcomes | TV {fmt(v["partitions"]["tv"][j])} | J {fmt(v["partitions"]["gain"][j])} |')
        lines.append(f'| {label} | Unrestricted retained pair | Pulse + float64 phase | TV_full {fmt(v["partitions"]["full_tv"])} | Lower bound from width1; incidence upper bound |')
    lines += ['', 'The extension coarse row uses the new family and nested-law constraints; the',
        'historical row uses #13 unchanged. Different confidence bounds are not changes',
        'in the underlying coarse estimand. A refined true TV cannot decrease, while its',
        'estimated upper bound may widen. Tree consistency and the event-incidence cap',
        'prevent adding hundreds of empty-cell penalties into a vacuous TV<=1 result.', '',
        '## What remains unresolved', '',
        'A positive J lower bound would show dependence lost by the specified coarse map.',
        'It would not identify its causal origin. Current-intervention interpretations',
        'require the declared per-history joint-setting probabilities and trial correspondence.',
        'The grid preserves pulse identity and one-tag phase cells, but not all decoder',
        'phase digits. K_full interpolation and the full-pair TV cap are conservative',
        'enclosures; they are not full-resolution likelihood fits or analog-delay limits.',
        'All upper bounds apply to time-averaged laws, not instantaneous effects.', '',
        'See [recovery report](recovery-report.md) for paired empirical-background power,',
        '[conclusions](conclusions.md) for the five research answers, and',
        '[statistics](statistics.md) for coverage, consistency and uncertainty assumptions.', '',
        '![CDF bands](figures/cdf-bands.svg)', '',
        '![Resolution comparison](figures/resolution.svg)', '']
    return '\n'.join(lines)


def main():
    p=argparse.ArgumentParser(); p.add_argument('--check',action='store_true'); a=p.parse_args()
    r=run()
    # Keep each sensitivity case on one line while preserving full bands.
    content=artifact_json(r)
    for name,text in [('results.json',content),('report.md',report(r))]:
        if a.check:
            if (HERE/name).read_text()!=text: raise SystemExit(name+' differs')
        else: (HERE/name).write_text(text)
    print('fine analysis '+('verified' if a.check else 'written'))
if __name__=='__main__': main()
