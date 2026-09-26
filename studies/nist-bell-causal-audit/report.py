"""Render the checked numerical report and a deterministic scientific figure."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

HERE=Path(__file__).resolve().parent


def render(output):
    result=json.loads((HERE/'results/analysis.json').read_text())
    counts=json.loads((HERE/'results/counts.json').read_text())
    raw=json.loads((HERE/'results/reconstruction.json').read_text())
    primary=[r for r in result['rows'] if r['phase']=='nominal' and r['slots']==[4,5,6] and r['block']=='all']
    lines=['# Recorded marginals and conditional sensitivity','',
           'The source contains 107,109,596 aligned stored rows. Both raw setting streams',
           'and the archived click-update behavior reconcile exactly over that interval.',
           'The analysis below uses bitwise-OR reconstructed clicks, includes no-click rows,',
           'and keeps ambiguous settings in the population. No spacelike causal claim is enabled.','',
           '## Primary three-pulse record','',
           'The marginal differences describe valid-setting rows. The intervals instead use',
           'the full-row score, include worst-case completion of 19,968 uncertain rows,',
           'and require the assignment, record and causal premises in [statistics.md](../statistics.md).',
           'Here joint-setting TV allowance and ordinary-leakage allowance are **assumed zero**.',
           'Bounds are on absolute average signed effects, not average absolute influences.','',
           '| Direction | Receiver setting | Descriptive difference (per million) | Conditional absolute upper limit |',
           '|---|---:|---:|---:|']
    for r in primary:
        lines.append(f"| {r['direction'].replace('_to_',' → ')} | {r['receiver_setting']} | {r['descriptive_gap']*1e6:.3f} | {100*r['ideal_assignment_interval']['absolute_upper']:.4f}% |")
    lines+=['','The simultaneous interval construction allocates total error 0.01 across the',
            'declared family and every contiguous interval; it tolerates arbitrary device',
            'memory under the stated conditional assignment premise. Every reported interval',
            'contains zero. This does not validate the premise or certify no signaling.','',
            '## Assumed calibration sensitivity','',
            'Largest primary upper limit across both directions and receiver settings:','',
            '| Assumed joint-setting TV | Assumed leakage gap | Conditional upper limit |',
            '|---:|---:|---:|']
    for eps in counts['protocol']['joint_setting_tv_scenarios']:
        for leakage in counts['protocol']['ordinary_leakage_gap_scenarios']:
            u=max(r['interval']['absolute_upper'] for r in result['primary_sensitivity'] if r['assumed_joint_setting_tv']==eps and r['assumed_leakage_gap']==leakage)
            lines.append(f'| {eps:g} | {leakage:g} | {100*u:.4f}% |')
    lines+=['','These allowances are sensitivity parameters, not measured calibration results.',
            'Lack of a valid history-conditional calibration prevents promoting these',
            'numbers to an apparatus-certified causal bound.','',
            '## Reconstruction and window sensitivity','',
            '| Side | Nominal words corrected | Narrow-radius words changed | Wide-radius words changed |',
            '|---|---:|---:|---:|']
    for side,r in raw.items():
        v=r['changed_words'];lines.append(f"| {side} | {v['nominal']} | {v['narrow']} | {v['wide']} |")
    lines+=['','## Selection diagnostic','',
            'Conditioning on a sender click yields the following receiver differences.',
            'These are selected correlations, not causal effects; the primary analysis',
            'does not make this selection. Their size illustrates why complete trials matter.','',
            '| Direction | Receiver setting | Sender-click-selected difference |',
            '|---|---:|---:|']
    for r in result['selection_diagnostics']:
        lines.append(f"| {r['direction'].replace('_to_',' → ')} | {r['receiver_setting']} | {r['descriptive_selected_gap']:.4f} |")
    old=np.array(counts['stored_primary']);new=np.array(counts['tables'])[0,0].sum(axis=0)
    lines+=['','Correcting repeated-index click loss adds the following primary-window',
            'receiver clicks across all setting codes (counts, not inferred effects):','',
            f"- Alice: {int(new[:,: ,1,:].sum()-old[:,:,1,:].sum())}.",
            f"- Bob: {int(new[:,:,:,1].sum()-old[:,:,:,1].sum())}.",
            '', 'The original diagnostic table is reproduced after assigning the two A=3',
            'no-click rows to A=2. The causal audit does not adopt that ambiguous-setting',
            'assignment. This establishes count reconciliation, not a reproduction of the',
            'published Bell-test p-value or a challenge to that result.','',
            '![Drift and conditional sensitivity](diagnostics.svg)','',
            'The left panel shows descriptive contrasts for ten contiguous index blocks.',
            'The right panel shows full-run conditional limits under ideal joint assignment',
            'and zero ordinary leakage for all pulse groups and phase radii. Narrow/wide',
            'change each phase radius by one timetagger bin (78.125 ps), keeping its center',
            'fixed. Pulse grouping is distinct from phase-radius variation. Neither is a',
            'new spacelike-separation calibration. No best window is selected.','']
    output.mkdir(parents=True,exist_ok=True)
    (output/'report.md').write_text('\n'.join(lines))
    plt.rcParams.update({'svg.hashsalt':'nist-causal-audit-v1','font.size':9,'svg.fonttype':'none'})
    fig,axes=plt.subplots(1,2,figsize=(11,4),layout='constrained')
    colors=['#0072B2','#D55E00','#009E73','#CC79A7']
    for color,p in zip(colors,primary):
        rows=[r for r in result['rows'] if r['phase']=='nominal' and r['slots']==[4,5,6] and r['block']!='all' and r['direction']==p['direction'] and r['receiver_setting']==p['receiver_setting']]
        axes[0].plot([r['block']+1 for r in rows],[r['descriptive_gap']*1e6 for r in rows],'.-',color=color,label=f"{p['direction'].replace('_to_','→')}, receiver {p['receiver_setting']}")
    axes[0].axhline(0,color='gray',lw=.7);axes[0].set(xlabel='Contiguous row block',ylabel='Descriptive probability gap × 10⁶',title='Recorded marginal drift (no confidence claim)')
    axes[0].legend(fontsize=8)
    for vi,phase in enumerate(['nominal','narrow','wide']):
        ys=[]
        for size in (1,3,5):
            ys.append(max(r['ideal_assignment_interval']['absolute_upper']*100 for r in result['rows'] if r['phase']==phase and len(r['slots'])==size and r['block']=='all'))
        axes[1].plot([1,3,5],ys,'o-',label=phase)
    axes[1].set(xlabel='Number of pulses in record',ylabel='Largest conditional upper limit (%)',title='Assumed ideal assignment; assumed leakage = 0')
    axes[1].legend();axes[1].set_xticks([1,3,5])
    fig.savefig(output/'diagnostics.svg',metadata={'Date':None})
    svg=output/'diagnostics.svg'
    svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')
    fig.savefig(output/'diagnostics.png', dpi=120)
    plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output-dir',type=Path,required=True)
    a=p.parse_args();render(a.output_dir)
