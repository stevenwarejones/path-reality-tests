#!/usr/bin/env python3
"""Deterministic scientific figures from committed results, never from new fits."""
import argparse
import io
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from common import HERE,LEVELS,PHASE_CELLS
plt.rcParams.update({'svg.hashsalt':'fine-timing-fixed-v1','svg.fonttype':'none',
                     'font.family':'DejaVu Sans','font.size':10})


def save(fig,name,check):
    buf=io.StringIO(); fig.savefig(buf,format='svg',metadata={'Date':None}); plt.close(fig)
    path=HERE/'figures'/name; text='\n'.join(line.rstrip() for line in buf.getvalue().splitlines())+'\n'
    if check:
        if path.read_text()!=text: raise SystemExit(name+' differs')
    else: path.parent.mkdir(exist_ok=True); path.write_text(text)


def run(check=False):
    r=json.loads((HERE/'results.json').read_text()); recovery=json.loads((HERE/'recovery-results.json').read_text())
    chosen=[v for v in r['cases'] if v['epsilon']==0 and v['envelope']=='event_supported']
    fig,axes=plt.subplots(2,2,figsize=(12,7),layout='constrained')
    for ax,v in zip(axes.flat,chosen):
        band=np.array(v['cdf']['band'])*1e6; x=np.arange(len(band))
        ax.fill_between(x,band[:,0],band[:,1],step='post',color='#2376aa',alpha=.23,label='Simultaneous contrast band')
        ax.plot(x,band[:,0],color='#2376aa',linewidth=.6); ax.plot(x,band[:,1],color='#2376aa',linewidth=.6)
        ax.axhline(0,color='black',linewidth=.7)
        for boundary in (162,324): ax.axvline(boundary,color='gray',ls=':',linewidth=.7)
        ax.set_xticks([0,90,161,162+90,323,324+90,485],['4:1','4:91','4:162','5:91','5:162','6:91','6:162'])
        ax.tick_params(axis='x',labelsize=8)
        ax.set_title(f'{v["receiver"].title()}, local setting {v["local_setting"]}')
        ax.set_xlabel('Pulse bit : phase upper boundary (strict cut)'); ax.set_ylabel('CDF contrast (full-trial ppm)')
    fig.suptitle('486 fixed cumulative cuts; ideal assignment, event-supported completion\nBlue regions are inferential bands, not descriptive histograms',fontsize=12)
    save(fig,'cdf-bands.svg',check)
    fig,axes=plt.subplots(1,2,figsize=(12,4.8),layout='constrained')
    for v in chosen:
        label=f'{v["receiver"]}/{v["local_setting"]}'
        axes[0].plot(LEVELS,np.array(v['partitions']['tv'])[:,1]*1e6,'o-',label=label)
        axes[1].plot(LEVELS[1:],np.array(v['partitions']['gain'])[1:,1]*1e6,'o-',label=label)
    axes[0].set_title('TV upper bound by retained resolution'); axes[1].set_title('Upper bound on added TV beyond coarse')
    for ax in axes:
        ax.set_ylim(bottom=0); ax.set_ylabel('Full-trial ppm'); ax.set_xlabel('Frozen partition'); ax.grid(alpha=.2); ax.legend()
    fig.suptitle('All measured lower bounds are zero; finer resolution increases uncertainty',fontsize=12)
    save(fig,'resolution.svg',check)
    fig,axes=plt.subplots(2,2,figsize=(12,7),layout='constrained')
    for ax,actual in zip(axes.flat,chosen):
        base=[v for v in recovery['cases'] if v['receiver']==actual['receiver'] and v['local_setting']==actual['local_setting'] and v['mode']=='within' and v['scale']==1 and v['envelope']=='event_supported']
        for label,kind,key,color in [('Coarse','tv','coarse','#777777'),('CDF','cdf',None,'#087e8b'),('Fine TV','tv','width1','#6f42a8'),('Resolved J','gain','width1','#c45718')]:
            values=[v['outcomes'][kind] if key is None else v['outcomes'][kind][key] for v in base]
            rates=np.array([a['rate'] for a in values]); intervals=np.array([a['interval95'] for a in values])
            ax.errorbar([v['truth']['tv'][-1]*1e6 for v in base],rates,yerr=np.array([rates-intervals[:,0],intervals[:,1]-rates]),label=label,color=color,marker='o',capsize=2)
        ax.set_title(f'{actual["receiver"].title()}, local setting {actual["local_setting"]}')
        ax.set_ylim(-.03,1.03); ax.set_xlabel('Injected within-core TV (full-trial ppm)'); ax.set_ylabel('Recovery fraction'); ax.grid(alpha=.2)
    axes[0,0].legend(fontsize=9)
    fig.suptitle('Paired simulations at archive size with event-supported uncertainty\n400 repetitions per point; pointwise exact 95% Monte Carlo intervals',fontsize=12)
    save(fig,'recovery.svg',check)
    # Sample-size comparison complements the archive-size paired curves.
    fig,axes=plt.subplots(2,2,figsize=(12,7),layout='constrained')
    for ax,actual in zip(axes.flat,chosen):
        for scale in (.25,1.,4.):
            vs=[v for v in recovery['cases'] if v['receiver']==actual['receiver'] and v['local_setting']==actual['local_setting'] and v['mode']=='within' and v['scale']==scale and v['envelope']=='event_supported']
            ax.plot([v['truth']['tv'][-1]*1e6 for v in vs],[v['outcomes']['tv']['width1']['rate'] for v in vs],'o-',label=f'{scale:g} × rows and uncertainty')
        ax.set_title(f'{actual["receiver"].title()}, local setting {actual["local_setting"]}'); ax.set_xlabel('Injected TV (full-trial ppm)'); ax.set_ylabel('width1 recovery fraction'); ax.set_ylim(-.03,1.03); ax.grid(alpha=.2); ax.legend(fontsize=8)
    fig.suptitle('Size sensitivity in the declared independent-block model',fontsize=12)
    save(fig,'sample-size.svg',check)

if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('--check',action='store_true'); a=p.parse_args(); run(a.check)
    print('fine figures '+('verified' if a.check else 'written'))
