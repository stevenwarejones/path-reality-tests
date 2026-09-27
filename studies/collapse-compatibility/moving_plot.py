"""Deterministic visual summary of the moving-field response calculation."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from moving_response import lateral_frozen_factor
from combinations import cube_response


def generate(data,output):
    matplotlib.rcParams.update({'svg.fonttype':'none','svg.hashsalt':'collapse-motion-v1','font.size':10})
    fig,axes=plt.subplots(1,3,figsize=(15,4),layout='constrained')
    a=data['benchmark']['a_lambda_over_tau_s_minus2'];q=cube_response(1e-7)[1]
    speeds=np.logspace(-4,5,350)
    response=np.array([a*q*lateral_frozen_factor(2*np.pi*.003,v,.046,1e-7)/5.7e-34 for v in speeds])
    axes[0].loglog(speeds,response,label='Axis-aligned frozen response')
    axes[0].axhline(1,color='k',ls='--',label='Published torque envelope')
    vmin,vmax=data['velocity_audit']['lpf_speed_minmax_m_s']
    axes[0].axvspan(vmin-2000,vmax+2000,color='orange',alpha=.3,label='Nominal LPF ±2 km/s')
    axes[0].set(xlabel='Relative speed (m/s)',ylabel='Torque / published envelope',title='Motion data restrict a hidden variable',ylim=(1e-8,1e8))
    axes[0].legend(fontsize=8)
    rows=data['comparisons'];x=np.arange(len(rows));width=.35
    axes[1].bar(x-width/2,[r['force_fraction_upper'] for r in rows],width,label='Assembly + sensor window')
    axes[1].bar(x+width/2,[r['torque_fraction_upper'] for r in rows],width,label='LPF, all orientations')
    axes[1].axhline(1,color='k',ls='--');axes[1].set_yscale('log')
    axes[1].set(xticks=x,xticklabels=['0','2','10'],xlabel='Assumed LPF velocity error radius (km/s)',ylabel='Response upper bound / envelope',title='Same-frame witness remains below bounds',ylim=(1e-3,2))
    axes[1].legend(fontsize=8)
    rows=data['sodium']['rows'];delta=[r['max_change_over_measured_visibility_se'] for r in rows]
    axes[2].plot(np.arange(1,len(rows)+1),delta,'.',color='tab:green')
    axes[2].set(xlabel='Sodium scan (archive order)',ylabel='Maximum visibility change / measured SE',title='Signed mass/velocity mixture')
    axes[2].ticklabel_format(axis='y',style='sci',scilimits=(0,0))
    fig.savefig(output,metadata={'Date':None})
    plt.close(fig)
    if Path(output).suffix=='.svg':
        p=Path(output);p.write_text('\n'.join(x.rstrip() for x in p.read_text().splitlines())+'\n')
