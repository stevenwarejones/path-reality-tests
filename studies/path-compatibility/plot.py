#!/usr/bin/env python3
"""Regenerate the two scientific figures from checked numerical results."""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import boundary as b
from fractions import Fraction as F

ROOT=Path(__file__).resolve().parent

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output-dir',type=Path,default=ROOT/'figures');args=ap.parse_args()
    args.output_dir.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    q=np.linspace(float(b.G),1,1000)
    d=np.array([float(b.required_disturbance(F(str(v)))) for v in q])
    fig,ax=plt.subplots(figsize=(8,4.7),layout='constrained')
    ax.fill_between(q,d,1,color='#cbe8d5',label='Compatible: explicit two-state models')
    ax.plot(q,d,color='#185b39',lw=2,label='Sharp full-table boundary')
    single=np.maximum(0,(float(b.A)-q*float(b.F0))/(1-float(b.F0)))
    ax.plot(q,single,'--',color='#9e5a12',label='Negative-success inequality alone')
    ax.scatter([.64],[float(F(297,1225))],color='#185b39',zorder=4)
    ax.annotate('q = 0.64: d min = 0.24245',(.64,float(F(297,1225))),xytext=(.68,.31),arrowprops={'arrowstyle':'-'})
    ax.axvspan(.45,float(b.G),color='#ececec',label='Negative marginal already exceeds cap')
    ax.set(xlim=(.45,1),ylim=(0,.5),xlabel='Ontic negative-response cap q',ylabel='Identity-mixture disturbance d',title='Reference table: exact compatibility within the declared model')
    ax.legend(loc='upper right',fontsize=9)
    fig.savefig(args.output_dir/'compatibility-boundary.png',dpi=160);plt.close(fig)
    r=json.loads((ROOT/'results.json').read_text())['wen'];x=np.array(r['x_coordinates']);target=np.array(r['deposited_complex_slice'])
    fig,axes=plt.subplots(1,2,figsize=(11,4.2),layout='constrained',sharex=True)
    for part,ax in enumerate(axes):
        ax.plot(x,target[:,part],'ko',ms=4,label='Deposited repeated-measurement mean')
        for name,label in [('all_y','All y'),('strip_y_1','Strip |y| <= 1'),('central_y','Central y=0')]:
            item=r['roi_results'][name];a=np.array(item['reconstructed_before_scale']);z=a[:,0]+1j*a[:,1]
            c=complex(*item['comparison']['fitted_complex_scale']);v=c*z
            ax.plot(x,v.real if part==0 else v.imag,label=label,alpha=.85)
        ax.set(xlabel='Output position / delta x',ylabel='Arbitrary units',title=['Real component','Imaginary component'][part])
    axes[0].legend(fontsize=8)
    fig.suptitle('One image set versus deposited slice; each ROI fitted by one common complex scale',fontsize=11)
    fig.savefig(args.output_dir/'wen-example-reconstruction.png',dpi=160);plt.close(fig)

if __name__=='__main__':main()
