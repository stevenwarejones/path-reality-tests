"""Figure of certificate redundancy; explicitly not an exclusion plot."""
import argparse
import io
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

HERE=Path(__file__).resolve().parent


def figure(png=None):
    r=json.loads((HERE/'results/boundary.json').read_text());s=r['source_removal']
    plt.rcParams.update({'svg.hashsalt':'moving-noise-boundary-v1','font.size':10})
    fig,axes=plt.subplots(1,2,figsize=(10,4.3),layout='constrained')
    labels=['Sodium only','Sensor only','Both']
    vals=[s['sodium_only_certificate_variance'],s['sensor_only_certificate_variance'],s['both_certificate_variance']]
    axes[0].scatter(vals,[2,1,0],s=85,color=['#17657d','#ac6c2d','#17657d'])
    axes[0].set_yticks([2,1,0],labels);axes[0].set_xscale('log')
    axes[0].axvline(r['benchmark']['variance_s_minus2'],color='#657066',ls='--',label='Constructed variance')
    axes[0].set_xlabel('Sufficient certificate ceiling (s⁻²)')
    axes[0].set_ylim(-.6,2.6);axes[0].set_title('Adding the sensor gives no gain')
    axes[0].legend(loc='upper right',fontsize=8);axes[0].grid(axis='x',alpha=.2)
    tail=np.array([0,1e-6,1e-5,1e-4,3.4055187e-4,1e-3])
    from noise_kernels import quantum_tv_bound,harmonic_visibility_bound
    src=json.loads((HERE/'results/sources.json').read_text())
    v=r['benchmark']['variance_s_minus2'];ratios=[]
    for eta in tail:
        eps=quantum_tv_bound(v,195500,.025,eta)
        ratios.append(max(harmonic_visibility_bound(eps,x['optical_mean_probability'],x['optical_visibility'],x['ols_row_l1'])/x['visibility_se'] for x in src['sodium']['scans']))
    axes[1].plot(tail,ratios,'o-',color='#17657d');axes[1].axhline(1,color='#ac6c2d',ls='--',label='1 measured SE')
    axes[1].set_xscale('symlog',linthresh=1e-6);axes[1].set_xlim(0,1e-3);axes[1].set_yscale('log')
    axes[1].set_xlabel('Unbounded mass/exposure population weight')
    axes[1].set_ylabel('Maximum visibility bound / measured SE')
    axes[1].set_title('Unmeasured tails limit the certificate')
    axes[1].legend(fontsize=8);axes[1].grid(alpha=.2)
    fig.suptitle('Conditional response bounds — not empirical exclusions',fontsize=13)
    stream=io.StringIO();fig.savefig(stream,format='svg',metadata={'Date':None})
    if png: fig.savefig(png,dpi=130)
    plt.close(fig)
    return '\n'.join(line.rstrip() for line in stream.getvalue().splitlines())+'\n'


def main():
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');p.add_argument('--png',type=Path)
    args=p.parse_args();content=figure(args.png);path=HERE/'results/boundary.svg'
    if args.check:
        if content!=path.read_text():raise SystemExit('Figure differs')
    else:path.write_text(content)
    print('Boundary figure reproduces' if args.check else 'Boundary figure generated')

if __name__=='__main__':main()
