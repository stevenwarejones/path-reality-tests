#!/usr/bin/env python3
"""Generate scientific figures from the committed numerical results."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from common import HERE, LAGS


def plot(directory):
    result = json.loads((HERE/'results.json').read_text())
    directory.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
    selected = [v for v in result['cancellation'] if v['period'] == 'full' and v['epsilon'] == 0 and v['envelope'] == 'event_supported']
    fig, ax = plt.subplots(figsize=(9, 5.5), layout='constrained')
    labels = [f'{v["receiver"].title()} · local {v["local_setting"]} · {v["condition"]}' for v in selected]
    for i, v in enumerate(selected):
        lo, hi = np.array(v['hidden_tv'])*1e6
        ax.plot([lo, hi], [i, i], color='#236878', linewidth=4, solid_capstyle='butt')
        ax.scatter([lo, hi], [i, i], s=20, color='#236878')
        ax.text(hi+2, i, f'{hi:.1f}', va='center', fontsize=9)
    ax.set(yticks=range(len(labels)), yticklabels=labels, xlabel='Hidden TV beyond pooled distribution (ppm of all interior rows)',
           xlim=(-3, 140), title='Cancellation remains unresolved in the two fixed local-history conditions')
    ax.invert_yaxis()
    ax.grid(axis='x', alpha=.2)
    fig.suptitle('Simultaneous conditional bounds · ideal assignment · event-supported completion', fontsize=10)
    fig.savefig(directory/'cancellation.png', dpi=170)
    plt.close(fig)
    fig, axes = plt.subplots(2, 2, figsize=(11, 7), sharex=True, sharey=True, layout='constrained')
    for row, side in enumerate(('alice', 'bob')):
        for col, y in enumerate((0, 1)):
            ax = axes[row, col]
            records = [v for v in result['lag_scores'] if v['receiver'] == side and v['local_setting'] == y
                       and v['period'] == 'full' and v['feature'] == 'any' and v['envelope'] == 'event_supported']
            records.sort(key=lambda v: v['lag'])
            central = np.array([v['trusted_score'] for v in records])*1e6
            bounds = np.array([v['completion_range'] for v in records])*1e6
            ax.errorbar(np.arange(len(records)), central, yerr=[central-bounds[:, 0], bounds[:, 1]-central],
                        fmt='o', color='#236878', capsize=3, linewidth=1.5)
            ax.axhline(0, color='gray', linewidth=1)
            ax.axvspan(4.5, 8.5, color='#d8eae7', alpha=.35)
            ax.set(title=f'{side.title()} · local setting {y}', xticks=np.arange(len(LAGS)),
                   xticklabels=[f'{lag:+d}' for lag in LAGS], xlim=(-.5, 8.5))
            ax.grid(axis='y', alpha=.2)
            if row == 1:
                ax.set_xlabel('Lag: receiver Oᵢ against remote Xᵢ₊ₗ')
            if col == 0:
                ax.set_ylabel('Signed count score (ppm)')
    fig.suptitle('Temporal fingerprints: descriptive scores\nBars show missing-record completions, NOT statistical confidence; shaded lags use later settings', fontsize=12)
    fig.savefig(directory/'lags.png', dpi=170)
    plt.close(fig)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output-dir', type=Path, required=True)
    plot(p.parse_args().output_dir)
