"""Offline recomputation of scores, physical witnesses, source-removal and narrative."""
import argparse,json
import numpy as np
from sqd_sources import HERE,split
from sqd_models import certificate,parameters,density_probability,dormant_flag_probability
from sqd_inference import deviance,simultaneous_intervals,context_certificate
from sqd_equivalence import fiber_interval,fiber_point,invariants,synthetic_certificate,leakage_population
from removal import acceptance
from certificates import close


def load(name):return json.loads((HERE/'results'/name).read_text())


def verify_offline():
    rows=load('observations.json');fits=load('fits.json');bank=load('removal.json');eq=load('equivalence.json')
    train=split(rows);family=np.array([r['family'] for r in rows]);k=np.array([r['k'] for r in rows]);n=np.array([r['n'] for r in rows])
    if len(rows)!=6661 or train.sum()!=3785:raise ValueError('unexpected source geometry')
    for label,f in fits.items():
        p=np.array(f['probabilities'])
        if p.shape!=(len(rows),) or np.any(~np.isfinite(p)) or np.any(p<0) or np.any(p>1):raise ValueError('invalid saved predictions')
        c=certificate(f['x'],f['model'])
        if c['trace_preservation_max_error']>1e-10 or c['choi_min_eigenvalue']< -1e-10:raise ValueError('nonphysical gate certificate')
        for fam in ('GST','RB'):
            for part,tt in [('train',train),('validation',~train)]:
                mask=tt&(family==fam);saved=f['scores'][fam+'_'+part]
                actual={'rows':int(mask.sum()),'shots':int(n[mask].sum()),'deviance':float(deviance(k[mask],n[mask],p[mask]).sum()),'rms_probability_error':float(np.sqrt(np.mean((k[mask]/n[mask]-p[mask])**2)))}
                close(actual,saved)
    for label,b in bank['bank'].items():
        certificate(b['x'],b['model']);p=np.array(b['probabilities'])
        for source in ('GST','RB','joint'):
            mask=train if source=='joint' else train&(family==source)
            close(acceptance([r for r,m in zip(rows,mask) if m],p[mask]),b['acceptance'][source])
    for source in ('GST','RB','joint'):
        accepted=[(name,b) for name,b in bank['bank'].items() if b['acceptance'][source]['accepted']]
        spread=np.ptp(np.array([b['probabilities'] for _,b in accepted]),axis=0);v=bank['ranges'][source]
        if v['accepted']!=len(accepted) or set(v['labels'])!={a for a,b in accepted}:raise ValueError('bank membership mismatch')
        for fam in ('GST','RB'):
            close(float(spread[(~train)&(family==fam)].max()),v[fam+'_validation_max_attainable_spread'])
    close(synthetic_certificate(),eq['synthetic'])
    for source,v in eq['measured_fit_fibers'].items():
        x=fits['leakage_'+source]['x'];close(list(fiber_interval(x)),v['loss_interval']);close(invariants(x),v['invariants'])
        a,b=v['parameters_a'],v['parameters_b']
        for w in [(),(1,2,0),(1,1,2,0)*10]:
            pa=density_probability(w,a,'leakage');pb=density_probability(w,b,'leakage');pc=dormant_flag_probability(w,a)
            if max(abs(pa-pb),abs(pa-pc))>2e-9:raise ValueError('physical equivalence failure')
        close(leakage_population(1000,a),v['population_at_1000_a']);close(leakage_population(1000,b),v['population_at_1000_b'])
    close(json.loads(json.dumps(context_certificate(rows))),load('audit.json')['context_certificate'])
    controls=load('controls.json');cal=load('calibration.json')
    if len(controls['null_scores'])!=39 or len(cal['cases'])!=5:raise ValueError('incomplete controls')
    for name,c in controls['controls'].items():
        if len(c['qubit_scores'])!=20 or len(c['generating_class_scores'])!=20:raise ValueError('incomplete control')
        if int(np.sum(np.array(c['qubit_scores'])>controls['empirical_95_threshold']))!=c['narrow_qubit_rejections']:raise ValueError('control count mismatch')
        cc=cal['cases'][name]
        if len(cc['null_scores'])!=39 or cc['threshold']!=sorted(cc['null_scores'])[37]:raise ValueError('calibration mismatch')
        if int(np.sum(np.array(c['qubit_scores'])>cc['threshold']))!=cc['qubit_exceedances']:raise ValueError('adapted count mismatch')
    if len(load('dependence.json')['qubit_scores'])!=20:raise ValueError('incomplete dependence control')
    return rows,fits,bank,eq


def make_report():
    rows,fits,bank,eq=verify_offline();audit=load('audit.json');gate=load('gate-a.json');control=load('controls.json');cal=load('calibration.json');dep=load('dependence.json')
    lines=['# Shared quantum dynamics: measured combination and resource obstruction','',
      '**Question:** Can one ordinary set of operations predict both acquisition families, and what physical resources do their records identify?','',
      '**Result:** The interleaved GST/RB combination removes explicit ordinary-model ambiguities that either family leaves open. It does not identify leakage versus a bounded dormant classical flag. Within the declared incoherent leakage family, an exact all-word transformation trades loss, return, readout response and depolarization without changing any terminal probability. The sharp conditional compatibility interval is derived analytically.','',
      '**Prediction limitation:** None of the saved joint restricted-model fits passes every long-GST simultaneous cell interval. RB predictions are much better. This is not a complete temporal explanation, a rejection of all stationary qubit channels, or a necessary nonzero leakage/memory cost.','',
      'The endpoint is a restricted operational obstruction with explicit measured witnesses. Publication-level novelty is unestablished; the original paper already studied alternating errors. See [proof and model scope](../theory.md), [source audit](../sources.md), and [prior art](../prior-art.md).','',
      '## Claim status','', '| Claim | Status |','|---|---|',
      '| Counts, words, denominators, preserved row order | Measured; hash-pinned and independently parsed |',
      '| Same acquisition: GST/RB interleaving | Primary-paper assertion; timestamps absent from condensed files |',
      '| Shared-gate physical fits and source-removal examples | Conditional numerical compatibility examples; restricted models, local optimizers |',
      '| Exact leakage fiber and dormant-flag equality | Analytically certified within stated channel class; independently checked Kraus realizations |',
      '| Fitted leakage populations are actual device populations | Unavailable; no operational subspace monitor |',
      '| Globally necessary leakage or finite memory | Not established; compatible qubit training examples remain |',
      '| Correlated drift and injected-error power | Simulated ordinary controls at actual exposures; not experimental discovery |',
      '| Broader theory violation, novel-priority claim | Not claimed |','',
      '## Matched measured geometry','', '| Family | Rows | Shots | Train rows | Validation rows |','|---|---:|---:|---:|---:|']
    for fam in ('GST','RB'):
        f=fits['qubit_joint']['scores'];a=audit['families'][fam]
        lines.append(f"| {fam} | {a['rows']} | {a['total_shots']} | {f[fam+'_train']['rows']} | {f[fam+'_validation']['rows']} |")
    lines+=['','Training is GST primitive length ≤1,024 and RB length ≤1,000. Longer sequences are held out for frozen prediction. No random split of duplicate circuits is used. Both files use Gi/Gx/Gy, one preparation, and binary terminal readout. Gates, preparation and readout are shared across families within the declared stationarity assumption. This does not connect Nairobi epochs.','',
      'The independent Gate A example is Gy⁴: GST records 0/50 plus counts, RB 3/285 in its first matching row. Local weighted Jacobian diagnostics at the same generic 18-parameter leakage point are:','', '| Design | Thresholded rank | Parameter count |','|---|---:|---:|']
    for f,v in gate['diagnostics'].items():lines.append(f"| {f} | {v['rank_relative_1e-6']} | {v['parameter_count']} |")
    lines+=['','The relative cutoff is 1e-6 of the largest singular value. Units and finite differences matter. These are diagnostics, not global rank certificates. The exact all-word fiber is the structural result.','',
      '## Frozen source-removal predictions','', '| Model | Fit source | GST validation deviance | RB validation deviance | GST / RB simultaneous-cell violations |','|---|---|---:|---:|---:|']
    train=split(rows);fam=np.array([r['family'] for r in rows]);k=np.array([r['k'] for r in rows]);n=np.array([r['n'] for r in rows])
    for name,f in fits.items():
        p=np.array(f['probabilities']);viol=[]
        for ff in ('GST','RB'):
            mask=(~train)&(fam==ff);lo,hi=simultaneous_intervals(k[mask],n[mask]);viol.append(int(np.sum((p[mask]<lo)|(p[mask]>hi))))
        lines.append(f"| {f['model']} | {f['source']} | {f['scores']['GST_validation']['deviance']:.2f} | {f['scores']['RB_validation']['deviance']:.2f} | {viol[0]} / {viol[1]} |")
    lines+=['','Deviances are descriptive, not Wilks-calibrated global exclusions. Each displayed validation-family cell check is separately 95% simultaneous under iid Bernoulli shots per circuit (no independence between circuits is needed for its union bound). It excludes that fixed prediction vector when violated, not the full model class. Training/validation conditional independence would be needed to treat the fitted predictor as fixed for a prospective validation test. These data are retrospective and lack chronology.','',
      '## Identification gain beyond a shifted optimizer','',
      'An exploratory bank of 30 fully physical shared-operation models contains the twelve primary fits and fixed alternating-amplitude profiles. Nuisances are refitted using training counts only. A single 95% joint count region uses 3,785 exact binomial cell intervals and eight Hoeffding length-bin averages, splitting alpha equally. Source removal drops constraints from that same region; it never changes the physical model.','',
      '| Available source constraints | Bank models retained | Maximum attainable long-GST probability spread | Maximum attainable long-RB probability spread |','|---|---:|---:|---:|']
    for source,v in bank['ranges'].items():lines.append(f"| {source} | {v['accepted']} | {v['GST_validation_max_attainable_spread']:.6f} | {v['RB_validation_max_attainable_spread']:.6f} |")
    gst_only=[name for name,b in bank['bank'].items() if b['acceptance']['GST']['accepted'] and not b['acceptance']['joint']['accepted']]
    rb_only=[name for name,b in bank['bank'].items() if b['acceptance']['RB']['accepted'] and not b['acceptance']['joint']['accepted']]
    lines+=['',f'Adding RB removes {len(gst_only)} bank examples permitted by GST; adding GST removes {len(rb_only)} permitted by RB. These are explicit physical witnesses of complementary constraints, not just two separate best fits.', '',
      '**Scope:** These spreads describe an attainable subset of the full prediction region. They are lower bounds on its diameter, not outer confidence bounds, and their reduction does not prove that all physical predictions lie in the displayed joint spread. The sampling bounds can exclude the individual saved prediction vectors even when the bank was chosen adaptively, because the count region is simultaneous. They do not make the selected bank exhaustive.','',
      '## Sharp conditional leakage fibers','',
      'For each fitted leakage observable vector, the following interval is analytically necessary and sufficient along its declared fixed-coordinate fiber. The parameter box is imposed for numerical sensitivity, not measured hardware calibration. Distinct points keep all measured-word probabilities equal and also agree on every unmeasured word.','',
      '| Fit source | Loss-rate interval per primitive | Population after 1,000 gates, example A / B | Max probability difference over all 6,661 rows |','|---|---|---|---:|']
    for source,v in eq['measured_fit_fibers'].items():
        lo,hi=v['loss_interval'];lines.append(f"| {source} | [{lo:.8g}, {hi:.8g}] | {v['population_at_1000_a']:.7f} / {v['population_at_1000_b']:.7f} | < 1e-9 |")
    lines+=['','These are **not confidence intervals on actual leakage**. They condition on one local fitted observable vector; other fits and model classes remain possible. The joint example is a physical training-compatible model whose long-GST predictions fail. Its equivalence persists whether or not that particular vector fits the hardware. A separate two-state classical dormant flag reproduces the same channel recursion and all terminal probabilities. No arbitrary circuit-label memory is used.','',
      'A calibrated leakage projector after one primitive measures the loss rate directly and breaks this fiber. A known-contrast monitor also works with propagated calibration uncertainty. More terminal binary words, another device’s leakage data, or an uncalibrated extra detector output do not supply that information.','',
      '## Ordinary adversaries and calibration failure','',
      'All controls generate counts at the 3,785 actual training-row exposures. Every qubit fit refits all qubit nuisance coordinates; generating-class fits also refit nuisances. Each mechanism has 20 simulations. A 39-replicate null calibration at the measured joint qubit fit is deliberately compared with a separate-pilot, nuisance-adapted 39-replicate calibration. One start is used for these exploratory simulations; local failure is never an exclusion certificate.','',
      '| Injected ordinary mechanism | Fixed-reference exceedances / 20 | Pilot-adapted exceedances / 20 |','|---|---:|---:|']
    for name,c in control['controls'].items():lines.append(f"| {name} | {c['narrow_qubit_rejections']} | {cal['cases'][name]['qubit_exceedances']} |")
    lines+=['','The fixed-reference threshold is not portable across SPAM nuisance values: it produces false alarms even for imperfect reset and readout that are inside the fitted qubit class. The pilot-adapted column is a conditional empirical sensitivity estimate, not a uniformly calibrated test. With only 20 trials, 0/20 or 20/20 has a two-sided exact 95% interval approximately [0,.168] or [.832,1]; intermediate rates have similarly substantial uncertainty. Generating-class fit scores and convergence flags are retained, so residuals are not classified as new physics.','',
      f"A further 20 simulations share a coherent sign across hypothetical blocks of 64 file rows. Their qubit deviance range is {min(dep['qubit_scores']):.1f}–{max(dep['qubit_scores']):.1f}. These blocks are an adversary, not reconstructed experimental jobs. Marginal drift-fit scores are retained; their independent-shot likelihood is not a correct joint likelihood for these correlated blocks.", '',
      'The weakest structural direction is exact: along the leakage fiber the count distributions are identical under the same sampling law, so the power of any terminal-count test is at most its size at every exposure. A synthetic pair with different loss rates and a calibrated monitor explicitly verifies recovery of the loss rate from the added observation. No known-nuisance-only recovery is offered as evidence of terminal identification.','',
      '## Dependence and the stopping boundary','',
      f"There are {len(audit['context_certificate']['words'])} words measured in both files. A simultaneous shared-word context-cost certificate gives lower bound {audit['context_certificate']['max_half_context_gap_lower']:.1f} under its iid-within-pool assumption. It does not force context dependence, leakage, or memory.", '',
      'No timestamps, independent job blocks, or reset logs are deposited in these selected records. Huge total exposure does not establish independence. With arbitrary dependence within the acquisition, no nontrivial confidence statement is claimed; the certified statistical lower bound is zero. The deterministic channel equalities and the conditional physical interval survive because they do not use sampling independence.','',
      '**What is proved impossible:** distinguishing the declared leakage-fiber points or their two-state dormant-flag encodings using any sequence of the specified gates, reset, and binary terminal readout alone. **What is merely not found:** an adequate frozen long-sequence fit in the four restricted numerical families, a globally necessary resource cost, and a matched leakage monitor for this acquisition. General CPTP fits, coherent leakage, richer bounded memory, and real acquisition drift are not ruled out.','',
      'A second exact ambiguity equates the alternating bit with an external two-phase slot clock under an assumed one-slot-per-primitive schedule. The records do not establish that schedule. A calibrated identity wait advancing the clock but not the gate-count flag separates the constructed pair; an unknown idle would not.','',
      '## Reproduction limits','',
      'Offline checks independently verify Kraus/Choi physics, synthetic and fitted equivalence witnesses, every saved count score, source-removal membership, and this report. They use committed reduced counts and saved predictions. The separate full-source route re-downloads hash-pinned bytes, reconstructs words/denominators, recomputes all primary and bank predictions, verifies independent density-matrix propagation, and regenerates source geometry. Numerical refits and simulation calibrations are explicit additional commands and are not promised bitwise unique.','']
    return '\n'.join(lines)


if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--check',action='store_true');args=a.parse_args();text=make_report();p=HERE/'results/report.md'
    if args.check:
        if p.read_text()!=text:raise SystemExit('report differs')
        print('Offline physical certificates, scores, removal sets and report reproduce')
    else:p.write_text(text)
