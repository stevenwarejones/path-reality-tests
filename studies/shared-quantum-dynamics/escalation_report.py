"""Offline escalation verification and generated claim-limited research report."""
import argparse,json,hashlib,itertools,math
import numpy as np
from sqd_sources import HERE,split
from sqd_general import certificate,density_probability
from sqd_repeat import BINS,GRID,exact_upper,rational_repeat
from sqd_inference import deviance,simultaneous_intervals
from transfer_audit import setup,probabilities


def load(name):return json.loads((HERE/'results'/name).read_text())

def make_report():
    rows=load('observations.json');train=split(rows);fam=np.array([r['family'] for r in rows]);k=np.array([r['k'] for r in rows]);n=np.array([r['n'] for r in rows]);fits=load('general-fits.json')
    if (HERE/'results/general-polish-initial.json').exists():fits['joint_continuation_120']=load('general-polish-initial.json')
    if (HERE/'results/general-polish.json').exists():fits['joint_continuation']=load('general-polish.json')
    lines=['# General CPTP escalation: measured bounds and a broader obstruction','',
        'Can shared ordinary dynamics predict the acquisition, and what costs do the combined records identify? The answer remains split: we now have general-channel local fits, an optimizer-independent outer certificate, and a broader exact leakage obstruction. An all-length stationary qubit construction is compatible with a declared conservative joint count region. We do **not** have a successful held-out stationary-qubit predictor or a global exclusion.','',
        '## General CPTP fits','',
        'All three gates are arbitrary stationary qubit CPTP maps, with arbitrary physical preparation/readout up to a common unitary gauge. The original GST≤1024 / RB≤1000 training split is retained. The all-length fit is an explicitly in-sample diagnostic. Optimizer success is a termination flag, not a global optimality certificate.','',
        '| Fit source | Training/all-fit deviance | Evaluations | Terminated | Optimality diagnostic | Long GST violations | Long RB violations |',
        '|---|---:|---:|---|---:|---:|---:|']
    for name,r in fits.items():
        p=np.array(r['probabilities']);assert len(p)==len(rows)
        cert=certificate(r['x']);assert cert['trace_preservation_max_error']<1e-10 and cert['choi_min_eigenvalue']>=-1e-10
        for family in ['GST','RB']:
            for part,tt in [('train',train),('validation',~train)]:
                mask=tt&(fam==family);lo,hi=simultaneous_intervals(k[mask],n[mask]);s=r['scores'][family+'_'+part]
                assert abs(s['deviance']-float(deviance(k[mask],n[mask],p[mask]).sum()))<1e-6
                assert s['cell_violations']==int(np.sum((p[mask]<lo)|(p[mask]>hi)))
        lines.append(f"| {name} | {r['deviance']:.6f} | {r['nfev']} | {r['success']} | {r['optimality']:.4g} | {r['scores']['GST_validation']['cell_violations']} | {r['scores']['RB_validation']['cell_violations']} |")
    geometry=load('general-geometry.json')
    lines += ['', 'The observable Jacobian retains gauge redundancies: 52 factor coordinates describe at most 31 generic observable degrees of freedom. Relative 10⁻⁶ diagnostic ranks are GST '+str(geometry['sources']['GST']['rank_relative_1e-6'])+', RB '+str(geometry['sources']['RB']['rank_relative_1e-6'])+', joint '+str(geometry['sources']['joint']['rank_relative_1e-6'])+'. Each threshold scales with its own largest singular value, so these numbers are not a monotone dimension test or global identification theorem. In particular, no missing structural degree of freedom is demonstrated by this local rank comparison.']
    from general_compatibility import calculate as compatibility
    from certificates import close
    comp=compatibility();close(comp,load('general-compatibility.json'))
    assert comp['compatible']
    lines += ['', '**Constructive full-acquisition compatibility:** the all-length model has zero cell violations in a simultaneous region allocating α=0.025 to all 6,661 cells and α=0.025 to ten family/length means. Minimum cell slack is '+f"{comp['minimum_cell_slack']:.8f}"+'; minimum slack across all constraints is '+f"{comp['minimum_all_constraint_slack']:.8f}"+', compared with a 2×10⁻⁸ numerical verification tolerance. Cell endpoints use exact outward binomial-tail arithmetic. This explicit stationary qubit model has no leakage level and no nontrivial classical memory. **No strictly positive leakage or memory cost follows from this region.** This is an in-sample, numerically verified physical compatibility example, not held-out validation, a globally optimal likelihood fit, or adequacy under every stronger statistic. The region/mean audit is an explicitly retrospective diagnostic; it does not replace the frozen predictive test.']
    lines += ['', 'Violations use the earlier simultaneous Clopper–Pearson intervals separately within each evaluation family. For all_joint these are in-sample diagnostics, **not held-out predictions**. Budget-limited or poorly stationary local solutions cannot settle whether a better general model exists. Even disappearance of cell violations would not prove a full sampling-model goodness of fit.','',
        '## Full-class outer prediction certificate','',
        'For a whole word channel Φ_w, common reset and binary effect, every stationary qubit CPTP model obeys p(ww) ≤ min(1,(2√p(w)+√p(empty))²). See the independent proof in [escalation-theory.md](../escalation-theory.md). Six exact outward binomial bounds, with Bonferroni tail 1/120, give simultaneous 95% coverage under within-row binomial sampling. There is no intermediate reset in ww. These doubled sequences are predictions, not newly acquired observations.','',
        '| RB row (zero-based) | Word length | Plus / shots | Doubled length | Joint outer upper bound | Either source removed: this relaxation |',
        '|---|---:|---:|---:|---:|---:|']
    bounds=load('repeat-bounds.json');a0=exact_upper(bounds['empty']['k'],bounds['empty']['n'],120)
    assert a0==bounds['empty_upper_numerator']
    lookup={(r['family'],r['row']):r for r in rows};hashes={r['word_sha256'] for r in rows}
    for r,(low,high) in zip(bounds['targets'],BINS):
        original=lookup[('RB',r['row'])]
        for key in ['k','n','word_sha256','length']:assert r[key]==original[key]
        assert hashlib.sha256(bytes(r['word'])).hexdigest()==r['word_sha256']
        eligible=[v for v in rows if v['family']=='RB' and low<=v['length']<=high]
        assert min(eligible,key=lambda v:(v['word_sha256'],v['row']))['row']==r['row']
        a1=exact_upper(r['k'],r['n'],120);upper=rational_repeat(a0,a1)
        assert a1==r['rb_upper_numerator'] and upper==r['joint_upper_numerator']
        assert r['joint_upper']==upper/GRID
        assert r['target_already_recorded']==(hashlib.sha256(bytes(r['word'])*2).hexdigest() in hashes)
        lines.append(f"| {r['row']} | {r['length']} | {r['k']}/{r['n']} | {r['doubled_length']} | {r['joint_upper']:.10f} | 1 |")
    lines += ['',f"GST empty circuit: {bounds['empty']['k']}/{bounds['empty']['n']}; outward upper bound {a0/GRID:.10f}. Exact integer binomial tails and integer square-root rounding certify the printed bounds at denominator 10¹⁰.",'',
        '**What removal establishes:** both sources are required by this particular certificate. **What it does not establish:** that either full single-source feasible set reaches one, or that these bounds improve on every implicit constraint of GST alone. The earlier 30-model attainable spread is retained only as historical exploratory evidence; it is not promoted to a full-class region. With unrestricted within-row dependence the reported statistical upper bound is one.','',
        '## Why an observable-basis relaxation was insufficient','']
    obs=load('observable-audit.json');A=np.array(obs['segment_first_numerators'],dtype=object);B=np.array(obs['segment_second_numerators'],dtype=object);scale=obs['determinant_segment_denominator']
    det=lambda M:sum((-1)**sum(p[i]>p[j] for i in range(4) for j in range(i+1,4))*math.prod(int(M[i,p[i]]) for i in range(4)) for p in itertools.permutations(range(4)))
    assert det(A)==int(obs['first_exact_determinant']) and det(B)==int(obs['second_exact_determinant']) and det(A)*det(B)<0
    for row in obs['rows']:
        original=lookup[(row['family'],row['row'])]
        for key in ['k','n','word_sha256','length']:assert row[key]==original[key]
    ol,oh=simultaneous_intervals(np.array([r['k'] for r in obs['rows']]),np.array([r['n'] for r in obs['rows']]),alpha=.025*16/3785)
    np.testing.assert_allclose(ol.reshape(4,4),obs['lower'],atol=1e-12,rtol=0)
    np.testing.assert_allclose(oh.reshape(4,4),obs['upper'],atol=1e-12,rtol=0)
    for M in [A,B]:
        assert np.all(np.array(M,float)/scale>=np.array(obs['lower'])) and np.all(np.array(M,float)/scale<=np.array(obs['upper']))
        seen={}
        for value,row in zip(M.ravel(),obs['rows']):
            key=row['word_sha256']
            if key in seen:assert seen[key]==value
            seen[key]=value
    lines += ['A nominally complete four-fiducial GST matrix has nonzero point-estimate singular values, but the entrywise simultaneous count box contains a singular matrix. Two saved rational matrices inside that box have opposite **exact integer determinant signs**, so their segment crosses determinant zero. Repeated circuit entries are shared. Uniform inversion of this box is therefore impossible. This certifies failure of this relaxation, not impossibility of stronger CPTP-aware bounds: the singular witness need not have a common physical realization or satisfy other records.','',
        '## Gate-dependent, nonunital leakage obstruction','',
        'The exact similarity now permits gate-dependent loss/return, arbitrary computational CPTP channels, anisotropy and nonunital translation. Unknown physical return-state polarization absorbs the translation. Positivity and readout restrict the similarity parameter κ to an explicit interval. Every word remains indistinguishable, while leakage populations and loss rates scale by κ. A two-state classical flag also realizes this entire block model. This is a theorem about a model class, not a fitted device population.','']
    t=load('transfer-audit.json');old,new,locked,state,E,E2,kk,kp,tp,lower=setup()
    assert abs(max(lower)-t['common_kappa_lower'])<1e-12
    selected=t['locked_return_break_examples'];p=probabilities(selected,old,state,E);q=probabilities(selected,new,state,E2);r=probabilities(selected,locked,state,E2)
    for i,row in enumerate(selected):
        assert abs(p[i]-row['original_probability'])<1e-12 and abs(q[i]-row['free_return_probability'])<1e-12 and abs(r[i]-row['locked_return_probability'])<1e-12
    lines += [f"Synthetic control on all 6,661 recorded words: κ={t['kappa']}, common conditional interval lower endpoint {t['common_kappa_lower']:.6f}, maximum all-word numerical difference {t['all_recorded_words_max_error']:.3g}. Independent qutrit Kraus TP error {t['trace_preservation_max_error']:.3g}; minimum Choi eigenvalue {t['choi_min_eigenvalue']:.3g}; minimum transformed return-state eigenvalue {t['transformed_return_state_min_eigenvalue']:.6f}.",'',
        '| Actual recorded geometry | Maximum change when return state is incorrectly locked to I/2 |','|---|---:|']
    for row in selected:lines.append(f"| {row['family']} row {row['row']}, length {row['length']} | {row['absolute_break']:.9f} |")
    power=t['locked_return_detectability']
    lines += ['', 'The locked-return perturbation is weak at the deposited exposures even with known nuisance parameters. Under independent binomial sampling, D(P_original || P_locked) is '+f"{power['joint']['independent_binomial_KL_original_vs_locked']:.8f}"+' jointly. Pinsker gives power ≤ α+√(D/2): a level-0.05 test has power at most '+f"{power['joint']['level_005_power_upper_by_pinsker']:.6f}"+' for this synthetic pair (GST alone '+f"{power['GST']['level_005_power_upper_by_pinsker']:.6f}"+', RB alone '+f"{power['RB']['level_005_power_upper_by_pinsker']:.6f}"+'). These are numerical evaluations of an analytic power ceiling, not an experimental fit. Nuisance refitting cannot improve a uniformly valid test past the simple-pair ceiling. Along the exact free-return orbit, power equals size for every test because the complete count laws coincide.']
    from lift_certificate import calculate as lift_calculate
    lift=load('lift-certificate.json');offline=load('lift-certificate.json')
    for m in offline['models']:m.pop('full_source_max_probability_error',None)
    close(lift_calculate(),offline)
    witness=lift['recorded_example']
    lines += ['', '## Exact lifts of the compatible qubit model','',
        'The stronger construction starts from the compatible all-length general qubit map Ψ_g. For τ=I/2 it chooses Φ_g=(Ψ_g−λ_g R_τ)/(1−λ_g), return state τ_g=τ+(Ψ_g(τ)−τ)/η_g, and leakage effect z=Tr(Eτ). Complete positivity holds exactly when λ_g≤2λ_min(J_g); return positivity requires η_g≥||d_g||. The collapse Q(ρ_C,L)=ρ_C+Lτ intertwines the lifted channel with Ψ_g. Therefore the qubit, leakage and two-state-flag realizations give **identical probabilities for every word**, not only this deposit. These are sharp conditions for this construction with the qubit representation fixed, not gauge-invariant device bounds.','',
        'Unlike the original restricted-fit examples, these lifts share the prediction vector inside the full-acquisition joint count region. They require the particular unknown leakage response z=Tr(Eτ) and unknown return polarization; calibrated apparatus information could exclude them.','',
        '| Ordinary realization | Conditional leakage after recorded Gx^8192 | Maximum difference over all 6,661 probabilities |','|---|---:|---:|']
    for i,m in enumerate(lift['models']):
        assert m['intertwining_max_error']<1e-12 and m['trace_preservation_max_error']<1e-12 and m['choi_min_eigenvalue']>=-1e-12
        lines.append(f"| Lift fraction {m['fraction_of_selected_loss']} | {100*witness['conditional_leakage_populations'][i]:.6f}% | {m['full_source_max_probability_error']:.3g} |")
    lines += ['',f"The recorded GST row is {witness['row']} (zero-based), with {witness['k']}/{witness['n']} plus counts; all models predict {witness['common_probability']:.9f}. The population differences are stipulated model resources, **not observed populations**. The all-word theorem, not more sampling of the existing binary readout, supplies the nonidentification conclusion. A calibrated leakage-sensitive monitor after one primitive gate distinguishes zero from its positive loss rate."]
    lines += ['', 'The GST example is `GxGxGx(Gi)^8192`. These are simulated differences on actual words, not observed residual attribution or powered detection claims. Freeing return polarization restores exact equivalence. A calibrated subspace monitor after one Gi would distinguish 0.001 from 0.0006 in this example. More of the same unknown binary readout cannot.','',
        '## Claim status after escalation','',
        '| Claim | Status |','|---|---|',
        '| Arbitrary stationary qubit CPTP fits implemented | Physical local numerical examples; optimizer limits retained |',
        '| Five full-class joint outer prediction bounds | Certified analytically and by exact count-tail arithmetic; conditional sampling |',
        '| Strict contraction of complete GST-only feasible prediction region by RB | Not established; certificate removal is weaker than full-model removal |',
        '| Nonunital/gate-dependent all-word leakage equivalence | Analytic similarity plus exact lifts of the compatible qubit fit; conditional CP budgets and Kraus controls |',
        '| Full-acquisition joint count-region compatibility | Explicit stationary qubit point; in-sample, numerical physical certificate |',
        '| Necessary device leakage or finite-memory cost | No positive lower bound follows from this declared region; actual populations not measured |',
        '| All stationary qubit models excluded by acquisition | Not certified |',
        '| Major novelty beyond original GST/RB and leakage-identification work | Not established; precise comparison in prior-art.md |','',
        'The original paper already compared GST with RB, used general Markovian gate-set estimates, and studied alternating errors. This escalation closes a model-generality gap in our implementation, supplies a valid outer bound in place of the finite-bank interpretation, and exposes return-state calibration as another anchor needed to identify leakage. It does not claim that reanalysis of this known pairing alone earns the intended major-result bar.','']
    return '\n'.join(lines)


if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--check',action='store_true');args=a.parse_args();text=make_report();path=HERE/'results/escalation.md'
    if args.check:
        if path.read_text()!=text:raise ValueError('escalation report differs')
    else:path.write_text(text)
    print('Escalation physical, count, bound, source-selection and report checks passed')
