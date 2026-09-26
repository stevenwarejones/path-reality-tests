"""Generated portfolio findings and candidate-removal comparisons."""
def render(audit,result):
    lines=['# Expanded dataset-combination findings','',
        'Main result: a normalized single-qubit probability deformation is exactly indistinguishable from one ordinary control-axis warp throughout the audited integer gate-depth circuit family. Combining depths cannot identify it. A complementary-preparation witness can separate the pair ideally, but contrast/calibration uncertainty removes that prospective gain. No measured Born-rule exclusion is obtained. The waveguide runner-up supplies a separate exact optical ambiguity and phase-cycling design.','',
        '## Public source audit','',
        '| Archive | Evidence | Backend | Jobs | Circuit rows | Largest first-harmonic residual |','|---|---|---|---:|---:|---:|']
    for name,d in audit['ibm'].items():
        residual=max([v['residual_max_abs'] for v in d['depths'].values()],default=None)
        val='not mapped' if residual is None else f'{residual:.7f}'
        lines.append(f"| {name} | {d['evidence_type']} | {d['backend']} | {d['jobs']} | {d['circuit_rows']} | {val} |")
    lines+=['','Hardware residuals are descriptive, with job-level errors in [portfolio-audit.json](portfolio-audit.json). They are not a theta fit or discovery significance. The simulator is not independent hardware evidence. Missing benchmark gate labels are not inferred from order.','',
        '| Photon archive | Outcome matrix | Available event settings | Max event/P discrepancy | Max POVM row-sum error |','|---|---|---:|---:|---:|']
    for name,d in audit['photon'].items():
        lines.append(f"| {name} | {d['outcome_matrix_shape']} | {len(d['event_reconciliation'])} | {max(r['max_frequency_discrepancy'] for r in d['event_reconciliation']):.7f} | {d['response_row_sum_max_error']:.7f} |")
    lines+=['','Zero outcomes are retained. The supplied detector POVMs are derived under assumed coherent-state statistics; they are not independent Born-free calibration. PNR classified events and fitted probability products differ and cannot be pooled as independent replications.','',
        '## Exact single-versus-joint result','',
        'For fθ(p)=p+θp(1−p)(2p−1), hθ(x)=2fθ((1+x)/2)−1, choose cos(gθ(φ))=hθ(cosφ). Since sin(nπ/2) is 0 or ±1, fθ(p_n(φ))=p_n(gθ(φ)) for every integer n. Independent complex circuit multiplication verifies the convention. Both datasets individually and jointly remain ambiguous with the same one-parameter warp. This is an exact conditional observation-map theorem, not an assertion that all real residuals follow this model.','',
        '| θ | Largest numerical equality error over checked depths | Ideal complementary squared radius |','|---:|---:|---:|']
    for r in result['qubit_equivalence']:
        lines.append(f"| {r['theta']} | {r['max_exact_family_residual']:.2g} | {r['quadrature_squared_radius']:.7f} |")
    p=result['qubit_prospective']
    lines+=['','## Raw-count feasibility and removal','',
        f"Each of four prospective preparations uses {p['shots_per_preparation']:,} shots, matching the hardware pooled per-angle depth-1 scale. These new preparations were **not** measured in the public archive. Simultaneous exact binomial intervals use family alpha 0.01.",'',
        '| Scenario | Rejections / 1000 | Rejections with radius 1.01 | MC 95% interval |','|---|---:|---:|---|']
    for r in p['rows']:
        lo,hi=r['mc_95_interval']
        lines.append(f"| {r['scenario']} | {r['rejected']} | {r['relaxed_calibration_rejected']} | [{lo:.4f}, {hi:.4f}] |")
    lines+=['','Removing the complementary preparation restores exact ambiguity at any shot count. Relaxing the preparation/measurement consistency removes the unit-disk certificate; radius 1.01 is a sensitivity scenario, not a measured bound. Loss of contrast and negative theta can remain inside the ordinary disk, so the witness is not a general identifiability result.','',
        'The photon-calibration negative control separately verifies (P T)R_d=P(T R_d) with independent responses for two detectors; removing either detector does not remove this ordinary stochastic gauge. Treating T as a Born-rule modification would be a category error.','',
        '## What is and is not new','',
        'Known ingredients include quantum control gauges, phase cycling, PSD coherence and concentration bounds. This PR provides an explicit all-depth counterexample, public byte-level audit, complementary-preparation cost/sensitivity, and an optical runner-up. Physics novelty is not established. The expanded [candidate inventory and bridges](../combinations.md) gives four proposals and precise reasons the available combinations do not yet support a probability-rule exclusion.','',
        'Highest-value next input: independently certified complementary ±X/±Y preparations, with common channel/readout and contrast bounds. The optical alternative is a calibrated antipodal all-open phase exposure. No author contact, new cloud experiment, or raw redistribution was performed.','']
    return '\n'.join(lines)
