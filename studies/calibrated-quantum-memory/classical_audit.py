#!/usr/bin/env python3
"""Verify unrestricted classical fits before checking optional quantum counterparts."""
import argparse,json
from pathlib import Path
import numpy as np
import identifiability as physical
import classical_models as models
import joint_statistics as statistics


def assess(model,counts,labels):
    models.check_model(model)
    nums,dens=models.probabilities(model,labels)
    aggregate=statistics.pearson_region(counts,nums,dens)
    cells=statistics.cell_region(counts,nums,dens)
    decisions=[aggregate['aggregate_region_membership'],cells['cell_region_membership']]
    joint='inside' if decisions==['inside','inside'] else 'excluded' if 'excluded' in decisions else 'rounding_unresolved'
    p=np.asarray(nums/dens[...,None],float);expected=counts.sum(-1)[...,None]*p
    deviance=float(2*np.where(counts>0,counts*np.log(np.maximum(counts,1)/expected),0).sum())
    return dict(**aggregate,**cells,joint_region_membership=joint,exact_classical_physicality=True,
                multinomial_deviance_display=round(deviance,6))


def analyze(source_dir,candidate_file):
    catalog=json.loads(candidate_file.read_text());reports=[]
    for model in catalog['models']:
        counts,labels,delays=physical.load_counts(source_dir,model['mapping'])
        if delays!=model['delays']:raise ValueError('delay order mismatch')
        # This result cannot depend on any quantum-counterpart certificate.
        classical=assess(model,counts,labels)
        report=dict(mapping=model['mapping'],family=model['family'],
                    internal_classical_record_sizes=[len(p['before']['choi']) for p in model['processes']],
                    quantum_constraint_used_in_classical_discovery=False,unrestricted_classical_candidate=classical)
        if 'counterpart_certificate' in model:
            if classical['joint_region_membership']!='inside':raise ValueError('assess a viable classical fit first')
            mixed=models.mix_instruments(model,model['counterpart_post_fit_mixture'])
            paired=assess(mixed,counts,labels)
            if paired['joint_region_membership']!='inside':raise ValueError('paired model fails complete region')
            certificate=model['counterpart_certificate'];value=models.npt_expectation(mixed,certificate)
            _,endpoint=models.quantum_instrument(mixed,certificate['memory_weight'])
            paired.update(exact_quantum_instrument_physicality=True,exact_full_probability_identity=True,
                          quantum_memory_weight=certificate['memory_weight'],
                          quantum_memory_weight_decimal=float(physical.F(*certificate['memory_weight'])),
                          fixed_instrument_endpoint=[endpoint.numerator,endpoint.denominator],
                          npt_delay=delays[certificate['delay_index']],
                          exact_npt_expectation=[value.numerator,value.denominator],npt_expectation_decimal=float(value),
                          isolated_probe=models.isolated_probe(mixed,certificate))
            report['paired_post_fit_candidate']=paired
        reports.append(report)
    return dict(kind='complete_fixed_region_classical_and_paired_model_certificate',familywise_alpha=.05,
                empirical_quantum_memory_certified=False,complementary_measured_source_gain_established=False,
                historical_mapping_verified=False,independent_fixed_exposure_acquisition_verified=False,
                sampling_assumptions='Conditional independent fixed-exposure multinomial rows, stable within rows; both flags retained.',
                statistical_region_changed_during_search=False,reports=reports)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source-dir',type=Path,required=True)
    p.add_argument('--candidate-file',type=Path,default=physical.HERE/'results/ibm-classical-candidates.json')
    p.add_argument('--output',type=Path,default=physical.HERE/'results/ibm-classical-audit.json')
    p.add_argument('--check',action='store_true');args=p.parse_args()
    text=json.dumps(analyze(args.source_dir,args.candidate_file),indent=2,sort_keys=True)+'\n'
    if args.check:
        if args.output.read_text()!=text:raise SystemExit('classical audit differs')
        print('Classical physicality, both statistical components and separate quantum counterparts verified')
    else:args.output.write_text(text);print(args.output)


if __name__=='__main__':main()
