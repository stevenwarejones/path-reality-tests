"""Exact witnesses: information lost by quantization, aggregation and causal ambiguity."""
from fractions import Fraction as F
from math import floor


def total_variation(a, b):
    return sum((abs(a.get(k, F(0))-b.get(k, F(0))) for k in a.keys() | b.keys()), F(0))/2


def pushforward(law, function):
    out = {}
    for point, probability in law.items():
        value = function(point)
        out[value] = out.get(value, F(0))+probability
    return out


def quantization(p=F(3, 10000), offset=F(49, 100)):
    """Identical full recorded laws for any common distribution of integer tags."""
    if not 0 <= p <= 1 or not 0 < offset < F(1, 2):
        raise ValueError('invalid quantization witness')
    # A finite nontrivial tag law; the proof applies atom by atom to any tag law.
    tags = {17: F(1, 3), 104: F(2, 3)}
    laws = [{None: 1-p, **{F(k)+s*offset: p*w for k, w in tags.items()}} for s in (-1, 1)]
    quantizer = lambda t: None if t is None else floor(t+F(1, 2))
    recorded = [pushforward(law, quantizer) for law in laws]
    assert recorded[0] == recorded[1]
    latent_tv = total_variation(*laws)
    assert latent_tv == p
    return dict(recorded_tv=str(total_variation(*recorded)), latent_tv=str(latent_tv),
                conditional_event_tv='1', mean_event_shift_bins=str(2*offset),
                illustrative_shift_ps=str(2*offset*F(625, 8)),
                universal_claim='same quantized law for every shared integer-tag distribution')


def coarse_graining(p=F(3, 10000)):
    a, b = {None: 1-p, F(-2): p}, {None: 1-p, F(-1): p}
    category = lambda t: 'none' if t is None else ('early' if t < 0 else 'core' if t < 2 else 'late')
    assert pushforward(a, category) == pushforward(b, category)
    return dict(fine_record_tv=str(total_variation(a, b)), coarse_category_tv='0',
                distinction='different digitized tags inside the same category')


def cancellation():
    p, d = F(3, 10000), F(1, 10000)
    probabilities = [[p+(2*x-1)*(2*g-1)*d for x in (0, 1)] for g in (0, 1)]
    pooled = [sum((probabilities[g][x] for g in (0, 1)), F(0))/2 for x in (0, 1)]
    weighted_tv = sum((abs(v[1]-v[0])/2 for v in probabilities), F(0))
    assert pooled[0] == pooled[1] and weighted_tv == 2*d
    return dict(pooled_tv='0', state_averaged_tv=str(weighted_tv), hidden_tv=str(weighted_tv),
                state_contrasts=[str(v[1]-v[0]) for v in probabilities])


def temporal_aliasing():
    # Stationary binary symmetric Markov X, response only to X_(i-1).
    rho, d = F(4, 5), F(1, 10000)
    contrasts = {str(l): str(2*d*rho**abs(l+1)) for l in (-64, -16, -4, -1, 0, 1, 4, 16, 64)}
    assert F(contrasts['-1']) == 2*d and F(contrasts['1']) > 0
    # Direct causal and latent common-cause SCMs have the same P(X,D).
    p = F(3, 10000)
    direct = {(x, o): F(1, 2)*(p+(2*x-1)*d if o else 1-p-(2*x-1)*d)
              for x in (0, 1) for o in (0, 1)}
    common = {(x, o): sum((F(1, 2)*(p+(2*z-1)*d if o else 1-p-(2*z-1)*d)
                           for z in (0, 1) if x == z), F(0)) for x in (0, 1) for o in (0, 1)}
    assert direct == common
    return dict(markov_stay_probability='9/10', ordinary_response_delay_rows=1,
                lag_contrasts=contrasts, observational_tv_between_causal_models='0',
                direct_model_interventional_gap=str(2*d), common_cause_interventional_gap='0',
                missing_assumption='exogeneity relative to latent state')


def selection_collider():
    table = {(x, o): F(1, 4) for x in (0, 1) for o in (0, 1)}
    selected = {k: v for k, v in table.items() if k[0] == k[1]}
    assert len(selected) == 2 and sum(selected.values()) == F(1, 2)
    return dict(pooled_gap='0', selected_gap='1',
                rule='select O=X, which uses the current setting and outcome',
                prevention='only predictable pretrial states in cancellation analysis')


def arbitrary_record_lifts():
    """Distinct observed arm laws still leave their sub-bin mean gap unidentified."""
    recorded = [{-2: F(1, 2), 3: F(1, 2)}, {-1: F(1, 3), 2: F(2, 3)}]
    means = []
    for amplitude in (F(0), F(49, 100)):
        lifted = [{F(k)+(2*x-1)*amplitude: p for k, p in law.items()}
                  for x, law in enumerate(recorded)]
        for x in (0, 1):
            assert pushforward(lifted[x], lambda t: floor(t+F(1, 2))) == recorded[x]
        means.append(sum((k*p for k, p in lifted[1].items()), F(0))-
                     sum((k*p for k, p in lifted[0].items()), F(0)))
    assert means[1]-means[0] == F(49, 50)
    return dict(recorded_law_change_in_either_arm='0',
                latent_event_mean_gaps=[str(v) for v in means],
                unidentified_gap_difference_bins=str(means[1]-means[0]))


def all_witnesses():
    return dict(quantization=quantization(), arbitrary_record_lifts=arbitrary_record_lifts(), coarse_graining=coarse_graining(),
                cancellation=cancellation(), temporal_aliasing=temporal_aliasing(),
                selection_collider=selection_collider())
