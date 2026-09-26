"""Finite-menu design and exact non-identifiability certificates.

The restricted phase families and calibration bounds below are prospective
assumptions, not confidence regions extracted from Wen's deposited images.
"""
from fractions import Fraction as F
import math


def main_trial_budget(margin,alpha=0.01,beta=0.1):
    margin=float(margin)
    if not 0<margin<=1 or not 0<alpha<1 or not 0<beta<1:
        raise ValueError('positive margin and risks in (0,1) required')
    # Z is -1,+1 for the two output ports and 0 for no detection.
    # Under H+ E Z>=margin; under H- E Z<=-margin. Hoeffding, range width 2.
    n=math.ceil(2*math.log(1/min(alpha,beta))/margin**2)
    return {'eligible_main_trials':n,'type_I_and_II_upper_bound':math.exp(-n*margin**2/2),
            'calibration_trials':None,'total_eligible_trials':None}


def analyze():
    eta=F(4,5); visibility=F(9,10)
    # H+: slope/pi in [1/20,1/16]; H-: its reflection. Unknown common
    # reference offset/pi in [-1/40,1/40]. Positions in deposited delta-x units.
    # x=4: theta/pi in [7/40,11/40] -> sin(theta)>=1/2.
    # x=8: theta/pi in [3/8,21/40] -> sin(theta)>=9/10.
    # These are analytic conservative bounds; no sampled-grid minimization.
    candidates=[]
    for name,x,sine_bound in [('reference_Y_x4',4,F(1,2)),('reference_Y_x8',8,F(9,10))]:
        margin=eta*visibility*sine_bound
        candidates.append({'measurement':name,'position_delta_x':x,
                           'certified_sine_lower_bound':str(sine_bound),
                           'signed_score_margin':str(margin),
                           **main_trial_budget(margin)})
    for name,certificate in [
        ('reference_X_x8','cos(theta)=cos(-theta); choose the allowed offset zero in both families'),
        ('reference_intensity','abs(exp(i theta)*psi)^2=abs(psi)^2'),
        ('repeat_local_pointer_images','conj(exp(i theta)*psi)*(exp(i theta)*K)=conj(psi)*K')]:
        candidates.append({'measurement':name,'worst_case_separation':0,
                           'eligible_main_trials':None,'certificate':certificate})
    return {'kind':'conditional prospective design; no new intervention data',
            'hypotheses':{'positive_slope_pi_per_delta_x':['1/20','1/16'],
                          'negative_slope_pi_per_delta_x':['-1/16','-1/20'],
                          'reference_offset_pi':['-1/40','1/40']},
            'assumed_calibration':{'minimum_efficiency':str(eta),'minimum_visibility':str(visibility),
               'stationary_independent_eligible_trials':True,
               'warning':'These bounds are assumptions requiring separate certification; raw grayscale is not a trial record.'},
            'alpha':0.01,'beta':0.1,'candidates':candidates,
            'recommended_within_declared_menu':'reference_Y_x8',
            'recommendation_scope':'smallest certified main-trial budget among the specified candidates, not global apparatus optimization',
            'unrestricted_completion_result':{'worst_case_separation':0,
              'certificate':'If both slope families include zero, the same completion belongs to both and has identical predictions for every measurement. No uniform discrimination test can have alpha+beta<1.',
              'required_extra_assumption':'Disjoint restricted phase families and an independently bounded reference offset, efficiency and visibility.'},
            'calibration_cost_status':'Unavailable from deposit. No end-to-end total trial budget is claimed.'}
