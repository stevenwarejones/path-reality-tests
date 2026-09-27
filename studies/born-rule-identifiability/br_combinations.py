"""Single-qubit probability/control equivalence and complementary preparation design."""
import math
import numpy as np
from scipy.stats import beta


def probability(p, theta):
    """Binary projective diagnostic; ensembles are averaged after this pure-state rule."""
    p=np.asarray(p,float)
    if not np.isfinite(p).all() or np.any((p<0)|(p>1)) or not -1<=theta<=1:
        raise ValueError('Require p in [0,1], theta in [-1,1]')
    return p+theta*p*(1-p)*(2*p-1)


def h(x, theta):
    return 2*probability((np.asarray(x)+1)/2,theta)-1


def ideal(phi, depth):
    # Exactly {0,+1,0,-1}, without accumulated trigonometric error for large depths.
    s=(0,1,0,-1)[int(depth)%4]
    return (1+s*np.cos(phi))/2


def warped_angle(phi, theta):
    """One common control function, no independently fitted per-angle parameters."""
    phi=np.asarray(phi)
    return np.copysign(np.arccos(np.clip(h(np.cos(phi),theta),-1,1)),np.sin(phi))


def interval(k,n,alpha):
    if not 0<=k<=n or n<=0 or not 0<alpha<1: raise ValueError('Invalid binomial data')
    return (0. if k==0 else float(beta.ppf(alpha/2,k,n-k+1)),
            1. if k==n else float(beta.ppf(1-alpha/2,k+1,n-k)))


def quadrature_decision(counts, shots, alpha=.01, radius=1.):
    """Four independent prespecified binomial proportions, simultaneous CP intervals.

    Differences for +/-X and +/-Y eliminate common readout offset. Radius is an
    externally justified ordinary contrast bound, not estimated from these counts.
    """
    if len(counts)!=4 or radius<0: raise ValueError('Invalid design')
    bounds=[interval(int(k),shots,alpha/4) for k in counts]
    diff=[(bounds[0][0]-bounds[1][1],bounds[0][1]-bounds[1][0]),
          (bounds[2][0]-bounds[3][1],bounds[2][1]-bounds[3][0])]
    lower=sum((0. if lo<=0<=hi else min(abs(lo),abs(hi)))**2 for lo,hi in diff)
    return lower>radius**2,lower


def models(theta=.01):
    phi=np.linspace(-math.pi,math.pi,257)
    warp=warped_angle(phi,theta)
    errors={str(n):float(np.max(abs(probability(ideal(phi,n),theta)-ideal(warp,n)))) for n in [1,2,3,5,9,70]}
    x=1/math.sqrt(2)
    contrast=float(h(x,theta))
    return {'theta':theta,'max_exact_family_residual':max(errors.values()),
        'depth_residuals':errors,'quadrature_squared_radius':2*contrast**2,
        'maximum_axis_warp_rad':float(np.max(abs(warp-phi))),
        'depth1_only_identifiable':False,'depth5_only_identifiable':False,
        'joint_depths_identifiable_with_common_warp':False,
        'added_complementary_preparations_radius':2*contrast**2}


def feasibility(shots, observed_contrast):
    """Raw binomial controls at an archive-derived shot scale; no measured exclusion."""
    rng=np.random.default_rng(7538941);reps=1000;theta=.01;x=1/math.sqrt(2)
    cases=[('ordinary',0.,False,1.),('probability_deformation',theta,False,1.),
           ('ordinary_common_axis_warp',theta,True,1.),
           ('deformation_99_percent_contrast',theta,False,.99),
           ('deformation_hardware_reference_contrast',theta,False,observed_contrast),
           ('negative_deformation',-theta,False,1.)]
    rows=[]
    for name,th,warp,visibility in cases:
        if warp:
            angle=float(warped_angle(math.pi/4,th)); cx,cy=math.cos(angle),math.sin(angle)
        else: cx=cy=float(h(x,th))
        cx*=visibility;cy*=visibility
        p=np.array([(1+cx)/2,(1-cx)/2,(1+cy)/2,(1-cy)/2])
        rejected=0;relaxed=0
        for counts in rng.binomial(shots,p,size=(reps,4)):
            rejected+=quadrature_decision(counts,shots)[0]
            # Remove independent orthogonality/calibration: conservative enlarged radius.
            relaxed+=quadrature_decision(counts,shots,radius=1.01)[0]
        rows.append({'scenario':name,'rejected':int(rejected),'relaxed_calibration_rejected':int(relaxed),
                     'visibility':visibility,
                     'replications':reps,'true_squared_radius':cx*cx+cy*cy,
                     'mc_95_interval':list(interval(rejected,reps,.05))})
    # The new preparation is indispensable: every original-depth binomial probability
    # is exactly the same for the alternative and the ordinary warped model.
    return {'evidence_type':'simulation_and_exact_observation_map','seed':7538941,
        'shots_per_preparation':shots,'source_of_shot_scale':'median pooled per-angle hardware count for depth 1; not an executed quadrature experiment',
        'replications':reps,'alpha':.01,'rows':rows,
        'removal':{'remove_complementary_preparation':'exact indistinguishability for all integer depths',
                   'remove_orthogonality_certificate':'no unit-radius witness; enlarged radius tested'},
        'scope':'known pure antipodal orthogonal preparations, common quantum channel/POVM, stable readout, independent binomial trials'}
