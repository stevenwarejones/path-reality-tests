"""Moving Gaussian random-potential responses; SI, ensemble dephasing only.

No outcome-selection dynamics is asserted. Bounds are for stated rigid source
functions and normalized finite-observation kernels, not arbitrary instruments.
"""
from functools import lru_cache
import numpy as np
from scipy import constants as co, special as sp
from scipy.integrate import quad
from numpy.polynomial.legendre import leggauss
from response import M0, ou_time


def top_hat_covariance(length, displacement, rc):
    """Integral exp(-rc² k²)|2 sin(k L/2)/k|² exp(i k x) dk."""
    x = np.asarray(displacement)
    def primitive(y):
        return y*sp.erf(y/(2*rc))+2*rc/np.sqrt(np.pi)*np.exp(-y*y/(4*rc*rc))
    return np.maximum(0.,np.pi*(primitive(x+length)+primitive(x-length)-2*primitive(x)))


def lateral_ou_factor(omega, tau, speed, length, rc):
    """Exact factor for translation along a separable undifferentiated side.

    Applicable to the layer's x/y translation and a cube translating along its
    torque axis. The response is single-body S1=2 lambda Q0 * this factor.
    Integration uses separation, not a grid that misses a narrowing Lorentzian.
    """
    if min(tau,speed)<0 or min(length,rc)<=0: raise ValueError('Invalid parameters')
    if tau==0: return 1.
    if speed==0: return 1/(1+(omega*tau)**2)
    if speed*tau/rc<1e-3:
        return quad(lambda z: np.exp(-z)*np.cos(omega*tau*z)*
                    top_hat_covariance(length,speed*tau*z,rc)/top_hat_covariance(length,0,rc),0,32,epsabs=1e-11)[0]
    integral=quad(lambda x: np.exp(-x/(speed*tau))*np.cos(omega*x/speed)*
                  top_hat_covariance(length,x,rc)/top_hat_covariance(length,0,rc),
                  0,length+12*rc,points=[length],epsabs=length*1e-10,limit=200)[0]
    return integral/(speed*tau)


def lateral_frozen_factor(omega,speed,length,rc):
    """Limit tau*lateral_ou_factor = pi W_parallel(omega/v)/(v I0)."""
    if speed<=0: raise ValueError('Treat zero speed separately')
    k=omega/speed
    return np.pi*length**2*np.sinc(k*length/(2*np.pi))**2*np.exp(-(rc*k)**2)/(speed*top_hat_covariance(length,0,rc))


@lru_cache(maxsize=32)
def sphere_nodes(radius,rc,n=80):
    """Radial overlap-volume representation of a Gaussian-smoothed sphere."""
    z,w=leggauss(n); r=radius*(z+1); w=w*radius
    volume=4*np.pi*radius**3/3
    overlap=np.pi*(4*radius+r)*(2*radius-r)**2/12
    return r,w*4*np.pi*r*r*overlap/volume**2


def sphere_overlap(distance,radius,density,rc):
    """Mass-normalized spatial covariance g(s), without temporal factor."""
    s=np.atleast_1d(np.asarray(distance,float))
    r,w=sphere_nodes(radius,rc)
    q=s[:,None]*r[None,:]/rc**2
    angular=np.ones_like(q)
    np.divide(-np.expm1(-q),q,out=angular,where=q!=0)
    # exp(-(r-s)^2/4rc²) (1-exp(-rs/rc²))/(rs/rc²)
    value=(np.exp(-(s[:,None]-r)**2/(4*rc**2))*angular)@w
    mass=4*np.pi*radius**3*density/3
    result=(mass/M0)**2*value
    return float(result[0]) if np.ndim(distance)==0 else result


def moving_sphere_exponent_per_a(time,tau,speed,radius,density,separation,rc,cos_angle=0.):
    """Exponent divided by a=lambda/tau, exact prescribed rigid paths.

    tau=np.inf is the static temporal spectral atom. Separation is fixed;
    both paths translate at speed, angle between separation and velocity fixed.
    """
    if time<=0 or tau<=0 or speed<0 or abs(cos_angle)>1: raise ValueError('Invalid benchmark')
    g=sphere_overlap(0,radius,density,rc)-sphere_overlap(separation,radius,density,rc)
    if speed==0:
        return g*(time*time/2 if np.isinf(tau) else tau*ou_time(time,tau))
    cutoff=min(speed*time/rc,(2*radius+separation)/rc+12)
    def integrand(z):
        x=z*rc; t=x/speed
        plus=np.sqrt(x*x+separation**2+2*x*separation*cos_angle)
        minus=np.sqrt(x*x+separation**2-2*x*separation*cos_angle)
        diff=sphere_overlap(x,radius,density,rc)-.5*(sphere_overlap(plus,radius,density,rc)+sphere_overlap(minus,radius,density,rc))
        return (time-t)*np.exp(-t/tau)*diff*rc/speed
    scale=g*time*min(time,rc/speed)
    return quad(integrand,0,cutoff,epsabs=max(scale*1e-10,1e-300),limit=200)[0]


def covariance_length_bound(q0,total_variation,diameter,rc):
    """Integral |Q(vt)|dt <= Q0*Leff/v, uniformly over directions.

    For Q= hbar²/m0² int ds ds' g(s)g(s') exp(-|x+s-s'|²/4rc²),
    with signed generalized source g of finite total variation supported in a
    set of diameter D. Cauchy-Schwarz gives Q0, support gives the Gaussian tail.
    """
    if min(q0,total_variation,diameter,rc)<=0: raise ValueError('Invalid source bound')
    b=(co.hbar/M0*total_variation)**2
    if b<q0*(1-1e-10): raise ValueError('Inconsistent source norm')
    z=np.sqrt(max(0.,np.log(b/q0)))
    return diameter+2*rc*z+rc*np.sqrt(np.pi)*sp.erfcx(z)


def moving_psd_upper_per_a(q0,total_variation,diameter,rc,speed):
    """Uniform one-sided PSD upper / a, any tau>0; not a lower bound."""
    if speed<=0: raise ValueError('Zero speed needs finite-observation response')
    return 2*q0*covariance_length_bound(q0,total_variation,diameter,rc)/speed


def periodogram_weights(n,dt,frequency,window='blackmanharris',detrend=False):
    """One-sided periodogram |b*x|², with sum window² normalization.

    Gaussian covariance propagation is exact even with detrending. Overlap and
    inter-bin dependence are NOT treated as independent gamma observations.
    """
    from scipy.signal.windows import blackmanharris
    w=blackmanharris(n,sym=False) if window=='blackmanharris' else np.ones(n)
    t=np.arange(n)*dt
    b=np.sqrt(2*dt/(w@w))*w*np.exp(-2j*np.pi*frequency*t)
    if detrend:
        design=np.column_stack([np.ones(n),np.linspace(-1,1,n)])
        b=b-(b@design)@np.linalg.solve(design.T@design,design.T)
    return b


def frozen_periodogram_per_a(q0,weights):
    return float(q0/2*abs(np.sum(weights))**2)


def relative_pair_factor(k_dot_baseline):
    """PSD of (F1-F2)/2 relative to one body's PSD, shared frozen field."""
    return np.sin(np.asarray(k_dot_baseline)/2)**2


def transverse_ou_loss(separation,rc,speed,tau):
    """Integrated local point-path loss, for velocity perpendicular to arms."""
    if speed<=0 or tau<=0: raise ValueError('Positive speed/memory required')
    factor=np.sqrt(np.pi)*rc/(speed*tau)*sp.erfcx(rc/(speed*tau))
    return factor*(-np.expm1(-np.asarray(separation)**2/(4*rc*rc)))


def rigid_box_force_q0(length,width,height,density,rc):
    i2=2*np.sqrt(np.pi)/rc*(-np.expm1(-height*height/(4*rc*rc)))
    return co.hbar**2*rc**3/(np.pi**1.5*M0**2)*density**2*top_hat_covariance(length,0,rc)*top_hat_covariance(width,0,rc)*i2


def sphere_force_q0(radius,density,rc):
    mass=4*np.pi*radius**3*density/3
    a=radius/rc
    if a>20:
        # Exact radial integral, expressed using exp(-a²); stable for large a.
        shape=3/a**4*(1-2/a**2+(1+2/a**2)*np.exp(-a*a))
        return (co.hbar*mass/M0/rc)**2*shape
    val=quad(lambda z:z**4*np.exp(-z*z)*(3*sp.spherical_jn(1,a*z)/(a*z))**2,0,10,epsabs=1e-13)[0]
    return 4/(3*np.sqrt(np.pi))*(co.hbar*mass/M0/rc)**2*val


def assembly_psd_upper_per_a(rc,speed):
    """All cross terms retained through a triangle inequality on spectral norms.

    Flexure 0<=mode(x)<=1 along z: its Q0 is bounded by the rigid Si beam.
    No exact undocumented component placement is fabricated. Glue uses a mass
    bound (1 pL, density <=4000 kg/m³), not a guessed shape.
    """
    from response import layered_force_coefficient
    load_q=layered_force_coefficient(rc)
    load_tv=113e-6*82e-6*(2*7170+46*(7170-2200))
    load=moving_psd_upper_per_a(load_q,load_tv,np.linalg.norm([113e-6,82e-6,47*370e-9]),rc,speed)
    l,w,h,rho=450e-6,57e-6,2.5e-6,2330.
    beam_q=rigid_box_force_q0(l,w,h,rho,rc)
    beam=moving_psd_upper_per_a(beam_q,2*rho*l*w,np.linalg.norm([l,w,h]),rc,speed)
    r,rho_m=15.5e-6,7430.
    magnet_q=sphere_force_q0(r,rho_m,rc)
    magnet=moving_psd_upper_per_a(magnet_q,2*np.pi*rho_m*r*r,2*r,rc,speed)
    glue_mass=4e-12
    glue=(co.hbar*glue_mass/M0)**2*np.sqrt(np.pi)/(rc*speed)
    parts={'load':load,'beam':beam,'magnet':magnet,'glue_mass_bound':glue}
    return {'parts':parts,'all_cross_terms_upper':sum(np.sqrt(list(parts.values())))**2,
            'q0_upper':(np.sqrt(load_q)+np.sqrt(beam_q)+np.sqrt(magnet_q)+co.hbar*glue_mass/(np.sqrt(2)*M0*rc))**2}


def blackman_dc_gain(n,dt,frequency,symmetric=False):
    """Exact squared transfer of a constant to a sampled Blackman PSD bin."""
    divisor=n-1 if symmetric else n
    a=[.42,-.5,.08]
    def d(cycles):
        # Stable complex finite geometric sum; reduce mod 1 first.
        x=(cycles+.5)%1-.5
        return n*np.sinc(n*x)/np.sinc(x)*np.exp(-1j*np.pi*(n-1)*x)
    cycles=frequency*dt
    total=a[0]*d(cycles)
    for j in [1,2]: total+=a[j]/2*(d(cycles-j/divisor)+d(cycles+j/divisor))
    # Exact window norm from geometric sums, including symmetric convention.
    coeff={0:a[0],1:a[1]/2,-1:a[1]/2,2:a[2]/2,-2:a[2]/2}
    norm=sum(c*cc*d((j+k)/divisor).real for j,c in coeff.items() for k,cc in coeff.items())
    return float(2*dt*abs(total)**2/norm)


def triangle_path_exponent_per_a(mass_u,beam_speed,wind_speed,rc,tau,n=80):
    """Independent two-time Gaussian phase integral for prescribed paraxial arms.

    Wind is perpendicular to the splitting axis. This checks the local limit;
    it is not identified with the full quantum grating experiment.
    """
    from response import L,D
    flight=L/beam_speed;duration=2*flight
    dmax=co.h*flight/(mass_u*M0*D)
    nodes,weights=leggauss(n)
    times=np.r_[flight*(nodes+1)/2,flight+flight*(nodes+1)/2]
    tw=np.tile(weights*flight/2,2)
    def arm(t):return dmax*np.minimum(t/flight,2-t/flight)
    if wind_speed==0:
        aa=arm(times)
        gram=np.exp(-(aa[:,None]-aa[None,:])**2/(16*rc**2))-np.exp(-(aa[:,None]+aa[None,:])**2/(16*rc**2))
        return .5*mass_u**2*np.sum(tw[:,None]*tw[None,:]*np.exp(-abs(times[:,None]-times[None,:])/tau)*gram)
    scale=min(rc/wind_speed,tau)
    total=0.
    for t,w in zip(times,tw):
        end=min(t/scale,32.)
        lag=scale*end*(nodes+1)/2
        d1,d2=arm(t),arm(t-lag)
        kernel=np.exp(-lag/tau-(wind_speed*lag/(2*rc))**2)*(np.exp(-(d1-d2)**2/(16*rc**2))-np.exp(-(d1+d2)**2/(16*rc**2)))
        total+=w*np.dot(weights,kernel)*scale*end/2
    return mass_u**2*total
