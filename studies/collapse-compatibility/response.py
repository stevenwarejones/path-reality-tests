"""Gaussian mass-density responses. All conventions are explicit in derivation.md."""
import numpy as np
from scipy import constants as co, special as sp
from scipy.integrate import quad

M0 = co.m_u
D = 133e-9
L = .98375


def velocity_nodes(n=15, mean=158., sigma=9.):
    z = np.sort(np.r_[sp.ndtri(np.linspace(1/(n+1), n/(n+1), n)), [-2.5,-2,-1.5,1.5,2,2.5]])
    v = mean+sigma*z
    if np.any(v <= 0): raise ValueError("Nonpositive velocity")
    widths = np.r_[np.diff(v)[0]/2, (v[2:]-v[:-2])/2, np.diff(v)[-1]/2]
    w = widths*np.exp(-z*z/2)
    return v, w/w.sum()


def optical_coefficients(mass_u, velocity, powers_w, classical=False):
    """Independent vectorized signed B_-1 B_2 B_-1 and B_0^3, perfect loss gratings."""
    m,v = np.broadcast_arrays(np.asarray(mass_u)*M0, np.asarray(velocity))
    alpha = -4.45e-30*m/M0/22.9897
    absorption = (.5837*(m/M0/1000)**(2/3)-1.5)*4.24e-20
    n = np.array([16*D*absorption*p/((2*np.pi)**1.5*co.hbar*co.c*waist*v)
                  for p,waist in zip(powers_w,[560e-6,540e-6,575e-6])])
    phi = 8*np.sqrt(2*np.pi)*alpha*powers_w[1]/(co.hbar*co.c*540e-6*v)
    xi = (L/v)/(D*D*m/co.h)
    a = phi*(np.pi*xi if classical else np.sin(np.pi*xi))
    b = -n[1]/2 if classical else n[1]*(np.sin(np.pi*xi/2)**2-.5)
    den = a-b
    ratio = np.divide(a+b,den,out=np.zeros_like(a),where=abs(den)>1e-20)
    middle = np.exp(-n[1]/2)*sp.jv(2,np.sqrt((a*a-b*b).astype(complex)))*ratio
    middle = np.where(abs(den)<=1e-20,np.exp(-n[1]/2)*(a+b)**2/8,middle)
    s1 = sp.ive(1,n[0]/2)*middle*sp.ive(1,n[2]/2)
    if np.max(abs(s1.imag))>1e-10: raise ArithmeticError("Complex harmonic")
    return np.prod(sp.ive(0,n/2),axis=0),s1.real


def triangular_white_coefficient(mass_u, velocity, rc):
    """Point-particle log damping / lambda, seconds, two equal free flights."""
    if rc<=0: raise ValueError("Invalid rc")
    m,v = np.broadcast_arrays(np.asarray(mass_u),np.asarray(velocity))
    t = L/v
    x = co.h*t/(m*M0*D)/(2*rc)
    loss = np.where(x<1e-3,x*x/3-x**4/10+x**6/42,1-np.sqrt(np.pi)*sp.erf(x)/(2*x))
    return 2*t*m*m*loss


def author_sphere_tau_exponent(mass_u,velocity,hbar_over_sigma_q=1e-8):
    """Reproduce author's 200-node sphere quadrature: -log(R)*tau_e."""
    m,v = np.broadcast_arrays(np.asarray(mass_u)*M0,np.asarray(velocity))
    t = L/v
    rho = (3*m/(4*np.pi*971))**(1/3)/hbar_over_sigma_q
    a = co.h*t/(D*m*hbar_over_sigma_q)
    z = np.linspace(1e-10,8,200)
    az = a[...,None]*z
    integral = np.trapezoid(np.exp(-z*z/2)*sp.spherical_jn(1,rho[...,None]*z)**2*(1-sp.sici(az)[0]/az),z,axis=-1)
    return 2*t*(m/co.m_e)**2*9*np.sqrt(2/np.pi)/rho**2*integral


def sphere_static_coefficient(radius,density,separation,rc):
    """Homogeneous sphere log damping/(lambda*time), dimensionless."""
    if min(radius,density,separation,rc)<=0: raise ValueError("Invalid sphere")
    mass = 4*np.pi/3*radius**3*density
    def integrand(z):
        u = z*radius/rc
        form = 1. if abs(u)<1e-7 else 3*sp.spherical_jn(1,u)/u
        x = z*separation/rc
        loss = x*x/6-x**4/120 if abs(x)<1e-3 else 1-np.sinc(x/np.pi)
        return z*z*np.exp(-z*z)*form*form*loss
    value,_ = quad(integrand,0,10,epsabs=1e-11,limit=500)
    return 4/np.sqrt(np.pi)*(mass/M0)**2*value


def ou_time(time,tau):
    """Double temporal integral for exp(-|t-s|/tau)/(2 tau)."""
    if time<0 or tau<0: raise ValueError("Negative time")
    if tau==0: return time
    x = time/tau
    if x<1e-3: return time*(x/2-x*x/6+x**3/24-x**4/120)
    return time+tau*np.expm1(-x)


def ou_spectrum(omega,tau):
    if tau<0 or not np.isfinite(tau): raise ValueError("Invalid tau")
    return 1/(1+(np.asarray(omega)*tau)**2)


def moving_local_ou_factor(velocity,rc,tau):
    """Local straight-path overlap integral only; NOT a colored Talbot-Lau map."""
    if velocity<=0 or rc<=0 or tau<0: raise ValueError("Invalid parameters")
    if tau==0: return 1.
    x = rc/(velocity*tau)
    return np.sqrt(np.pi)*x*sp.erfcx(x)


def layered_force_coefficient(rc,thickness=370e-9,reference_mass=co.m_u):
    """Isolated 47-layer load TWO-SIDED force PSD/lambda; excludes rest of assembly."""
    if min(rc,thickness,reference_mass)<=0: raise ValueError("Invalid geometry")
    rho = np.array([7170. if k%2==0 else 2200. for k in range(47)])
    jumps = np.diff(np.r_[0.,rho,0.])
    edges = np.arange(48)*thickness
    iz = np.sqrt(np.pi)/rc*np.sum(jumps[:,None]*jumps[None,:]*np.exp(-(edges[:,None]-edges[None,:])**2/(4*rc*rc)))
    def lateral(length):
        x = length/(2*rc)
        return 2*np.pi*(length*sp.erf(x)+2*rc/np.sqrt(np.pi)*np.expm1(-x*x))
    return co.hbar**2*rc**3/(np.pi**1.5*reference_mass**2)*lateral(113e-6)*lateral(82e-6)*iz


def prescribed_path_kernel(times,weights,path_a,path_b,omega,rc,mass_ratio=1.):
    """Quadrature of prescribed-path Gaussian response; not a quantum path integral."""
    t,w,a,b = map(np.asarray,(times,weights,path_a,path_b))
    if rc<=0 or mass_ratio<=0 or a.shape!=b.shape or len(a)!=len(t) or len(w)!=len(t): raise ValueError("Invalid paths")
    if not all(np.isfinite(x).all() for x in [t,w,a,b]) or np.any(w<=0) or np.any(np.diff(t)<=0): raise ValueError("Invalid nodes")
    def overlap(x,y):
        return np.exp(-np.sum((x[:,None]-y[None,:])**2,axis=-1)/(4*rc*rc))
    gram = overlap(a,a)+overlap(b,b)-overlap(a,b)-overlap(b,a)
    return .5*mass_ratio**2*np.sum(w[:,None]*w[None,:]*np.cos(omega*(t[:,None]-t[None,:]))*gram)

