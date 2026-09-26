"""Exact full-table compatibility for the declared reference ontic interface.

Fractions are proof certificates checked independently of any optimizer.
No empirical calibration residual is interpreted as an ontic distance.
"""
from fractions import Fraction as F

A=F(1369,15625)
B=F(144,15625)
F0=F(49,625)
G=F(337,625)
TARGET=((A,F(7056,15625)),(B,F(7056,15625)))
Q_MAX=1-B/F0
CROSS=((1-F0)*Q_MAX-A)/(1-2*F0)


def required_disturbance(q):
    """Sharp minimum for the FULL table, arbitrary finite ontic cardinality.

The proof is in derivation.md: three universal inequalities and a two-state
attainer. None exists when q < the measured negative marginal G.
"""
    q=F(q)
    if not 0<=q<=1:
        raise ValueError('q must lie in [0,1]')
    if q<G:
        return None
    return max(F(1,50),(A-q*F0)/(1-F0),1-q-B/F0)


def attaining_model(q,d=None):
    q=F(q)
    minimum=required_disturbance(q)
    if minimum is None:
        raise ValueError('negative marginal exceeds response cap')
    d=minimum if d is None else F(d)
    if not minimum<=d<=1:
        raise ValueError('d must be at least the sharp minimum and at most 1')
    u=min(q,Q_MAX)
    x=(A-F0*u)/(1-F0)
    # Indexing: state S,F; pointer minus,plus; destination S,F.
    kernel=(((u,F(0)),(B/F0,1-u-B/F0)),
            ((x,F(49,100)),(F(0),1-x-F(49,100))))
    disturbance=tuple(tuple((sum(kernel[l][m][j] for m in range(2))
                            -(1-d)*(l==j))/d for j in range(2)) for l in range(2))
    model={'preparation':(F0,1-F0),'response':(F(1),F(0)),
           'kernel':kernel,'disturbance_kernel':disturbance,'q':q,'d':d}
    verify_model(model)
    return model


def verify_model(model):
    mu=model['preparation']; r=model['response']; k=model['kernel']
    D=model['disturbance_kernel']; q=model['q']; d=model['d']
    assert sum(mu)==1 and min(mu)>=0
    for l in range(2):
        assert sum(map(sum,k[l]))==1 and min(v for row in k[l] for v in row)>=0
        assert sum(k[l][0])<=q
        assert sum(D[l])==1 and min(D[l])>=0
        for j in range(2):
            assert sum(k[l][m][j] for m in range(2)) == (1-d)*(l==j)+d*D[l][j]
    table=tuple(tuple(sum(mu[l]*k[l][m][j]*(r[j] if f==0 else 1-r[j])
                             for l in range(2) for j in range(2))
                       for f in range(2)) for m in range(2))
    assert table==TARGET
    assert sum(mu[l]*r[l] for l in range(2))==F0


def positive_floor(q,d,f_interval,b_error=0,prep_tv=0,readout_error=0):
    """Conditional floor on observed b, not a confidence interval.

f_interval bounds measured bypass f before mismatch allowances. prep_tv and
readout_error bound the difference between the model's bypass and that value;
b_error bounds discrepancy of the observed positive-success probability.
All are supplied assumptions, not derived from SD columns or tomography residuals.
"""
    lo,hi=map(F,f_interval)
    errors=tuple(map(F,(b_error,prep_tv,readout_error)))
    if not 0<=lo<=hi<=1 or any(not 0<=e<=1 for e in errors):
        raise ValueError('invalid interval or error budget')
    be,pe,re=errors
    lo=max(F(0),lo-pe-re); hi=min(F(1),hi+pe+re)
    c=1-F(q)-F(d)
    return max(F(0),min(c*lo,c*hi)-be)
