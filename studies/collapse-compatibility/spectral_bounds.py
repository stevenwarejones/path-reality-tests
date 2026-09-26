"""Finite LP diagnostics and exact certificates for explicitly piecewise-affine kernels."""
from fractions import Fraction as F
import numpy as np
from scipy.optimize import linprog


def finite_bound(response,bounds,objective):
    a,b,c = map(lambda v: np.asarray(v,float),(response,bounds,objective))
    if a.ndim!=2 or b.shape!=(a.shape[0],) or c.shape!=(a.shape[1],): raise ValueError("Invalid dimensions")
    if not all(np.isfinite(v).all() for v in [a,b,c]) or any(np.any(v<0) for v in [a,b,c]): raise ValueError("Invalid responses")
    gap = np.flatnonzero((a==0).all(axis=0)&(c>0))
    if len(gap): return {"status":"unbounded","escape_column":int(gap[0]),"scope":"finite supplied responses"}
    fit = linprog(-c,A_ub=a,b_ub=b,bounds=(0,None),method="highs")
    if not fit.success: raise ArithmeticError(fit.message)
    y = -fit.ineqlin.marginals
    pv,dv,g = float(max(0,np.max(a@fit.x-b))),float(max(0,np.max(c-y@a))),float(y@b-c@fit.x)
    if pv>1e-8 or dv>1e-8 or abs(g)>1e-8 or min(y)<-1e-10: raise ArithmeticError("LP residual failure")
    return {"status":"bounded","value":float(c@fit.x),"primal":fit.x.tolist(),"dual":y.tolist(),
            "primal_violation":pv,"dual_violation":dv,"duality_gap":g,"scope":"finite numerical result, not continuum"}


def affine_certificate(knots,kernels,macro,tail_slopes,macro_tail_slope,weights,bounds):
    """Exact [0,infinity) domination including tail. This validates declared functions, not physical kernels."""
    def rational(x):
        if not isinstance(x,(str,int)): raise ValueError("Use exact rational strings/integers")
        return F(x)
    if tail_slopes is None or macro_tail_slope is None: raise ValueError("Explicit tail required")
    k = list(map(rational,knots)); a = [list(map(rational,row)) for row in kernels]
    c = list(map(rational,macro)); slopes = list(map(rational,tail_slopes)); cs = rational(macro_tail_slope)
    y,b = list(map(rational,weights)),list(map(rational,bounds))
    if not k or k[0]!=0 or any(x>=z for x,z in zip(k,k[1:])): raise ValueError("Invalid knots")
    if len(c)!=len(k) or not a or any(len(row)!=len(k) for row in a) or not len(a)==len(y)==len(b)==len(slopes): raise ValueError("Dimensions disagree")
    if any(v<0 for row in a for v in row) or any(v<0 for v in c+y+b+slopes+[cs]): raise ValueError("Nonnegative functions required")
    residuals = [sum(yi*row[j] for yi,row in zip(y,a))-c[j] for j in range(len(k))]
    tail = sum(yi*s for yi,s in zip(y,slopes))-cs
    if min(residuals)<0 or tail<0: raise ValueError("Continuum domination failed")
    return {"upper_bound":str(sum(yi*bi for yi,bi in zip(y,b))),"knot_slacks":list(map(str,residuals)),
            "tail_slope_slack":str(tail),"evidence":"formal_conditional",
            "scope":"exact rational affine domination on [0,infinity); not Lean-verified"}

