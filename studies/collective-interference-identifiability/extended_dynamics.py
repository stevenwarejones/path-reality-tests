#!/usr/bin/env python3
"""Date-separated confinement/alignment/NNN diagnostic; no class exclusion.

Monte Carlo assesses each deposited point under its own multinomial law.
It does not turn local optimizer failure into infeasibility of a family.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,json
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import minimize
from scipy.special import xlogy
from scipy.stats import beta
from numpy.polynomial.hermite import hermgauss
from dynamics import close
HERE=Path(__file__).resolve().parent
SCALE=np.array([120,110,8,7,20,20,3,3,.4,.08,.96])
BOUNDS=[(60,200),(60,200),(-20,20),(-20,20),(-100,100),(-100,100),(0,25),(0,25),(-.49,.2),(0,.3),(.85,.999)]


def july7_records():
    cs=json.loads((HERE/'results/matched-regions.json').read_text())['cases'][:15]
    return [dict(time_ms=c['time_ms'],input_y=j-1,y_start=-5,x_start=-5,width=12,shots=r['shots'],counts=r['counts'],intervals=r['intervals']) for c in cs for j,r in enumerate(c['singletons'])]


def july2_records(cache):
    import xarray as xr
    from control_audit import audit
    audit(cache)
    d=xr.open_dataarray(cache/'selected/1D_singles.nc').sel(key=30)
    a=d.values;a=a[:,~np.all(np.isnan(a),axis=(0,2,3))];mask=~np.isnan(a[1,0]);ys,xs=np.where(mask)
    a=np.nan_to_num(a);ct=a[1,:,ys.min():ys.max()+1,xs.min():xs.max()+1].sum(axis=0).astype(int).ravel()
    assert a.shape[1]==931 and np.all(a[1].sum(axis=(1,2))<=1)
    result=[dict(time_ms=2.45,input_y=0,y_start=int(d.y.values[ys.min()])-25,x_start=int(d.x.values[xs.min()])-24,width=12,shots=931,counts=ct.tolist()+[931-int(ct.sum())])]
    c=json.loads((HERE/'results/matched-regions.json').read_text())['cases'][-1];assert c['time_ms']==4.65
    for j,r in enumerate(c['singletons']):
        result.append(dict(time_ms=4.65,input_y=j-1,y_start=-9,x_start=-8,width=18,shots=r['shots'],counts=r['counts'],intervals=r['intervals']))
    return result


def prediction(p,records,order=31,half=18):
    jx,jy,j2x,j2y,fx,fy,kx,ky,off,sig,eta=p
    coord=np.arange(-half,half+1);n=len(coord);z,w=hermgauss(order);w/=np.sqrt(np.pi)
    times=sorted(set(r['time_ms'] for r in records));t=(np.array(times)[:,None]+off)*(1+sig*np.sqrt(2)*z[None,:])/1000
    def axis(j,j2,f,k,inputs):
        H=np.diag(f*coord+k*coord**2)-j*(np.eye(n,k=1)+np.eye(n,k=-1))-j2*(np.eye(n,k=2)+np.eye(n,k=-2))
        en,V=eigh(H)
        return np.einsum('oe,tze,je->tzjo',V,np.exp(-2j*np.pi*t[:,:,None]*en),V[np.array(inputs)+half],optimize=True)
    px=abs(axis(jx,j2x,fx,kx,[0])[:,:,0])**2;py=abs(axis(jy,j2y,fy,ky,[-1,0,1]))**2
    result=[]
    for r in records:
        i=times.index(r['time_ms']);j=r['input_y']+1;W=r['width'];ys=r['y_start']+half;xs=r['x_start']+half
        assert 0<=ys<=n-W and 0<=xs<=n-W
        q=eta*np.einsum('zy,zx,z->yx',py[i,:,j,ys:ys+W],px[i,:,xs:xs+W],w).ravel()
        q=np.r_[q,1-q.sum()];assert np.all(q>=0)
        result.append(q)
    return result


def fit(records):
    bounds=list(BOUNDS);bounds[8]=(max(-1.5,-min(r['time_ms'] for r in records)+.01),.2)
    def loss(v):return -sum(np.sum(xlogy(r['counts'],np.maximum(q,1e-200))) for r,q in zip(records,prediction(v*SCALE,records,11)))
    candidates=[]
    for seed in range(3):
        start=np.array([122,107,7,6,0,0,2,2,-.374,.08,.963])
        if seed:start[4:8]=[5*seed,-7*seed,seed,seed*2]
        f=minimize(loss,start/SCALE,method='L-BFGS-B',bounds=[(a/s,b/s) for (a,b),s in zip(bounds,SCALE)],options={'maxiter':600,'ftol':1e-11,'gtol':1e-5,'maxls':30})
        candidates.append(dict(parameters=(f.x*SCALE).tolist(),nll=float(f.fun),solver_success=bool(f.success)))
    return min(candidates,key=lambda r:r['nll'])


def evaluate(records,parameters,replicates=20000,predictor=None,orders=(31,41),halves=(18,24)):
    predictor=prediction if predictor is None else predictor
    qs=predictor(parameters,records,orders[0],halves[0]);fine=predictor(parameters,records,orders[1],halves[1])
    discrepancy=max(float(np.max(abs(a-b))) for a,b in zip(qs,fine));assert discrepancy<1e-8
    statistic=0.;failures=[]
    for i,(r,q) in enumerate(zip(records,qs)):
        ct=np.array(r['counts']);statistic+=2*np.sum(xlogy(ct,ct/np.maximum(r['shots']*q,1e-250)))
        if 'intervals' in r:
            for k,(n,p) in enumerate(zip(ct,q)):
                lo,hi=map(lambda x:float(F(x)),r['intervals'][str(n)])
                if not lo<=p<=hi:failures.append(dict(record=i,time_ms=r['time_ms'],input_y=r['input_y'],output=k,count=int(n),probability=float(p),lower=lo,upper=hi))
        else:
            assert r['time_ms']==2.45 and r['input_y']==0 and r['width']==12
            cert=json.loads((HERE/'results/cell-certificate.json').read_text());grid=q[:-1].reshape(12,12)
            checks=[('row:'+str(y),float(grid[y].sum()),cert['row_intervals'][y]) for y in range(12)]
            checks.append(('loss',float(q[-1]),cert['row_intervals'][-1]))
            checks.extend(('conditional:'+str(y)+','+str(x),float(grid[y,x]/grid[y].sum()),cert['conditional_intervals'][y][x]) for y in range(12) for x in range(12))
            for label,p,interval in checks:
                lo,hi=map(lambda x:float(F(x)),interval)
                if not lo<=p<=hi:failures.append(dict(record=i,time_ms=r['time_ms'],input_y=0,output=label,probability=p,lower=lo,upper=hi))
    # A point-null simulation diagnostic, not guessed fitted degrees of freedom.
    rng=np.random.default_rng(927);sim=np.zeros(replicates)
    for r,q in zip(records,qs):
        draws=rng.multinomial(r['shots'],q,size=replicates)
        sim+=2*np.sum(xlogy(draws,draws/np.maximum(r['shots']*q,1e-250)),axis=1)
    exceed=int(np.sum(sim>=statistic))
    upper=1. if exceed==replicates else float(beta.isf(.01,exceed+1,replicates-exceed))
    return dict(shots=sum(r['shots'] for r in records),independent_records=len(records),parameters=parameters,
        multinomial_deviance=float(statistic),retained_cell_interval_failures=failures,
        numerical_domain_and_quadrature_difference=discrepancy,
        monte_carlo=dict(replicates=replicates,seed=927,exceedances=exceed,tail_probability_upper_99=upper),
        status='Deposited point diagnostic only; neither a global fit optimum nor a physical-family exclusion.')


def main():
    p=argparse.ArgumentParser();p.add_argument('--cache',type=Path);p.add_argument('--fit',action='store_true');p.add_argument('--check',action='store_true');a=p.parse_args()
    source=HERE/'results/july2-dynamics-records.json'
    if a.cache:
        records=july2_records(a.cache)
        if a.check:assert records==json.loads(source.read_text())
        else:source.write_text(json.dumps(records,sort_keys=True,separators=(',',':'))+'\n')
    records={'July2':json.loads(source.read_text()),'July7':july7_records()}
    proposal=HERE/'results/extended-dynamics-parameters.json'
    if a.fit:
        params={date:fit(rs) for date,rs in records.items()};proposal.write_text(json.dumps(params,sort_keys=True,indent=2)+'\n')
    params=json.loads(proposal.read_text());result={date:evaluate(rs,params[date]['parameters']) for date,rs in records.items()}
    out=HERE/'results/extended-dynamics.json'
    if a.check:close(result,json.loads(out.read_text()))
    else:out.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
