"""Portable numerical utilities. No network, installation, or implicit source paths."""
import numpy as np
import pandas as pd
from scipy import stats

def bh(values):
    p=np.asarray(values,float);out=np.full(p.shape,np.nan);ok=np.flatnonzero(np.isfinite(p));order=ok[np.argsort(p[ok],kind='stable')]
    if len(order):out[order]=np.minimum.accumulate((p[order]*len(order)/np.arange(1,len(order)+1))[::-1])[::-1].clip(0,1)
    return out

def corr(x,y,cov=None,rank_covariates=True):
    x=np.asarray(x,float);y=np.asarray(y,float);c=np.empty((len(x),0)) if cov is None else np.asarray(cov,float)
    if c.ndim==1:c=c[:,None]
    ok=np.isfinite(x)&np.isfinite(y)&np.isfinite(c).all(1);x=x[ok];y=y[ok];c=c[ok];n=len(x);k=c.shape[1]
    if n<20 or np.std(x)==0 or np.std(y)==0:return {'n':n,'rho':np.nan,'p':np.nan}
    rx=stats.rankdata(x);ry=stats.rankdata(y)
    if k:
        if rank_covariates:c=pd.DataFrame(c).rank().to_numpy()
        C=np.column_stack([np.ones(n),c])
        if np.linalg.matrix_rank(C)!=C.shape[1]:raise ValueError('Rank deficient correlation covariates')
        rx=rx-C@np.linalg.lstsq(C,rx,rcond=None)[0];ry=ry-C@np.linalg.lstsq(C,ry,rcond=None)[0]
    r=float(np.corrcoef(rx,ry)[0,1]);df=n-k-2;p=float(2*stats.t.sf(abs(r)*np.sqrt(df/max(1e-16,1-r*r)),df))
    return {'n':n,'rho':r,'p':p}

def ols(y,x,cov):
    y=np.asarray(y,float);x=np.asarray(x,float);c=np.asarray(cov,float)
    if c.ndim==1:c=c[:,None]
    ok=np.isfinite(y)&np.isfinite(x)&np.isfinite(c).all(1);y=y[ok];x=x[ok];c=c[ok];n=len(y)
    if n<30 or np.std(y)==0 or np.std(x)==0:return {'n':n,'beta':np.nan,'p':np.nan}
    X=np.column_stack([np.ones(n),x,c]);rank=np.linalg.matrix_rank(X)
    if rank!=X.shape[1]:raise ValueError('Rank deficient OLS design')
    b=np.linalg.lstsq(X,y,rcond=None)[0];e=y-X@b;df=n-rank;se=float(np.sqrt((e@e/df)*np.linalg.inv(X.T@X)[1,1]));p=float(2*stats.t.sf(abs(b[1]/se),df));q=stats.t.ppf(.975,df)
    C=np.column_stack([np.ones(n),c]);be=y-C@np.linalg.lstsq(C,y,rcond=None)[0];sst=((y-y.mean())**2).sum()
    return {'n':n,'beta':float(b[1]),'se':se,'lower':float(b[1]-q*se),'upper':float(b[1]+q*se),'p':p,'std_beta':float(b[1]*np.std(x)/np.std(y)),'R2_base':float(1-be@be/sst),'R2_full':float(1-e@e/sst),'delta_R2':float((be@be-e@e)/sst),'condition_number':float(np.linalg.cond(X))}

def zrows(frame):return frame.sub(frame.mean(axis=1),axis=0).div(frame.std(axis=1,ddof=0).replace(0,np.nan),axis=0)

def frozen_score(frame,weights,negative=False):
    """Original scoring rule: gene z scores, global weights, missing z contribution zero."""
    weights=weights/weights.sum();z=zrows(frame.reindex(weights.index));out=z.mul(weights,axis=0).sum(axis=0,min_count=1)
    return -out if negative else out
