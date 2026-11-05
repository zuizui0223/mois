"""Broad documented-species assemblage associations; no absence/causal estimand."""
from pathlib import Path
import argparse,json
import numpy as np
import pandas as pd
from scipy.stats import t

def absorb(a,groups,weights,tol=1e-10):
    a=np.asarray(a,float).copy(); codes=[pd.factorize(g)[0] for g in groups]
    for iteration in range(10000):
        old=a.copy()
        for code in codes:
            den=np.bincount(code,weights=weights)
            for j in range(a.shape[1]):
                means=np.bincount(code,weights=weights*a[:,j])/den
                a[:,j]-=means[code]
        if np.max(np.abs(a-old))<tol:return a
    raise ValueError('Fixed effects failed to converge')

def fit(y,x,groups,clusters,weights=None):
    n,k=x.shape;w=np.ones(n) if weights is None else np.asarray(weights,float)
    z=absorb(np.column_stack([y,x]),groups,w)
    xx=z[:,1:]*np.sqrt(w[:,None]); yy=z[:,0]*np.sqrt(w)
    if np.linalg.matrix_rank(xx)<k:raise ValueError('Unidentified model')
    beta=np.linalg.lstsq(xx,yy,rcond=None)[0];e=yy-xx@beta
    bread=np.linalg.inv(xx.T@xx);score=xx*e[:,None]
    # Conservative upper bound on nuisance rank for disconnected FE graphs.
    df=n-k-1-sum(pd.Series(g).nunique()-1 for g in groups)
    if df<=0:raise ValueError('No residual degrees of freedom')
    def meat(g):
        code,labels=pd.factorize(g);ng=len(labels)
        if ng<3:raise ValueError('Insufficient independent clusters')
        sums=np.stack([np.bincount(code,weights=score[:,j]) for j in range(k)],axis=1)
        return sums.T@sums*(ng/(ng-1))*((n-1)/df),ng
    m,ng=meat(clusters[0]);dft=ng-1
    if len(clusters)==2:
        second,n2=meat(clusters[1]);joint,_=meat(pd.Series(list(zip(clusters[0],clusters[1]))))
        m=m+second-joint;dft=min(ng,n2)-1
    cov=bread@m@bread
    if np.any(np.diag(cov)<0):raise ValueError('Nonpositive clustered variance')
    se=np.sqrt(np.diag(cov));critical=t.ppf(.975,dft)
    return beta,se,beta-critical*se,beta+critical*se

def run(source,environment,climate,out):
    source=Path(source);out=Path(out);out.mkdir(parents=True,exist_ok=True)
    d=pd.read_csv(source/'primary_cell_year_source.csv.gz')
    e=pd.read_csv(environment);c=pd.read_csv(climate)
    keys=['cell025','year']
    if e.duplicated(keys).any() or c.duplicated(keys).any():raise ValueError('Duplicate environmental keys')
    d=d.merge(e,on=keys,how='left',validate='one_to_one').merge(c,on=keys,how='left',validate='one_to_one')
    cols=['radiance','population_density','temperature_c','precipitation_mm']
    complete=d[cols].notna().all(axis=1)&d.population_density.ge(0)&d.precipitation_mm.ge(0)
    d['environment_complete']=complete
    d.to_csv(out/'all_primary_rows_with_environment_and_missingness.csv.gz',index=False)
    receipt={'source_rows':len(d),'complete_rows':int(complete.sum()),'excluded_year_counts':d.loc[~complete].year.value_counts().sort_index().to_dict(),'estimand':'percentage of documented species that have single-state white flowers; association conditional on sampled frame','inference':'spatial clusters; genus sensitivity additionally clusters genera; not a complete phylogenetic model'}
    d=d.loc[complete].copy();d['log_alan']=np.log1p(d.radiance.clip(lower=0));d['log_population']=np.log1p(d.population_density)
    d['log_precipitation']=np.log1p(d.precipitation_mm);d['log_effort']=np.log1p(d.n_precise_observations);d['log_richness']=np.log(d.n_species)
    predictors=['log_alan','log_population','temperature_c','log_precipitation','log_effort','log_richness']
    rows=[];unsupported=[]
    def estimate(label,frame,response,fe,clusters,weights=None):
        try:b,s,lo,hi=fit(response,frame[predictors].to_numpy(),[frame[g] for g in fe],[frame[g] for g in clusters],weights)
        except ValueError as error:
            if label=='primary_region_year':raise
            unsupported.append({'model':label,'reason':str(error),'n_rows':len(frame)})
            (out/'unsupported_sensitivities.json').write_text(json.dumps(unsupported,indent=2)+'\n')
            return
        for name,beta,se,l,h in zip(predictors,b,s,lo,hi):rows.append({'model':label,'term':name,'estimate':beta,'se':se,'ci95_low':l,'ci95_high':h,'n_rows':len(frame),'n_cells':frame.cell025.nunique(),'n_blocks5':frame.block5.nunique()})
        pd.DataFrame(rows).to_csv(out/'coefficients.csv',index=False)
    estimate('primary_region_year',d,100*d.white_species_fraction,['region10','year'],['block5'])
    for hemi in ['north','south']:
        v=d.loc[d.latitude.ge(0) if hemi=='north' else d.latitude.lt(0)]
        estimate(hemi,v,100*v.white_species_fraction,['region10','year'],['block5'])
    repeated=d.groupby('cell025').year.transform('nunique').ge(3);v=d.loc[repeated]
    estimate('within_cell_year',v,100*v.white_species_fraction,['cell025','year'],['block5'])
    g=pd.read_csv(source/'all_documented_genus_cell_years.csv.gz')
    g=g.rename(columns={'n_species':'genus_n_species','n_white_species':'genus_n_white_species'}).merge(d,on=keys,validate='many_to_one')
    # Each cell-year has total weight 1, irrespective of genus richness.
    weights=1/g.groupby(keys).genus.transform('size')
    estimate('within_genus_region_year',g,100*g.genus_n_white_species/g.genus_n_species,['genus','region10','year'],['block5','genus'],weights)
    pd.DataFrame(rows).to_csv(out/'coefficients.csv',index=False)
    receipt['complete_year_counts']=d.year.value_counts().sort_index().to_dict()
    receipt['unsupported_sensitivities']=unsupported
    receipt['all_sensitivities_estimated']=not unsupported
    receipt['complete_cells']=d.cell025.nunique();receipt['complete_blocks']=d.block5.nunique()
    (out/'model_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',required=True);p.add_argument('--environment',required=True);p.add_argument('--climate',required=True);p.add_argument('--out',required=True);a=p.parse_args();run(a.source,a.environment,a.climate,a.out)
