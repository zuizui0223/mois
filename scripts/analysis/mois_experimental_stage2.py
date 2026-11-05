"""Control-matched diel definitions and covariance-aware species profiles."""
from pathlib import Path
import argparse,json,hashlib
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scripts.analysis.mois_macro_fit import fit

def controlled_ratios(night,day,closed):
    denominator=(night+day).where((night>=0)&(day>=0)&(closed>=0)&((night+day)>0))
    return night/denominator,(night-closed)/denominator

def run(source,published,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True);d=pd.read_csv(source);p=pd.read_csv(published)
    assert hashlib.sha256(Path(published).read_bytes().replace(b'\r\n',b'\n')).hexdigest()=='77ac95c42ad876143aa2c7bd8cadc799f7b5987e31b589a55ae0aa0c340480d3'
    p=p.loc[p.treatment.eq('day_night')&p.accepted_name.isin(d.species)].copy()
    p['raw_control_matched'],p['net_night_relative_outcome']=controlled_ratios(p['effectiveness_value_night.pollination'],p['effectiveness_value_day.pollination'],p['effectiveness_value_complete.exclusion'])
    p.to_csv(out/'all_experimental_comparisons_with_control_status.csv',index=False)
    means=p.groupby(['accepted_name','study_ID'])[['raw_control_matched','net_night_relative_outcome']].mean().groupby('accepted_name').mean().reset_index().rename(columns={'accepted_name':'species'})
    phy=p[['accepted_name','phylo']].drop_duplicates().rename(columns={'accepted_name':'species'});assert phy.species.is_unique and phy.phylo.is_unique
    traits=d[['species','white','night_share']].drop_duplicates();assert traits.species.is_unique
    traits=traits.merge(means,on='species',validate='one_to_one').merge(phy,on='species',validate='one_to_one')
    traits.to_csv(out/'experimental_control_and_phylogeny_ledger.csv',index=False)
    matched=set(traits.loc[traits.net_night_relative_outcome.notna(),'species']);rows=[]
    variants=[('original_all','night_share',set(traits.species)),('original_control_species','night_share',matched),('raw_control_matched','raw_control_matched',matched),('net_control_matched','net_night_relative_outcome',matched)]
    for label,column,species in variants:
        v=d.loc[d.species.isin(species)].drop(columns=['night_centered','white_night'],errors='ignore').copy()
        mapping=traits.set_index('species')[column];v['night_centered']=v.species.map(mapping)-.5;v['white_night']=v.white*v.night_centered
        names=['white','night_centered','white_night','log_population','temperature_c','log_precipitation','log_effort']
        b,se,lo,hi=fit(np.log1p(v.radiance.clip(lower=0)),v[names].to_numpy(),[v.region10,v.year],[v.block5,v.species],1/v.groupby('species').species.transform('size'))
        for name,bb,ss,ll,hh in zip(names,b,se,lo,hi):rows.append({'analysis':label,'term':name,'estimate':bb,'se':ss,'ci95_low':ll,'ci95_high':hh,'n_species':v.species.nunique(),'n_rows':len(v)})
    pd.DataFrame(rows).to_csv(out/'control_definition_coefficients.csv',index=False)
    # Species intercepts adjusted for identical region/year/environment controls.
    # Their full spatially clustered covariance is passed to phylogenetic models.
    w=1/d.groupby('species').species.transform('size');species=sorted(d.species.unique())
    species_design=pd.get_dummies(pd.Categorical(d.species,categories=species),dtype=float)
    nuis=pd.get_dummies(d[['region10','year']].astype(str),drop_first=True,dtype=float).reset_index(drop=True)
    cols=['log_population','temperature_c','log_precipitation','log_effort'];controls=d[cols].reset_index(drop=True)
    centered=controls-np.average(controls,axis=0,weights=w)
    x=np.column_stack([species_design,nuis,centered]);y=np.log1p(d.radiance.clip(lower=0))
    model=sm.WLS(y.to_numpy(),x,weights=w.to_numpy()).fit(cov_type='cluster',cov_kwds={'groups':d.block5.to_numpy()})
    assert model.model.rank==x.shape[1]
    n=len(species);cov=np.asarray(model.cov_params())[:n,:n]
    if np.linalg.eigvalsh(cov).min()<=0:raise ValueError('Species sampling covariance is not positive definite')
    pd.DataFrame(cov,index=species,columns=species).to_csv(out/'adjusted_species_full_sampling_covariance.csv')
    profile=pd.DataFrame({'species':species,'adjusted_log_alan':np.asarray(model.params)[:n]}).merge(traits,on='species',validate='one_to_one')
    profile.to_csv(out/'adjusted_species_profiles.csv',index=False)
    receipt={'species':len(species),'control_matched_species':len(matched),'first_stage_rows':len(d),'first_stage_rank':int(model.model.rank),'sampling_covariance_min_eigenvalue':float(np.linalg.eigvalsh(cov).min()),'sampling_covariance':'full spatial-cluster covariance, not independent species SEs','net_definition':'(night-complete_exclusion)/(night+day); negative net contrasts retained; not a fractional nocturnal dependence probability','first_stage':'species fixed effects plus region/year/environment, centered covariates; second stage is a separate comparative sensitivity','trait_ratio_uncertainty':'not yet propagated as measurement error'}
    (out/'experimental_stage2_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',required=True);p.add_argument('--published',required=True);p.add_argument('--out',required=True);a=p.parse_args();run(a.source,a.published,a.out)
