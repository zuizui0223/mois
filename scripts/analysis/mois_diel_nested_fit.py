"""Species-balanced sampled ALAN association with colour and diel treatment ratio."""
from pathlib import Path
import argparse,json
import numpy as np
import pandas as pd
from scripts.analysis.mois_macro_fit import fit

def run(source,environment,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    d=pd.read_csv(source);e=pd.read_csv(environment)
    keys=['cell025','year'];cols=keys+['radiance','population_density','temperature_c','precipitation_mm','n_precise_observations','block5','region10']
    if e.duplicated(keys).any():raise ValueError('Duplicate environmental keys')
    d=d.merge(e[cols],on=keys,how='left',validate='many_to_one')
    raw=len(d)
    required=cols+['white','night_share']
    d['environment_complete']=d[required].notna().all(axis=1)&d.population_density.ge(0)&d.precipitation_mm.ge(0)
    d['missing_fields']=d[required].isna().apply(lambda row:'|'.join(row.index[row]),axis=1)
    d.to_csv(out/'experimental_candidate_join_and_exclusions.csv.gz',index=False)
    support=d.groupby('species').agg(candidate_rows=('year','size'),retained_rows=('environment_complete','sum'),candidate_cells=('cell025','nunique')).reset_index()
    support['excluded_rows']=support.candidate_rows-support.retained_rows
    support.to_csv(out/'experimental_environment_support_by_species.csv',index=False)
    d=d.loc[d.environment_complete].copy()
    d['log_population']=np.log1p(d.population_density);d['log_precipitation']=np.log1p(d.precipitation_mm)
    d['log_effort']=np.log1p(d.n_precise_observations);d['night_centered']=d.night_share-.5;d['white_night']=d.white*d.night_centered
    names=['white','night_centered','white_night','log_population','temperature_c','log_precipitation','log_effort']
    # Every retained experimental species contributes the same total weight.
    w=1/d.groupby('species').species.transform('size')
    b,se,lo,hi=fit(np.log1p(d.radiance.clip(lower=0)),d[names].to_numpy(),[d.region10,d.year],[d.block5,d.species],w)
    pd.DataFrame({'term':names,'estimate':b,'se':se,'ci95_low':lo,'ci95_high':hi}).to_csv(out/'experimental_nested_coefficients.csv',index=False)
    d.to_csv(out/'experimental_nested_analysis_rows.csv.gz',index=False)
    r={'candidate_species_cell_years':raw,'complete_species_cell_years':len(d),'complete_species':d.species.nunique(),'complete_blocks5':d.block5.nunique(),'response':'log1p sampled annual VIIRS radiance at recorded species cells','weight':'equal total weight per experimental species','fixed_effects':['region10','year'],'clusters':['block5','species'],'night_share':'relative night/(night+day) treatment reproductive outcome, centered on 0.5; not exclusion-corrected pollination dependence','claim_ceiling':'association of sampled habitat light with literature colour and treatment ratio; neither local nocturnal pollination nor causal effects; no full phylogenetic correction'}
    (out/'experimental_nested_model_receipt.json').write_text(json.dumps(r,indent=2)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',required=True);p.add_argument('--environment',required=True);p.add_argument('--out',required=True);a=p.parse_args();run(a.source,a.environment,a.out)
