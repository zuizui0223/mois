"""Keep the original ROI9 version as a distinct, same-year image supplement."""
from pathlib import Path
import argparse,json,hashlib
import numpy as np
import pandas as pd
from scripts.analysis.mois_macro_fit import fit

def aggregate(d):
    keys=['cell025','year']
    c=d.groupby(keys).agg(n_photos=('observation_id','size'),n_valid_roi=('strict9','sum'),roi_white_photo_weighted=('strict9_white','mean'),n_nominal_taxa=('inat_taxon_id','nunique')).reset_index()
    within=d.loc[d.strict9].groupby(keys+['inat_taxon_id']).strict9_white.mean().reset_index()
    balanced=within.groupby(keys).agg(roi_white_fraction=('strict9_white','mean'),n_roi_taxa=('inat_taxon_id','nunique')).reset_index()
    c=c.merge(balanced,on=keys,how='left',validate='one_to_one');c['n_roi_taxa']=c.n_roi_taxa.fillna(0).astype(int)
    return c

def assemble(source,broad,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    assert hashlib.sha256(Path(source).read_bytes()).hexdigest()=='f033661d76e07898a5c13199a55c5d0ce6e55353daa48b0058c456979e90644c'
    cols=['observation_id','inat_taxon_id','source_subarm','observed_on','latitude','longitude','strict9','strict9_white']
    d=pd.read_csv(source,usecols=cols);assert d.observation_id.is_unique
    d['year']=pd.to_numeric(d.observed_on.str[:4],errors='coerce')
    d=d.loc[d.source_subarm.eq('historical_main')&d.year.between(2015,2024)&d.latitude.between(-60,60)&d.longitude.ge(-180)&d.longitude.lt(180)].copy()
    d['year']=d.year.astype(int);d['cell025']=np.floor((d.longitude+180)/.25).astype(int).astype(str)+'_'+np.floor((d.latitude+90)/.25).astype(int).astype(str)
    assert d.loc[d.strict9,'strict9_white'].between(0,1).all()
    d.loc[~d.strict9,'strict9_white']=np.nan
    keys=['cell025','year'];b=pd.read_csv(broad,usecols=keys)
    d=d.merge(b,on=keys,validate='many_to_one')
    c=aggregate(d)
    c['roi_success_fraction']=c.n_valid_roi/c.n_photos
    c['white_eligible']=c.n_valid_roi.ge(5)&c.n_roi_taxa.ge(3)
    c['success_eligible']=c.n_photos.ge(10)&c.n_nominal_taxa.ge(3)
    c.to_csv(out/'roi_cell_year_source.csv',index=False)
    r={'version':'original strict9 white-area fraction','source_arm':'historical_main only','source_observations_in_broad_frame':len(d),'source_nominal_taxa':d.inat_taxon_id.nunique(),'valid_roi_observations':int(d.strict9.sum()),'cell_years':len(c),'white_eligible_cell_years':int(c.white_eligible.sum()),'success_eligible_cell_years':int(c.success_eligible.sum()),'coordinate_accuracy':'not available in canonical image table; this supplement does not inherit the precise-coordinate primary cohort claim','unknown_roi_is_not_nonwhite':True,'interpretation':'ROI white area is not botanical species colour, flowering confirmation, or abundance'}
    (out/'roi_source_receipt.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r),flush=True)

def model(source,environment,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True);d=pd.read_csv(source);e=pd.read_csv(environment)
    keys=['cell025','year'];cols=keys+['radiance','population_density','temperature_c','precipitation_mm','n_precise_observations','block5','region10']
    d=d.merge(e[cols],on=keys,how='left',validate='one_to_one')
    d['environment_complete']=d[cols].notna().all(axis=1)
    d.to_csv(out/'roi_candidates_with_environment.csv',index=False)
    d=d.loc[d.environment_complete].copy();d['log_alan']=np.log1p(d.radiance.clip(lower=0));d['log_population']=np.log1p(d.population_density)
    d['log_precipitation']=np.log1p(d.precipitation_mm);d['log_effort']=np.log1p(d.n_precise_observations);d['log_taxa']=np.log(d.n_nominal_taxa)
    predictors=['log_alan','log_population','temperature_c','log_precipitation','log_effort','log_taxa'];rows=[]
    for name,outcome,gate in [('original_ROI9_taxon_balanced','roi_white_fraction','white_eligible'),('ROI9_photo_weighted_sensitivity','roi_white_photo_weighted','white_eligible'),('technical_ROI9_success','roi_success_fraction','success_eligible')]:
        v=d.loc[d[gate]]
        b,se,lo,hi=fit(v[outcome].to_numpy()*100,v[predictors].to_numpy(),[v.region10,v.year],[v.block5])
        for term,bb,ss,ll,hh in zip(predictors,b,se,lo,hi):rows.append({'model':name,'term':term,'estimate':bb,'se':ss,'ci95_low':ll,'ci95_high':hh,'n_cell_years':len(v),'n_blocks5':v.block5.nunique()})
    pd.DataFrame(rows).to_csv(out/'roi_sameyear_coefficients.csv',index=False)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',required=True);p.add_argument('--broad');p.add_argument('--environment');p.add_argument('--out',required=True);a=p.parse_args()
    model(a.source,a.environment,a.out) if a.environment else assemble(a.source,a.broad,a.out)
