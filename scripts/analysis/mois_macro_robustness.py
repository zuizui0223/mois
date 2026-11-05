"""Post-primary exhaustive influence checks; no effect-selected subsetting."""
from pathlib import Path
import argparse,json
import numpy as np
import pandas as pd
from scripts.analysis.mois_macro_fit import fit

PREDICTORS=['log_alan','log_population','temperature_c','log_precipitation','log_effort','log_richness']

def transform(d):
    d=d.copy();d['log_alan']=np.log1p(d.radiance.clip(lower=0));d['log_population']=np.log1p(d.population_density)
    d['log_precipitation']=np.log1p(d.precipitation_mm);d['log_effort']=np.log1p(d.n_precise_observations);d['log_richness']=np.log(d.n_species)
    return d

def remove_genus(d,g,genus):
    c=g.loc[g.genus.eq(genus),['cell025','year','n_species','n_white_species']].rename(columns={'n_species':'removed_species','n_white_species':'removed_white'})
    x=d.merge(c,on=['cell025','year'],how='left',validate='one_to_one')
    present=x.removed_species.notna()
    x['n_species']-=x.removed_species.fillna(0);x['n_white_species']-=x.removed_white.fillna(0);x['n_genera']-=present.astype(int)
    x=x.loc[x.n_species.ge(10)&x.n_genera.ge(5)].copy()
    assert x.n_white_species.between(0,x.n_species).all()
    x['white_species_fraction']=x.n_white_species/x.n_species
    return transform(x)

def run(source,genera,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    d=pd.read_csv(source);d=transform(d.loc[d.environment_complete]);s=d.loc[d.latitude.lt(0)].copy()
    s['sector']=np.select([s.longitude.lt(-30),s.longitude.lt(60)],['west_of_minus30','minus30_to_60'],default='east_of_60')
    support=s.groupby('sector').agg(n_rows=('year','size'),n_cells=('cell025','nunique'),n_blocks=('block5','nunique')).reset_index();support.to_csv(out/'southern_sector_support.csv',index=False)
    rows=[];failures=[]
    def estimate(name,v,predictors=PREDICTORS,groups=None):
        try:
            groups=groups or [v.region10,v.year]
            b,se,lo,hi=fit(100*v.white_species_fraction,v[predictors].to_numpy(),groups,[v.block5])
            for term,bb,ss,ll,hh in zip(predictors,b,se,lo,hi):rows.append({'analysis':name,'term':term,'estimate':bb,'se':ss,'ci95_low':ll,'ci95_high':hh,'n_rows':len(v),'n_cells':v.cell025.nunique(),'n_blocks':v.block5.nunique()})
        except ValueError as e:failures.append({'analysis':name,'n_rows':len(v),'error':str(e)})
        pd.DataFrame(rows).to_csv(out/'all_sensitivity_coefficients.csv',index=False)
        (out/'unsupported_models.json').write_text(json.dumps(failures,indent=2)+'\n')
    estimate('southern_baseline',s)
    d['south']=d.latitude.lt(0).astype(int);interaction=[]
    for p in PREDICTORS:
        name=p+'_south';d[name]=d[p]*d.south;interaction.append(name)
    estimate('hemisphere_slope_interactions',d,PREDICTORS+interaction,[d.region10,d.year.astype(str)+'_'+d.south.astype(str)])
    for sector in sorted(s.sector.unique()):
        estimate('only_sector_'+sector,s.loc[s.sector.eq(sector)])
        estimate('omit_sector_'+sector,s.loc[~s.sector.eq(sector)])
    for year in range(2015,2025):estimate('omit_year_'+str(year),s.loc[s.year.ne(year)])
    g=pd.read_csv(genera);g=g.merge(s[['cell025','year']],on=['cell025','year'],validate='many_to_one')
    ranking=g.groupby('genus').size().sort_values(ascending=False,kind='stable').rename('southern_cell_year_presence').reset_index()
    ranking.to_csv(out/'genus_presence_ranking_before_influence.csv',index=False)
    for genus in ranking.head(10).genus:estimate('omit_genus_'+genus,remove_genus(s,g,genus))
    r={'baseline_rows':len(d),'southern_rows':len(s),'attempted_models':len(set(x['analysis'] for x in rows))+len(failures),'successful_models':len(set(x['analysis'] for x in rows)),'unsupported_models':failures,'selection':'all3 longitude sectors, all10 years, top10 genera by support; no effect-driven choice','stage':'post-primary sensitivity, not preregistered confirmatory tests','multiple_comparisons':'individual confidence intervals; do not select significant rows as new primary'}
    (out/'robustness_receipt.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',required=True);p.add_argument('--genera',required=True);p.add_argument('--out',required=True);a=p.parse_args();run(a.source,a.genera,a.out)
