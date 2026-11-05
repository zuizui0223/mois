"""Versioned experimental diel evidence and same-year macro source contracts."""
from pathlib import Path
import json,hashlib
import numpy as np
import pandas as pd

def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def name_lookup(names,accepted):
    n=names.loc[names.wcvp_match_state.eq('unique_exact_name_to_accepted_id')]
    pairs=pd.concat([n[['source_species_exact_binomial','wcvp_accepted_plant_name_id']].rename(columns={'source_species_exact_binomial':'name'}),accepted[['wcvp_accepted_name','wcvp_accepted_plant_name_id']].rename(columns={'wcvp_accepted_name':'name'})]).dropna().drop_duplicates()
    groups=pairs.groupby('name').wcvp_accepted_plant_name_id.agg(set)
    return {name:next(iter(ids)) for name,ids in groups.items() if len(ids)==1}

def resolve_species(aliases,lookup):
    ids={lookup[x] for x in aliases if x in lookup}
    return next(iter(ids)) if len(ids)==1 else None

def night_share(night,day):
    total=night+day
    return night.div(total.where((total>0)&(night>=0)&(day>=0)))

def join_exposure(x,e):
    if e.duplicated(['site_id','year']).any():raise ValueError('Duplicate same-year exposure')
    return x.merge(e,on=['site_id','year'],how='left',validate='many_to_one')

def species_weights(x):return 1/x.groupby('species').species.transform('size')

def build_species(root,out):
    root=Path(root);out=Path(out);out.mkdir(parents=True,exist_ok=True)
    names_path=root/'signed/taxonomy/mois350399_wcvp_source_taxon_exact_name_resolution.csv.gz'
    accepted_path=root/'signed/taxonomy/mois350399_wcvp_accepted_native_distribution_TDWG3.csv.gz'
    assert digest(names_path)=='dc5bb0e3686fde0cfd327e4336245d0ceae3c7657b3bd12841bb814cc94d1eb2'
    assert digest(accepted_path)=='9cd6aa7ca29f00c64283bcc08dcb778e212b644b837e64070051339b11d8b7b3'
    names=pd.read_csv(names_path,dtype=str);accepted=pd.read_csv(accepted_path,dtype=str)
    lookup=name_lookup(names,accepted)
    path=root/'signed/diel_figshare_v1/published_analysis_export.csv'
    assert digest(root/'signed/diel_figshare_v1/DvN_analysis_dataframe final.RData')=='b0aa1fa67fc484fa0000ed5f08119733d710d1c21e41c7ac0a018c1c50f6345a'
    # R text connections use platform line endings; bind identical CSV content.
    assert hashlib.sha256(path.read_bytes().replace(b'\r\n',b'\n')).hexdigest()=='77ac95c42ad876143aa2c7bd8cadc799f7b5987e31b589a55ae0aa0c340480d3', 'Frozen published export differs'
    d=pd.read_csv(path);d=d.loc[d.treatment.eq('day_night')].copy()
    d['night_share']=night_share(d['effectiveness_value_night.pollination'],d['effectiveness_value_day.pollination'])
    rows=[]
    for sp,x in d.groupby('accepted_name'):
        aliases=set(x.plant_species.dropna())|set(x.old_name.dropna())|{sp}
        wid=resolve_species(aliases,lookup)
        perstudy=x.groupby('study_ID').night_share.mean()
        colors=x.color.dropna().unique()
        white=float(colors[0]=='white') if len(colors)==1 else np.nan
        rows.append(dict(species=sp,wcvp_id=wid,aliases='|'.join(sorted(aliases)),family=x.family.iloc[0],genus=x.genus.iloc[0],flower_colour=colors[0] if len(colors)==1 else 'unknown',white=white,night_share=float(perstudy.mean()),night_share_between_study_sd=float(perstudy.std()) if len(perstudy)>1 else np.nan,n_studies=x.study_ID.nunique(),n_comparisons=len(x),night_positive_comparisons=int(x['effectiveness_value_night.pollination'].gt(0).sum()),studies='|'.join(sorted(x.study_ID.unique())),resolution='unique_exact_source_alias' if wid else 'unresolved_or_conflicting',source='10.6084/m9.figshare.26004538.v1',actual_pollinator_taxon='not_recorded_by_meta_analysis'))
    ledger=pd.DataFrame(rows)
    # Several published entities mapped to one accepted ID cannot become duplicate species.
    dup=ledger.wcvp_id.notna()&ledger.wcvp_id.duplicated(keep=False)
    ledger.loc[dup,'resolution']='multiple_published_species_same_accepted_id'
    ledger.to_csv(out/'species_evidence_selection_ledger.csv',index=False)
    d.to_csv(out/'published_day_night_comparisons.csv',index=False)
    valid=ledger.loc[ledger.resolution.eq('unique_exact_source_alias')&ledger.white.notna()&ledger.night_share.notna()].copy()
    valid.to_csv(out/'eligible_experimental_species.csv',index=False)
    print(json.dumps(dict(published_species=len(ledger),comparisons=len(d),studies=d.study_ID.nunique(),exact_eligible_species=len(valid),white_species=int(valid.white.sum()),nonwhite_species=int((1-valid.white).sum()))),flush=True)
    return valid,names

def select_source(root,out):
    root=Path(root);out=Path(out)
    valid,names=build_species(root,out)
    path=root/'signed/canonical/canonical_observations_climate.csv.gz'
    assert digest(path)=='f033661d76e07898a5c13199a55c5d0ce6e55353daa48b0058c456979e90644c'
    bridge=names.loc[names.wcvp_match_state.eq('unique_exact_name_to_accepted_id'),['inat_taxon_id','wcvp_accepted_plant_name_id']].rename(columns={'wcvp_accepted_plant_name_id':'wcvp_id'})
    assert bridge.inat_taxon_id.is_unique
    cols=['photo_id','observation_id','inat_taxon_id','source_subarm','strict9','strict9_white','observed_on','latitude','longitude','wc2_bio05_10m','wc2_bio12_10m','wc2_elev_2p5m']
    x=pd.read_csv(path,usecols=cols,dtype={'photo_id':str,'observation_id':str,'inat_taxon_id':str},low_memory=False)
    assert x.photo_id.is_unique and x.observation_id.is_unique
    x['year']=pd.to_numeric(x.observed_on.str[:4],errors='coerce')
    x=x.loc[x.source_subarm.eq('historical_main')&x.year.between(2014,2020)&x.latitude.between(-60,60)&x.longitude.between(-180,180)].copy()
    x=x.merge(bridge,on='inat_taxon_id',validate='many_to_one').merge(valid,on='wcvp_id',validate='many_to_one')
    x['year']=x.year.astype(int)
    x['gx']=np.floor((x.longitude+180)/.05).astype(int);x['gy']=np.floor((x.latitude+90)/.05).astype(int)
    x['site_id']=x.gx.astype(str)+'_'+x.gy.astype(str)
    support=x.groupby('species').agg(n_observations=('observation_id','nunique'),n_sites=('site_id','nunique'),n_years=('year','nunique'),n_strict9=('strict9','sum')).reset_index()
    support['eligible']=(support.n_observations>=20)&(support.n_sites>=10)&(support.n_years>=2)
    support.to_csv(out/'species_photo_support_before_exposure.csv',index=False)
    x=x.loc[x.species.isin(support.loc[support.eligible,'species'])].copy()
    x.to_csv(out/'selected_original_observations.csv.gz',index=False)
    # Collapse repeated photographs to one species/site/year, no ROI-success gate.
    cells=x.groupby(['species','site_id','year'],as_index=False).agg(n_observations=('observation_id','size'),longitude=('gx','first'),latitude=('gy','first'),bio5=('wc2_bio05_10m','mean'),rain=('wc2_bio12_10m','mean'),elev=('wc2_elev_2p5m','mean'),roi_white=('strict9_white','mean'),roi_success=('strict9','sum'))
    cells['longitude']=(cells.longitude+.5)*.05-180;cells['latitude']=(cells.latitude+.5)*.05-90
    cells=cells.merge(valid,on='species',validate='many_to_one')
    cells.to_csv(out/'species_site_year_source.csv.gz',index=False)
    sites=cells[['site_id','year','longitude','latitude']].drop_duplicates()
    assert not sites.duplicated(['site_id','year']).any()
    sites.to_csv(out/'exposure_queue.csv',index=False)
    receipt=dict(selected_original_observations=len(x),selected_species=x.species.nunique(),species_site_year_rows=len(cells),unique_site_year_requests=len(sites),year_counts=x.year.value_counts().sort_index().to_dict(),selection_rule='>=20 observations, >=10 distinct 0.05 degree sites, >=2 years; no ALAN values or ROI success used',cohort_years=[2014,2020],canonical_sha256=digest(path),data_frame_not_population_census=True)
    (out/'source_selection_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt),flush=True)

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,default=Path('.'));p.add_argument('--out',type=Path,required=True);p.add_argument('--species-only',action='store_true');a=p.parse_args()
    (build_species if a.species_only else select_source)(a.root,a.out)
