#!/usr/bin/env python3
"""2020 MOIS original ROI9 x VIIRS quarterly radiance x WorldPop density.

Real, source-audited numeric records. No new photographs, no flower labels,
no MNAR scenarios. A cross-sectional photographic assemblage is not a census
of wild native plants or evidence of species disappearance.
"""
from __future__ import annotations
import argparse
import calendar
import concurrent.futures
import datetime as dt
import hashlib
import json
import math
import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests

CANONICAL_SHA = 'f033661d76e07898a5c13199a55c5d0ce6e55353daa48b0058c456979e90644c'
VIIRS='https://di-lawimagery1.img.arcgis.com/arcgis/rest/services/NighttimeLightsMDNB/ImageServer'
POP='https://worldpop.arcgis.com/arcgis/rest/services/WorldPop_Population_Density_1km/ImageServer'
YEAR=2020
MONTHS=(1,4,7,10)
GRID=.05


def check(cond,msg):
    if not bool(cond):raise ValueError('MOIS_ALAN_2020: '+msg)


def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for part in iter(lambda:f.read(4<<20),b''):h.update(part)
    return h.hexdigest()


def month_range_ms(y,m):
    a=dt.datetime(y,m,1,tzinfo=dt.timezone.utc)
    b=dt.datetime(y,m,calendar.monthrange(y,m)[1],23,59,59,tzinfo=dt.timezone.utc)
    return f'{int(a.timestamp()*1000)},{int(b.timestamp()*1000)}'


def select_cells(canonical,grid=GRID):
    x=canonical.copy()
    needed={'photo_id','observation_id','inat_taxon_id','source_subarm','strict9','strict9_white','observed_on','latitude','longitude','climate_complete','wc2_bio05_10m','wc2_bio12_10m','wc2_elev_2p5m'}
    check(needed<=set(x),'Canonical source columns missing')
    check(x.photo_id.is_unique and x.observation_id.is_unique,'Original photo / observation duplicates')
    strict=x.strict9.astype(str).str.lower().map({'true':True,'false':False})
    check(strict.notna().all(),'Invalid ROI status')
    x['strict9']=strict
    check(x.loc[~strict,'strict9_white'].isna().all(),'Unknown ROI9 must NOT have a colour')
    x['year']=pd.to_numeric(x.observed_on.astype(str).str.slice(0,4),errors='coerce')
    x=x.loc[x.source_subarm.eq('historical_main') & x.year.eq(YEAR) & x.climate_complete &
            x.latitude.between(-60,60) & x.longitude.between(-180,180)].copy()
    x=x.loc[np.isfinite(x[['wc2_bio05_10m','wc2_bio12_10m','wc2_elev_2p5m']].to_numpy(float)).all(axis=1)]
    check(len(x)>10000 and x.strict9.sum()>10000,'Insufficient frozen main original source')
    x['gx']=np.floor((x.longitude+180)/grid).astype(int)
    x['gy']=np.floor((x.latitude+90)/grid).astype(int)
    x['cell_id']=x.gx.astype(str)+'_'+x.gy.astype(str)
    old=x.groupby('cell_id',sort=True).agg(
        n_attempted=('photo_id','size'),n_attempted_taxa=('inat_taxon_id','nunique'),
        n_strict9=('strict9','sum'),
        bio5=('wc2_bio05_10m','mean'),bio12=('wc2_bio12_10m','mean'),
        elev=('wc2_elev_2p5m','mean'),
        gx=('gx','first'),gy=('gy','first')).reset_index()
    good=x.loc[x.strict9].copy()
    check(good.strict9_white.between(0,1).all(),'True ROI9 has invalid white fraction')
    by_taxon=good.groupby(['cell_id','inat_taxon_id'],observed=True).strict9_white.mean().reset_index()
    by_cell=by_taxon.groupby('cell_id',observed=True).agg(white_taxon_equal=('strict9_white','mean'),
                n_measured_taxa=('inat_taxon_id','nunique')).reset_index()
    old=old.merge(by_cell,on='cell_id',how='left',validate='one_to_one')
    # An observation-biased, ROI-success-selected cell audit. Not floristic absence.
    ok=(old.n_attempted>=5)&(old.n_attempted_taxa>=3)&(old.n_strict9>=3)&(old.n_measured_taxa>=3)
    out=old.loc[ok].copy()
    out['longitude']=(out.gx+0.5)*grid-180
    out['latitude']=(out.gy+0.5)*grid-90
    out['block_lon']=np.floor(out.longitude/5).astype(int)
    out['block_lat']=np.floor(out.latitude/5).astype(int)
    out['block_id']=out.block_lon.astype(str)+'_'+out.block_lat.astype(str)
    out['strict9_yield']=out.n_strict9/out.n_attempted
    check(out.cell_id.is_unique and len(out)>=200,'Analysis cells vanished')
    return out.sort_values('cell_id').reset_index(drop=True),{
        'photo_attempts_eligible_before_cell_gate':len(x),
        'strict9_measurements_before_cell_gate':int(x.strict9.sum()),
        'source_subarm':'historical_main','year':YEAR,'cell_size_degrees':grid,
        'all_2020_attempted_cells':int(len(old)),
        'selected_cells':int(len(out)),
        'selected_photo_attempts':int(out.n_attempted.sum()),
        'selected_strict9_success':int(out.n_strict9.sum()),
        'selected_nominal_taxa_in_cell_sum':int(out.n_measured_taxa.sum()),
        'cell_selection_conditions':['attempted_photos>=5','attempted_taxa>=3','strict9_photos>=3','strict9_taxa>=3'],
        'sampling_is_historical_opportunistic_photographed_subset':True,
        'global_wild_native_flora_inference':False}


def _parse_viirs_response(obj,y,m):
    if obj.get('error'):return None,None,'API_ERROR'
    value=str(obj.get('value') or '').strip()
    if value in ('','NoData','None','null'):
        return None,None,'NO_DATA'
    try:v=float(value)
    except ValueError:return None,None,'NOT_NUMERIC'
    items=obj.get('catalogItems') or {}
    features=items.get('features') or []
    if len(features)!=1:return None,None,'SOURCE_MONTH_NOT_SINGLE'
    a=features[0].get('attributes') or {}
    if int(a.get('year',-1))!=y or int(a.get('month',-1))!=m:
        return None,None,'SOURCE_MONTH_DRIFT'
    vals=(obj.get('properties') or {}).get('Values') or []
    if not vals:return None,None,'NO_CLOUD_FREE_COUNT'
    token=str(vals[0]).split()
    if len(token)<2:return None,None,'INVALID_CLOUD_COUNT'
    try:cloud=float(token[1]);radiance_from_props=float(token[0])
    except ValueError:return None,None,'INVALID_CLOUD_COUNT'
    if not np.isfinite(v) or not np.isfinite(radiance_from_props) or not np.isfinite(cloud):
        return None,None,'NONFINITE'
    if abs(v-radiance_from_props)>1e-3:return None,None,'PROPERTY_VALUE_MISMATCH'
    if cloud<2:return None,cloud,'INSUFFICIENT_CLOUD_FREE'
    return v,cloud,'OK'


def identify(base,lon,lat,year,month=None, *,timeout=20,session=None):
    session=session or requests
    q={'f':'json','geometry':json.dumps({'x':lon,'y':lat,'spatialReference':{'wkid':4326}}),
       'geometryType':'esriGeometryPoint','returnGeometry':'false','returnCatalogItems':'true'}
    if month is not None:
        q['time']=month_range_ms(year,month)
        q['renderingRule']=json.dumps({'rasterFunction':'Average Monthly Radiance (Raw Values)'})
    else:q['time']=month_range_ms(year,1).split(',')[0]+','+month_range_ms(year,12).split(',')[1]
    for attempt in range(3):
        try:
            req=session.get(base+'/identify',params=q,timeout=timeout)
            if req.status_code in (429,500,502,503,504):
                if attempt<2:time.sleep(.4*(attempt+1));continue
                return {'source_error':'http_'+str(req.status_code)}
            req.raise_for_status()
            return req.json()
        except (requests.Timeout,requests.ConnectionError,requests.HTTPError,ValueError) as exc:
            if attempt<2:time.sleep(.4*(attempt+1));continue
            return {'source_error':str(exc)[:150]}
    return {'source_error':'MAX_RETRIES'}


def retrieve_cell(row):
    lon=float(row['longitude']);lat=float(row['latitude']);ret={'cell_id':str(row['cell_id'])}
    pop=identify(POP,lon,lat,YEAR)
    if 'source_error' in pop or pop.get('error'):
        ret['population_density']=np.nan;ret['population_status']='API_ERROR'
    else:
        try:
            p=float(pop.get('value'))
            ret['population_density']=p if np.isfinite(p) and p>=0 else np.nan
            ret['population_status']='OK' if np.isfinite(ret['population_density']) else 'NO_DATA'
        except (ValueError,TypeError):
            ret['population_density']=np.nan;ret['population_status']='NO_DATA'
    radiances=[];cvs=[]
    for m in MONTHS:
        obj=identify(VIIRS,lon,lat,YEAR,month=m)
        if 'source_error' in obj:
            v,cv,status=None,None,'API_ERROR'
        else:v,cv,status=_parse_viirs_response(obj,YEAR,m)
        ret[f'radiance_{m:02d}']=np.nan if v is None else v
        ret[f'cloud_free_{m:02d}']=np.nan if cv is None else cv
        ret[f'viirs_status_{m:02d}']=status
        if status=='OK':radiances.append(v);cvs.append(cv)
    ret['viirs_valid_quarter_months']=len(radiances)
    ret['viirs_quarterly_median']=float(np.median(radiances)) if len(radiances)>=3 else np.nan
    ret['viirs_quarterly_mean']=float(np.mean(radiances)) if len(radiances)>=3 else np.nan
    ret['viirs_summed_cloudfree_quarters']=float(np.sum(cvs)) if len(radiances)>=3 else np.nan
    # Radiances < 0 due noise are preserved in raw and clipped only for log transform.
    ret['source_complete']=bool(len(radiances)>=3 and ret['population_status']=='OK')
    return ret


def fetch_exposure(cells,workers=10):
    check(1<=workers<=16,'Source server concurrency must be bounded')
    jobs=cells[['cell_id','longitude','latitude']].to_dict('records')
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        records=list(pool.map(retrieve_cell,jobs))
    q=pd.DataFrame(records)
    check(q.cell_id.is_unique and len(q)==len(cells),'Source results did not match source cells')
    return cells.merge(q,on='cell_id',validate='one_to_one')


def match_cells(cells,alancaliper=.65,popcaliper=.6,max_dist_deg=3.0):
    """Deterministic unique 1:1 local matches, outcomes are not used for matching."""
    a=cells.copy()
    valid=(a.source_complete & a.viirs_quarterly_median.notna() &
           a.population_density.notna() & a.white_taxon_equal.notna())
    a=a.loc[valid].copy()
    if not len(a):return pd.DataFrame(),{'eligible_cells':0,'matched_pairs':0}
    a['log_v']=np.log1p(np.maximum(0,a.viirs_quarterly_median.to_numpy(float)))
    a['log_pop']=np.log1p(a.population_density.to_numpy(float))
    a['log_rain']=np.log1p(np.maximum(0,a.bio12.to_numpy(float)))
    pairs=[];seen=set();candidates=0
    for group,df in a.groupby('block_id',sort=True):
        if len(df)<2:continue
        records=df.sort_values('cell_id').to_dict('records')
        for i,x in enumerate(records):
            for y in records[i+1:]:
                dv=abs(x['log_v']-y['log_v'])
                if dv<alancaliper:continue
                if abs(x['log_pop']-y['log_pop'])>.6:continue
                if abs(x['bio5']-y['bio5'])>3.0:continue
                if abs(x['log_rain']-y['log_rain'])>.5:continue
                if abs(x['elev']-y['elev'])>300:continue
                distance=math.hypot((x['latitude']-y['latitude']),
                                    (x['longitude']-y['longitude'])*math.cos(math.radians((x['latitude']+y['latitude'])/2)))
                if distance>max_dist_deg:continue
                hi,lo=(x,y) if x['log_v']>y['log_v'] else (y,x)
                cost=abs(x['log_pop']-y['log_pop'])/.6+abs(x['bio5']-y['bio5'])/3+abs(x['log_rain']-y['log_rain'])/.5+abs(x['elev']-y['elev'])/300+distance/max_dist_deg
                candidates+=1
                pairs.append((cost,hi['cell_id'],lo['cell_id'],hi,lo,distance))
    pairs.sort(key=lambda z:(z[0],z[1],z[2]))
    rows=[]
    for cost,hid,lid,hi,lo,dist in pairs:
        if hid in seen or lid in seen:continue
        seen.add(hid);seen.add(lid)
        rows.append({'bright_cell_id':hid,'dark_cell_id':lid,'match_block':hi['block_id'],
            'match_cost':cost,'distance_degrees_lat_adjusted':dist,
            'bright_viirs_nW_cm2_sr':hi['viirs_quarterly_median'],
            'dark_viirs_nW_cm2_sr':lo['viirs_quarterly_median'],
            'bright_population_per_km2':hi['population_density'],
            'dark_population_per_km2':lo['population_density'],
            'bright_white_taxon_equal':hi['white_taxon_equal'],
            'dark_white_taxon_equal':lo['white_taxon_equal'],
            'bright_strict_yield':hi['strict9_yield'],
            'dark_strict_yield':lo['strict9_yield'],
            'white_diff_bright_minus_dark_pp':100*(hi['white_taxon_equal']-lo['white_taxon_equal']),
            'strict_yield_diff_pp':100*(hi['strict9_yield']-lo['strict9_yield']),
            'bright_measured_taxa':int(hi['n_measured_taxa']),
            'dark_measured_taxa':int(lo['n_measured_taxa'])})
    return pd.DataFrame(rows), {'eligible_cells':len(a),'candidate_pairs_before_nonreuse':candidates,'matched_pairs':len(rows),
        'matcher':'same 5° block, log1p density <=0.6, BIO5 <=3°C, log rainfall<=.5, elevation<=300m, max center separation 3°, log1p(VIIRS)diff>=.65',
        'uses_outcome_for_matching':False}


def outcome_report(cells,pairs,source_meta):
    valid=cells[cells.source_complete]
    report={'schema':'mois_2020_quarterly_VIIRS_pop_density_original_strict9_cell_pilot_v1',
        'source':source_meta,'source_exposure_endpoints':{'viirs_monthly':VIIRS,'worldpop_density':POP},
        'quarter_months':list(MONTHS),
        'radiance_units':'nW/cm2/sr',
        'population_units':'people/km2',
        'population_year':2020,
        'four_months_NOT_true_annual_VIIRS_composite':True,
        'spatial_point_is_five_km_cell_center_NOT_instantaneous_flower_exposure':True,
        'existing_original_photo_images_downloaded':0,
        'new_flower_ROI_inference':0,
        'HuiWalter_or_MNAR_sensitivity_not_run':True,
        'photographed_taxon_colour_composition_NOT_species_absence':True,
        'independent_moth_pollination_trait_not_established':True,
        'not_causal_ALAN_effect':True,
        'cells_exposed':int(len(cells)),
        'VIIRS_cloud_pop_complete':int(len(valid)),
        'photographs_in_exposure_complete_cells':int(valid.n_attempted.sum()),
        'strict9_success_in_exposure_complete_cells':int(valid.n_strict9.sum()),
        'pairs':int(len(pairs)),
        'matched_cell_white_difference_mean_pp':float(pairs.white_diff_bright_minus_dark_pp.mean()) if len(pairs)>=20 else None,
        'matched_cell_white_difference_median_pp':float(pairs.white_diff_bright_minus_dark_pp.median()) if len(pairs)>=20 else None,
        'matched_roi_yield_difference_mean_pp':float(pairs.strict_yield_diff_pp.mean()) if len(pairs)>=20 else None,
        'white_assoc_reportable_if_at_least_20_pairs': bool(len(pairs)>=20),
    }
    if len(pairs):
        report['n_bright_lower_white_mass']=int((pairs.white_diff_bright_minus_dark_pp<0).sum())
        report['n_bright_higher_white_mass']=int((pairs.white_diff_bright_minus_dark_pp>0).sum())
        report['median_pop_density_relative_gap']=float((abs(pairs.bright_population_per_km2-pairs.dark_population_per_km2)/
                              (1+pairs.bright_population_per_km2+pairs.dark_population_per_km2)).median())
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--canonical',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--workers',type=int,default=10)
    args=p.parse_args()
    check(not args.out.exists(),'Do not overwrite signed prior result')
    check(sha(args.canonical)==CANONICAL_SHA,'Original canonical observation input hash changed')
    data=pd.read_csv(args.canonical,compression='gzip',low_memory=False)
    x,meta=select_cells(data)
    args.out.mkdir(parents=True)
    x.to_csv(args.out/'candidate_source_cells_2020.csv',index=False)
    print('MOIS_ALAN_REAL_SELECTED_CELLS',json.dumps(meta,sort_keys=True),flush=True)
    data=None
    enriched=fetch_exposure(x,workers=args.workers)
    enriched.to_csv(args.out/'viirs_population_source_matched_cells_2020.csv',index=False)
    pairs,match=match_cells(enriched)
    pairs.to_csv(args.out/'population_matched_bright_dark_pairs_2020.csv',index=False)
    report=outcome_report(enriched,pairs,meta)
    report['matching']=match
    report['input_sha256']={'canonical':CANONICAL_SHA}
    report['output_sha256']={p.name:sha(p) for p in args.out.iterdir()}
    (args.out/'receipt.json').write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+'\n')
    print('MOIS_ALAN_VIIRS_POP2020_REAL_SOURCE_MATCH',json.dumps({k:report[k] for k in ('cells_exposed','VIIRS_cloud_pop_complete','pairs','matched_cell_white_difference_mean_pp','matched_roi_yield_difference_mean_pp')},sort_keys=True),flush=True)

if __name__=='__main__':main()