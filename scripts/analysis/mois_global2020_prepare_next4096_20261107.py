"""Reconstruct the exact MOIS 2020 all-land source, select 4096 next pending.

Must match the independently reconstructed local 4096-source manifest SHA.
No photo downloads; unknowns and all source arms are retained.
"""
from __future__ import annotations
import argparse,hashlib,json,time
from pathlib import Path
import pandas as pd,numpy as np,geopandas as gpd
from scripts.analysis.mois_global_2020_GSHHG_scope_20261107 import (
    prepare_land,assign_points,aggregate_targets)

SHAS={
 'canonical':'f033661d76e07898a5c13199a55c5d0ce6e55353daa48b0058c456979e90644c',
 'legacy':'3562db1e5783e77d3bb408bd33817b4d7cccf31ad4d4ad42586469708fd407e3',
 'previous256':'303424afe40b41e7bda2ec0ae456dd27820a50e17f7ca23204c81036adf21466',
 'previous1024':'7aa0193178d664818bca777749074fc88fe0b9e594e83e5eb34d6f72110efce9',
 'previous64':'a513250fa9911a35c7f40380684ae97f52d7cdfea7bf828cdcd7cbe1be34192c'}
MANIFEST='fa06170ce59e41b6d7e486c236acf7d004091a035f7cbddbb3632bcd788ffa65'

def check(ok,msg):
    if not bool(ok):raise ValueError("MOIS_2020_4096_SOURCE: "+msg)

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(4<<20),b''):h.update(b)
    return h.hexdigest()

def select(q,existing,legacy_cells,expected_pending=34726,selected_n=4096):
    """Target selection independent of strict9, flower colour, and geography."""
    no=set(existing)
    same=q.exposure_site_kind.eq('same_land_grid_centre')
    no.update(q.loc[same & q.cell_id.isin(legacy_cells),'target_id'])
    d=q.loc[q.land_class.isin(['island','mainland']) & ~q.target_id.isin(no)].copy()
    check(len(d)==expected_pending,'previous source population not reconciled')
    d['blind_hash']=d.target_id.map(
        lambda x:hashlib.sha256(('MOIS-global-v1|'+x).encode()).hexdigest())
    d=d.sort_values('blind_hash').head(selected_n)
    check(len(d)==selected_n and d.target_id.is_unique,"duplicated or insufficient target cohort")
    check(not d[['query_longitude','query_latitude']].duplicated().any(),
          'duplicate raster coordinate inside target sample')
    return d

def run(canonical,shoreline,legacy,previous256,previous1024,previous64,out,previous4096=None):
    files={'canonical':canonical,'legacy':legacy,'previous256':previous256,
           'previous1024':previous1024,'previous64':previous64}
    for name,p in files.items():check(sha(p)==SHAS[name],'source SHA changed '+name)
    shore=Path(shoreline)
    coast=json.loads((shore/'source_receipt.json').read_text())
    check(coast.get('original_archive_sha256')==
        '8dbbe7e071e77e9e75f2d639239099ebca8d5c16d6a07df8169729d49f15cf41',
        'coastline source version changed')
    for file,dig in coast['extracted_sha256'].items():
        check(sha(shore/file)==dig,'GSHHG source component changed '+file)
    out=Path(out);check(not out.exists(),'refuse to overwrite previous manifest')
    t0=time.monotonic()
    d=pd.read_csv(canonical,low_memory=False,dtype={'photo_id':str,'observation_id':str},
        usecols=['photo_id','observation_id','source_subarm','strict9',
                 'observed_on','latitude','longitude'])
    check(len(d)==1032882 and d.photo_id.is_unique and d.observation_id.is_unique,
          'original photo/observation units changed')
    d['year']=pd.to_numeric(d.observed_on.astype(str).str[:4],errors='coerce')
    d=d.loc[d.year.eq(2020)].copy()
    check(len(d)==74023,'original photo year cohort changed')
    d['unit_id']=d.photo_id
    land=prepare_land(gpd.read_file(shore/'GSHHS_h_L1.shp'))
    lakes=gpd.read_file(shore/'GSHHS_h_L2.shp')
    point=assign_points(d,land,lakes)
    d=d.merge(point,on='unit_id',validate='one_to_one')
    q=aggregate_targets(d)
    check(len(q)==38403,'original 2020 GSHHG cell-landmass count changed')
    ij=q.cell_id.str.split('_',expand=True).astype(int)
    q['longitude']=(ij[0]+.5)*.05-180
    q['latitude']=(ij[1]+.5)*.05-90
    centres=q[['cell_id','longitude','latitude']].drop_duplicates('cell_id').rename(columns={'cell_id':'unit_id'})
    ctr=assign_points(centres,land,lakes).rename(columns={
         'unit_id':'cell_id','land_id':'centre_land_id',
         'land_class':'centre_land_class','area_km2':'centre_area_km2'})
    q=q.merge(ctr,on='cell_id',validate='many_to_one')
    same=q.land_id.notna()&q.land_id.eq(q.centre_land_id)&q.land_class.isin(['island','mainland'])
    check(int(same.sum())==35275,'source grid-center same-land count changed')
    q['exposure_site_kind']=np.where(same,'same_land_grid_centre','original_photo_anchor')
    q['query_longitude']=np.where(same,q.longitude,q.anchor_longitude)
    q['query_latitude']=np.where(same,q.latitude,q.anchor_latitude)
    q['land_area_km2']=q.land_id.map(land.set_index('land_id').area_km2)
    legacy_df=pd.read_csv(legacy,dtype={'cell_id':str})
    check(len(legacy_df)==1059 and int(legacy_df.source_complete.sum())==1044,
          '2020 source VIIRS availability changed')
    old_cells=set(legacy_df.loc[legacy_df.source_complete,'cell_id'])
    check(int((same&q.cell_id.isin(old_cells)).sum())==1014,
          'old 2020 GSHHG target-availability relation changed')
    old=[]
    for name,p in [('previous256',previous256),('previous1024',previous1024),('previous64',previous64)]:
        r=pd.read_csv(p,dtype={'target_id':str},low_memory=False)
        check(len(r)==int(name.replace('previous','')) and r.target_id.is_unique,
              'existing prior source result lost '+name)
        old.extend(r.target_id)
    if previous4096 is not None:
        check(sha(previous4096)==MANIFEST,'first 4096 source selection manifest changed')
        extra=pd.read_csv(previous4096,dtype={'target_id':str})
        check(len(extra)==4096 and extra.target_id.is_unique,'previous new 4096 IDs changed')
        old.extend(extra.target_id)
    expected_old=5440 if previous4096 is not None else 1344
    check(len(old)==len(set(old))==expected_old and set(old)<=set(q.target_id),
          'prior source target ID collision')
    all_remaining=previous4096 is not None
    expected_n=30630 if all_remaining else 4096
    source_sha='f9805adb2533c6a26a632e6d011b2468637e7603bc7b1afb209b06ab38a2fd35' if all_remaining else MANIFEST
    sampled=select(q,old,old_cells,expected_pending=expected_n,selected_n=expected_n)
    fields=['target_id','land_class','query_longitude','query_latitude',
            'n_attempted','n_strict9','source_subarms']
    out.mkdir(parents=True)
    target=out/('source_locked_2020_remaining30630.csv' if all_remaining else 'source_locked_2020_new4096.csv')
    sampled[fields].to_csv(target,index=False,float_format='%.12g',lineterminator='\n')
    check(sha(target)==source_sha,'source manifest differs from independently computed exact global queue')
    report={'schema':'mois_source_locked_2020_after4096_remaining30630_noarea_v1' if all_remaining else 'mois_source_locked_2020_next4096_noarea_v1',
      'canonical_sha256':SHAS['canonical'],
      'coastline_archive_sha256':coast['original_archive_sha256'],
      'source_SHA256':SHAS,
      'original_2020_source_photo_count':len(d),
      'original_2020_source_grid_landmass_targets':len(q),
      'original_2020_same_landmass_grid_centres':35275,
      'previous_legacy_2020_source_sites':1014,
      'previous_new_requests_count':len(old),
      'pending_known_land_target_pool':expected_n,
      'new_selected_targets':len(sampled),
      'island':int(sampled.land_class.eq('island').sum()),
      'mainland':int(sampled.land_class.eq('mainland').sum()),
      'island_below5km2':int((sampled.land_class.eq('island')&sampled.land_area_km2.lt(5)).sum()),
      'outside_plusminus60':int(sampled.query_latitude.abs().gt(60).sum()),
      'photo_attempts_referenced':int(sampled.n_attempted.sum()),
      'previous_strict9_success_referenced':int(sampled.n_strict9.sum()),
      'manifest_sha256':sha(target),
      'selection_independent_of_flower_colour_ROI_and_photo_success':True,
      'island_area_cutoff':None,'latitude_cutoff':None,
      'new_photo_downloads_or_ROI_inference':0,
      'elapsed_to_reconstruct_seconds':round(time.monotonic()-t0,1)}
    (out/'source_receipt.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print('MOIS_NEXT4096_RECONSTRUCTED_AND_SOURCE_LOCKED',
          json.dumps(report,sort_keys=True),flush=True)
    return report

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('canonical','shoreline','legacy','previous256','previous1024','previous64','out'):
        p.add_argument('--'+n,required=True,type=Path)
    p.add_argument('--previous4096',type=Path)
    a=p.parse_args()
    run(a.canonical,a.shoreline,a.legacy,a.previous256,a.previous1024,a.previous64,a.out,a.previous4096)
