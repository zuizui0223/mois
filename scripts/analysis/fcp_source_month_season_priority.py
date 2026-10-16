#!/usr/bin/env python3
"""Season-balanced opportunity ranking from the signed original iNaturalist source calendar.

Never interprets source photo opportunity as blooming, never selects based on
phenotype, and never changes the 1,000-taxon RGFCA random-screen cohort.
"""
from __future__ import annotations
import argparse, csv, hashlib, json, math
from collections import defaultdict
from pathlib import Path

SOURCE_DATE='2026-09-27'
SOURCE_TAXA=1000
SOURCE_PHOTOS=3430502
SOURCE_STATUS='COMPLETE_ORIGINAL_SOURCE_CALENDAR_CENSUS_NO_FLOWER_PIXELS'


def checksum(p:Path)->str:
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()


def local_season(month:int,latitude_cell:int)->str:
    lat=-90+(latitude_cell+.5)/4
    if abs(lat)<23.5:return 'tropical_no_four_season_label'
    season=('winter' if month in (12,1,2) else 'spring' if month in (3,4,5)
            else 'summer' if month in (6,7,8) else 'autumn')
    if lat<0:season={'winter':'summer','spring':'autumn','summer':'winter','autumn':'spring'}[season]
    return season


def work(site:Path,receipt:Path,out:Path,*,min_obs=24,latest_climate_year=2025,
         test_only=False)->dict:
    if out.exists() and any(out.iterdir()):raise ValueError('Do not overwrite signed pre-colour season screen')
    doc=json.loads(receipt.read_text(encoding='utf-8'))
    if not (doc.get('status')==SOURCE_STATUS and
            doc.get('original_source_site_month_file_sha256')==checksum(site) and
            doc.get('source_photo_opportunity_is_not_flowering_detection') is True and
            doc.get('new_flower_colour_or_white_frequency_estimated') is False and
            doc.get('n_taxa')==(2 if test_only else SOURCE_TAXA) and
            doc.get('n_original_photo_ids')==(16 if test_only else SOURCE_PHOTOS)):
        raise ValueError('Input must be SHA-signed, all-original-date photo opportunity census')
    if min_obs<2 or latest_climate_year>2025 or latest_climate_year<1958:
        raise ValueError('Unsupported observation depth or climate availability')
    summary=defaultdict(lambda:[0,0,0])
    names={}
    with site.open(newline='',encoding='utf-8') as f:
        rows=csv.DictReader(f)
        mandatory={'inat_taxon_id','species','year_month','lat_0p25_cell','lon_0p25_cell','all_source_photo_ids_with_5km_accuracy','distinct_source_observations','distinct_source_observers'}
        if not mandatory.issubset(rows.fieldnames or []):raise ValueError('Source-only site-month schema changed')
        for x in rows:
            tax=x['inat_taxon_id'];species=x['species']
            if tax in names and names[tax]!=species:raise ValueError('Taxon species identity drift')
            names[tax]=species
            date=x['year_month'];year=int(date[:4]);month=int(date[5:7])
            if len(date)!=7 or date[4]!='-' or not 1<=month<=12 or not 1958<=year<=2026:
                continue
            ilat=int(x['lat_0p25_cell']);ilon=int(x['lon_0p25_cell'])
            if not 0<=ilat<=719 or not 0<=ilon<=1439:raise ValueError('Source original 0.25-degree geography invalid')
            season=local_season(month,ilat)
            if season.startswith('tropical') or year>latest_climate_year:continue
            photos=int(x['all_source_photo_ids_with_5km_accuracy'])
            obs=int(x['distinct_source_observations']);reporters=int(x['distinct_source_observers'])
            if photos<obs or obs<reporters or min(photos,obs,reporters)<0:raise ValueError('Inconsistent original independent observations')
            if obs<2 or reporters<2:continue
            sy=year+int(month==12)
            if sy>latest_climate_year:continue
            agg=summary[(tax,species,ilat,ilon,sy,season)]
            agg[0]+=photos;agg[1]+=obs;agg[2]+=1
    blocks=defaultdict(dict)
    for (tax,name,ilat,ilon,year,season),counts in summary.items():
        if counts[1]>=min_obs and counts[0]>=min_obs:
            blocks[(tax,name,ilat,ilon,year)][season]=counts
    pairs=[]
    for (tax,name,ilat,ilon,year),grp in blocks.items():
        if len(grp)<2:continue
        ordered=sorted(grp.items(),key=lambda x:(-x[1][1],-x[1][0],x[0]))
        best=[]
        for i,(a,A) in enumerate(ordered):
            for b,B in ordered[i+1:]:
                k=(min(A[1],B[1]),min(A[0],B[0]),A[1]+B[1])
                best.append((k,a,b,A,B))
        (depth,a,b,A,B)=max(best,key=lambda x:x[0])
        pairs.append({'inat_taxon_id':tax,'species':name,'site_cell_lat_0p25':ilat,
          'site_cell_lon_0p25':ilon,'season_year':year,'local_season_a':a,
          'local_season_b':b,'source_photos_a':A[0],'source_photos_b':B[0],
          'independent_source_observations_a':A[1],
          'independent_source_observations_b':B[1],
          'min_obs_per_season':depth[0],'min_photos_per_season':depth[1],
          'each_season_has_at_least_two_observers_per_original_month':True,
          'independent_observers_ACROSS_seasons_verified':False,
          'is_flowering_photo_or_colour':False})
    pairs.sort(key=lambda r:(-r['min_obs_per_season'],-r['min_photos_per_season'],int(r['inat_taxon_id']),r['site_cell_lat_0p25'],r['site_cell_lon_0p25'],r['season_year']))
    out.mkdir(parents=True,exist_ok=True)
    file=out/'season_depth_24_source_only_geocell_year_priority.csv'
    fields=list(pairs[0]) if pairs else ['inat_taxon_id','species','site_cell_lat_0p25','site_cell_lon_0p25','season_year','local_season_a','local_season_b','source_photos_a','source_photos_b','independent_source_observations_a','independent_source_observations_b','min_obs_per_season','min_photos_per_season','each_season_has_at_least_two_observers_per_original_month','independent_observers_ACROSS_seasons_verified','is_flowering_photo_or_colour']
    with file.open('x',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=fields,lineterminator='\n');w.writeheader();w.writerows(pairs)
    report={
       'schema':'fcp_original_outcome_blind_temporal_geocell_priority_v1',
       'status':'SYNTHETIC_NOT_REAL' if test_only else 'REAL_FULL_ORIGINAL_SOURCE_TEMPORAL_OPPORTUNITY_NOT_FLOWER_COLOUR',
       'min_independent_original_observations_per_season':min_obs,
       'last_verified_climate_year':latest_climate_year,
       'n_eligible_geocell_years':len(pairs),
       'n_eligible_taxa':len({p['inat_taxon_id'] for p in pairs}),
       'no_2026_missing_climate_replaced_by_2025_values':True,
       'no_flower_or_colour_measurements_reused':True,
       'observer_disjoint_between_seasons_not_yet_audited':True,
       'input_site_source_sha256':checksum(site),
       'input_complete_calendar_receipt_sha256':checksum(receipt),
       'output_priority_sha256':checksum(file),
       'full_42111_species_93m_original_photo_census_not_implied':True,
       'candidate_photos_NOT_proven_flowering_or_white':True}
    (out/'season_depth_24_source_only_receipt.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    return report


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--original-site-months',type=Path,required=True)
    p.add_argument('--source-calendar-receipt',type=Path,required=True)
    p.add_argument('--outdir',type=Path,required=True)
    p.add_argument('--min-observations-per-season',type=int,default=24)
    a=p.parse_args()
    print(json.dumps(work(a.original_site_months,a.source_calendar_receipt,a.outdir,
          min_obs=a.min_observations_per_season),sort_keys=True))
if __name__=='__main__':main()
