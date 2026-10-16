#!/usr/bin/env python3
"""Source-gated WorldClim monthly-climatology and seasonal photographed-white support.

Monthly precipitation normals, not observed weather. WET/DRY refer to the
original site's own mean monthly rainfall cycle; small-variation localities
remain ASEASONAL. The source is a sealed Mois ROI sample (not 93M pixels).
No 13-variable model is fitted or biological selection inferred.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
import os
from collections import defaultdict, Counter
from pathlib import Path

NORMALS = "WorldClim 2.1 1970-2000 monthly precipitation normals 10m"
MANIFEST_SCHEMA = "mois_worldclim21_10m_monthly_precipitation_12_geotiffs_v1"
RESULT_SCHEMA = "mois_original_site_monthly_rainfall_white_support_v1"
RESULT_STATUS = "SOURCE_JOINED_MONTHLY_NORMALS_NO_HISTORICAL_WEATHER_OR_13ENV_REGRESSION"
LABELS = ("WET", "DRY", "TRANSITION", "ASEASONAL", "UNKNOWN_MONTH", "CLIMATE_UNAVAILABLE")
SOURCE_SEASON_SCHEMA = "mois_private_photo_seasonal_trait_support_v1"
SOURCE_PUBLIC_SCHEMA = "mois_photo_season_opportunity_v1"
PHENO_HUE = {"white", "yellow_orange", "red_pink", "blue_purple"}
MIN_CLASSIFIED = 3
MIN_OBSERVERS = 2
BAND_SIZE_MAX = 12_000_000
FIELDS = [
    "measurement_id", "inat_taxon_id", "genus", "observation_id", "observer_id",
    "latitude", "longitude", "bio12_mm", "month", "domain", "local_season",
    "annotation_status", "region_stratum", "precip_stratum", "sample_observation_weight",
    "sample_photo_weight", "morph", "measurement_status"
]
MONTHLY_EXTRA = ["monthly_precipitation_mm", "site_mean_monthly_mm",
                 "monthly_minus_site_mean_mm", "site_monthly_cv",
                 "site_precip_profile_annual_total_mm", "local_month_class",
                 "site_cell_row", "site_cell_col"]


def sha256_file(path: Path) -> str:
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(4*1024*1024), b''):
            h.update(chunk)
    return h.hexdigest()


def rainfall_label(monthly:list[float]|None, month:int|None)->dict:
    """Compare a month's normal with the SAME site's twelve-month cycle."""
    if monthly is None:
        return {"label":"CLIMATE_UNAVAILABLE", "cv":None,"mean":None,
                "monthly_mm":None,"anomaly_mm":None,"annual_mm":None}
    if (len(monthly)!=12 or any(not math.isfinite(z) or z<0 or z>10000 for z in monthly)):
        raise ValueError("Monthly precipitation profile is not twelve finite nonnegative mm values")
    total=sum(monthly)
    mean=total/12
    cv=(sum((z-mean)**2 for z in monthly)/12)**.5/mean if mean>0 else 0.
    if month is None:
        return {"label":"UNKNOWN_MONTH","cv":cv,"mean":mean,
                "monthly_mm":None,"anomaly_mm":None,"annual_mm":total}
    if not 1<=month<=12:
        raise ValueError("Original observation month must be 1..12")
    val=monthly[month-1]
    if cv<0.15 or max(monthly)-min(monthly)<30:
        label="ASEASONAL"
    else:
        # Interpolated quartiles. Ties are included together and never
        # arbitrarily split into wet/dry using the calendar month number.
        ordered=sorted(monthly)
        def quantile(q):
            idx=11*q;lo=math.floor(idx);hi=math.ceil(idx)
            return ordered[lo]+(idx-lo)*(ordered[hi]-ordered[lo])
        q25,q75=quantile(.25),quantile(.75)
        label=("WET" if val>=q75 and val>mean else
               "DRY" if val<=q25 and val<mean else "TRANSITION")
    return {"label":label,"cv":cv,"mean":mean,"monthly_mm":val,
            "anomaly_mm":val-mean,"annual_mm":total}


class MonthlyNormals:
    """Byte-check 12 co-registered WorldClim 2.1 precipitation GeoTIFFs; no web."""

    def __init__(self, manifest_path:Path, expected_hash:str, *, test_small=False):
        if len(expected_hash)!=64 or sha256_file(manifest_path)!=expected_hash:
            raise ValueError("Twelve-month rainfall manifest SHA256 mismatch")
        meta=json.loads(manifest_path.read_text(encoding='utf8'))
        if (meta.get('schema')!=MANIFEST_SCHEMA or meta.get('source')!=NORMALS
            or meta.get('units')!='mm/month' or meta.get('resolution')!='10m'
            or set(meta.get('months',{}))!={f'{m:02}' for m in range(1,13)}):
            raise ValueError("Not the complete pinned WorldClim 2.1 monthly rainfall normals")
        import rasterio
        import numpy as np
        band=[];reference=None
        self.raster_hashes={}
        for month in range(1,13):
            spec=meta['months'][f'{month:02}']
            if not isinstance(spec,dict) or set(spec)!={'path','sha256'}:
                raise ValueError("Missing unique raster path and exact SHA256")
            path=Path(spec['path'])
            if not path.is_absolute():
                path=manifest_path.parent/path
            actual=sha256_file(path)
            if len(spec['sha256'])!=64 or actual!=spec['sha256']:
                raise ValueError(f"WorldClim month {month:02} GeoTIFF SHA256 mismatch")
            with rasterio.open(path) as src:
                trans=src.transform
                crs=src.crs
                if (src.count!=1 or crs is None or crs.to_epsg()!=4326
                    or abs(trans.b)>1e-12 or abs(trans.d)>1e-12
                    or trans.a<=0 or trans.e>=0
                    or src.width*src.height>BAND_SIZE_MAX
                    or not test_small and (
                        abs(trans.a-1/6)>1e-5 or abs(trans.e+1/6)>1e-5
                        or src.width<2000 or src.height<800
                        or abs(src.bounds.left+180)>1 or abs(src.bounds.right-180)>1
                        or abs(src.bounds.top-90)>1 or not -91<=src.bounds.bottom<=-55)):
                    raise ValueError("Month raster not aligned global north-up EPSG4326 10m precip")
                key=(src.width,src.height,tuple(trans)[:6],crs.to_string())
                if reference is None:
                    reference=key
                elif key!=reference:
                    raise ValueError("Twelve monthly rainfall rasters are not co-registered")
                if src.scales[0]!=1. or src.offsets[0]!=0.:
                    raise ValueError("Unexpected monthly GeoTIFF unit scale, mm/month not guaranteed")
                val=src.read(1,masked=True)
                band.append(np.ma.filled(val.astype(float),np.nan))
            self.raster_hashes[f'{month:02}']=actual
        self.monthly=np.stack(band,axis=0)
        self.transform=rasterio.Affine(*reference[2])
        self.width,self.height=reference[:2]
        self.manifest_sha256=expected_hash
        self.cache={}
        self.test_small=test_small

    def cell(self,lat_raw:str,lon_raw:str):
        if not lat_raw or not lon_raw:
            return None
        lat,lon=float(lat_raw),float(lon_raw)
        if not math.isfinite(lat) or not math.isfinite(lon) or not -90<=lat<=90 or not -180<=lon<=180:
            raise ValueError("Invalid original source coordinate")
        cf,rf=(~self.transform)*(lon,lat)
        row,col=math.floor(rf),math.floor(cf)
        if not 0<=col<self.width or not 0<=row<self.height:
            return None
        return row,col

    def extract(self,lat_raw:str,lon_raw:str):
        cell=self.cell(lat_raw,lon_raw)
        if cell is None:
            return None,None
        if cell not in self.cache:
            row,col=cell
            vals=self.monthly[:,row,col]
            self.cache[cell]=None if not all(
                math.isfinite(float(x)) and 0<=float(x)<=10000 for x in vals
            ) else [float(x) for x in vals]
        return self.cache[cell],cell


def _summary(stat:dict)->dict:
    def ratio(n,d):return n/d if d else None
    return {"n_observations":stat['n'], "n_classifiable":stat['class'],
            "n_white":stat['white'],
            "observation_weighted_classifiability":ratio(stat['classw'],stat['allw']),
            "observation_weighted_white_given_classifiable":ratio(stat['whitew'],stat['classw']),
            "distinct_observers":len(stat['observers']),
            "distinct_taxa":len(stat['taxa'])}


def _new():
    return {"n":0,"class":0,"white":0,"allw":0.,"classw":0.,"whitew":0.,
            "observers":set(),"taxa":set()}


def audit(private_season_csv:Path,private_season_report:Path,public_season_report:Path,
          source_collection_receipt:Path,normals:MonthlyNormals,
          output_private_csv:Path,output_private_json:Path,output_public_json:Path,
          *,min_classifiable=MIN_CLASSIFIED,min_observers=MIN_OBSERVERS,
          min_public=10,synthetic=False)->dict:
    if min_classifiable<1 or min_observers<1 or min_public<1:
        raise ValueError("Matched seasonal support and privacy limits must be positive")
    priv=json.loads(private_season_report.read_text(encoding='utf8'))
    pub=json.loads(public_season_report.read_text(encoding='utf8'))
    if (priv.get('schema')!=SOURCE_SEASON_SCHEMA
        or priv.get('private_derived_photo_rows_sha256')!=sha256_file(private_season_csv)
        or priv.get('source_sample_receipt_sha256')!=sha256_file(source_collection_receipt)
        or pub.get('schema')!=SOURCE_PUBLIC_SCHEMA
        or pub.get('source_sealed_roi_join_sha256')!=priv.get('source_join_sha256')
        or pub.get('synthetic_test_only') is not synthetic
        or pub.get('old_genus_13_environment_estimates_updated') is not False
        or pub.get('all_13_photo_site_environment_predictors_available_here') is not False):
        raise ValueError("Season input not a complete source-attested sampled photo dataset")
    for p in (output_private_csv,output_private_json,output_public_json):
        if p.exists() or p.with_name(p.name+'.partial').exists():
            raise ValueError("Refuse overwrite or interrupted monthly climate/flower output")
        p.parent.mkdir(parents=True,exist_ok=True)
    stage=output_private_csv.with_name(output_private_csv.name+'.partial')
    by_label=defaultdict(_new)
    by_species_label=defaultdict(_new)
    by_site_label=defaultdict(_new)
    count=Counter();duplicates=set();annual_compared=0
    total=classified=missing_date=0
    try:
        with private_season_csv.open(newline='',encoding='utf8') as inp,stage.open('w',newline='',encoding='utf8') as out:
            reader=csv.DictReader(inp)
            if not set(FIELDS).issubset(reader.fieldnames or []):
                raise ValueError("Monthly source missing original selected photo observation and season fields")
            writer=csv.DictWriter(out,fieldnames=FIELDS+MONTHLY_EXTRA)
            writer.writeheader()
            for row in reader:
                ident=row['measurement_id']
                if not ident or ident in duplicates:
                    raise ValueError("Missing or duplicated original photo identity")
                duplicates.add(ident)
                if not row['inat_taxon_id'].isdigit() or not row['genus']:
                    raise ValueError("Missing original species or genus")
                month_str=row['month']
                if month_str=='UNKNOWN_MONTH':
                    month=None
                elif month_str.isdecimal() and 1<=int(month_str)<=12:
                    month=int(month_str)
                else:
                    raise ValueError("Observed month is not the frozen original seasonal audit label")
                lat,lon=row['latitude'],row['longitude']
                if bool(lat)!=bool(lon):
                    raise ValueError("Half-missing original coordinates")
                profile,cell=normals.extract(lat,lon)
                z=rainfall_label(profile,month)
                if profile is not None and row['bio12_mm'].strip():
                    bio12=float(row['bio12_mm'])
                    if not math.isfinite(bio12) or bio12<0 or abs(bio12-z['annual_mm'])>12:
                        raise ValueError("Monthly rain normals do not reconstruct original site BIO12 within 12 mm; raster version/resolution/unit mismatch")
                    annual_compared+=1
                if (row['measurement_status']=='classified_four_state_morph')!=(row['morph'] in PHENO_HUE):
                    raise ValueError("Original photo ROI status does not match white/four-colour outcome")
                ow=float(row['sample_observation_weight']);pw=float(row['sample_photo_weight'])
                if not (math.isfinite(ow) and ow>=1 and math.isfinite(pw) and pw>=1):
                    raise ValueError("Selection weights missing from original sample")
                label=z['label']
                if label not in LABELS:
                    raise ValueError("Unrecognized monthly rainfall group")
                fields={**{key:row[key] for key in FIELDS},
                        'monthly_precipitation_mm':z['monthly_mm'],
                        'site_mean_monthly_mm':z['mean'],
                        'monthly_minus_site_mean_mm':z['anomaly_mm'],
                        'site_monthly_cv':z['cv'],
                        'site_precip_profile_annual_total_mm':z['annual_mm'],
                        'local_month_class':label,
                        'site_cell_row':cell[0] if cell else '',
                        'site_cell_col':cell[1] if cell else ''}
                writer.writerow(fields)
                for key,store in (
                    ((row['domain'],label),by_label),
                    ((row['inat_taxon_id'],row['domain'],label),by_species_label),
                    ((row['inat_taxon_id'],row['domain'],cell,label),by_site_label),
                ):
                    d=store[key]
                    good=row['measurement_status']=='classified_four_state_morph'
                    white=good and row['morph']=='white'
                    d['n']+=1;d['class']+=int(good);d['white']+=int(white)
                    d['allw']+=ow;d['classw']+=ow*good;d['whitew']+=ow*white
                    if row['observer_id'].strip():
                        d['observers'].add(row['observer_id'].strip())
                    d['taxa'].add(row['inat_taxon_id'])
                count[label]+=1
                total+=1;classified+=int(row['measurement_status']=='classified_four_state_morph')
                missing_date+=int(month is None)
        if total!=pub['selected_source_observations'] or classified!=pub['source_classifiable_photographs']:
            raise ValueError("Monthly rainfall source differs from original sealed whole-sample seasonal totals")
        if not annual_compared and not synthetic:
            raise ValueError("No original annual BIO12 site raster extraction calibration was possible")
        matched=[]
        for tag,store in [('same_species',by_species_label),('same_species_same_worldclim_cell',by_site_label)]:
            grouped=defaultdict(dict)
            for key,s in store.items():
                if (key[-1] in ('DRY','WET') and s['class']>=min_classifiable
                    and len(s['observers'])>=min_observers and s['classw']>0):
                    grouped[key[:-1]][key[-1]]=s
            group_summaries=defaultdict(list)
            for identity,phases in grouped.items():
                if 'DRY' in phases and 'WET' in phases:
                    dry=phases['DRY'];wet=phases['WET']
                    if tag=='same_species_same_worldclim_cell' and identity[-1] is None:
                        continue
                    change=wet['whitew']/wet['classw']-dry['whitew']/dry['classw']
                    overlap=len(dry['observers']&wet['observers'])
                    group_summaries[identity[1]].append((change,overlap))
            for domain,vals in group_summaries.items():
                matched.append({
                    'comparison':tag,'domain':domain,'n_matched_units':len(vals),
                    'equal_unit_mean_white_difference_wet_minus_dry':
                        sum(v[0] for v in vals)/len(vals),
                    'n_units_with_observer_overlap':sum(v[1]>0 for v in vals),
                    'n_positive':sum(v[0]>0 for v in vals),
                    'n_negative':sum(v[0]<0 for v in vals),
                    'n_zero':sum(v[0]==0 for v in vals),
                })
        public_matches=[]
        for item in sorted(matched,key=lambda x:(x['comparison'],x['domain'])):
            if item['n_matched_units']>=min_public:
                public_matches.append(item)
            else:
                public_matches.append({
                    'comparison':item['comparison'],'domain':item['domain'],
                    'privacy_suppressed_matching_units_below':min_public
                })
        cells=[]
        for (domain,label),s in sorted(by_label.items()):
            if s['n']>=min_public:
                cells.append({'domain':domain,'local_rainfall_label':label,**_summary(s)})
            else:
                cells.append({'domain':domain,'local_rainfall_label':label,
                              'suppressed_below':min_public})
        result={
            'schema':RESULT_SCHEMA, 'status':RESULT_STATUS,
            'source_season_table_sha256':sha256_file(private_season_csv),
            'source_sample_collection_receipt_sha256':sha256_file(source_collection_receipt),
            'precip_monthly_normals_manifest_sha256':normals.manifest_sha256,
            'normals_climatology_period':'1970-2000',
            'normals_are_not_observed_weather_at_photo_date':True,
            'n_source_selected_original_observations':total,
            'n_classifiable_source_photos':classified,
            'n_missing_source_observation_month':missing_date,
            'n_sites_with_comparable_annual_bio12':annual_compared,
            'source_bio12_annual_reconstructed_within_mm':12,
            'local_label_totals':dict(sorted(count.items())),
            'domain_monthly_rainfall_conditioned_photo_white':cells,
            'matched_same_species_and_site_wet_dry':public_matches,
            'minimum_classified_source_observations_per_wet_dry_phase':min_classifiable,
            'minimum_distinct_observers_per_phase':min_observers,
            'minimum_public_matched_units':min_public,
            'seasonality_aseasonal_if_cv_below':0.15,
            'seasonality_aseasonal_if_monthly_max_min_span_below_mm':30,
            'same_species_colour_change_does_not_identify_individual_plant_change':True,
            'source_observation_annotation_bias_not_resolved':True,
            'global_13_environment_regression_not_run':True,
            'synthetic_test_only':bool(synthetic),
            'note':'Wet/dry month labels use local 1970-2000 monthly climatology only, with unknown/aseasonal distinct; observational colour outcomes are conditional on ROI classifiability.'
        }
        priv={
            'schema':RESULT_SCHEMA+'_private',
            'public_source_season_table_sha256':sha256_file(private_season_csv),
            'joined_private_monthly_photo_table_sha256':sha256_file(stage),
            'unsuppressed_matched_units':matched,
            'private_n_matched_units_by_comparison_and_domain':{
                item['comparison']+'/'+item['domain']:item['n_matched_units'] for item in matched
            },
            'source_original_monthly_profile_file_hashes':normals.raster_hashes,
            'privacy_risk_contains_raster_cell_coordinates':True,
            'status':'RESTRICTED_ONLY'
        }
        temp_priv=output_private_json.with_name(output_private_json.name+'.partial')
        temp_pub=output_public_json.with_name(output_public_json.name+'.partial')
        temp_priv.write_text(json.dumps(priv,indent=2,sort_keys=True)+'\n',encoding='utf8')
        temp_pub.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n',encoding='utf8')
        os.replace(stage,output_private_csv)
        os.replace(temp_priv,output_private_json)
        os.replace(temp_pub,output_public_json)
        return result
    except Exception:
        for p in (stage,output_private_json.with_name(output_private_json.name+'.partial'),
                  output_public_json.with_name(output_public_json.name+'.partial')):
            if p.exists():p.unlink()
        raise


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--private-seasonal-photo-table',type=Path,required=True)
    p.add_argument('--private-seasonal-report',type=Path,required=True)
    p.add_argument('--source-seasonal-public-report',type=Path,required=True)
    p.add_argument('--sample-collection-receipt',type=Path,required=True)
    p.add_argument('--monthly-rainfall-manifest',type=Path,required=True)
    p.add_argument('--monthly-rainfall-manifest-sha256',required=True)
    p.add_argument('--out-private-monthly-photo-table',type=Path,required=True)
    p.add_argument('--out-private-monthly-support',type=Path,required=True)
    p.add_argument('--out-public-summary',type=Path,required=True)
    p.add_argument('--min-classifiable',type=int,default=MIN_CLASSIFIED)
    p.add_argument('--min-observers',type=int,default=MIN_OBSERVERS)
    p.add_argument('--min-public',type=int,default=10)
    args=p.parse_args()
    normals=MonthlyNormals(args.monthly_rainfall_manifest,args.monthly_rainfall_manifest_sha256)
    summary=audit(args.private_seasonal_photo_table,args.private_seasonal_report,
                  args.source_seasonal_public_report,args.sample_collection_receipt,
                  normals,args.out_private_monthly_photo_table,args.out_private_monthly_support,
                  args.out_public_summary,min_classifiable=args.min_classifiable,
                  min_observers=args.min_observers,min_public=args.min_public)
    print(json.dumps({k:summary[k] for k in (
        'status','n_source_selected_original_observations',
        'n_classifiable_source_photos','n_sites_with_comparable_annual_bio12')},sort_keys=True))


if __name__=='__main__':
    main()
