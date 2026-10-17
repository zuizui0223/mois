#!/usr/bin/env python3
"""MOIS 247032 classified-photo source: separate site-static and time-variable support.

This does not reanalyze FCP's within-species D/W/spatial diversity, assume that
RGB-labelled images are genetically maintained floral morphs, or turn WorldClim
BIO normals into observation-month weather.
"""
from __future__ import annotations
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd

HUES=("white","yellow_orange","red_pink","blue_purple")
STATIC=("latitude","longitude","wc_elevation_m","soil_pH","soil_SOC","soil_N",
        "soil_clay","soil_available_water_proxy")
CLIMATOLOGY=("wc_bio1","wc_bio5","wc_bio12","wc_bio15")
COLUMNS=("photo_id","inat_taxon_id","species","morph","old_mois_or_new_fcp",
         "observation_id","observed_on","time_observed_at",*STATIC,*CLIMATOLOGY)


def require(ok,message):
    if not ok: raise ValueError(message)


def hash_file(p):
    h=hashlib.sha256()
    with Path(p).open("rb") as f:
        for b in iter(lambda:f.read(1<<20),b""):h.update(b)
    return h.hexdigest()


def local_season(month,lat):
    """Month labels are by hemisphere, and tropical months stay ungrouped."""
    if not np.isfinite(lat) or not 1<=month<=12:
        return "unknown"
    if abs(lat)<23.5:
        return "tropical_month_only"
    if month in (12,1,2):
        season="winter"
    elif month in (3,4,5):
        season="spring"
    elif month in (6,7,8):
        season="summer"
    else:
        season="autumn"
    if lat<0:
        season={"winter":"summer","spring":"autumn",
                "summer":"winter","autumn":"spring"}[season]
    return season


def classify_support(d):
    if d.photo_id.isna().any() or d.photo_id.duplicated().any():
        raise ValueError("Duplicate or unknown original photo")
    known=d.morph.isin(HUES)
    d["is_classified"]=known
    d["white_binary"]=d.morph.eq("white").astype(np.int8)
    for col in STATIC+CLIMATOLOGY:
        d[col]=pd.to_numeric(d[col],errors="coerce")
    valid_geo=d.latitude.between(-90,90)&d.longitude.between(-180,180)
    year_date=pd.to_datetime(d.observed_on,errors="coerce")
    d["calendar_date"]=year_date
    d["month"]=year_date.dt.month
    d["year"]=year_date.dt.year
    d["local_season"]=[local_season(m,x) for m,x in zip(
        d.month.fillna(0).astype(int),d.latitude.fillna(np.nan))]
    d["hemisphere"]=np.select(
        [d.latitude.ge(23.5),d.latitude.le(-23.5),d.latitude.between(-23.5,23.5)],
        ["north_temperate","south_temperate","tropics"],default="unknown")
    d["geocell025_lat"]=np.floor((d.latitude+90)*4)
    d["geocell025_lon"]=np.floor((d.longitude+180)*4)
    valid_clock=d.time_observed_at.astype("string").notna()&(
        d.time_observed_at.astype("string").str.strip()!="")
    require(not valid_clock.any(),"Unexpected original time-of-day requires timezone/hour validation before any daypart labels")
    per_source={}
    for source,part in d.groupby("old_mois_or_new_fcp",dropna=False):
        n=len(part)
        per_source[str(source)]={
          "all_photo_ids":n,
          "hue_classified":int(part.is_classified.sum()),
          "georeferenced":int((valid_geo&d.index.isin(part.index)).sum()),
          "observation_calendar_day_available":int(part.calendar_date.notna().sum()),
          "observation_clock_time_validated":int(valid_clock.loc[part.index].sum()),
          "static_source_coverage_photo_ids":{v:int(part[v].notna().sum()) for v in STATIC},
          "long_term_climate_normal_source_coverage":{v:int(part[v].notna().sum()) for v in CLIMATOLOGY},
        }
    return d,per_source


def describe(d,source,hue,field,bins=4):
    tmp=d[(d.old_mois_or_new_fcp==source)&d.is_classified&
          d[field].notna()].copy()
    if len(tmp)<bins*15:
        return []
    try:
        tmp["bin"]=pd.qcut(tmp[field],q=bins,duplicates="drop")
    except ValueError:
        return []
    result=[]
    for s,sub in tmp.groupby("bin",observed=True):
        result.append({"original_photo_source":source,"variable":field,
                       "range":str(s),"photo_count":len(sub),
                       "unique_taxa":sub.inat_taxon_id.nunique(),
                       "white_photo_fraction":round(float(sub.white_binary.mean()),6),
                       "sampling_is_not_species_balanced":True,
                       "not_a_causal_or_independent_effect":True})
    return result


def source_season_summary(d):
    # Distinguish photograph acquisition opportunity and labelled observation.
    dated=d[d.calendar_date.notna()&d.latitude.notna()&d.longitude.notna()].copy()
    monthly=[]
    for (src,year,month,hem),q in dated.groupby(
        ["old_mois_or_new_fcp","year","month","hemisphere"],observed=True):
        if not len(q):continue
        n=int(q.is_classified.sum())
        monthly.append({"original_source":src,
                        "year":int(year),"month":int(month),
                        "hemisphere":hem,"source_photo_ids":len(q),
                        "classified_hues":n,
                        "white_photo_count":int(q.white_binary.sum()),
                        "white_fraction_among_classified_photos":
                             round(int(q.white_binary.sum())/n,6) if n else None,
                        "unique_taxa":int(q.inat_taxon_id.nunique())})
    return dated,monthly


def site_year_support(d):
    f=d[d.is_classified & d.calendar_date.notna() &
        d.latitude.notna() & d.longitude.notna()].copy()
    f=f.loc[f.local_season.isin(["spring","summer","autumn","winter"])]
    f=f.loc[f.geocell025_lat.between(0,719)&f.geocell025_lon.between(0,1439)]
    f["latcell"]=f.geocell025_lat.astype(int)
    f["loncell"]=f.geocell025_lon.astype(int)
    by_year=f.groupby(["inat_taxon_id","year","local_season"],
                      observed=True).size().unstack(fill_value=0)
    by_site_year=f.groupby(["inat_taxon_id","latcell","loncell","year","local_season"],
                           observed=True).size().unstack(fill_value=0)
    by_site=f.groupby(["inat_taxon_id","latcell","loncell","local_season"],
                      observed=True).size().unstack(fill_value=0)
    thresholds={}
    for n in (1,2,4,8):
        thresholds[str(n)]={
            "species_year_with_2_local_seasons":int(((by_year>=n).sum(axis=1)>=2).sum()),
            "species_0p25cell_year_with_2_local_seasons":int(((by_site_year>=n).sum(axis=1)>=2).sum()),
            "species_0p25cell_any_year_with_2_local_seasons":int(((by_site>=n).sum(axis=1)>=2).sum()),
        }
    return thresholds


def run(photos_gz:Path,union_receipt:Path,outdir:Path):
    require(not outdir.exists() or not any(outdir.iterdir()),"No rewriting a frozen result tree")
    receipt=json.loads(union_receipt.read_text(encoding="utf-8"))
    require(hash_file(photos_gz)==receipt["photo_hue_union_csv_sha256"],
            "Real photographic union SHA256 drift")
    require(receipt["combined_unique_photo_ids"]==247032 and
            receipt["combined_classified_hue_photo_ids"]==120310,
            "Original measured-photo counts changed")
    d=pd.read_csv(photos_gz,compression="gzip",usecols=list(COLUMNS),
                  dtype={"photo_id":"string","inat_taxon_id":"string",
                         "observation_id":"string"},low_memory=False)
    require(len(d)==247032,"All combined source photo rows required")
    d,per_source=classify_support(d)
    dated,months=source_season_summary(d)
    static=[]
    for v in ("latitude","longitude","wc_elevation_m",
              "soil_pH","soil_SOC","soil_N","soil_clay","soil_available_water_proxy",
              "wc_bio1","wc_bio5","wc_bio12","wc_bio15"):
        static.extend(describe(d,"historical",HUES,v))
    outdir.mkdir(parents=True,exist_ok=True)
    files={}
    for name,records in (("original_year_month_hemi_photo_support.csv",months),
                         ("original_site_static_climatology_hue_descriptives.csv",static)):
        p=outdir/name
        pd.DataFrame(records).to_csv(p,index=False,lineterminator="\n")
        files[name]=hash_file(p)
    thresholds=site_year_support(d)
    coverage={
      "schema":"mois_247032_classified_photo_static_vs_real_observation_time_support_v1",
      "source_union_sha256":receipt["photo_hue_union_csv_sha256"],
      "original_source_photo_ids":int(len(d)),
      "original_classified_photo_ids":int(d.is_classified.sum()),
      "source_stratified_support":per_source,
      "static_site_variables_with_existing_source_values":list(STATIC),
      "longterm_climatic_normals_NOT_monthly_weather":list(CLIMATOLOGY),
      "actual_monthly_weather_variables_NOT_YET_JOINED":["tmin","tmax","ppt","srad","def","vpd","soil"],
      "calendar_observation_day_present":int(d.calendar_date.notna().sum()),
      "actual_UTC_or_local_observed_time_PRESENT_and_timezone_verified":False,
      "all_photo_daypart_labels_UNKNOWN":True,
      "at_least_two_local_seasons_by_sampling_support":thresholds,
      "monthly_descriptive_table_rows":len(months),
      "static_site_hue_quartile_table_rows":len(static),
      "output_sha256":files,
      "FCP_within_species_D_W_and_50km_local_depletion_not_analyzed":True,
      "no_missing_value_imputation":True,
      "biological_scope":"MOIS global photo-hue composition, geographic/edaphic baseline and original sampled calendar-month support, NOT adaptive genetic polymorphism or photographic flower-stage verification",
      "interpretation_limit":"Source composition and site/season turnover may confound photo hues; no month-specific real weather coefficient or causal within-site adaptation was estimated",
    }
    (outdir/"mois_247k_static_vs_time_support_receipt.json").write_text(
        json.dumps(coverage,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("MOIS_247K_PHOTO_STATIC_AND_SEASON_SUPPORT",json.dumps({
       "photo_ids":len(d),"classifiable":int(d.is_classified.sum()),
       "calendar_day":coverage["calendar_observation_day_present"],
       "matched":thresholds,
    },sort_keys=True))
    return coverage


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--photo-hue-union",type=Path,required=True)
    p.add_argument("--photo-hue-union-receipt",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    run(a.photo_hue_union,a.photo_hue_union_receipt,a.output)

if __name__=="__main__":main()
