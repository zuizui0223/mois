#!/usr/bin/env python3
"""MOIS source-photo community season contrast and geographic-scale sensitivity.

Community means all photographed taxa within one original coarse lat/lon cell
and meteorological season-year, NOT individual-plant or same-taxon plasticity.
Real original photo ROI labels are not botanically confirmed open petals.
The 267425 genuine measured pixel data and the 9597 actual verified original
historical public dates are SHA/ID authenticated before any statistical test.
"""
from __future__ import annotations
import argparse
import hashlib
import io
import json
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

from scripts.analysis.mois_seasonal_267k_climate_normals import (
    HUES, MOIS_ROWS, MOIS_HUES, PHOTO_DATA, PHOTO_RECEIPT, read_signed, require,
)
from scripts.analysis.mois_historical9600_taxonomic_breadth_vs_local_season_support import (
    make_panel, DATE_INPUT, DATE_RECEIPT, PHOTO_SHA, DATE_SHA,
)

SCALES=(0.25,0.5,1.0,2.0,5.0)
THRESHOLDS=(3,5,10)
ZONES=("north_extratropical","south_extratropical")
SEASONS=("summer","winter")
N_BOOT=999


def matched_community(p: pd.DataFrame, *, width:float, min_obs:int,
                      min_taxa:int=2) -> pd.DataFrame:
    """One row for each source cell-season-year with BOTH summer and winter.

    At least min_obs distinct original iNaturalist observation IDs and
    min_taxa original named taxa are required in EACH season independently.
    Source single-species and photographic effort are not field abundances.
    """
    require(width>0 and min_obs>=1 and min_taxa>=1,"Invalid locality or sample minimum")
    d=p.loc[p.classified&p.date.notna()&
            p.zone.isin(ZONES)&p.season.isin(SEASONS)&
            p.latitude.between(-90,90)&p.longitude.between(-180,180)].copy()
    d["latitude_cell"]=np.floor((d.latitude+90.0)/width).astype(int)
    d["longitude_cell"]=np.floor((d.longitude+180.0)/width).astype(int)
    keys=["zone","latitude_cell","longitude_cell","season_year"]
    tab=d.groupby(keys+["season"],observed=True).agg(
        n_photos=("photo_id","size"),
        n_observations=("observation_id","nunique"),
        n_taxa=("inat_taxon_id","nunique"),
        n_white=("white","sum")).unstack(fill_value=0)
    required=[("n_observations",s) for s in SEASONS]+[("n_taxa",s) for s in SEASONS]
    if tab.empty or any(c not in tab.columns for c in required):
        return pd.DataFrame(columns=keys+[
            "summer_photos","winter_photos","summer_observations","winter_observations",
            "summer_taxa","winter_taxa","summer_white","winter_white",
            "summer_minus_winter_white_pp"])
    mask=(tab[("n_observations","summer")]>=min_obs)&(
          tab[("n_observations","winter")]>=min_obs)&(
          tab[("n_taxa","summer")]>=min_taxa)&(
          tab[("n_taxa","winter")]>=min_taxa)
    ok=tab.loc[mask]
    out=pd.DataFrame(index=ok.index)
    for season in SEASONS:
        for feature,name in (
           ("n_photos","photos"),("n_observations","observations"),
           ("n_taxa","taxa"),("n_white","white")):
            out[season+"_"+name]=ok[(feature,season)].astype(int)
    out["summer_minus_winter_white_pp"]=100*(
        out.summer_white/out.summer_photos -
        out.winter_white/out.winter_photos)
    return out.reset_index()


def cell_cluster_interval(table:pd.DataFrame, *, draws:int=N_BOOT, seed:int=20261017):
    """Bootstrap 0.25–5 degree geographic CELLS (not site-years or photos).

    One coordinate cell can be visited in multiple years; give it one sampled
    cluster to avoid claiming its years are independent new locations.
    Descriptive percentile intervals are not full multiple-testing maxT CIs.
    """
    if table.empty:return None
    point=float(table.summer_minus_winter_white_pp.mean())
    tab=table.groupby(["latitude_cell","longitude_cell"],observed=True)[
        "summer_minus_winter_white_pp"].agg(["sum","count"])
    n_cells=len(tab)
    if n_cells<10:
        return {"equal_cell_year_white_difference_pp":round(point,6),
                "bootstrap_95pp":None,"n_unique_geocells":n_cells,
                "uncertainty_status":"INSUFFICIENT_INDEPENDENT_CELLS"}
    rng=np.random.default_rng(seed)
    sums=tab["sum"].to_numpy(float)
    ns=tab["count"].to_numpy(float)
    sampled=rng.integers(0,n_cells,size=(draws,n_cells))
    effects=sums[sampled].sum(axis=1)/ns[sampled].sum(axis=1)
    bounds=np.quantile(effects,[0.025,0.975])
    return {"equal_cell_year_white_difference_pp":round(point,6),
            "bootstrap_95pp":[round(float(x),6) for x in bounds],
            "n_unique_geocells":n_cells,
            "uncertainty_status":"EXPLORATORY_GEOCELL_CLUSTER_PERCENTILE_NOT_MAXT"}


def diagnostics(p:pd.DataFrame) -> tuple[pd.DataFrame,dict]:
    records=[]
    for label,subset in (("before_historical_date_recovery",
                          p[p.source_calendar_type.eq("original_only")]),
                         ("after_original9600_date_recovery",p)):
        for width in SCALES:
            for n in THRESHOLDS:
                frame=matched_community(subset,width=width,min_obs=n,min_taxa=2)
                for zone in ZONES:
                    sub=frame[frame.zone.eq(zone)]
                    summary=cell_cluster_interval(
                        sub,draws=N_BOOT,seed=20261017+int(100*width)+n)
                    records.append({
                        "photo_date_source":label,"geocell_width_degrees":width,
                        "minimum_distinct_observation_IDs_each_season":n,
                        "minimum_named_taxa_each_season":2,
                        "hemisphere":zone,"coarse_cell_season_years":int(len(sub)),
                        "distinct_source_geocells":int(sub[
                            ["latitude_cell","longitude_cell"]].drop_duplicates().shape[0]),
                        "photographs_in_both_seasons":int(
                            sub.summer_photos.sum()+sub.winter_photos.sum()) if len(sub) else 0,
                        "n_summer_photos":int(sub.summer_photos.sum()) if len(sub) else 0,
                        "n_winter_photos":int(sub.winter_photos.sum()) if len(sub) else 0,
                        "equal_cell_year_white_summer_winter_pp":
                            None if summary is None else summary[
                                "equal_cell_year_white_difference_pp"],
                        "site_cluster_percentile95_low":
                            None if summary is None or summary["bootstrap_95pp"] is None
                            else summary["bootstrap_95pp"][0],
                        "site_cluster_percentile95_high":
                            None if summary is None or summary["bootstrap_95pp"] is None
                            else summary["bootstrap_95pp"][1],
                        "uncertainty_status":"NO_COMPARABLE_PHOTO_LOCALITIES" if summary is None
                            else summary["uncertainty_status"],
                    })
    out=pd.DataFrame(records)
    return out,{
        "sample_sizes_and_site_year_photo_differences":len(out),
        "bootstrap_clustered_by_coarse_geocell_not_original_photographs":True,
        "unadjusted_globally_pooled_northern_or_southern_seasonal_ratios_not_equivalent_to_paired_locality_year":True,
        "at_least_two_original_taxa_both_seasons_required":True,
        "min_original_distinct_observation_ID_not_independent_photographers":True,
    }


def run(photo_zip:Path,date_zip:Path,output:Path):
    require(not output.exists() or not any(output.iterdir()),
            "Do not overwrite historical photo/source outputs")
    original,receipt,source_sha=read_signed(photo_zip,PHOTO_DATA,PHOTO_RECEIPT,
        "source_photo_level_csv_sha256",
        ["photo_id","observation_id","inat_taxon_id","morph","latitude",
         "longitude","observed_on","old_mois_or_new_fcp"])
    require(source_sha==PHOTO_SHA and len(original)==MOIS_ROWS and
            original.morph.isin(HUES).sum()==MOIS_HUES,
            "Actual 267425 original source photo/hue bytes drifted")
    with zipfile.ZipFile(date_zip) as z:
        dated_source=z.read(DATE_INPUT)
        date_receipt=json.loads(z.read(DATE_RECEIPT))
    require(hashlib.sha256(dated_source).hexdigest()==DATE_SHA and
            date_receipt["new_original_photo_date_CSV_sha256"]==DATE_SHA and
            date_receipt["exact_original_photo_observation_and_current_date_matches"]==9597,
            "Historical original 9600 photo-ID source date provenance drifted")
    d=pd.read_csv(io.BytesIO(dated_source),dtype=str)
    p=make_panel(original,d,True)
    before=make_panel(original,d,False)
    require(int((p.classified&p.date.notna()).sum())==86172 and
            int((before.classified&before.date.notna()).sum())==81867,
            "Current photo/observed date hue support changed")
    p["source_calendar_type"]=np.where(
        p.photo_id.isin(d.original_photo_id),
        "original_plus_public_verified" ,"original_only")
    original_date_ids=set(original.loc[original.observed_on.notna(),"photo_id"])
    p["source_calendar_type"]=np.where(
        p.photo_id.isin(original_date_ids),"original_only",
        "original_plus_public_verified")
    rows,protocol=diagnostics(p)
    output.mkdir(parents=True,exist_ok=True)
    csv=output/"mois_real_86172_photo_local_season_community_white_site_year_scaling.csv"
    rows.to_csv(csv,index=False,lineterminator="\n",float_format="%.6f")
    key=lambda cohort,width,n,zone:rows[
       (rows.photo_date_source.eq(cohort))&
       (rows.geocell_width_degrees.eq(width))&
       (rows.minimum_distinct_observation_IDs_each_season.eq(n))&
       (rows.hemisphere.eq(zone))].iloc[0].to_dict()
    result={
        "schema":"mois_86172_authentic_photographic_hue_spatial_sampling_scale_for_seasonal_community_v1",
        "photo_original_source_sha256":source_sha,
        "date_9600_original_source_sha256":DATE_SHA,
        "source_photo_ids":MOIS_ROWS,"classified_original_photo_ids":MOIS_HUES,
        "initial_dated_classified_original_photos":81867,
        "after_historical_dates_dated_classified_original_photos":86172,
        "minimum_distinct_taxa_in_each_season":2,
        "fixed_local_meteorological_season_year":True,
        "spatial_scale_degrees":list(SCALES),
        "observer_effort_not_standardized":True,
        "same_0p25_degree_site_year_min3_observations":{
          zone:{
           "before":key("before_historical_date_recovery",0.25,3,zone),
           "after":key("after_original9600_date_recovery",0.25,3,zone)}
          for zone in ZONES},
        "same_0p25_degree_site_year_min5_observations":{
          zone:{
           "before":key("before_historical_date_recovery",0.25,5,zone),
           "after":key("after_original9600_date_recovery",0.25,5,zone)}
          for zone in ZONES},
        "after_min5_multiscale":{
          str(w):{zone:key("after_original9600_date_recovery",w,5,zone)
               for zone in ZONES} for w in SCALES},
        "sampling_uncertainty":protocol,
        "results_CSV_sha256":hashlib.sha256(csv.read_bytes()).hexdigest(),
        "interpretive_limit":[
         "A source photocommunity pools different taxa and observers; not same-taxon pigment plasticity.",
         "A 0.25 degree geocell is tens of kilometres wide and not a biological population.",
         "Photos are automated hue ROIs, no independent open-flower organ confirmation.",
         "No independent photographer effort or count of absent flowering was observed.",
         "1999-site/group photo confidence intervals are descriptive source cluster bootstrap, not simultaneous confidence intervals.",
         "No actual observation-year monthly temperature, water deficit or solar mechanism tested.",
         "FCP within species W/D local phenotype-space inferences and historical genus moisture maxT unchanged.",
        ],
    }
    target=output/"mois_real86172_photo_community_seasonal_geocell_scale_decision.json"
    target.write_text(json.dumps(result,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    print("MOIS_REAL_PHOTOGRAPHIC_COMMUNITY_LOCALITY_YEAR_SCALE",
          json.dumps({zone:{k:result["same_0p25_degree_site_year_min5_observations"][zone]["after"][k]
                            for k in ("coarse_cell_season_years","photographs_in_both_seasons")}
                     for zone in ZONES},sort_keys=True),flush=True)
    return result


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--original-photo-zip",required=True,type=Path)
    ap.add_argument("--original-date-9600-zip",required=True,type=Path)
    ap.add_argument("--output",required=True,type=Path)
    args=ap.parse_args()
    run(args.original_photo_zip,args.original_date_9600_zip,args.output)
