#!/usr/bin/env python3
"""MOIS genuine photo-date sampling: seasonally changing four-hue classifiability.

There are THREE distinct rates: classifiable photos / dated original photo
attempts (c), white labels / classified photos (p), and white labelled photos
/ all dated attempts (y = c*p). The symmetric identity
  delta(y) = mean(c)*delta(p) + mean(p)*delta(c)
shows source photographic eligibility, NOT a causal decomposition of flowering.

Original 267425 photo IDs, 129795 four-hue ROI classifications and the signed
historical S3 + independently numeric-observation/ID verified fallback dates
are source-locked. No unclassifiable photo is coded nonwhite or nonflowering.
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
    PHOTO_DATA, PHOTO_RECEIPT, HUES, read_signed, labelled_calendar, require,
)

ORIGINAL_PHOTOS = 267425
ORIGINAL_HUES = 129795
DATED_PHOTOS = 261578
DATED_HUES = 127504
HISTORICAL_PHOTOS = 100543
SOURCE_SHA = "fbb54b58dc91e40e2e166081212588acc18adb19af4bdbcd1c7fe1778bee28d9"
FINAL_HISTORIC_SHA = "ac5056e25cb9322618dda3b4f2e4e2a97907dfb16749d84c65ee9f3149143564"
HISTORICAL_NAME = "mois_original100543_s3_or_verified_API_dates_photo_source.csv.gz"
FINAL_SOURCE_PROOF = "mois_real127504_original_photo_S3_plus_verified_API_local_season_result.json"
WIDTHS = (0.25, 0.5, 1.0)
GEOCELL_BOOTSTRAPS = 999


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def date_and_original_hues(photo_zip: Path, final_zip: Path) -> pd.DataFrame:
    d, receipt, source_hash = read_signed(
        photo_zip, PHOTO_DATA, PHOTO_RECEIPT, "source_photo_level_csv_sha256",
        ["photo_id", "observation_id", "inat_taxon_id", "morph", "observed_on",
         "latitude", "longitude", "old_mois_or_new_fcp"])
    require(len(d)==ORIGINAL_PHOTOS and d.photo_id.is_unique and
            d.morph.isin(HUES).sum()==ORIGINAL_HUES and source_hash==SOURCE_SHA,
            "Original real photographed pixel/ROI source identity changed")
    with zipfile.ZipFile(final_zip) as z:
        raw = z.read(HISTORICAL_NAME)
        proof = json.loads(z.read(FINAL_SOURCE_PROOF))
    require(sha256(raw)==FINAL_HISTORIC_SHA and
            proof["full_historical_join"]["historical_original_photo_ids"]==HISTORICAL_PHOTOS and
            proof["full_historical_join"]["real_original_photo_hue_and_date_after_all_verified_rejoins"]==DATED_HUES and
            proof["historical_original_photo_classifications_frozen"]==45953,
            "Signed final historical date source or classifier source changed")
    h = pd.read_csv(io.BytesIO(raw), compression="gzip", low_memory=False,
                    dtype={"photo_id": str, "original_numeric_observation_id": str,
                           "original_taxon_id": str},
                    usecols=["photo_id", "original_numeric_observation_id",
                             "original_taxon_id", "original_frozen_ROI_morph",
                             "qualified_observed_on"])
    require(len(h)==HISTORICAL_PHOTOS and h.photo_id.is_unique,
            "Historical source photo/date identity duplicate or missing")
    merged = d.merge(h, on="photo_id", how="left",validate="one_to_one")
    old = merged.old_mois_or_new_fcp.eq("historical")
    require(old.sum()==HISTORICAL_PHOTOS and
            merged.loc[old,"original_numeric_observation_id"].notna().all() and
            merged.loc[old,"observation_id"].astype(str).eq(
                merged.loc[old,"original_numeric_observation_id"].astype(str)).all() and
            merged.loc[old,"inat_taxon_id"].astype(str).eq(
                merged.loc[old,"original_taxon_id"].astype(str)).all() and
            merged.loc[old,"morph"].fillna("").eq(
                merged.loc[old,"original_frozen_ROI_morph"].fillna("")).all(),
            "An original historical photograph, observation, taxon or hue changed")
    require(not (merged.observed_on.notna()&
                 merged.qualified_observed_on.notna()).any(),
            "An historic public-date match would overwrite original observed_on")
    merged["observed_on"] = merged.observed_on.fillna(
        merged.qualified_observed_on)
    keep = merged[["photo_id","observation_id","inat_taxon_id",
                   "morph","observed_on","latitude","longitude",
                   "old_mois_or_new_fcp"]].copy()
    q = labelled_calendar(keep)
    require(q.date.notna().sum()==DATED_PHOTOS and
            (q.date.notna() & q.classified).sum()==DATED_HUES and
            q.photo_id.is_unique,
            "Original dated classifiable ROI and attempted source denominator drift")
    return q


def symmetric_white_detection(c_s:float,c_w:float,p_s:float,p_w:float) -> dict:
    """Exact symmetric decomposition in percentage points of attempts."""
    require(all(0<=x<=1 for x in (c_s,c_w,p_s,p_w)),
            "Photo classification and white-label rates must be probabilities")
    rate = 100*(c_s*p_s-c_w*p_w)
    hue = 100*(c_s+c_w)/2*(p_s-p_w)
    classifiability = 100*(p_s+p_w)/2*(c_s-c_w)
    require(abs(rate-hue-classifiability)<1e-10,
            "Original photo-classifiability decomposition not additive")
    return {
        "summer_minus_winter_white_label_yield_per_all_attempts_pp":round(rate,6),
        "source_white_hue_share_conditional_on_classified_component_pp":round(hue,6),
        "source_four_hue_classifiability_component_pp":round(classifiability,6),
    }


def globe_by_zone(d:pd.DataFrame) -> list[dict]:
    q = d[d.date.notna() & d.zone.isin(
        ["north_extratropical","south_extratropical"]) &
        d.season.isin(["summer","winter"]) &
        d.latitude.between(-90,90)&d.longitude.between(-180,180)]
    rec=[]
    for zone,group in q.groupby("zone",observed=True):
        s=group[group.season.eq("summer")]
        w=group[group.season.eq("winter")]
        require(len(s)>0 and len(w)>0 and s.classified.sum()>0 and
                w.classified.sum()>0, "No real source classifiable season")
        cs=float(s.classified.mean())
        cw=float(w.classified.mean())
        ps=float(s.white.sum()/s.classified.sum())
        pw=float(w.white.sum()/w.classified.sum())
        rec.append({
            "hemisphere":zone,
            "summer_actual_original_photo_attempts":len(s),
            "winter_actual_original_photo_attempts":len(w),
            "summer_four_hue_classified_original_photos":int(s.classified.sum()),
            "winter_four_hue_classified_original_photos":int(w.classified.sum()),
            "summer_white_labeled_original_photos":int(s.white.sum()),
            "winter_white_labeled_original_photos":int(w.white.sum()),
            "summer_classifiability_fraction":round(cs,8),
            "winter_classifiability_fraction":round(cw,8),
            "summer_white_fraction_among_classified":round(ps,8),
            "winter_white_fraction_among_classified":round(pw,8),
            "summer_minus_winter_classifiability_pp":round(100*(cs-cw),6),
            "summer_minus_winter_white_conditional_pp":round(100*(ps-pw),6),
            **symmetric_white_detection(cs,cw,ps,pw),
        })
    require(len(rec)==2, "Both hemispheres required")
    return rec


def geocell_bootstrap(table:pd.DataFrame, metric:str, *,
                      draws:int=GEOCELL_BOOTSTRAPS, seed:int=20261017):
    if table.empty:
        return {"estimate_pp":None,"cluster_bootstrap95_pp":None,
                "status":"NO_PAIRED_CLASSIFIABLE_PHOTO_SOURCE"}
    estimate=float(table[metric].mean())
    # All years in the same approximate source cell are sampled together.
    block=table.groupby(["latitude_cell","longitude_cell"])[metric].agg(
        ["sum","count"])
    if len(block)<10:
        return {"estimate_pp":round(estimate,6),"cluster_bootstrap95_pp":None,
                "status":"HOLD_FEWER_THAN_10_GEOCELLS"}
    a=block["sum"].to_numpy(float)
    n=block["count"].to_numpy(float)
    rng=np.random.default_rng(seed)
    ix=rng.integers(0,len(a),size=(draws,len(a)))
    sim=a[ix].sum(axis=1)/n[ix].sum(axis=1)
    return {
        "estimate_pp":round(estimate,6),
        "cluster_bootstrap95_pp":[round(float(z),6)
                                  for z in np.quantile(sim,[0.025,0.975])],
        "status":"EXPLORATORY_GEOCELL_CLUSTER_PERCENTILE_NOT_MAXT",
    }


def site_year_selection(d:pd.DataFrame,width:float,min_obs:int=5,
                        min_taxa:int=2) -> list[dict]:
    """Align unknown+known ROI photo attempts to cells selected by known hues.

    Locality membership is selected using original classified photo-observation
    counts to be EXACTLY the same site×season-years as the MOIS local hue
    comparison, not an entirely new unmatched site panel.
    """
    q=d[d.date.notna()&d.zone.isin(
        ["north_extratropical","south_extratropical"]) &
        d.season.isin(["summer","winter"]) &
        d.latitude.between(-90,90)&d.longitude.between(-180,180)].copy()
    q["latitude_cell"]=np.floor((q.latitude+90)/width).astype(int)
    q["longitude_cell"]=np.floor((q.longitude+180)/width).astype(int)
    keys=["zone","latitude_cell","longitude_cell","season_year","season"]
    all_photo=q.groupby(keys,observed=True).agg(
        n_attempts=("photo_id","size"),n_classified=("classified","sum"),
        n_white=("white","sum")).reset_index()
    known=q[q.classified].groupby(keys,observed=True).agg(
        known_observation_IDs=("observation_id","nunique"),
        known_source_taxa=("inat_taxon_id","nunique")).reset_index()
    joined=all_photo.merge(known,on=keys,how="left",validate="one_to_one")
    index=["zone","latitude_cell","longitude_cell","season_year"]
    values=["n_attempts","n_classified","n_white",
            "known_observation_IDs","known_source_taxa"]
    tab=joined.pivot(index=index,columns="season",values=values).fillna(0)
    for name in ("known_observation_IDs","known_source_taxa"):
        for season in ("summer","winter"):
            if (name,season) not in tab.columns:
                tab[(name,season)]=0
    eligible=(
        (tab[("known_observation_IDs","summer")]>=min_obs)&
        (tab[("known_observation_IDs","winter")]>=min_obs)&
        (tab[("known_source_taxa","summer")]>=min_taxa)&
        (tab[("known_source_taxa","winter")]>=min_taxa))
    tab=tab.loc[eligible]
    values=[]
    for key,row in tab.iterrows():
        zone,latitude_cell,longitude_cell,year=key
        s_attempt=int(row[("n_attempts","summer")])
        w_attempt=int(row[("n_attempts","winter")])
        s_class=int(row[("n_classified","summer")])
        w_class=int(row[("n_classified","winter")])
        s_white=int(row[("n_white","summer")])
        w_white=int(row[("n_white","winter")])
        assert s_class>=min_obs and w_class>=min_obs
        values.append({
          "zone":zone,"latitude_cell":int(latitude_cell),
          "longitude_cell":int(longitude_cell),"season_year":int(year),
          "summer_all_photo_attempts":s_attempt,
          "winter_all_photo_attempts":w_attempt,
          "summer_classifiable":s_class,"winter_classifiable":w_class,
          "summer_white_labels":s_white,"winter_white_labels":w_white,
          "source_classifiability_season_difference_pp":
              100*(s_class/s_attempt-w_class/w_attempt),
          "source_white_share_season_difference_pp":
              100*(s_white/s_class-w_white/w_class),
        })
    pair=pd.DataFrame(values)
    out=[]
    for zone in ("north_extratropical","south_extratropical"):
        site=pair.loc[pair.zone.eq(zone)].copy() if len(pair) else pd.DataFrame()
        cl=geocell_bootstrap(site,"source_classifiability_season_difference_pp",
                             seed=20261017)
        hu=geocell_bootstrap(site,"source_white_share_season_difference_pp",
                             seed=20261017)
        out.append({
            "hemisphere":zone,"site_width_degrees":width,
            "minimum_original_classifiable_observation_IDs_per_season":min_obs,
            "minimum_classifiable_named_taxa_per_season":min_taxa,
            "n_selected_geocell_season_years":int(len(site)),
            "n_unique_geocells":int(site[[
                "latitude_cell","longitude_cell"]].drop_duplicates().shape[0])
                if len(site) else 0,
            "classified_and_missing_photo_attempts_in_paired_seasons":
                int(site.summer_all_photo_attempts.sum()+
                    site.winter_all_photo_attempts.sum()) if len(site) else 0,
            "classified_photo_records_in_paired_seasons":
                int(site.summer_classifiable.sum()+site.winter_classifiable.sum())
                if len(site) else 0,
            "classifiability_season_source_diagnostic":cl,
            "conditional_white_share_source_diagnostic":hu,
        })
    return out


def analyze(photos:Path,history:Path,output:Path) -> dict:
    require(not output.exists() or not any(output.iterdir()),
            "Never replace prior source photo/date/season results")
    d=date_and_original_hues(photos,history)
    global_rates=globe_by_zone(d)
    spatial=[]
    for width in WIDTHS:
        spatial.extend(site_year_selection(d,width))
    output.mkdir(parents=True,exist_ok=True)
    table=output/"mois_127504_source_photo_hue_classifiability_season_geocell.csv"
    pd.DataFrame(spatial).to_csv(table,index=False,lineterminator="\n")
    outcome={
        "schema":"mois_original127504_dated_hue_photo_season_classifiability_source_QC_v1",
        "original_source_photo_sha256":SOURCE_SHA,
        "final_historical_source_photo_date_sha256":FINAL_HISTORIC_SHA,
        "total_actual_original_photo_attempts":ORIGINAL_PHOTOS,
        "total_actual_original_classifiable_hues":ORIGINAL_HUES,
        "dated_source_photo_attempts":DATED_PHOTOS,
        "dated_and_classifiable_original_photo_ROIs":DATED_HUES,
        "source_original_photo_attempts_without_colour_but_with_date":
            DATED_PHOTOS-DATED_HUES,
        "global_photo_attempt_selection_identity":global_rates,
        "matched_geocell_season_year_source_classifiability":spatial,
        "geocell_source_table_sha256":sha256(table.read_bytes()),
        "photo_source_scope":[
           "The source classifiability / white detection identity is algebraic, NOT pigment biology.",
           "Original ROI technical missing is neither white nor coloured nor a known absent flower.",
           "Matched locality/year subset selected on classifiable source observations; this selection must not be treated as randomized.",
           "Coarse 0.25-degree cells remain tens of km across, not breeding populations.",
           "Photo-specific botany, independent observer effort and real observed-month climate not validated.",
           "Exploratory geographically clustered percentile intervals not multicomparison maxT.",
           "Existing MOIS genus/moisture and FCP intraspecific W/D estimands unchanged.",
        ]
    }
    (output/"mois_127504_source_photo_classifiability_season_QC_receipt.json").write_text(
        json.dumps(outcome,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    print("MOIS_OBSERVED_PHOTO_CLASSIFIABILITY_BY_SEASON",
          json.dumps({"full_photo_ids":ORIGINAL_PHOTOS,
            "dated_colored_ROIs":DATED_HUES,"geocell025":spatial[:2],
            "hemisphere_globally":global_rates},sort_keys=True),flush=True)
    return outcome


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--original-photo-zip",required=True,type=Path)
    p.add_argument("--final-historical-date-zip",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    q=p.parse_args()
    analyze(q.original_photo_zip,q.final_historical_date_zip,q.output)
