#!/usr/bin/env python3
"""MOIS: final original 267425-photo date-only union using official S3 + checked API.

Full S3 2026-09-27 photo_id→observation_uuid→observed_on is canonical for
99,732/100,543 historical original photos. For S3-missing original photo IDs,
an exact photo_id PLUS original numeric observation_id separately verified in
the *prior finished* iNaturalist API 9,600 pilot is permitted as a source-
tagged fallback. Do NOT replace even one S3 date with newer API metadata;
record dated-photo conflicts and keep unresolved photographs as missing.
All original photo pixels, hues, taxa, geography and annual moisture models
are byte/SHA frozen. These are automatic photo ROIs, NOT independently
botanically confirmed open flowers or pigment plasticity.
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
    HUES, MOIS_ROWS, MOIS_HUES, PHOTO_DATA, PHOTO_RECEIPT,
    read_signed, require, labelled_calendar, photo_season_table,
)
from scripts.analysis.mois_taxon_season_composition_decomposition import (
    symmetric_photo_decomposition,
)
from scripts.analysis.mois_original_photo_community_season_geoscale import (
    matched_community, cell_cluster_interval,
)
from scripts.analysis.mois_seasonal_spatial_identifiability_20261017 import (
    eligible_site_photo_support,
)
from scripts.analysis.mois_historical9600_taxonomic_breadth_vs_local_season_support import (
    DATE_INPUT, DATE_RECEIPT, DATE_SHA, PHOTO_SHA,
)

HISTORICAL_N=100543
S3_FILENAME="mois_historical_100543_original_photo_id_s3_uuid_date_real_snapshot.csv.gz"
S3_RECEIPT="mois_full_historical_photo_dates_official_monthly_s3_source_receipt.json"
S3_TABLE_SHA="d9f599d4ff6d479d8714344b014e604ade1af90a3e8dd8640a19d1be48407566"
HISTORIC_OUTPUT="mois_original100543_s3_or_verified_API_dates_photo_source.csv.gz"
COMMUNITY_CSV="mois_real127504_hues_community_summer_winter_source_geocell_year.csv"
SEASON_CSV="mois_real127504_hues_summer_winter_photo_source_season_rates.csv"
OUTPUT_JSON="mois_real127504_original_photo_S3_plus_verified_API_local_season_result.json"
SCALES=(0.25,0.5,1.0,2.0,5.0)
NAMES=("north_extratropical","south_extratropical")


def checksum(data:bytes)->str:
    return hashlib.sha256(data).hexdigest()


def source_tables(photo_zip:Path,s3_zip:Path,api_zip:Path):
    photo,r,sha=read_signed(photo_zip,PHOTO_DATA,PHOTO_RECEIPT,
        "source_photo_level_csv_sha256",
        ["photo_id","observation_id","inat_taxon_id","morph",
         "observed_on","old_mois_or_new_fcp","latitude","longitude"])
    require(len(photo)==MOIS_ROWS and sha==PHOTO_SHA and
            int(photo.morph.isin(HUES).sum())==MOIS_HUES,
            "Source photo pixels/nominal hue labels changed")
    with zipfile.ZipFile(s3_zip) as z:
        raw=z.read(S3_FILENAME)
        s3proof=json.loads(z.read(S3_RECEIPT))
    require(checksum(raw)==S3_TABLE_SHA and
            s3proof["recovered_historical_photo_table_SHA256"]==S3_TABLE_SHA and
            s3proof["photo_ids_with_verified_snapshot_observation_date"]==99732 and
            s3proof["full_historical_attempted_source_photo_count"]==HISTORICAL_N and
            s3proof["source_photo_sha256"]==PHOTO_SHA,
            "Actual 99732-photograph official original S3 snapshot drift")
    s3=pd.read_csv(io.BytesIO(raw),compression="gzip",
                   dtype={"original_photo_id":str,
                          "original_numeric_observation_id":str,
                          "original_taxon_id":str},low_memory=False)
    with zipfile.ZipFile(api_zip) as z:
        api_csv=z.read(DATE_INPUT)
        api_proof=json.loads(z.read(DATE_RECEIPT))
    require(checksum(api_csv)==DATE_SHA and
            api_proof["new_original_photo_date_CSV_sha256"]==DATE_SHA and
            api_proof["exact_original_photo_observation_and_current_date_matches"]==9597 and
            api_proof["original_full_source_photo_sha256"]==PHOTO_SHA,
            "Current public photo and numeric observation ID backup evidence drift")
    api=pd.read_csv(io.BytesIO(api_csv),dtype=str)
    return photo,s3,api,sha,api_proof,s3proof


def rejoin(photo:pd.DataFrame,s3:pd.DataFrame,api:pd.DataFrame):
    """Preserve one strict photo source row; keep both API/S3 date columns."""
    hist=photo[photo.old_mois_or_new_fcp.eq("historical")].copy()
    require(len(hist)==HISTORICAL_N and hist.observed_on.isna().all() and
            hist.photo_id.is_unique and hist.observation_id.is_unique,
            "Original historical IDs or blank original date source changed")
    require(len(s3)==HISTORICAL_N and s3.original_photo_id.is_unique,
            "Original full S3 photo-index coverage changed")
    joined=hist.merge(s3,left_on="photo_id",right_on="original_photo_id",
                      validate="one_to_one",how="left")
    require(len(joined)==HISTORICAL_N and
            joined.original_photo_id.notna().all(),
            "Original 100543 photographed IDs not all present in snapshot rejoin")
    for col,right in (("photo_id","original_photo_id"),
                      ("observation_id","original_numeric_observation_id"),
                      ("inat_taxon_id","original_taxon_id")):
        require(joined[col].astype(str).eq(joined[right].astype(str)).all(),
                "Historical photo/taxon/observation ID shifted between frozen sources")
    require(joined.morph.fillna("").eq(
        joined.original_frozen_ROI_morph.fillna("")).all(),
        "Original historical hue predictions have changed")
    require(len(api)==9600 and api.original_photo_id.is_unique,
            "Original 9600 pilot ID source changed")
    api_rows=api[["original_photo_id","original_observation_id","source_taxon_id",
                  "provenance_status","retrieved_observed_on"]].copy()
    api_rows.columns=["photo_id","api_original_observation_id",
                      "api_source_taxon_id","api_provenance_status",
                      "api_public_observed_on"]
    merged=joined.merge(api_rows,on="photo_id",how="left",validate="one_to_one")
    old_api=merged.api_original_observation_id.notna()
    require(int(old_api.sum())==9600 and
            merged.loc[old_api,"observation_id"].astype(str).eq(
                merged.loc[old_api,"api_original_observation_id"].astype(str)).all() and
            merged.loc[old_api,"inat_taxon_id"].astype(str).eq(
                merged.loc[old_api,"api_source_taxon_id"].astype(str)).all(),
            "A completed API source photo/observation/taxon identity was changed")
    s3_date=merged.snapshot_observed_on.notna()
    valid_api=merged.api_provenance_status.eq(
        "source_photo_and_observation_verified_with_current_public_date") & (
        merged.api_public_observed_on.notna())
    both=s3_date & valid_api
    disagreement=both & merged.snapshot_observed_on.ne(merged.api_public_observed_on)
    api_fallback=~s3_date & valid_api
    require(int(s3_date.sum())==99732 and int(api_fallback.sum())==85 and
            int(both.sum())==9512 and int(disagreement.sum())==1,
            "Source original API/S3 photo-date overlap or conflict changed")
    # Snapshot date retains priority on the one public metadata disagreement.
    merged["qualified_date_source"]=np.select(
        [s3_date,api_fallback],["s3_20260927","iNat_current_API_20261017"],
        default="unresolved_original_photo_date")
    merged["qualified_observed_on"]=merged.snapshot_observed_on.where(
        s3_date,merged.api_public_observed_on.where(api_fallback))
    approved=merged.loc[merged.qualified_observed_on.notna(),
                        ["photo_id","qualified_observed_on"]].copy()
    source=photo.merge(approved,on="photo_id",how="left",validate="one_to_one")
    require(len(source)==MOIS_ROWS and not (source.observed_on.notna()&
         source.qualified_observed_on.notna()).any(),
         "Recovered date would overwrite an immutable original observed_on")
    source["effective_observed_on"]=source.observed_on.where(
        source.observed_on.notna(),source.qualified_observed_on)
    as_date=source[["photo_id","observation_id","inat_taxon_id","morph",
        "latitude","longitude","old_mois_or_new_fcp",
        "effective_observed_on"]].rename(
          columns={"effective_observed_on":"observed_on"})
    p=labelled_calendar(as_date)
    require(int(p.classified.sum())==MOIS_HUES and
            int((p.classified&p.date.notna()).sum())==127504 and
            int(p.date.notna().sum())==261578,
            "Original photo/date four hue denominators changed")
    return merged,p,{
        "original_historical_photo_ids":HISTORICAL_N,
        "historical_date_S3_source_verified":int(s3_date.sum()),
        "historic_public_API_photo_id_and_numeric_observation_id_verified_fallback":int(
            api_fallback.sum()),
        "history_original_dates_total_with_both_source_precedence":int(
            merged.qualified_observed_on.notna().sum()),
        "still_unresolved_original_historical_dates":int(
            merged.qualified_observed_on.isna().sum()),
        "current_API_and_monthly_S3_dates_both_valid":int(both.sum()),
        "current_API_vs_official_S3_observed_on_disagreements":int(disagreement.sum()),
        "S3_primary_date_retained_on_conflicting_API_observation":True,
        "newly_source_date_recovered_historical_original_four_hue_photos":int(
            (merged.qualified_observed_on.notna() & merged.morph.isin(HUES)).sum()),
        "original_photos_with_observation_date_after_all_verified_rejoins":int(
            p.date.notna().sum()),
        "real_original_photo_hue_and_date_after_all_verified_rejoins":int(
            (p.classified&p.date.notna()).sum()),
        "dated_and_classified_original_source_nominal_taxa":int(
            p.loc[p.classified&p.date.notna(),"inat_taxon_id"].nunique()),
    }


def seasonal_taxon_decomposition(p:pd.DataFrame,zone:str):
    q=p[p.classified&p.date.notna()&p.zone.eq(zone)&
        p.season.isin(["summer","winter"])&
        p.latitude.between(-90,90)&p.longitude.between(-180,180)]
    S=q[q.season.eq("summer")]
    W=q[q.season.eq("winter")]
    st=set(S.inat_taxon_id)
    wt=set(W.inat_taxon_id)
    shared=st & wt
    if len(shared)<2:return None
    tab=q[q.inat_taxon_id.isin(shared)].groupby(
        ["inat_taxon_id","season"],observed=True).white.agg(
        ["size","sum"]).unstack(fill_value=0)
    a=np.column_stack((
        tab["size"]["summer"].to_numpy(float),
        tab["size"]["winter"].to_numpy(float),
        tab["sum"]["summer"].to_numpy(float),
        tab["sum"]["winter"].to_numpy(float)))
    within,mix,matched_total,taxon_equal=symmetric_photo_decomposition(a)
    common_S=S[S.inat_taxon_id.isin(shared)]
    common_W=W[W.inat_taxon_id.isin(shared)]
    absent_S=S[~S.inat_taxon_id.isin(shared)]
    absent_W=W[~W.inat_taxon_id.isin(shared)]
    correction_S=100*(len(absent_S)/len(S))*(
        absent_S.white.mean()-common_S.white.mean())
    correction_W=-100*(len(absent_W)/len(W))*(
        absent_W.white.mean()-common_W.white.mean())
    total=100*(S.white.mean()-W.white.mean())
    require(abs(total-(within+mix+correction_S+correction_W))<1e-8,
            "Full photo observed season colour composition identity not additive")
    return {
        "summer_original_classified_photos":int(len(S)),
        "winter_original_classified_photos":int(len(W)),
        "summer_original_white_photo_ROIs":int(S.white.sum()),
        "winter_original_white_photo_ROIs":int(W.white.sum()),
        "summer_photographed_named_taxa":len(st),
        "winter_photographed_named_taxa":len(wt),
        "common_taxa_both_seasons":len(shared),
        "total_photo_taxa_union":len(st|wt),
        "shared_taxa_photo_Jaccard":round(len(shared)/len(st|wt),6),
        "original_source_photo_white_summer_minus_winter_pp":round(total,6),
        "common_taxon_photo_white_rate_difference_pp":round(float(within),6),
        "common_taxon_relative_photo_frequency_shift_pp":round(float(mix),6),
        "summer_only_taxon_photo_support_correction_pp":round(float(correction_S),6),
        "winter_only_taxon_photo_support_correction_pp":round(float(correction_W),6),
        "common_taxa_equal_taxon_photo_white_difference_pp":round(float(taxon_equal),6),
        "is_exact_additive_source_photo_algebra_NOT_causal_biological_mechanism":True,
    }


def local_support(p:pd.DataFrame):
    q=p[p.classified&p.date.notna()&
        p.zone.isin(["north_extratropical","south_extratropical"])&
        p.season.isin(["summer","winter"])&
        p.latitude.between(-90,90)&p.longitude.between(-180,180)]
    result=[]
    for width in (0.25,0.5,1.,2.,5.):
        for minimum in (3,5,10):
            communities=matched_community(q,width=width,min_obs=minimum,min_taxa=2)
            for zone in ("north_extratropical","south_extratropical"):
                local=communities[communities.zone.eq(zone)]
                ci=cell_cluster_interval(local,draws=999,
                       seed=20261017+int(100*width)+minimum)
                result.append({
                    "geographic_cell_degrees":width,
                    "minimum_distinct_original_observation_ID_per_summer_and_winter":minimum,
                    "minimum_named_original_taxa_per_season":2,
                    "hemisphere":zone,
                    "same_cell_same_meteorological_season_year_comparisons":int(len(local)),
                    "n_source_classified_photos_in_both_seasons":int(
                        (local.summer_photos.sum()+local.winter_photos.sum()) if len(local) else 0),
                    "n_independent_geocells":int(local[
                      ["latitude_cell","longitude_cell"]].drop_duplicates().shape[0]),
                    "equal_site_year_summer_white_photo_difference_pp":
                        None if ci is None else ci["equal_cell_year_white_difference_pp"],
                    "locality_cluster_percentile_bootstrap95_pp":
                        None if ci is None else ci["bootstrap_95pp"],
                    "inference_status":"NO_REPLICATES" if ci is None
                        else ci["uncertainty_status"],
                })
    for zone in ("north_extratropical","south_extratropical"):
        sub=q[q.zone.eq(zone)]
        z=eligible_site_photo_support(sub,.25,True,2)
        result.append({"geographic_cell_degrees":.25,
            "minimum_distinct_original_observation_ID_per_summer_and_winter":2,
            "minimum_named_original_taxa_per_season":1,
            "hemisphere":zone,
            "comparison_unit":"same_named_taxon_WITHIN_same_cell_and_seasonyear",
            "same_cell_same_meteorological_season_year_comparisons":int(z["matched_groups"]),
            "n_source_classified_photos_in_both_seasons":int(z["matched_photos"])})
    return result


def run(photo_archive:Path,s3_archive:Path,api_archive:Path,outdir:Path):
    require(not outdir.exists() or not any(outdir.iterdir()),
            "Do not overwrite a complete signed ecological source result")
    photos,s3,api,sha,apip,s3p=source_tables(
        photo_archive,s3_archive,api_archive)
    history,original,source=rejoin(photos,s3,api)
    by_zone={z:seasonal_taxon_decomposition(original,z)
             for z in ("north_extratropical","south_extratropical")}
    community=local_support(original)
    outdir.mkdir(parents=True,exist_ok=True)
    path=outdir/HISTORIC_OUTPUT
    history.to_csv(path,index=False,
        compression={"method":"gzip","mtime":0},
        lineterminator="\n")
    quarter=photo_season_table(original)
    quarter["source"]="date-only-original-full-S3-and-strict-current-API-fallback"
    season=outdir/SEASON_CSV
    quarter.to_csv(season,index=False,lineterminator="\n")
    local=outdir/COMMUNITY_CSV
    pd.DataFrame(community).to_csv(local,index=False,lineterminator="\n")
    result={
        "schema":"mois_original267425_photos_full100543_S3_verified_API_dates_local_seasons_v1",
        "photo_original_full_source_SHA256":sha,
        "official_monthly_20260927_full_S3_99732_historic_dates_archive_SHA256":S3_TABLE_SHA,
        "independent_original_9600_current_iNat_API_dates_archive_SHA256":DATE_SHA,
        "independent_current_API_date_sample_count":9600,
        "original_100543_historical_photo_IDs_exactly_preserved":True,
        "full_original_measured_photo_ids":MOIS_ROWS,
        "full_original_genuine_four_hue_classified_photo_ids":MOIS_HUES,
        "historical_original_photo_classifications_frozen":45953,
        "original_historical_and_FCP_source_photo_IDs_never_double_counted":True,
        "original_source_date_present_classified_before_recovery":81867,
        "provider_observed_on_is_current_snapshot_or_API_metadata_not_original_EXIF":True,
        "one_conflicting_date_preserved_as_S3_20260927_canonical":True,
        "full_historical_join":source,
        "zone_season_original_photo_taxonomic_composition_decomposition":by_zone,
        "site_photo_community_and_taxon_local_source_season_support":community,
        "source_historical_photo_rows_csv_SHA256":checksum(path.read_bytes()),
        "all_source_photo_season_counts_csv_SHA256":checksum(season.read_bytes()),
        "matched_locality_photo_community_sensitivity_CSV_SHA256":checksum(local.read_bytes()),
        "any_photo_specific_open_flower_organs_botanically_confirmed_by_this_work":False,
        "original_colour_roi_D_W_or_50km_FCP_polymorphism_estimator_imported":False,
        "historical_MOIS_annual_genus_rainfall_vapour_maxT_result_overwritten":False,
        "observation_year_monthly_temperature_solar_precipitation_VPD_or_deficit_modelled":False,
        "noncausal_limits":[
            "Different taxa and photographer effort contribute to any matched photographed COMMUNITY summer-winter difference.",
            "0.25 degree coordinates cover tens of km; not field-level populations.",
            "Photo-specific open petal anthesis not botanically verified by classification or public observation-level tags.",
            "999 geographic-cell cluster percentile bootstraps are exploratory and not maxT familywise intervals.",
            "Current public iNaturalist observed_on is mutable and S3/API snapshots are versioned separately.",
            "Unresolved original photographic dates are never inferred from nearby climate or other taxa.",
        ],
    }
    (outdir/OUTPUT_JSON).write_text(json.dumps(
        result,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    print("MOIS_FINAL_HISTORICAL_DATES_PLUS_LOCAL_PHOTO_SEASON",
      json.dumps({"dated_classified":source[
          "real_original_photo_hue_and_date_after_all_verified_rejoins"],
          "taxa":source["dated_and_classified_original_source_nominal_taxa"],
          "historic_unresolved":source["still_unresolved_original_historical_dates"],
          "local_0p25_n5":{z:[r for r in community if
                r["geographic_cell_degrees"]==.25 and
                r["minimum_distinct_original_observation_ID_per_summer_and_winter"]==5 and
                r["hemisphere"]==z][0][
                "same_cell_same_meteorological_season_year_comparisons"]
                for z in by_zone}},sort_keys=True),flush=True)
    return result


if __name__=="__main__":
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument("--original-photo-zip",type=Path,required=True)
    a.add_argument("--full-official-s3-zip",type=Path,required=True)
    a.add_argument("--old9600-api-zip",type=Path,required=True)
    a.add_argument("--output",type=Path,required=True)
    q=a.parse_args()
    run(q.original_photo_zip,q.full_official_s3_zip,
        q.old9600_api_zip,q.output)
