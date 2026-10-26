#!/usr/bin/env python3
"""Audit the *composition*, not merely size, of MOIS 1,043,930 real photos.

An exact original-photo-ID and original observation-ID audit of frozen
1,043,930 pixel-attempt source and 390,541 strict continuous nine-fraction
measured flower ROIs. Membership in FCP's 447,500-photo fixed 100/taxon
cohort comes from its original 447,500 photo ID CSV, NOT from a guessed
historical source label. Use frozen 42,111-taxon, 93,209,351-photo
2026-09-27 full-index species opportunity counts to compare taxon support
but NOT to infer photo-level latitude/season inclusion probabilities.

Output aggregate metrics only, not individual images, photo IDs or
personal observer IDs. No original-photo image downloads or inferences.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scripts.analysis.mois_real1043930_global_zonal_colour import latitude_zone, season_label

SCHEMA="mois_real1043930_original_photo_population_structure_and_missingness_v1"
TARGET=1_043_930
STRICT=390_541
AMBIGUOUS=36_224
FOUR_HUES=459_235
PHOTO_SHA="b20bbefe15020c21c3b5b0f7edf3f8e08148a4ee7dda8bf772a89bf6b126f447"
PHOTO_UNION_SHA="9a4f3dae6fbcee967b71cb29576e71d52c48e1f1f790b7d3ea4750c4bcb0e6d4"
FCP_SHA="8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6"
UPSTREAM_SHA="ef26c16252dd785f2d08bde3b18938c7624d73a54386fa3d3974dceb09ed97a9"
FRAME_TAXA=42_111
FRAME_TAXA_WITH_PHOTOS=42_095
FRAME_PHOTOS=93_209_351
PALETTE=("white","yellow","orange","red","pink","magenta","purple","blue","bronze")

def require(ok,message):
    if not ok:raise ValueError("MOIS_STRUCTURE: "+message)

def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda:f.read(4<<20),b""):h.update(block)
    return h.hexdigest()

def finite_pct(n,d):
    return 100.0*float(n)/d if d else None

def top_share(counts,k):
    x=np.asarray(counts,dtype=np.int64)
    return finite_pct(int(np.sort(x)[-k:].sum()),int(x.sum()))

def concentration(counts):
    x=np.asarray(counts,dtype=float)
    if not len(x) or x.sum()==0:return None
    p=x/x.sum()
    return float(1/np.sum(p*p))

def source_rows(frame,grouping,label="source_group"):
    rows=[]
    for source,d in frame.groupby(grouping,observed=True,dropna=False,sort=True):
        known=d.loc[d.continuous9_strict]
        ids=d.inat_taxon_id.astype(str)
        rows.append({
            label:str(source),
            "photo_attempts":int(len(d)),
            "pct_of_all_1043930_photos":finite_pct(len(d),TARGET),
            "old_four_hue_classified":int(d.original_four_hue.isin(
                ("white","yellow_orange","red_pink","blue_purple")).sum()),
            "strict_nine_component_success":int(len(known)),
            "strict9_pct_of_attempts":finite_pct(len(known),len(d)),
            "strict9_pct_of_all_successes":finite_pct(len(known),STRICT),
            "ambiguous_nine_component_sensitivity":int(
                d.continuous9_exploratory_ambiguous.sum()),
            "strict9_colour_missing":int(len(d)-len(known)),
            "distinct_nominal_taxa_attempted":int(ids.nunique()),
            "distinct_nominal_taxa_strict9":int(
                known.inat_taxon_id.astype(str).nunique()),
            "distinct_named_genera_strict9":int(known.genus.nunique()),
            "distinct_original_observations":int(d.observation_id.replace("",np.nan).nunique()),
            "source_observation_ID_missing":int(d.observation_id.isna().sum()+
                                                 d.observation_id.eq("").sum()),
            "median_photos_per_distinct_observation":(
                float(d.loc[d.observation_id.ne(""),"observation_id"].value_counts().median())
                if d.observation_id.ne("").any() else None),
        })
    return rows

def audit(photo,receipt,union,union_receipt,fcp,fcp_receipt,upstream,geography,geo_receipt,out):
    photo,receipt,union,union_receipt,fcp,fcp_receipt,upstream,geography,geo_receipt,out=map(
        Path,(photo,receipt,union,union_receipt,fcp,fcp_receipt,upstream,geography,geo_receipt,out))
    require(not out.exists(),"Refuse overwrite source population composition audit")
    a=json.loads(receipt.read_text(encoding="utf-8"))
    u=json.loads(union_receipt.read_text(encoding="utf-8"))
    f=json.loads(fcp_receipt.read_text(encoding="utf-8"))
    require(a.get("schema")=="mois_real1043930_original9_with_FCP447500_independent_source_reuse_v1"
            and a.get("latest_real_photo_ID_attempts")==TARGET
            and a.get("strict_original9_vectors_now")==STRICT
            and a.get("ambiguous_original9_vectors_now")==AMBIGUOUS
            and a.get("output_photo_lake_SHA256")==PHOTO_SHA==sha(photo),
            "Not exact original photo 9-composition full-source archive")
    require(u.get("source_pixel_attempts")==TARGET
            and u.get("photo_parquet_sha256")==PHOTO_UNION_SHA==sha(union),
            "Original source merged photograph observation/photo identity changed")
    require(f.get("schema")=="fcp_all4730_permanent_measured_original_photo_results_v1"
            and f.get("genuinely_measured_photo_rows")==447500
            and f.get("genuinely_four_state_classified_photos")==227107
            and f.get("source_archive_sha256")==FCP_SHA==sha(fcp),
            "FCP 100photo per 4475taxa raw original source no longer immutable")
    require(sha(upstream)==UPSTREAM_SHA,
            "Original 93m-source 42,111 species opportunity photo counts SHA changed")
    origin=pd.read_parquet(photo,columns=["photo_id","inat_taxon_id","species",
            "source_arm_original","original_four_hue",
            "continuous9_strict","continuous9_exploratory_ambiguous",
            "continuous9_white","latitude","longitude","date"])
    prior=pd.read_parquet(union,columns=["photo_id","observation_id","arm",
             "inat_taxon_id","latitude","longitude","date","measurement_status"])
    require(len(origin)==len(prior)==TARGET
            and origin.photo_id.astype(str).is_unique
            and prior.photo_id.astype(str).is_unique
            and origin.photo_id.astype(str).tolist()==prior.photo_id.astype(str).tolist(),
            "Photo union and colour trait map must be one-to-one in identical photo order")
    require(origin.inat_taxon_id.astype(str).tolist()==prior.inat_taxon_id.astype(str).tolist()
            and origin.source_arm_original.astype(str).tolist()==prior.arm.astype(str).tolist(),
            "Source photo ID attribution or taxon was changed by nine-vector restoration")
    d=origin.copy()
    d["photo_id"]=d.photo_id.astype(str)
    d["observation_id"]=prior.observation_id.fillna("").astype(str)
    d["genus"]=d.species.fillna("").astype(str).str.split().str[0]
    geom=json.loads(geo_receipt.read_text(encoding="utf-8"))
    require(geom.get("real_original_source_photo_attempts")==TARGET
            and geom.get("output_geography_photo_ID_parquet_SHA256")==sha(geography),
            "Actual source GSHHG landmass pixel/photo-ID SHA changed")
    land=pd.read_parquet(geography,columns=["photo_id","land_type","landmass"])
    require(len(land)==TARGET
            and land.photo_id.astype(str).tolist()==d.photo_id.tolist(),
            "GSHHG original photo geographic IDs not exactly aligned")
    d["land_type"]=land.land_type.fillna("unassigned").astype(str).to_numpy()
    d["landmass"]=land.landmass.fillna("").astype(str).to_numpy()
    d["original_four_hue"]=d.original_four_hue.fillna("").astype(str)
    d["inat_taxon_id"]=d.inat_taxon_id.fillna("").astype(str)
    d["strict"]=d.continuous9_strict.astype(bool)
    d["ambiguity"]=d.continuous9_exploratory_ambiguous.astype(bool)
    require(d.photo_id.str.fullmatch(r"\d+").all()
            and int(d.strict.sum())==STRICT
            and int(d.ambiguity.sum())==AMBIGUOUS
            and not (d.strict&d.ambiguity).any()
            and int(d.original_four_hue.isin(
                ("white","yellow_orange","red_pink","blue_purple")).sum())==FOUR_HUES,
            "Original source photo classifiability and photo IDs drifted")
    fcp_src=pd.read_csv(fcp,compression="gzip",
        usecols=["photo_id","inat_taxon_id","global_classifiable"],
        dtype={"photo_id":str,"inat_taxon_id":str},low_memory=False)
    require(len(fcp_src)==447500 and fcp_src.photo_id.is_unique,
            "The FCP fixed 100photo cohort IDs were not preserved")
    fcp_photos=set(fcp_src.photo_id)
    d["original_FCP447500_member"]=d.photo_id.isin(fcp_photos)
    require(int(d.original_FCP447500_member.sum())==447500
            and d.loc[d.original_FCP447500_member,"source_arm_original"].eq("historical_main").all(),
            "FCP source exact photo-ID membership is different from frozen main archive")
    fcp_subset=d.loc[d.original_FCP447500_member,["photo_id","inat_taxon_id",
                         "continuous9_strict"]].set_index("photo_id")
    sample=fcp_src.set_index("photo_id").join(fcp_subset,
                                              how="left",rsuffix="_MOIS")
    require(sample.continuous9_strict.notna().all(),
            "FCP exact source photo ID has no original MOIS photo identity")
    # Legacy MOIS and FCP may have different ACCEPTED TAXON assignments
    # for the same photo ID. Count discrepancies, do not silently overwrite
    # either classification or discard 447,500 genuine image observations.
    fcp_original_taxon_disagreements=int(
        (~sample.inat_taxon_id.eq(sample.inat_taxon_id_MOIS)).sum())
    require(int(sample.continuous9_strict.sum())==227107,
            "Source FCP photo 9-vector success is inconsistent")
    d["collection_source"]=np.where(d.original_FCP447500_member,
        "FCP_447500_fixed100taxon",
        np.where(d.source_arm_original.eq("historical_main"),
                 "historical_main_nonFCP",d.source_arm_original))
    # Reliable original photo dates and source acquisition provenance
    date=pd.to_datetime(d.date,errors="coerce")
    d["year"]=date.dt.year.astype("Int64")
    d["month"]=date.dt.month
    d["latitude"]=pd.to_numeric(d.latitude,errors="coerce")
    d["longitude"]=pd.to_numeric(d.longitude,errors="coerce")
    d["hemisphere"]=np.where(d.latitude>=0,"north",
                       np.where(d.latitude<0,"south","unknown"))
    d["lat_zone"]=[latitude_zone(v) for v in d.latitude.to_numpy(float)]
    d["season"]=[season_label(lat,month) for lat,month
                 in zip(d.latitude.to_numpy(float),d.month.to_numpy(dtype=float,na_value=np.nan))]
    d["lat5"]=np.floor((d.latitude+90)/5)
    d["lon5"]=np.floor((d.longitude+180)/5)
    d["lat5"]=d.lat5.fillna(-1).astype(int)
    d["lon5"]=d.lon5.fillna(-1).astype(int)
    # Colour success is missing not white or non-white.
    d["white_fraction"]=pd.to_numeric(d.continuous9_white,errors="coerce")
    require(d.loc[~(d.strict|d.ambiguity),"white_fraction"].isna().all(),
            "Unmeasured source photo silently received a white colour")
    # Frozen complete upstream source counts: exact species IDs and >=0 source photos
    counts=pd.read_csv(upstream,dtype={"inat_taxon_id":str},low_memory=False)
    require(len(counts)==FRAME_TAXA
            and counts.inat_taxon_id.is_unique
            and int(pd.to_numeric(counts.n_unique_source_photo_ids,errors="raise").sum())==FRAME_PHOTOS
            and int(counts.n_unique_source_photo_ids.gt(0).sum())==FRAME_TAXA_WITH_PHOTOS,
            "Original all-photo eligible 42,111 taxa source counts differ")
    source_photo_totals=counts.set_index("inat_taxon_id").n_unique_source_photo_ids.astype(int)
    d["upstream_frame_member"]=d.inat_taxon_id.isin(source_photo_totals.index)
    # This matters: an old source could genuinely include taxa outside the
    # frozen 42,111 eligible taxa, not a missing measurement.
    photo_attempts=d.groupby("inat_taxon_id").size().rename("photo_attempts")
    photo_strict=d[d.strict].groupby("inat_taxon_id").size().rename("strict9_photo_ROIs")
    photo_4hue=d[d.original_four_hue.isin((
           "white","yellow_orange","red_pink","blue_purple"))].groupby(
               "inat_taxon_id").size().rename("legacy4_hue_ROIs")
    taxa=counts[["inat_taxon_id","species","n_unique_source_photo_ids"]].copy()
    taxa=taxa.set_index("inat_taxon_id")
    taxa=taxa.join(photo_attempts).join(photo_strict).join(photo_4hue)
    for col in ("photo_attempts","strict9_photo_ROIs","legacy4_hue_ROIs"):
        taxa[col]=taxa[col].fillna(0).astype(int)
    valid=d[d.strict]
    photo_cohorts=source_rows(d,"collection_source")
    arm_cohorts=source_rows(d,"source_arm_original","original_source_arm")
    geographic=source_rows(d,"lat_zone","latitude_zone")
    time=source_rows(d,"season","meteorological_season")
    hemisphere=source_rows(d,"hemisphere","hemisphere")
    island=source_rows(d,"land_type","island_mainland")
    cohort_hemi=source_rows(d,
       d["collection_source"]+"|"+d["hemisphere"],"source_by_hemisphere")
    cohort_island=source_rows(d,
       d["collection_source"]+"|"+d["land_type"],"source_by_island")
    within_genus=d[d.strict].groupby("genus").size()
    per_taxon=valid.groupby("inat_taxon_id").size()
    all_obs=d.observation_id.replace("",np.nan)
    complete_observations=all_obs.dropna()
    observation_reuse=complete_observations.value_counts()
    # Stratify by source-available taxon opportunity to identify where
    # more photographs cannot unlock new taxa versus genuine zero coverage.
    binned=pd.cut(taxa.n_unique_source_photo_ids,
        bins=[-1,0,9,99,999,9999,10**18],
        labels=["0","1-9","10-99","100-999","1000-9999","10000+"])
    opportunity=[]
    for key,subset in taxa.groupby(binned,observed=True):
        opportunity.append({
            "original_93m_opportunity_bin":str(key),
            "total_frame_taxa":len(subset),
            "original_photo_opportunities":int(subset.n_unique_source_photo_ids.sum()),
            "taxa_with_any_processed_image":int(subset.photo_attempts.gt(0).sum()),
            "taxa_with_any_strict9_photo":int(subset.strict9_photo_ROIs.gt(0).sum()),
            "strict_nine_photo_ROIs":int(subset.strict9_photo_ROIs.sum()),
        })
    # Count controls by independent units of analysis rather than photos.
    matched={
      "strict_taxon_in_both_hemispheres_min2_each":0,
      "strict_taxon_in_both_island_and_mainland_min2_each":0,
      "strict_taxon_in_both_hemispheres_min5_each":0,
      "strict_taxon_region_5deg_source_arm_summer_winter_min2_each":0,
      "strict_source_arm_genus_5deg_multi_taxa_count":0,
    }
    hh=valid.groupby(["inat_taxon_id","hemisphere"]).size().unstack(fill_value=0)
    if {"north","south"}.issubset(hh.columns):
        for n in (2,5):
            matched[f"strict_taxon_in_both_hemispheres_min{n}_each"]=int(
                (hh.north.ge(n)&hh.south.ge(n)).sum())
    landpair=valid[valid.land_type.isin(("island","mainland"))]
    landtax=landpair.groupby(["inat_taxon_id","land_type"]).size().unstack(fill_value=0)
    if {"island","mainland"}.issubset(landtax.columns):
        matched["strict_taxon_in_both_island_and_mainland_min2_each"]=int(
            (landtax.island.ge(2)&landtax.mainland.ge(2)).sum())
    cond=valid[valid.season.isin(("summer","winter"))]
    level=["inat_taxon_id","source_arm_original","lat5","lon5"]
    ss=cond.groupby(level+["season"],dropna=False).size().unstack(fill_value=0)
    if {"summer","winter"}.issubset(ss.columns):
        matched["strict_taxon_region_5deg_source_arm_summer_winter_min2_each"]=int(
            (ss.summer.ge(2)&ss.winter.ge(2)).sum())
    multi=valid.groupby(["source_arm_original","genus","lat5","lon5"],
                        observed=True).inat_taxon_id.nunique()
    matched["strict_source_arm_genus_5deg_multi_taxa_count"]=int(multi.ge(2).sum())
    strict_zone=d.groupby("lat_zone",observed=True)[["strict","ambiguity"]].sum()
    n_z=d.groupby("lat_zone",observed=True).size()
    all_regions=int(n_z.size)
    year_table=[]
    for year,sub in d.groupby("year",dropna=False):
        year_table.append({
           "year":str(year),
           "original_photo_attempts":len(sub),
           "strict9_photo_ROIs":int(sub.strict.sum()),
           "strict9_photo_rate":finite_pct(int(sub.strict.sum()),len(sub)),
        })
    report={
       "schema":SCHEMA,
       "source_photo_grain":"one source photo_id per attempted image; original observation_id separately restored",
       "full_93m_frame_snapshot":"2026-09-27",
       "source_metadata_total_photoIDs":FRAME_PHOTOS,
       "source_metadata_distinct_taxa":FRAME_TAXA,
       "source_metadata_taxa_with_photos":FRAME_TAXA_WITH_PHOTOS,
       "original_image_photo_attempts":TARGET,
       "strict9_original_flower_ROI_photos":STRICT,
       "exploratory_ambiguous9_photo_ROIs":AMBIGUOUS,
       "four_hue_original_classifiable":FOUR_HUES,
       "nine_strict_photo_fraction_of_attempts":finite_pct(STRICT,TARGET),
       "nine_plus_ambiguous_fraction_of_attempts":finite_pct(STRICT+AMBIGUOUS,TARGET),
       "four_hue_fraction_of_attempts":finite_pct(FOUR_HUES,TARGET),
       "classified4_no_strict9_fraction_count":FOUR_HUES-STRICT,
       "no_nine_component_photo_IDs":TARGET-STRICT-AMBIGUOUS,
       "unmeasured_four_hue_photoIDs_without_nine":(
           TARGET-FOUR_HUES-AMBIGUOUS),
       "strict9_source_taxa":int(per_taxon.size),
       "strict9_named_genera":int(within_genus.size),
       "any_processed_photo_taxa":int(d.inat_taxon_id.nunique()),
       "upstream_frame_taxa_with_any_processed_photo":int(taxa.photo_attempts.gt(0).sum()),
       "upstream_frame_taxa_with_strict9_colour":int(taxa.strict9_photo_ROIs.gt(0).sum()),
       "upstream_frame_taxa_NO_original_images_attempted":int(taxa.photo_attempts.eq(0).sum()),
       "upstream_frame_taxa_NO_valid_original_nine_colour":int(taxa.strict9_photo_ROIs.eq(0).sum()),
       "upstream_total_photoIDs_in_taxa_with_no_strict9_colour":int(
             taxa.loc[taxa.strict9_photo_ROIs.eq(0),"n_unique_source_photo_ids"].sum()),
       "source_photo_attempts_OUTSIDE_42111_upstream_taxon_frame":int(
             (~d.upstream_frame_member).sum()),
       "FCP_447500_same_photo_accepted_taxon_ID_disagreements":
             fcp_original_taxon_disagreements,
       "source_strict9_OUTSIDE_upstream_taxon_frame":int(
             (d.strict&~d.upstream_frame_member).sum()),
       "taxa_with_exactly_one_strict9_photo":int(per_taxon.eq(1).sum()),
       "taxa_with_2_to_4_strict9_photos":int(per_taxon.between(2,4).sum()),
       "taxa_with_at_least_5_strict9_photos":int(per_taxon.ge(5).sum()),
       "taxa_with_at_least_10_strict9_photos":int(per_taxon.ge(10).sum()),
       "taxa_with_at_least_30_strict9_photos":int(per_taxon.ge(30).sum()),
       "taxa_with_at_least_100_strict9_photos":int(per_taxon.ge(100).sum()),
       "median_valid_original9_photos_per_taxon":float(per_taxon.median()),
       "p90_valid_original9_photos_per_taxon":float(per_taxon.quantile(.9)),
       "median_attempts_per_taxon_with_attempt":float(
             photo_attempts.median()),
       "effective_taxon_count_inverse_simpson_of_strict_photo_distribution":
             concentration(per_taxon),
       "effective_genus_count_inverse_simpson_of_strict_photo_distribution":
             concentration(within_genus),
       "top10_taxa_pct_of_strict9_photo_ROIs":top_share(per_taxon,10),
       "top100_taxa_pct_of_strict9_photo_ROIs":top_share(per_taxon,100),
       "top1000_taxa_pct_of_strict9_photo_ROIs":top_share(per_taxon,1000),
       "top1000_taxa_pct_of_UPSTREAM_photo_opportunities":top_share(
             taxa.n_unique_source_photo_ids,1000),
       "dates_nonempty_parseable_original_photos":int(d.year.notna().sum()),
       "photo_min_observation_year":str(d.year.min()),
       "photo_max_observation_year":str(d.year.max()),
       "original_observation_ID_populated":int(len(complete_observations)),
       "distinct_source_original_observation_IDs":int(len(observation_reuse)),
       "original_observations_with_multiple_image_photos":int(
             observation_reuse.ge(2).sum()),
       "original_photo_IDs_sharing_source_observation":int(
             observation_reuse[observation_reuse.ge(2)].sum()),
       "source_selection_inclusion_probability_available_on_final_lake":False,
       "source_observer_ID_available_on_final_lake":False,
       "derived_strata":matched,
       "original_source_arms":arm_cohorts,
       "true_source_cohort_separating_FCP_and_nonFCP":photo_cohorts,
       "original_latitude_zones":geographic,
       "original_photo_meteorological_seasons":time,
       "original_photo_hemispheres":hemisphere,
       "island_mainland_composition":island,
       "original_source_by_hemisphere":cohort_hemi,
       "original_source_by_island_mainland":cohort_island,
       "upstream_taxon_opportunity_bins":opportunity,
       "source_photo_observation_years":year_table,
       "nine_palette_and_photo_metadata_joins_conserve_1043930_exact_photo_IDs":True,
       "raw_flower_petal_ground_truth_available":False,
       "green_not_in_original_nine_palette":True,
       "sRGB_lab_only_real_251_photo_pilot_not_global_photo_coverage":True,
       "no_photo_RGB_download_or_colour_remeasurement":True,
       "source_SHA256":{
          "original_nine_source":sha(photo),"original_observation_union":sha(union),
          "FCP447500_photo_ID_source":sha(fcp),
          "complete_93m_upstream_taxon_counts":sha(upstream),
          "highres_GSHHG_landmass":sha(geography)},
    }
    out.mkdir(parents=True)
    bysource=out/"mois_1043930_photo_source_composition.csv"
    pd.DataFrame(photo_cohorts).to_csv(bysource,index=False)
    pd.DataFrame(arm_cohorts).to_csv(
       out/"mois_1043930_original_acquisition_arm_composition.csv",index=False)
    pd.DataFrame(geographic).to_csv(
       out/"mois_1043930_latitude_zone_photo_ROI_coverage.csv",index=False)
    pd.DataFrame(opportunity).to_csv(
       out/"mois_1043930_upstream_taxon_opportunity_coverage.csv",index=False)
    pd.DataFrame(year_table).to_csv(
       out/"mois_1043930_photo_capture_year_coverage.csv",index=False)
    pd.DataFrame(cohort_hemi).to_csv(
       out/"mois_1043930_source_by_hemisphere_stratified_ROI.csv",index=False)
    pd.DataFrame(cohort_island).to_csv(
       out/"mois_1043930_source_by_island_stratified_ROI.csv",index=False)
    pd.DataFrame(island).to_csv(
       out/"mois_1043930_island_mainland_photo_ROI.csv",index=False)
    (out/"mois_1043930_actual_data_structure_and_bias_audit.json").write_text(
        json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("MOIS_REAL1043930_AUDIT_SOURCE_COMPOSITION",json.dumps({
       k:report[k] for k in (
          "original_image_photo_attempts","strict9_original_flower_ROI_photos",
          "strict9_source_taxa","strict9_named_genera","any_processed_photo_taxa",
          "upstream_frame_taxa_NO_valid_original_nine_colour",
          "taxa_with_exactly_one_strict9_photo","taxa_with_at_least_10_strict9_photos",
          "effective_taxon_count_inverse_simpson_of_strict_photo_distribution",
          "top1000_taxa_pct_of_strict9_photo_ROIs",
          "top1000_taxa_pct_of_UPSTREAM_photo_opportunities",
          "original_observations_with_multiple_image_photos",
          "source_photo_IDs_sharing_source_observation" if False else
          "original_photo_IDs_sharing_source_observation",
          "source_photo_attempts_OUTSIDE_42111_upstream_taxon_frame")},
        "source_arms":photo_cohorts,
        "geographic_source":geographic,
        "source_by_hemisphere":cohort_hemi,
        "source_by_island":cohort_island,
        "source_cross_comparison":matched,
        },sort_keys=True,allow_nan=False),flush=True)
    return report

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("photo","receipt","union","union-receipt",
              "fcp","fcp-receipt","upstream","geography","geo-receipt","output"):
        p.add_argument("--"+k,type=Path,required=True)
    z=p.parse_args()
    audit(z.photo,z.receipt,z.union,z.union_receipt,
          z.fcp,z.fcp_receipt,z.upstream,z.geography,z.geo_receipt,z.output)

if __name__=="__main__":main()
