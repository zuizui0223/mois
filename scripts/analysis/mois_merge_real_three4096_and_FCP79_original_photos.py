#!/usr/bin/env python3
"""MOIS: add finished original second/third 4096 ROI waves, preserve date sources.

Inputs are three COMPLETED source-verified GitHub Actions artifacts. This script
does NOT download or reclassify any photo pixels, or invent hue values for
technically unclassifiable original photographs. Photo-ID union is definitive;
overlapping photo labels, taxon IDs and original observation IDs must agree.
The 100543 historic S3/API dates are a separate evidence column, not a rewrite
of the original observed_on field. FCP D/W within-species inference excluded.
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

HUES=("white","yellow_orange","red_pink","blue_purple")
FCP_BASE_N=286826
FCP_BASE_HUES=139405
MOIS_THREE_N=260600
MOIS_THREE_HUES=125517
MOIS_THREE_NEW_N=8191
MOIS_THREE_NEW_HUES=3133
UNION_N=295017
UNION_HUES=142538
DATE_SOURCE_N=100543
OLD_DATE_HUES=137114
NEW_DATE_HUES=140247
NEW_DATE_ATTEMPTS=289170
FCP_BASE_SHA="a6f4accd600ef793d803c7796b7a3d5345a8a4708af1ea6c52304310c1d4170d"
MOIS_THREE_SHA="8da4d41e2f0bb0a4992394171fe5bad6dd62ee2328136772fe7c5659b3a91b20"
HISTORIC_DATE_SOURCE_SHA="ac5056e25cb9322618dda3b4f2e4e2a97907dfb16749d84c65ee9f3149143564"
FCP_SOURCE="mois_fcp_wave4_PARTIAL_completed_pixel_original_photo_union.csv.gz"
FCP_RECEIPT="mois_fcp_wave4_PARTIAL_real_completed_original_ROI_receipt.json"
MOIS_SOURCE="mois_all_real_original_hue_measured_248312_plus_new_waves.csv.gz"
MOIS_RECEIPT="mois_real_additional_photo_waves_receipt.json"
HISTORY_SOURCE="mois_original100543_s3_or_verified_API_dates_photo_source.csv.gz"
HISTORY_RECEIPT="mois_real127504_original_photo_S3_plus_verified_API_local_season_result.json"
UNION_FILE="mois_real295017_original_photo_pixels_three_waves_FCP180.csv.gz"
DATE_FILE="mois_real295017_original_photo_hues_plus_verified_observed_dates.csv.gz"
SEASON_FILE="mois_real295017_dated_four_hue_season_photo_rates.csv"
SOURCE_RECEIPT="mois_real295017_signed_photo_union_date_and_colour_receipt.json"
OUTCOLS=("photo_id","inat_taxon_id","species","morph","source_original_measured",
         "old_mois_or_new_fcp","measurement_status","latitude","longitude",
         "photo_license","observation_id","observed_on","time_observed_at",
         "wc_elevation_m","soil_pH","soil_SOC","soil_N","soil_clay",
         "soil_available_water_proxy","wc_bio1","wc_bio5","wc_bio12","wc_bio15")


def require(ok:bool,reason:str):
    if not ok:raise ValueError(reason)


def sha(data:bytes)->str:
    return hashlib.sha256(data).hexdigest()


def load_zip(zpath:Path,data_name:str,receipt_name:str,expected_sha:str):
    with zipfile.ZipFile(zpath) as z:
        raw=z.read(data_name)
        receipt=json.loads(z.read(receipt_name))
    require(sha(raw)==expected_sha,"Original measured photo/date data SHA drift: "+data_name)
    d=pd.read_csv(io.BytesIO(raw),compression="gzip",dtype={
        "photo_id":"string","observation_id":"string","inat_taxon_id":"string"},
        low_memory=False)
    require(d.photo_id.is_unique,"Duplicate original photo-ID within source: "+data_name)
    return d,receipt


def photo_union(fcp:pd.DataFrame,mois:pd.DataFrame)->tuple[pd.DataFrame,dict]:
    require(tuple(fcp.columns)==OUTCOLS and tuple(mois.columns)==OUTCOLS,
            "Original measured photograph schemas differ")
    require(len(fcp)==FCP_BASE_N and len(mois)==MOIS_THREE_N and
            fcp.morph.isin(HUES).sum()==FCP_BASE_HUES and
            mois.morph.isin(HUES).sum()==MOIS_THREE_HUES,
            "Signed actual pixels and photograph hue outcomes changed")
    overlap=fcp.merge(mois,on="photo_id",how="inner",suffixes=("_fcp","_mois"),
                      validate="one_to_one")
    require(len(overlap)==252409,
            "The 252408 core plus one shared real pixel must have identical photo IDs")
    checks={}
    # FCP current original observation_id is numeric iNaturalist, while the
    # independent MOIS third-wave source sometimes retains an observation UUID.
    # One shared photo with exactly matching taxon, date and original ROI hue
    # may thus have a source-namespace discrepancy; never double-count it.
    for key in ("inat_taxon_id","morph","observed_on"):
        same=overlap[f"{key}_fcp"].fillna("").astype(str).eq(
             overlap[f"{key}_mois"].fillna("").astype(str))
        checks[key+"_source_conflicts"]=int((~same).sum())
        require(same.all(),"Shared original source photos changed "+key)
    obs_match=overlap.observation_id_fcp.fillna("").astype(str).eq(
        overlap.observation_id_mois.fillna("").astype(str))
    conflicts=overlap.loc[~obs_match]
    checks["observation_ID_source_namespace_mismatches"]=len(conflicts)
    require(len(conflicts)==1 and conflicts.photo_id.iloc[0]=="232092182"
            and conflicts.observation_id_fcp.iloc[0].isdigit()
            and "-" in conflicts.observation_id_mois.iloc[0],
            "More than one original observation-ID namespace discrepancy")
    checks["duplicate_photo_IDs_placed_in_result_twice"]=0
    new=mois.loc[~mois.photo_id.isin(fcp.photo_id)].copy()
    require(len(new)==MOIS_THREE_NEW_N and
            new.morph.isin(HUES).sum()==MOIS_THREE_NEW_HUES and
            new.observed_on.notna().all(),
            "MOIS 2 new finished 4096 image-pixel waves not fully authentic")
    # The completed MOIS 8192 original image set overlaps one newer FCP pixel
    # result; exact photo-ID exclusion leaves 8191 genuinely additional photos.
    require(new.old_mois_or_new_fcp.eq("mois_new_original_ROI").all(),
            "Unrelated original FCP wave4 pixels are not this new MOIS photo wave")
    union=pd.concat((fcp,new),ignore_index=True)
    require(len(union)==UNION_N and union.photo_id.is_unique and
            union.morph.isin(HUES).sum()==UNION_HUES,
            "Original 286826 + MOIS 2 x 4096 photo-ID union failed")
    return union,{"source_shared_photo_IDs":len(overlap),
                  "source_shared_conflicts":checks,
                  "new_original_photo_pixels":len(new),
                  "original_photo_source_ID_overlap_count":len(overlap),
                  "new_original_four_hue_photo_labels":int(new.morph.isin(HUES).sum()),
                  "new_original_source_observed_on_dates":int(new.observed_on.notna().sum())}


def restore_historic_dates(union:pd.DataFrame,hist:pd.DataFrame)->pd.DataFrame:
    require(len(hist)==DATE_SOURCE_N and
            hist.photo_id.is_unique and
            hist.qualified_observed_on.notna().sum()==99817,
            "Historic true observed-date photo-ID repair incomplete")
    required=("photo_id","original_numeric_observation_id","original_taxon_id",
              "original_frozen_ROI_morph","qualified_observed_on")
    require(all(x in hist for x in required),"Historical verified source columns changed")
    allow=hist[list(required)].copy()
    d=union.merge(allow,on="photo_id",how="left",validate="one_to_one")
    is_history=d.old_mois_or_new_fcp.eq("historical")
    require(is_history.sum()==DATE_SOURCE_N and
            d.loc[is_history,"original_numeric_observation_id"].notna().all(),
            "Actual originally measured historic photographs missing from new union")
    for left,right in (("observation_id","original_numeric_observation_id"),
                       ("inat_taxon_id","original_taxon_id"),
                       ("morph","original_frozen_ROI_morph")):
        require(d.loc[is_history,left].fillna("").astype(str).eq(
            d.loc[is_history,right].fillna("").astype(str)).all(),
            "Historical original photo/date/hue/taxon identity conflict: "+left)
    require(not (d.observed_on.notna()&
                 d.qualified_observed_on.notna()).any(),
            "Do not overwrite a genuine original observed_on with public metadata")
    d["qualified_date"] = d.observed_on.fillna(d.qualified_observed_on)
    d["original_source_date_provenance"] = np.select(
        [d.observed_on.notna(),d.qualified_observed_on.notna()],
        ["original_source_observed_on","historic_snapshot_verified_S3_or_original_API"],
        default="unknown_no_observed_on")
    require(d.qualified_date.notna().sum()==NEW_DATE_ATTEMPTS and
            (d.qualified_date.notna()&d.morph.isin(HUES)).sum()==NEW_DATE_HUES,
            "Original source + public historic observation dates/hues do not reconcile")
    result=d[list(OUTCOLS)+["qualified_date","original_source_date_provenance"]].copy()
    return result


def season_table(d:pd.DataFrame):
    date=pd.to_datetime(d.qualified_date,errors="coerce")
    latitude=pd.to_numeric(d.latitude,errors="coerce")
    longitude=pd.to_numeric(d.longitude,errors="coerce")
    classified=d.morph.isin(HUES)
    white=d.morph.eq("white")
    region=np.select([latitude.ge(23.5),latitude.le(-23.5)],
                     ["north_extratropical","south_extratropical"],
                     default="tropical_or_unknown")
    month=date.dt.month
    local_season=pd.Series(((month % 12)//3+np.where(latitude.lt(0),2,0))%4,
                           index=d.index).map({0:"winter",1:"spring",2:"summer",3:"autumn"})
    local_season.loc[region=="tropical_or_unknown"]="tropical_month_only"
    source=d.copy()
    source["date"]=date;source["zone"]=region;source["season"]=local_season
    source["classified"]=classified;source["white"]=white
    source=source[date.notna()&latitude.between(-90,90)&longitude.between(-180,180)]
    t=source.groupby(["zone","season"],observed=True).agg(
        dated_original_photo_attempts=("photo_id","size"),
        original_four_hue_classifiable=("classified","sum"),
        original_white_photo_labels=("white","sum"),
        named_original_photo_taxa=("inat_taxon_id","nunique")).reset_index()
    t["white_fraction_given_classified"]=(
        t.original_white_photo_labels/t.original_four_hue_classifiable).replace(
            [np.inf,-np.inf],np.nan)
    t["source_photo_hue_classifiability"]=t.original_four_hue_classifiable/t.dated_original_photo_attempts
    return t


def run(fcp_zip:Path,mois_zip:Path,historic_zip:Path,outdir:Path):
    require(not outdir.exists() or not any(outdir.iterdir()),
            "Never overwrite source-locked photographic source result")
    a,receipt_a=load_zip(fcp_zip,FCP_SOURCE,FCP_RECEIPT,FCP_BASE_SHA)
    b,receipt_b=load_zip(mois_zip,MOIS_SOURCE,MOIS_RECEIPT,MOIS_THREE_SHA)
    h,receipt_h=load_zip(historic_zip,HISTORY_SOURCE,HISTORY_RECEIPT,HISTORIC_DATE_SOURCE_SHA)
    require(receipt_a["combined_unique_photo_IDs"]==FCP_BASE_N and
            receipt_a["combined_four_hue_classified_photo_IDs"]==FCP_BASE_HUES and
            receipt_a["FCP_completed_ROI_partitions"]==180,
            "FCP 180 authentic completed photo partitions changed")
    require(receipt_b["all_original_unique_photos"]==MOIS_THREE_N and
            receipt_b["all_four_hue_classified_photos"]==MOIS_THREE_HUES and
            receipt_b["new_original_source_photos"]==12288,
            "MOIS three 4096 original photo worker waves were not measured")
    require(receipt_h["full_original_measured_photo_ids"]==267425 and
            receipt_h["full_original_genuine_four_hue_classified_photo_ids"]==129795 and
            receipt_h["full_historical_join"]["history_original_dates_total_with_both_source_precedence"]==99817,
            "Validated historical original dates not reproduced")
    union,stats=photo_union(a,b)
    extra=restore_historic_dates(union,h)
    seasons=season_table(extra)
    outdir.mkdir(parents=True,exist_ok=True)
    src=outdir/UNION_FILE
    union.to_csv(src,index=False,compression={"method":"gzip","mtime":0},
                 lineterminator="\n")
    derived=outdir/DATE_FILE
    extra.to_csv(derived,index=False,compression={"method":"gzip","mtime":0},
                 lineterminator="\n")
    seas=outdir/SEASON_FILE
    seasons.to_csv(seas,index=False,lineterminator="\n",float_format="%.8f")
    def file_hash(p):return hashlib.sha256(p.read_bytes()).hexdigest()
    report={
      "schema":"mois_295017_genuine_all_original_photos_three4096_waves_FCP180_date_enriched_v1",
      "original_fcp79_pixel_source_sha256":FCP_BASE_SHA,
      "original_mois_three_complete_waves_pixel_source_sha256":MOIS_THREE_SHA,
      "historical100543_independent_qualified_date_source_sha256":HISTORIC_DATE_SOURCE_SHA,
      "original_measured_photographs":len(union),
      "original_measured_photo_four_hue_labels":int(union.morph.isin(HUES).sum()),
      "photo_ids_with_verified_observed_dates":int(extra.qualified_date.notna().sum()),
      "photo_ids_with_both_date_and_automatic_four_hue":int(
          (extra.qualified_date.notna()&extra.morph.isin(HUES)).sum()),
      "date_bearing_original_photo_attempts_without_classifiable_hue":int(
          (extra.qualified_date.notna()&~extra.morph.isin(HUES)).sum()),
      "old_mois_267425_photo_date_hues":OLD_DATE_HUES,
      "photographed_named_taxa_with_date_and_hue":int(
          extra.loc[extra.qualified_date.notna()&
                    extra.morph.isin(HUES),"inat_taxon_id"].nunique()),
      **stats,
      "files_SHA256":{src.name:file_hash(src),derived.name:file_hash(derived),
                      seas.name:file_hash(seas)},
      "historic_observation_dates_are_public_metadata_not_image_EXIF":True,
      "original_photo_hues_unclassified_remain_unknown_not_nonflowering":True,
      "no_photo_pixels_reprocessed_or_source_taxa_reassigned":True,
      "full_FCP_49800_ROI_partitions_not_automatically_complete":True,
      "source_cohort_not_independent_provider_replication":True,
      "FCP_within_species_D_W_50km_inference_not_imported":True,
      "no_causal_monthly_solar_water_pigment_fitness_claim":True,
    }
    target=outdir/SOURCE_RECEIPT
    target.write_text(json.dumps(report,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    print("MOIS_REAL_295017_ALL_ACTUAL_PHOTOS",json.dumps({
        k:report[k] for k in ("original_measured_photographs",
            "original_measured_photo_four_hue_labels",
            "photo_ids_with_verified_observed_dates",
            "photo_ids_with_both_date_and_automatic_four_hue",
            "new_original_photo_pixels","new_original_four_hue_photo_labels")
        },sort_keys=True),flush=True)
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fcp79-real-original-zip",required=True,type=Path)
    parser.add_argument("--mois3-real-original-zip",required=True,type=Path)
    parser.add_argument("--historic-verified-public-dates-zip",required=True,type=Path)
    parser.add_argument("--output",required=True,type=Path)
    args=parser.parse_args()
    run(args.fcp79_real_original_zip,args.mois3_real_original_zip,
        args.historic_verified_public_dates_zip,args.output)

if __name__=="__main__":main()
