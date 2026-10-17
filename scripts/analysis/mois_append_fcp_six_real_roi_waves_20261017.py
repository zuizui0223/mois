#!/usr/bin/env python3
"""Append already measured FCP original ROI pixels to MOIS global hue composition.

FCP photos are permitted to overlap. Scientific FCP D/W/local-polymorphism
tests are not imported or recomputed. Actual additions are source-measured,
deduplicated by photo ID, and retain missing-photo hue outcomes.
"""
from __future__ import annotations
import argparse
import csv
import gzip
import hashlib
import io
import json
import zipfile
from collections import Counter
from pathlib import Path

PREVIOUS_PHOTO_COUNT=247032
PREVIOUS_CLASSIFIED=120310
NEW_PHOTO_COUNT=1280
NEW_CLASSIFIED=497
FINAL_PHOTO_COUNT=248312
FINAL_CLASSIFIED=120807
HUES={"white","yellow_orange","red_pink","blue_purple"}
ORIGINAL_UNION_SHA="c344e34be5508964e0de092dfcc7df63997332e2288bda257c3f561e937e10aa"
FIELDS=(
    "photo_id","inat_taxon_id","species","morph","source_original_measured",
    "old_mois_or_new_fcp","measurement_status","latitude","longitude",
    "photo_license","observation_id","observed_on","time_observed_at",
    "wc_elevation_m","soil_pH","soil_SOC","soil_N","soil_clay",
    "soil_available_water_proxy","wc_bio1","wc_bio5","wc_bio12","wc_bio15"
)
SOURCES={
 "rgfca_live_wave1":(
   "original_roi_results/rgfca_original_new_photo_six_taxa_live_API_technical.csv",
   "original_roi_results/rgfca_original_new_photo_six_taxa_live_API_technical_receipt.json",
   "new_terminal_photo_csv_sha256",552,202),
 "rgfca_live_wave2":(
   "original_roi_results/rgfca_original_new_photo_six_taxa_wave2_live_API_technical.csv",
   "original_roi_results/rgfca_original_new_photo_six_taxa_wave2_live_API_technical_receipt.json",
   "new_terminal_photo_csv_sha256",482,203),
 "rgfca_monthly_first64":(
   "measured/first_original_geographically_matched_photo_terminal.csv",
   "measured/first_geo_pair_real_pixel_receipt.json",
   "original_geographic_source_join_terminal_sha256",64,27),
 "rgfca_monthly_second64":(
   "measured/first_original_geographically_matched_photo_terminal.csv",
   "measured/first_geo_pair_real_pixel_receipt.json",
   "original_geographic_source_join_terminal_sha256",64,27),
 "rgfca_original_season63":(
   "measured/season63_measured_postblind.csv",
   "measured/season63_real_roi_receipt.json",
   "files_sha256:season63_measured_postblind.csv",63,26),
 "asclepias_original_2025_64":(
   "postblind/first_original_source_matched_season_64_photo_roi_QC.csv",
   "postblind/first_original_source_matched_season_64_roi_receipt.json",
   "postblind_join_sha256",64,18),
}
EXTRA_SOURCE_LICENSE_FILES={
 "rgfca_monthly_first64":"prepared/first_matched_pair_original_photo_source_acquisition.csv",
 "rgfca_monthly_second64":"new_source/prepared/first_matched_pair_original_photo_source_acquisition.csv",
 "rgfca_original_season63":"prepared/season63_acquisition.csv",
}
EXPECTED_ADDED_BY_SOURCE={
 "rgfca_live_wave1":(548,199),
 "rgfca_live_wave2":(477,200),
 "rgfca_monthly_first64":(64,27),
 "rgfca_monthly_second64":(64,27),
 "rgfca_original_season63":(63,26),
 "asclepias_original_2025_64":(64,18),
}


def checksum(data:bytes)->str:
    return hashlib.sha256(data).hexdigest()


def sha_file(p:Path)->str:
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(1<<20),b""):h.update(b)
    return h.hexdigest()


def require(ok:bool,msg:str):
    if not ok:raise ValueError(msg)


def read_source(path:Path,kind:str)->tuple[list[dict],dict]:
    columns=SOURCES[kind]
    with zipfile.ZipFile(path) as z:
        term,receipt,sha_key,n,accepted=columns
        content=z.read(term)
        source=json.loads(z.read(receipt))
        expected=(source.get("files_sha256",{}).get("season63_measured_postblind.csv")
                  if sha_key.startswith("files_sha256:")
                  else source.get(sha_key))
        require(checksum(content)==expected,
                f"FCP measured original image terminal SHA mismatch in {kind}")
        records=list(csv.DictReader(io.StringIO(content.decode("utf-8-sig"))))
        require(len(records)==n and
                sum(r["morph"] in HUES for r in records)==accepted,
                f"FCP measured original hue census changed in {kind}")
        require(len({r["photo_id"] for r in records})==n,
                f"Duplicate image ids within ROI measurement {kind}")
        licenses={}
        if kind in EXTRA_SOURCE_LICENSE_FILES:
            with io.StringIO(z.read(EXTRA_SOURCE_LICENSE_FILES[kind]).decode("utf-8")) as fin:
                rdr=csv.DictReader(fin)
                for row in rdr:
                    measurement=row["measurement_id"]
                    require(measurement not in licenses,
                            f"Duplicated ROI source acquisition in {kind}")
                    licenses[measurement]=row.get("photo_license","")
        for r in records:
            r["photo_license"]=r.get("photo_license") or licenses.get(r.get("measurement_id",""),"")
            require(r["photo_id"].isdigit() and r["inat_taxon_id"].isdigit()
                    and r["species"] and r["measurement_status"],
                    f"Invalid signed ROI photo identity {kind}")
            require(r["morph"] in HUES|{"mixed_uncertain"} and
                    (r["measurement_status"]=="classified_four_state_morph")==
                    (r["morph"] in HUES),f"Invalid classifiability in {kind}")
        return records,{
          "original_roi_rows":n,"original_four_hue_classified":accepted,
          "verified_original_roi_csv_sha256":checksum(content),
          "source_original_receipt_sha256":checksum(z.read(receipt)),
          "selected_from_live_API":kind.startswith("rgfca_live_"),
          "photo_source_snapshot_not_uniform_across_all_six":True,
        }


def normalize(r:dict,kind:str)->dict:
    o={k:"" for k in FIELDS}
    for k in FIELDS:
        if k in r and r[k] is not None:
            o[k]=str(r[k]).strip()
    o["source_original_measured"]=kind
    o["old_mois_or_new_fcp"]="fcp_additional_ROI"
    return o


def run(prior:Path, source_archives:dict[str,Path],outdir:Path):
    require(not outdir.exists() or not any(outdir.iterdir()),"Never overwrite source results")
    require(sha_file(prior)==ORIGINAL_UNION_SHA,"Original 247032-photo measured MOIS source changed")
    if set(source_archives)!=set(SOURCES):
        raise ValueError("All six originally measured ROI sources required")
    source_data={}
    for kind,path in source_archives.items():
        source_data[kind]=read_source(path,kind)
    outdir.mkdir(parents=True,exist_ok=True)
    union=outdir/"mois_original_photo_hues_248312_after_FCP_ROI_additions.csv.gz"
    seen={}
    count=classified=0
    conflicts=[]
    by_source={}
    with union.open("wb") as fo:
        with gzip.GzipFile(filename="",fileobj=fo,mode="wb",mtime=0) as gz:
            with io.TextIOWrapper(gz,encoding="utf-8",newline="") as fout:
                writer=csv.DictWriter(fout,fieldnames=FIELDS,lineterminator="\n")
                writer.writeheader()
                with gzip.open(prior,"rt",newline="",encoding="utf-8") as f:
                    reader=csv.DictReader(f)
                    require(tuple(reader.fieldnames or ())==FIELDS,
                            "Original 247k photo schema differs")
                    for r in reader:
                        pid=r["photo_id"]
                        require(pid not in seen,"Original prior union duplicate image")
                        seen[pid]=(r["morph"],r["inat_taxon_id"])
                        writer.writerow(r)
                        count+=1
                        classified+=r["morph"] in HUES
                require((count,classified)==(PREVIOUS_PHOTO_COUNT,PREVIOUS_CLASSIFIED),
                        "Frozen prior 247032 source-photo denominator changed")
                for kind in SOURCES:
                    rows,receipt=source_data[kind]
                    sub=Counter()
                    for r in rows:
                        pid=r["photo_id"]
                        sub["attempted"]+=1
                        if pid in seen:
                            sub["previously_used_photo_ids"]+=1
                            old_morph,old_tax=seen[pid]
                            if old_morph in HUES and r["morph"] in HUES and old_morph!=r["morph"]:
                                conflicts.append({"source":kind,"photo_id":pid,
                                    "prior_photo_hue":old_morph,"new_photo_hue":r["morph"]})
                            continue
                        writer.writerow(normalize(r,kind))
                        seen[pid]=(r["morph"],r["inat_taxon_id"])
                        sub["new_photo_ids"]+=1
                        sub["new_classified_hues"]+=r["morph"] in HUES
                        count+=1
                        classified+=r["morph"] in HUES
                    require((sub["new_photo_ids"],sub["new_classified_hues"])==
                            EXPECTED_ADDED_BY_SOURCE[kind],
                            f"Real original source-overlap ROI extension drift: {kind}")
                    by_source[kind]={**receipt,**dict(sub)}
    require(not conflicts,"Two original measured sources disagree on classified source image labels")
    require((count,classified)==(FINAL_PHOTO_COUNT,FINAL_CLASSIFIED),
            "Original 248312 total photos or 120807 classified did not materialize")
    result={
     "schema":"mois_global_measured_hue_data_248312_source_union_v1",
     "status":"REAL_248312_ORIGINAL_PHOTOS_WITH_120807_CLASSIFIED_FOUR_HUES",
     "first_parent_original_photo_ids":PREVIOUS_PHOTO_COUNT,
     "first_parent_classified_original_photo_ids":PREVIOUS_CLASSIFIED,
     "FCP_already_measured_original_additional_photo_ids":NEW_PHOTO_COUNT,
     "FCP_already_measured_additional_classified_hue_ids":NEW_CLASSIFIED,
     "final_photo_ids":count,"final_classified_hue_ids":classified,
     "final_unclassified_photo_ids":count-classified,
     "all_six_FCP_original_ROI_receipts_verified":True,
     "FCP_photo_sources_can_overlap":True,
     "FCP_within_species_D_W_or_local_diversity_estimand_reused":False,
     "original_historic_mois_macro_coefficients_overwritten":False,
     "genuinely_new_original_photo_pixels_opened_by_this_script":0,
     "source_identifier_overlaps":by_source,
     "original_247032_source_sha256":ORIGINAL_UNION_SHA,
     "result_source_photo_csv_sha256":sha_file(union),
     "source_data_is_not_independent_biological_confirmation":True,
     "photographic_colour_is_not_genotype_or_verified_flower_stage":True,
     "distinct_original_measurement_sources":"Historical MOIS, FCP three high-depth cohorts, FCP six technical ROI-pilot groups",
     "not_all_FCP_live_API_photo_identifiers_come_from_the_20260927_monthly_source":True,
     "next_environmental_model_requires_source_standardization":True,
    }
    (outdir/"mois_248312_frozen_measured_photo_addition_receipt.json").write_text(
        json.dumps(result,indent=2,sort_keys=True)+"\n")
    print("MOIS_248312_ACTUAL_MEASURED_PHOTOS",json.dumps({
        k:result[k] for k in ("final_photo_ids","final_classified_hue_ids",
             "FCP_already_measured_original_additional_photo_ids",
             "FCP_already_measured_additional_classified_hue_ids")},sort_keys=True))
    return result


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--prior-247032-original-photo-union",required=True,type=Path)
    for name in SOURCES:
        p.add_argument("--"+name,required=True,type=Path)
    p.add_argument("--outdir",required=True,type=Path)
    a=p.parse_args()
    run(a.prior_247032_original_photo_union,
        {name:getattr(a,name) for name in SOURCES},a.outdir)
