#!/usr/bin/env python3
"""Accept ONLY dated + actual classified bloom photos from FULL original FCP ROI.

This is the exact complete 42,111-taxon / 93,209,351-source-photo pathway,
NOT the old measured 100543 photos and NOT the <=40 source selection.

Full ROI reassembly outputs "newly_remeasured_photos.csv", which includes
the original snapshot observed_on and location and each original actual
photo ROI-v4+EfficientSAM nine-palette measurement. Here the safe
photographed-flower gate is: valid ORIGINAL observed_on day AND original
ROI's "classified_four_state_morph". A Flowers observation-level iNaturalist
annotation is neither required nor treated as a photo-level flowering label.

No photo with unclassified/missing flower ROI contributes a white denominator
or numerator. Original source image-classification failure and missing-date
denominators are retained independently, never called a nonwhite flower.
This requires the FULL sealed, remeasured source, and therefore cannot
produce a global measured result before the full source/ROI work is complete.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import sqlite3
from collections import Counter,defaultdict
from pathlib import Path

try:
    from scripts.acquisition.mois_full_source_dated_flower_candidates_20261016 import (
        source_observation_day,climatic_season,sha,FULL_TAXA,
        FULL_PHOTO_OPPORTUNITIES,SCHEMA as DATED_SCHEMA
    )
except ModuleNotFoundError:
    import sys
    sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
    from scripts.acquisition.mois_full_source_dated_flower_candidates_20261016 import (
        source_observation_day,climatic_season,sha,FULL_TAXA,
        FULL_PHOTO_OPPORTUNITIES,SCHEMA as DATED_SCHEMA
    )

SCHEMA="mois_original_full42k_photo_date_verified_flower_roi_hue_month_v1"
STATUS="COMPLETE_FCP_FULL_ORIGINAL_ROI_SNAPSHOT_DATE_AND_BLOOM_HUE_ACCEPTED"
SOURCE_STATUS="COMPLETE_FRESH_MEASURED_PHOTO_CENSUS"
MORPHS={"white","yellow_orange","red_pink","blue_purple"}
BIO=("white","yellow","orange","red","pink","magenta","purple","blue","bronze")
MODEL_COMMIT="9fae6ccdf684a46026f72ba12e98de2c5c54bf2a"
MEASUREMENT_BLOB="b1d7091db44682559bbe4b48a9f3b8db648501a6"
DETECTOR_SHA="f1aaeec4664fe2c178e5cf2bc1f508977bef3e4aa7b40613026cb8ae3de789d5"
SAM_ENCODER_SHA="84ed466ffcc5c1f8d08409bc34a23bb364ab2c15e402cb12d4335a42be0e0951"
SAM_DECODER_SHA="a62f8fa5ea080447c0689418d69e58f1e83e0b7adf9c142e2bd9bcc8045c0b11"
SOURCE_FIELDS={
    "inat_taxon_id","species","photo_id","observation_id","observer_id",
    "observed_on","latitude","longitude","morph","measurement_status",
    "global_classifiable",*(f"flower_fraction_{b}" for b in BIO)
}


def verify_complete(original_roi_csv:Path,original_receipt:Path,
                    method_attestation:Path,*,synthetic=False)->dict:
    doc=json.loads(original_receipt.read_text(encoding="utf8"))
    if (doc.get("schema")!="fcp_fresh_complete_remeasurement_v1"
        or doc.get("status")!=SOURCE_STATUS
        or doc.get("historical_measurements_reused") is not False
        or doc.get("n_terminal_measurements")!=doc.get("n_frozen_licensed_photo_ids")
        or doc.get("n_classifiable",-1)+doc.get("n_unclassifiable",-1)!=doc.get("n_terminal_measurements")
        or doc.get("new_photo_table_sha256")!=sha(original_roi_csv)
        or doc.get("method_attestation_sha256")!=sha(method_attestation)
        or not isinstance(doc.get("source_firewall_sha256"),str)
        or len(doc["source_firewall_sha256"])!=64):
        raise ValueError("Not a fully sealed, actual original FCP photo image-ROI remeasurement")
    method=json.loads(method_attestation.read_text(encoding="utf8"))
    if not synthetic and (
        doc.get("n_frozen_licensed_photo_ids")!=FULL_PHOTO_OPPORTUNITIES
        or doc.get("n_source_taxa")!=FULL_TAXA
        or doc.get("study_id")!="fcp_fresh_uncapped_20261014_v1"
        or method.get("source_commit")!=MODEL_COMMIT
        or method.get("measurement_git_blob")!=MEASUREMENT_BLOB
        or method.get("detector_sha256")!=DETECTOR_SHA
        or method.get("efficient_sam_encoder_sha256")!=SAM_ENCODER_SHA
        or method.get("efficient_sam_decoder_sha256")!=SAM_DECODER_SHA
        or method.get("location_blind") is not True
        or method.get("no_old_measurement_csv_inputs") is not True
        or method.get("persisted_image_pixels") is not False
    ):
        raise ValueError("Actual ROI-v4 detector/method or full uncapped 42111 photo frame not admitted")
    return doc


def analyze(original_roi_csv:Path,full_roi_receipt:Path,
            method_attestation:Path,private_blooms:Path,public_summary:Path,
            *,synthetic=False,min_public_cell=10)->dict:
    for path in (private_blooms,public_summary):
        if path.exists() or path.with_name(path.name+".INCOMPLETE").exists():
            raise ValueError("Cannot overwrite original full photographed flower calendar output")
    if min_public_cell<10:
        raise ValueError("Public source cell privacy minimum cannot be reduced")
    source=verify_complete(original_roi_csv,full_roi_receipt,method_attestation,
                           synthetic=synthetic)
    private_blooms.parent.mkdir(parents=True,exist_ok=True)
    public_summary.parent.mkdir(parents=True,exist_ok=True)
    part=private_blooms.with_name(private_blooms.name+".INCOMPLETE")
    public_part=public_summary.with_name(public_summary.name+".INCOMPLETE")
    scratch=private_blooms.with_name(private_blooms.name+".sqlite.INCOMPLETE")
    count=Counter()
    group_counts=defaultdict(Counter)
    conn=None
    try:
        conn=sqlite3.connect(scratch)
        conn.execute("PRAGMA journal_mode=DELETE")
        conn.execute("PRAGMA synchronous=NORMAL")
        conn.execute("""CREATE TABLE species_season(
            zone TEXT NOT NULL, month INTEGER NOT NULL, taxon INTEGER NOT NULL,
            photo_n INTEGER NOT NULL, white_n INTEGER NOT NULL,
            PRIMARY KEY(zone,month,taxon)) WITHOUT ROWID""")
        buffered=Counter()
        with original_roi_csv.open(newline="",encoding="utf8") as source_csv, \
             part.open("w",encoding="utf8") as output:
            reader=csv.DictReader(source_csv)
            if not SOURCE_FIELDS.issubset(reader.fieldnames or []):
                raise ValueError("Complete FCP ROI source lacks actual photo ID, date or nine original hue values")
            for r in reader:
                count["source_all_fresh_photo_outcomes"]+=1
                raw_status=r["measurement_status"]
                hue=r["morph"]
                is_flower_hue=raw_status=="classified_four_state_morph"
                if is_flower_hue!=(hue in MORPHS) or (
                    str(r["global_classifiable"]).lower()!=( "true" if is_flower_hue else "false")
                ):
                    raise ValueError("Original actual ROI flower classification outcome changed")
                count["source_classified_floral_photographs_regardless_of_date"]+=int(is_flower_hue)
                day,date_status=source_observation_day(r["observed_on"])
                count[date_status]+=1
                if day is None:continue
                lat=float(r["latitude"]);lon=float(r["longitude"])
                if not (math.isfinite(lat) and math.isfinite(lon)
                    and -90<=lat<=90 and -180<=lon<=180):
                    count["date_valid_but_missing_or_invalid_original_photo_coordinates"]+=1
                    continue
                zone,season=climatic_season(lat,day.month)
                g=group_counts[(zone,day.month)]
                g["all_date_and_geo_valid_source_photos"]+=1
                g["dated_original_four_hue_photo"]+=int(is_flower_hue)
                if not is_flower_hue:
                    g["dated_photo_roi_unclassifiable_no_negative_label"]+=1
                    continue
                if not str(r["photo_id"]).isdecimal() or not str(r["inat_taxon_id"]).isdecimal():
                    raise ValueError("Original flower positive photo ID or taxon changed")
                fractions=[float(r[f"flower_fraction_{b}"]) for b in BIO]
                if (any(not math.isfinite(v) or v<-1e-10 for v in fractions)
                    or abs(sum(fractions)-1)>1e-6):
                    raise ValueError("Real flower-colour ROI nine-palette lacks valid probability simplex")
                original_white=int(hue=="white")
                count["original_snapshot_dated_actual_roi_classified_flower_photos"]+=1
                count["original_snapshot_dated_classified_white_photos"]+=original_white
                g["dated_white_roi_classifiable_photos"]+=original_white
                tid=int(r["inat_taxon_id"])
                buffered[(zone,day.month,tid,original_white)]+=1
                if len(buffered)>=10000:
                    with conn:
                        conn.executemany("""INSERT INTO species_season(zone,month,taxon,photo_n,white_n)
                            VALUES(?,?,?,?,?) ON CONFLICT(zone,month,taxon) DO UPDATE SET
                            photo_n=photo_n+excluded.photo_n,white_n=white_n+excluded.white_n""",
                            [(z,m,t,n,n if w else 0)
                             for (z,m,t,w),n in buffered.items()])
                    buffered.clear()
                photo_row={
                    "original_photo_id":r["photo_id"],
                    "original_observation_uuid":r["observation_id"],
                    "original_taxon_id":tid,"original_observer":r["observer_id"],
                    "original_observed_on":day.isoformat(),"observed_month":day.month,
                    "latitude":lat,"longitude":lon,
                    "latitude_domain":zone,"local_season":season,
                    "original_ROI_v4_hue":hue,
                    "actual_photo_contains_classifiable_flower":True,
                    "original_snapshot_20260927_date_proven":True,
                }
                output.write(json.dumps(photo_row,sort_keys=True)+"\n")
        if buffered:
            with conn:
                conn.executemany("""INSERT INTO species_season(zone,month,taxon,photo_n,white_n)
                     VALUES(?,?,?,?,?) ON CONFLICT(zone,month,taxon) DO UPDATE SET
                     photo_n=photo_n+excluded.photo_n,white_n=white_n+excluded.white_n""",
                     [(z,m,t,n,n if w else 0) for (z,m,t,w),n in buffered.items()])
        if (count["source_all_fresh_photo_outcomes"]!=source["n_terminal_measurements"]
            or count["source_classified_floral_photographs_regardless_of_date"]!=
                source["n_classifiable"]):
            raise ValueError("Original source full photo/ROI denominators do not match FCP complete receipt")
        species_stats={}
        for zone,month,ns,mean_white in conn.execute("""SELECT zone,month,COUNT(*),
                AVG(1.0*white_n/photo_n) FROM species_season
                GROUP BY zone,month ORDER BY zone,month"""):
            species_stats[(zone,month)]=(ns,mean_white)
        public_month=[]
        for key,counts in sorted(group_counts.items()):
            zone,month=key
            base=counts["all_date_and_geo_valid_source_photos"]
            class_n=counts["dated_original_four_hue_photo"]
            white=counts["dated_white_roi_classifiable_photos"]
            if base<min_public_cell:
                public_month.append({
                    "latitude_domain":zone,"observation_month":month,
                    "public_suppressed_below":min_public_cell
                })
                continue
            n_taxa,species_white=species_stats.get(key,(0,None))
            public_month.append({
                "latitude_domain":zone,"observation_month":month,
                "n_original_photos_with_original_date_and_location":base,
                "n_original_photo_ROI_classifiable_bloom":class_n,
                "n_original_roi_white":white,
                "original_photo_ROI_classifiability_given_dated":class_n/base,
                "white_fraction_given_actual_photo_bloom_ROI":(
                    white/class_n if class_n else None),
                "species_equal_white_fraction_given_actual_photo_bloom_ROI":
                    species_white,
                "n_observed_named_taxa_with_classifiable_flower":n_taxa,
                "n_date_valid_unclassifiable_photos_not_assigned_hue":
                    counts["dated_photo_roi_unclassifiable_no_negative_label"],
            })
        report={
            "schema":SCHEMA,"status":STATUS,
            "original_completed_FCP_photo_source_receipt_sha256":sha(full_roi_receipt),
            "original_complete_real_photo_table_sha256":source["new_photo_table_sha256"],
            "original_frozen_method_attestation_sha256":sha(method_attestation),
            "n_original_FCP_all_source_photo_ROI_results":count["source_all_fresh_photo_outcomes"],
            "n_original_ROI_classifiable_flower_hues_all_dates":count[
                "source_classified_floral_photographs_regardless_of_date"],
            "n_original_date_and_site_known_ROI_bloom_hues":count[
                "original_snapshot_dated_actual_roi_classified_flower_photos"],
            "n_original_date_known_white_flower_photos":count[
                "original_snapshot_dated_classified_white_photos"],
            "original_observed_day_source_states":{
                k:v for k,v in count.items() if k.startswith("SOURCE_DATE_") or k=="ORIGINAL_SOURCE_DATE_KNOWN"
            },
            "species_equal_and_photo_weighted_calendar_month_source_white":public_month,
            "private_actual_bloom_photo_month_hue_manifest_sha256":sha(part),
            "full_42111_original_source_not_selectively_date_sampled":True,
            "per_species_source_photo_cap":None,
            "unannotated_original_photo_rejected_as_nonflower_only":False,
            "unclassifiable_photographs_assigned_fake_nonwhite":False,
            "this_is_frozen_snapshot_observed_on_not_current_API":True,
            "WorldClim_monthly_T_S_P_VPD_source_site_link_not_yet_run":True,
            "original_source_photos_hue_measured_anew":True,
            "data_source_is_not_random_wild_blooming_plant_census":True,
            "synthetic_test_only":bool(synthetic),
        }
        public_part.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n",encoding="utf8")
        conn.close();conn=None
        scratch.unlink()
        os.replace(part,private_blooms)
        os.replace(public_part,public_summary)
        return report
    except Exception:
        if conn is not None:conn.close()
        for p in (part,public_part,scratch):
            if p.exists():p.unlink()
        raise


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--real-full-FCP-remeasured-photo-csv",type=Path,required=True)
    p.add_argument("--real-FCP-full-complete-receipt",type=Path,required=True)
    p.add_argument("--original-ROI-method-attestation",type=Path,required=True)
    p.add_argument("--out-private-dated-bloom-photos",type=Path,required=True)
    p.add_argument("--out-public-month-white",type=Path,required=True)
    a=p.parse_args()
    result=analyze(a.real_full_FCP_remeasured_photo_csv,
        a.real_FCP_full_complete_receipt,a.original_ROI_method_attestation,
        a.out_private_dated_bloom_photos,a.out_public_month_white)
    print(json.dumps({
        "status":result["status"],
        "n_source_photo_ROI_results":result["n_original_FCP_all_source_photo_ROI_results"],
        "n_actual_bloom_ROI_hues_with_original_observed_on":result[
            "n_original_date_and_site_known_ROI_bloom_hues"],
        "n_white":result["n_original_date_known_white_flower_photos"],
    },sort_keys=True))


if __name__=="__main__":main()
