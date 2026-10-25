#!/usr/bin/env python3
"""SHA-authenticated 34,898 truly NEW full-resolution nested5 FCP image pixels.

A frozen precolour nested 5/N_taxon/M_photo work source is already globally
photo-ID de-duplicated against 1,009,032 completed original photo ROIs. Only
these untouched source photos may enter the unchanged YOLO11/EfficientSAM ROI-v4
evaluator. No source colour/flower-stage/native inference in the planner.
Unknowns and acquisition errors are kept, not converted to white flowers.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import json
import math
from pathlib import Path
import os
import shutil
import tempfile

from scripts.acquisition.mois_completed_pixel_checkpoint import (
    HUES, collect_wave, require, sha,
)
from scripts.acquisition.mois_fcp_u100_depth5_newonly_queue import (
    ARM, SCHEMA as ORIGINAL_SOURCE_SCHEMA,
)
from scripts.acquisition.mois_queue_from_full93m_target_tiers import (
    OUT_SCHEMA, WORKER_COLUMNS,
)

SCHEMA="mois_FCP_U100_nested5_34898_new_frozen_original_photo_roi_results_v1"
SOURCE_MANIFEST_SHA="89e6c924e172fcc8b862daf30ca88bade2d5f1e2141134950b03928d510e4931"
SOURCE_TAR_SHA="4c1105f4d348a090d46f5c55913867e0e8fc1f16a472aa387094d4c2034b6a9b"
NEW=34_898
WORK_UNITS=137
WORK_GROUPS=35
SOURCE_PHOTOS=101_210
PREVIOUS_ACTUAL=1_009_032
REUSED=66_312
ELIGIBLE_TAXA=20_242
CANDIDATES=80_968


def check_source(source, *, sha_expected=SOURCE_MANIFEST_SHA,
                 count_expected=NEW, units_expected=WORK_UNITS,
                 groups_expected=WORK_GROUPS,
                 prior_expected=PREVIOUS_ACTUAL,
                 reused_expected=REUSED,
                 taxa_expected=ELIGIBLE_TAXA,
                 candidate_expected=CANDIDATES,
                 all_first5_expected=SOURCE_PHOTOS):
    source=Path(source)
    manifest_file=source/"manifest.json"
    require(manifest_file.is_file() and sha(manifest_file)==sha_expected,
            "Original FCP U100 stage5 source manifest SHA differs from frozen release")
    doc=json.loads(manifest_file.read_text(encoding="utf-8"))
    require(
        doc.get("schema")==OUT_SCHEMA
        and doc.get("original_sampling_schema")==ORIGINAL_SOURCE_SCHEMA
        and doc.get("full_frame_photo_ids")==93_209_351
        and doc.get("source_arm")==ARM
        and doc.get("source_group")=="FCP_U100_nested5_rank2to5_only_NOMINAL_ELIGIBLE"
        and doc.get("colour_measured") is False
        and doc.get("new_photo_attempts_planned")==count_expected
        and doc.get("source_original_eligible_taxa")==taxa_expected
        and doc.get("original_incremental_rank2to5_candidate_photo_IDs")==candidate_expected
        and doc.get("original_first5_photo_source_IDs")==all_first5_expected
        and doc.get("previously_measured_historical_K1_FCP_1009032_photo_IDs_accounted")
            ==prior_expected
        and doc.get("previously_measurement_outcomes_reused_without_inference")
            ==reused_expected
        and doc.get("full_old_K1_FCP_reuse_photo_count_conservation")==all_first5_expected
        and doc.get("original_photo_technical_unsupported_new_photos_not_downloaded")==0
        and doc.get("source_photo_selection_probability_unchanged",True) is True
        and doc.get("individual_original_photo_nested_prefix_inclusion_probabilities_preserved") is True
        and doc.get("pending_image_work_is_only_source_photo_IDs_not_previously_inferred") is True
        and doc.get("old_unknown_ROI_not_treated_as_flower_colour") is True
        and doc.get("new_image_downloads_or_ROI_inferences_by_this_source_step")==0
        and doc.get("source_metadata_wild_native_bloom_not_independently_verified") is True
        and doc.get("technical_unknowns_preserved") is True
        and doc.get("original_photo_nested5_source_receipt_SHA256")
            =="037e74dd236d04c44f0d84c7b9d0b3e436e72fe6ddd470cc315c1dc56b0e197f"
        and doc.get("selection_probability") is None
        and doc.get("original_model_batch_size")==256
        and doc.get("units_per_group")==4
        and doc.get("work_units")==units_expected
        and doc.get("worker_ids")==list(range(units_expected))
        and doc.get("group_ids")==list(range(groups_expected)),
        "FCP nested5 original image work population/probability/dedup contract changed")
    counts=doc.get("photo_dispositions")
    require(isinstance(counts,dict)
            and counts.get("prior_original_photo_exact_ID_reuse")==reused_expected
            and counts.get("new_original_photo_ROI_to_infer")==count_expected
            and counts.get("original_photo_unsupported_extension_unknown",0)==0
            and sum(counts.values())==all_first5_expected
            and doc.get("full_old_K1_FCP_reuse_photo_count_conservation")==sum(counts.values()),
            "Source original measured/unknown/new photo partition counts changed")
    photos=set();obs=set();taxa=Counter();per_group=Counter();row_count=0
    manifest_shas=doc.get("photo_shard_sha256")
    require(isinstance(manifest_shas,dict) and
            set(manifest_shas)=={str(i) for i in range(units_expected)},
            "Missing SHA-attested frozen original photo worker shard")
    for unit in range(units_expected):
        path=source/"shards"/f"shard-{unit:03d}.csv"
        require(path.is_file() and sha(path)==manifest_shas[str(unit)],
                "Original FCP source photo shard SHA incorrect")
        with path.open(encoding="utf-8",newline="") as fh:
            reader=csv.DictReader(fh)
            require(tuple(reader.fieldnames or ())==WORKER_COLUMNS,
                    "Unexpected original photo work source column order")
            work=list(reader)
        expected_rows=min(256,count_expected-unit*256)
        require(len(work)==expected_rows and expected_rows>0,
                "An original model work unit is incomplete or overfilled")
        for row in work:
            pid=row["photo_id"]; uuid=row["observation_id"]
            taxon=row["inat_taxon_id"]; photo_url=row["photo_url_large"]
            require(pid.isdecimal() and taxon.isdecimal() and bool(uuid)
                    and pid not in photos and uuid not in obs
                    and row["source_arm"]==ARM and row["shard"]==str(unit)
                    and row["source_snapshot_date"]=="2026-09-27"
                    and row["measurement_id"].startswith("MOIS-")
                    and row["image_filename"]==row["measurement_id"]+".jpg",
                    "New FCP nested5 original photo source is repeated or replaced")
            try: prob=float(row["selection_probability"])
            except (TypeError,ValueError) as e:
                raise ValueError("Missing nested5 original per-photo design probability") from e
            require(math.isfinite(prob) and 0<prob<=.05,
                    "Invalid original 5/N_i/M_o nested photo probability")
            require(photo_url.startswith(
                f"https://inaturalist-open-data.s3.amazonaws.com/photos/{pid}/large.")
                and photo_url.rsplit(".",1)[-1] in ("jpg","jpeg","png"),
                "Signed original photo URL does not point to its real image ID")
            photos.add(pid);obs.add(uuid);taxa[taxon]+=1
            per_group[unit//4]+=1;row_count+=1
    require(len(photos)==len(obs)==row_count==count_expected
            and max(taxa.values())<=4
            and set(per_group)==set(range(groups_expected))
            and all(per_group[g]==min(1024,count_expected-g*1024)
                    for g in range(groups_expected)),
            "Original nested5 source-independent photos/group cardinality changed")
    result={
        "schema":SCHEMA+"_precolour_gate",
        "original_worldwide_photo_candidates_METADATA_ONLY":93_209_351,
        "independent_FCP_U100_eligible_nominal_taxa":taxa_expected,
        "nested_first5_original_photo_IDs":all_first5_expected,
        "rank2to5_original_source_candidates_before_prior_reuse":candidate_expected,
        "originally_measured_photo_IDs_reused":reused_expected,
        "new_distinct_original_photos_for_model":row_count,
        "original_model_work_units":units_expected,
        "original_model_worker_groups":groups_expected,
        "original_photo_design_probability_per_photo_validated":True,
        "source_image_manifest_SHA256":sha_expected,
        "photo_ID_and_observation_UUID_unique":True,
        "per_group_new_photo_counts":{str(g):per_group[g] for g in range(groups_expected)},
        "new_original_image_downloads_or_classifications_in_this_source_gate":0,
        "source_taxa_botanically_native_blooming_independently_verified":False,
    }
    print("MOIS_FCP_U100_NESTED5_TRUE_NEW_SOURCE_ORIGINAL_PHOTOS",json.dumps({
        "new_original_photo_ids":row_count,"units":units_expected,
        "groups":groups_expected,"old_photo_reuses_not_redownloaded":reused_expected,
        "model_calls":0},sort_keys=True),flush=True)
    return result


def collect(source,workers,output, *,
            manifest_sha=SOURCE_MANIFEST_SHA,
            expected_new=NEW, expected_units=WORK_UNITS,
            expected_groups=WORK_GROUPS,
            expected_prior=PREVIOUS_ACTUAL, expected_reuse=REUSED,
            expected_taxa=ELIGIBLE_TAXA, expected_candidates=CANDIDATES,
            expected_first5=SOURCE_PHOTOS):
    import pandas as pd
    from scripts.analysis.mois_reassemble_original_new4096_measurements import METHOD
    source,workers,output=Path(source),Path(workers),Path(output)
    require(not output.exists(),"Do not overwrite immutable measured original image ROIs")
    preflight=check_source(source,sha_expected=manifest_sha,
        count_expected=expected_new,units_expected=expected_units,
        groups_expected=expected_groups,prior_expected=expected_prior,
        reused_expected=expected_reuse,taxa_expected=expected_taxa,
        candidate_expected=expected_candidates,all_first5_expected=expected_first5)
    finished,audit=collect_wave(source,workers,METHOD)
    require(audit["wave_complete"] and audit["attempted"]==expected_new
            and len(finished)==expected_new
            and len(audit["completed_units"])==expected_units
            and audit["pending_photo_attempts"]==0
            and audit["classifiable"]+audit["acquisition_failed"]+
                audit["unclassifiable_acquired"]==expected_new,
            "Not all genuine model-completed nested5 original pixels were audited")
    all_ids={r["photo_id"] for r in finished}
    all_obs={r["observation_id"] for r in finished}
    require(len(all_ids)==len(all_obs)==expected_new
            and all(r["source_arm"]==ARM for r in finished),
            "Genuine source photo identity/observation group changed after ROI")
    f=pd.DataFrame(finished).sort_values(
        "photo_id",key=lambda ids:ids.map(int)).reset_index(drop=True)
    output.parent.mkdir(parents=True,exist_ok=True)
    stage=Path(tempfile.mkdtemp(prefix="mois-U100-5-genuine-ROI-",dir=output.parent))
    try:
        dest=stage/"new_FCP_depth5_original_photo_ROIs.parquet"
        f.to_parquet(dest,compression="zstd",index=False)
        report={
            "schema":SCHEMA,
            "new_original_photo_attempts_actual":expected_new,
            "actual_four_hue_ROIs":audit["classifiable"],
            "original_image_access_failures":audit["acquisition_failed"],
            "original_acquired_but_ROI_unclassifiable":audit["unclassifiable_acquired"],
            "original_source_SHA256":manifest_sha,
            "source_group_probability_and_ID_QA":preflight,
            "source_actual_completed_work_unit_SHA_receipts":audit["completed_unit_audit"],
            "genuine_completed_original_FCP_ROI_method":METHOD,
            "new_original_photo_result_parquet_SHA256":sha(dest),
            "prior_sampled_66312_photo_IDs_reused_outside_this_new_pixel_subset":expected_reuse,
            "source_34898_photo_work_units_audited":expected_units,
            "original_image_pixels_re_downloaded_during_this_collector":0,
            "true_botanical_wild_native_flowering_status_independently_reviewed":False,
            "probability_weighted_global_flora_white_rate_computed":False,
        }
        (stage/"receipt.json").write_text(
            json.dumps(report,indent=2,sort_keys=True,default=str)+"\n",encoding="utf-8")
        os.replace(stage,output)
    except BaseException:
        shutil.rmtree(stage,ignore_errors=True)
        raise
    print("MOIS_FCP_U100_NESTED5_TRUE_ORIGINAL_IMAGE_ROIS_COMPLETED",json.dumps({
        "image_attempts":expected_new,"classified_four_hue":audit["classifiable"],
        "acquisition_failed":audit["acquisition_failed"],
        "acquired_unclassifiable":audit["unclassifiable_acquired"],
        "no_double_fetch_or_reinference_in_collector":True},sort_keys=True),flush=True)
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest="operation",required=True)
    q=sub.add_parser("preflight")
    q.add_argument("--source",type=Path,required=True)
    q.add_argument("--output",type=Path,required=True)
    c=sub.add_parser("collect")
    for key in ("source","workers","output"):
        c.add_argument("--"+key,type=Path,required=True)
    args=parser.parse_args()
    if args.operation=="preflight":
        dest=args.output
        require(not dest.exists(),"Source gate already published")
        receipt=check_source(args.source)
        dest.parent.mkdir(parents=True,exist_ok=True)
        dest.write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    else:
        collect(args.source,args.workers,args.output)


if __name__=="__main__":
    main()
