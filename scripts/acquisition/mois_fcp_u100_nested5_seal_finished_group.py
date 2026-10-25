#!/usr/bin/env python3
"""Seal exactly one completed nested5 1..1024-photo original ROI model group.

Consume the already-running genuine frozen FCP YOLO11/EfficientSAM ROI-v4
workers without re-downloading any image or refitting the model. Verify all
source, per-photo terminal CSV, image SHA, acquisition failure, individual
nested prefix5 sample probabilities and frozen model setup. Keep missing
outcomes unknown and do NOT claim globally representative colour prevalence.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import json
import math
from pathlib import Path
import os
import shutil
import tempfile

from scripts.acquisition.mois_completed_pixel_checkpoint import HUES, collect_wave, require, sha
from scripts.acquisition.mois_fcp_u100_depth5_genuine_image_wave import (
    check_source, SOURCE_MANIFEST_SHA, NEW, WORK_UNITS, WORK_GROUPS,
    PREVIOUS_ACTUAL, REUSED, ELIGIBLE_TAXA, CANDIDATES, SOURCE_PHOTOS,
)
from scripts.acquisition.mois_fcp_u100_depth5_newonly_queue import ARM
from scripts.acquisition.mois_warm_original_roi import FROZEN_SOURCE_BLOB
from scripts.acquisition.mois_append_only_second8_original_worker_group import (
    original_photo_source_coverage,
)

SCHEMA="mois_FCP_U100_nested5_new34898_genuine_original_photo_completed_group_v1"
MODEL_RUN=40151706880


def seal(source, workers, group, output, *,
         expected_sha=SOURCE_MANIFEST_SHA, model_run=MODEL_RUN,
         planned=NEW, units=WORK_UNITS, groups=WORK_GROUPS,
         old_count=PREVIOUS_ACTUAL, reused=REUSED,
         eligible=ELIGIBLE_TAXA, candidates=CANDIDATES,
         all_first5=SOURCE_PHOTOS):
    import pandas as pd
    from scripts.analysis.mois_reassemble_original_new4096_measurements import METHOD

    source,workers,output=map(Path,(source,workers,output))
    require(not output.exists(),"Never replace immutable original-photo worker group")
    # Independently check all frozen original photo work source packets and
    # original per-photo selection probabilities before reading any colour.
    src=check_source(source,sha_expected=expected_sha,
        count_expected=planned,units_expected=units,groups_expected=groups,
        prior_expected=old_count,reused_expected=reused,taxa_expected=eligible,
        candidate_expected=candidates,all_first5_expected=all_first5)
    require(type(group) is int and group in range(groups),
            "Completed original FCP nested group does not exist in source")
    group_units=list(range(group*4,min((group+1)*4,units)))
    expected_photos=min(1024,planned-group*1024)
    require(expected_photos>0,
            "Source image group has no genuinely new photo IDs")
    actual,audit=collect_wave(source,workers,METHOD)
    require(audit["completed_units"]==group_units
            and audit["attempted"]==expected_photos==len(actual)
            and audit["planned_photo_attempts"]==planned
            and set(audit["missing_units"])==set(range(units))-set(group_units)
            and len(audit["completed_unit_audit"])==len(group_units),
            "Incomplete or duplicate frozen new FCP nested source work units")
    photographs=set();observations=set();hues=Counter()
    probabilities=set()
    for row in actual:
        pid=str(row["photo_id"])
        obs=str(row["observation_id"])
        state=row["morph"]
        status=row["measurement_status"]
        classified=state in HUES
        failed=status=="image_acquisition_failed"
        require(pid.isdecimal() and pid not in photographs
                and obs and obs not in observations
                and row["source_arm"]==ARM,
                "Original photo ID, observation or source colour arm repeated")
        require(state in (*HUES,"mixed_uncertain")
                and classified==(status=="classified_four_state_morph")
                and row["classified_four_hue"].lower()==str(classified).lower()
                and (failed or len(row.get("image_sha256",""))==64),
                "Real FCP original image technical unknown relabelled as flower colour")
        try: p=float(row["selection_probability"])
        except (ValueError,TypeError) as exc:
            raise ValueError("Missing original nested 5/N/M photo probability") from exc
        require(math.isfinite(p) and 0<p<=0.05,
                "Wrong FCP nested prefix5 per-photo source propensity")
        photographs.add(pid);observations.add(obs)
        hues[state]+=1;probabilities.add(p)
    classified=sum(hues[k] for k in HUES)
    require(classified==audit["classifiable"]
            and audit["acquisition_failed"]+audit["unclassifiable_acquired"]
                ==hues["mixed_uncertain"]
            and classified+audit["acquisition_failed"]
                +audit["unclassifiable_acquired"]==expected_photos,
            "Signed original-image ROI photo colour and unknown counts changed")
    diagpaths=list(workers.rglob(f"warm_roi_group_{group:03d}_diagnostic.json"))
    require(len(diagpaths)==1,
            "Missing or duplicated genuine warm YOLO/EfficientSAM group timing")
    diag=json.loads(diagpaths[0].read_text(encoding="utf-8"))
    pref=diag.get("photo_prefetch",{})
    require(diag.get("group")==group
            and diag.get("frozen_source_blob")==FROZEN_SOURCE_BLOB
            and diag.get("frozen_yolo_and_onnx_initializations")==1
            and diag.get("separate_original_main_calls")==len(group_units)
            and diag.get("new_photos_downloaded")==expected_photos
            and diag.get("method_or_thresholds_modified") is False
            and diag.get("old_photo_RGB_workers_modified") is False
            and diag.get("artificial_flower_labels_created") is False
            and pref.get("download_calls")==expected_photos
            and pref.get("consumed")==expected_photos
            and pref.get("permanent_failures")==audit["acquisition_failed"]
            and isinstance(diag.get("group_wall_seconds"),(int,float))
            and math.isfinite(diag["group_wall_seconds"])
            and diag["group_wall_seconds"]>0,
            "Completed genuine original ROI worker/provenance timing changed")
    output.parent.mkdir(parents=True,exist_ok=True)
    staged=Path(tempfile.mkdtemp(prefix="mois-FCP5-real-photo-group-",dir=output.parent))
    try:
        table=pd.DataFrame(actual)
        table["photo_id"]=table["photo_id"].astype(str)
        table=table.sort_values("photo_id",key=lambda z:z.map(int)).reset_index(drop=True)
        data=staged/"original_photo_ROIs.parquet"
        table.to_parquet(data,index=False,compression="zstd")
        coverage=original_photo_source_coverage(actual)
        coverage["schema"]="mois_FCP_U100_nested5_original_photo_taxon_date_lat_source_coverage_v1"
        coverage["botanically_verified_open_native_wild_flower_individuals"]=False
        coverage["original_photo_nested5_selection_probability_varies_per_photo"]=True
        coverage_file=staged/"original_photo_source_coverage.json"
        coverage_file.write_text(json.dumps(coverage,sort_keys=True,indent=2)+"\n")
        huefile=staged/"original_photo_hue_counts.csv"
        with huefile.open("w",newline="",encoding="utf-8") as stream:
            writer=csv.writer(stream,lineterminator="\n")
            writer.writerow(("morph","original_photo_attempts"))
            for name in (*HUES,"mixed_uncertain"):
                writer.writerow((name,hues[name]))
        receipt={
            "schema":SCHEMA,
            "frozen_original_model_run":model_run,
            "original_upstream_photo_IDs_METADATA_ONLY":93209351,
            "source_SHA256":expected_sha,
            "source_original_FCP_U100_first5_nominal_taxa":eligible,
            "source_101210_first5_original_photos":all_first5,
            "source_80968_new_nested_rank2to5_photo_candidates":candidates,
            "prior_complete_original_photo_IDs_reused_without_inference":reused,
            "prior_global_actual_photo_IDs_SHA_audited_before_new_selection":old_count,
            "all_source_genuinely_new_original_photo_attempts":planned,
            "original_worker_group":group,
            "original_completed_source_work_units":group_units,
            "original_terminal_source_receipts_SHA256":audit["completed_unit_audit"],
            "new_actual_original_photo_ID_attempts":expected_photos,
            "new_actual_four_hue_classified":classified,
            "new_original_image_access_failed":audit["acquisition_failed"],
            "new_original_acquired_hue_unclassifiable":audit["unclassifiable_acquired"],
            "new_distinct_original_observations":len(observations),
            "new_original_hue_counts":{k:hues[k] for k in (*HUES,"mixed_uncertain")},
            "distinct_original_per_photo_joint_probabilities":len(probabilities),
            "frozen_original_FCP_yolo_efficientSAM_roi_v4_method_attested":True,
            "per_original_photo_source_inclusion_probability_preserved":True,
            "old_photo_ID_reuses_not_relabelled":True,
            "unknown_technical_ROIs_never_cast_to_white":True,
            "original_image_bytes_downloaded_OR_colour_ROI_inferred_by_sealer":0,
            "botanical_wild_native_open_flower_status_independently_verified":False,
            "global_photo_observation_probability_complete_cohort":False,
            "worldwide_unbiased_white_flora_frequency_estimated":False,
            "model_group_wall_seconds":diag["group_wall_seconds"],
            "model_group_evaluator_seconds":diag["total_frozen_evaluator_seconds"],
            "original_photo_prefetch_stats":pref,
            "saved_file_SHA256":{
                data.name:sha(data),
                huefile.name:sha(huefile),
                coverage_file.name:sha(coverage_file),
                diagpaths[0].name:sha(diagpaths[0]),
            },
        }
        (staged/"original_group_receipt.json").write_text(
            json.dumps(receipt,sort_keys=True,indent=2)+"\n",encoding="utf-8")
        os.replace(staged,output)
        print("MOIS_FCP_U100_NESTED5_GENUINE_PHOTO_GROUP_SAVED",json.dumps({
            "original_group":group,
            "original_photo_attempts":expected_photos,
            "four_hue_classified":classified,
            "image_acquisition_failures":audit["acquisition_failed"],
            "acquired_unclassifiable":audit["unclassifiable_acquired"],
            "image_downloads_in_collector":0,
        },sort_keys=True),flush=True)
        return receipt
    except BaseException:
        shutil.rmtree(staged,ignore_errors=True)
        raise


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for arg in ("source","workers","output"):
        p.add_argument("--"+arg,type=Path,required=True)
    p.add_argument("--group",type=int,required=True)
    a=p.parse_args()
    seal(a.source,a.workers,a.group,a.output)


if __name__=="__main__":
    main()
