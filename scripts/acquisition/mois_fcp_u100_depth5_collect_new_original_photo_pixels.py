#!/usr/bin/env python3
"""Collect SHA-checked real original FCP U100 nested first5 newly measured pixels.

Only the exactly 34,898 genuinely never-before-measured source photo IDs from
the immutable million-photo-deduplicated FCP nested5 queue may be counted.
Individual inclusion probabilities are frozen at full nested prefix 5/N_i/M_o.
False flowering, wild/native labels and empty/unknown to white transformations
are not permitted.  This collector NEVER reacquires images or launches models.
"""
from __future__ import annotations

import argparse
from collections import Counter
import json
import math
from pathlib import Path
import shutil
import tempfile
import os

from scripts.acquisition.mois_completed_pixel_checkpoint import (
    HUES, collect_wave, require, sha,
)
from scripts.acquisition.mois_fcp_u100_depth5_newonly_queue import (
    SCHEMA as SOURCE_STAGE_SCHEMA,ARM,OUT_SCHEMA,
)
from scripts.analysis.mois_reassemble_original_new4096_measurements import METHOD

SCHEMA="mois_real_complete_FCP_U100_nested5_new34898_original_ROI_SHA_pixels_v1"
SOURCE_SHA="89e6c924e172fcc8b862daf30ca88bade2d5f1e2141134950b03928d510e4931"
PLANNED=34898
REUSED=66312
FIRST5=101210
DELTA=80968
UNITS=137
GROUPS=35


def require_manifest(path, *, expected_sha=SOURCE_SHA,
                     planned=PLANNED,reused=REUSED,
                     work_units=UNITS,groups=GROUPS,
                     total_first5=FIRST5,delta=DELTA):
    require(sha(path)==expected_sha,
            "Source-selected 34898 original-photo nested5 queue SHA changed")
    info=json.loads(Path(path).read_text(encoding="utf-8"))
    require(
        info.get("schema")==OUT_SCHEMA
        and info.get("original_sampling_schema")==SOURCE_STAGE_SCHEMA
        and info.get("full_frame_photo_ids")==93209351
        and info.get("source_arm")==ARM
        and info.get("source_group")=="FCP_U100_nested5_rank2to5_only_NOMINAL_ELIGIBLE"
        and info.get("colour_measured") is False
        and info.get("new_photo_attempts_planned")==planned
        and info.get("previously_measurement_outcomes_reused_without_inference")==reused
        and info.get("original_first5_photo_source_IDs")==total_first5
        and info.get("original_incremental_rank2to5_candidate_photo_IDs")==delta
        and info.get("source_original_eligible_taxa")==20242
        and info.get("original_photo_nested5_source_receipt_SHA256")==
            "037e74dd236d04c44f0d84c7b9d0b3e436e72fe6ddd470cc315c1dc56b0e197f"
        and info.get("previously_measured_historical_K1_FCP_1009032_photo_IDs_accounted")==1009032
        and info.get("selection_probability") is None
        and info.get("individual_original_photo_nested_prefix_inclusion_probabilities_preserved") is True
        and info.get("pending_image_work_is_only_source_photo_IDs_not_previously_inferred") is True
        and info.get("old_unknown_ROI_not_treated_as_flower_colour") is True
        and info.get("new_image_downloads_or_ROI_inferences_by_this_source_step")==0
        and info.get("source_metadata_wild_native_bloom_not_independently_verified") is True
        and info.get("technical_unknowns_preserved") is True
        and info.get("work_units")==work_units
        and info.get("worker_ids")==list(range(work_units))
        and info.get("group_ids")==list(range(groups))
        and info.get("units_per_group")==4
        and info.get("original_model_batch_size")==256
        and info.get("full_old_K1_FCP_reuse_photo_count_conservation")==total_first5
        and info.get("photo_dispositions",{}).get("new_original_photo_ROI_to_infer")==planned
        and info.get("photo_dispositions",{}).get("prior_original_photo_exact_ID_reuse")==reused
        and info.get("photo_dispositions",{}).get("original_photo_unsupported_extension_unknown")==0,
        "Actual original FCP nested depth5 selected photo source/probability contract changed")
    require(planned+reused==total_first5
            and planned>0 and work_units==(planned+255)//256,
            "Full photo-ID reuse and completed source work-unit counts not conserved")
    require(isinstance(info.get("photo_shard_sha256"),dict)
            and set(info["photo_shard_sha256"])=={str(n) for n in range(work_units)},
            "Source full original photo shard SHA inventory incomplete")
    return info


def verify_outcomes(rows, audit, *, expected_new=PLANNED):
    """Conserve signed per-photo terminal colours, unknowns and source probability."""
    require(isinstance(rows,list) and len(rows)==expected_new
            and audit.get("wave_complete") is True
            and audit.get("attempted")==expected_new
            and audit.get("planned_photo_attempts")==expected_new,
            "Unfinished newly selected original photo ROI work cannot be marked complete")
    identifiers=set();observations=set();counts=Counter()
    for row in rows:
        pid=row["photo_id"]
        obs=row["observation_id"]
        morph=row["morph"]
        status=row["measurement_status"]
        require(pid.isdecimal() and pid not in identifiers and bool(obs)
                and obs not in observations and row["source_arm"]==ARM,
                "Real nested depth5 original photo/observation or source-ID duplicate")
        identifiers.add(pid);observations.add(obs)
        p=float(row["selection_probability"])
        require(math.isfinite(p) and 0<p<=1,
                "Original full nested prefix5 photo-specific inclusion probability absent")
        cls=morph in HUES
        failed=status=="image_acquisition_failed"
        require(
            morph in (*HUES,"mixed_uncertain")
            and cls==(status=="classified_four_state_morph")
            and (failed or bool(row.get("image_sha256",""))),
            "Colour/classifiable and original image-SHA/nonresponse terminal contradiction")
        counts["attempts"]+=1
        counts["classified"]+=cls
        counts["acquisition_failed"]+=failed
        counts["acquired_unclassifiable"]+=not cls and not failed
        counts[morph]+=1
    require(counts["attempts"]==expected_new
            and counts["classified"]==audit["classifiable"]
            and counts["acquisition_failed"]==audit["acquisition_failed"]
            and counts["acquired_unclassifiable"]==audit["unclassifiable_acquired"]
            and counts["classified"]+counts["acquisition_failed"]
                +counts["acquired_unclassifiable"]==expected_new,
            "Original SHA-verified FCP nested5 actual technical unknowns lost")
    return dict(counts)


def run(source, workers, output, *, source_sha=SOURCE_SHA,
        expected_new=PLANNED,expected_reused=REUSED,
        work_units=UNITS,groups=GROUPS,total_first5=FIRST5,delta=DELTA):
    import pandas as pd
    source,workers,output=map(Path,(source,workers,output))
    require(not output.exists(),"Never overwrite original image ROI-v4 immutable release")
    info=require_manifest(source/"manifest.json",
        expected_sha=source_sha,planned=expected_new,reused=expected_reused,
        work_units=work_units,groups=groups,total_first5=total_first5,delta=delta)
    rows,audit=collect_wave(source,workers,METHOD)
    require(audit.get("wave_complete") and len(audit.get("completed_units",[]))==work_units,
            "Missing finished full-resolution original FCP nested5 pixel packets")
    counts=verify_outcomes(rows,audit,expected_new=expected_new)
    stage=Path(tempfile.mkdtemp(prefix="mois-FCP-new-nested5-ROI-",dir=output.parent))
    try:
        table=pd.DataFrame(rows)
        table["photo_id"]=table.photo_id.astype(str)
        table=table.sort_values("photo_id",key=lambda z:z.map(int)).reset_index(drop=True)
        parquet=stage/"new_FCP_U100_nested5_original_photo_ROIs.parquet"
        table.to_parquet(parquet,index=False,compression="zstd")
        receipt={
            "schema":SCHEMA,
            "original_K1_FCP_old_real_RGB_attempts_unchanged":1009032,
            "frozen_actual_pixel_source_run_and_workers_attested":True,
            "source_80968_nested_rank2to5_photo_candidates":delta,
            "source_first5_eligible_taxon_observation_photo_IDs":total_first5,
            "first5_source_taxa_nominal":20242,
            "previous_real_original_photo_ID_reuses_NOT_REINFERRED":expected_reused,
            "genuinely_new_original_photo_pixels_attempted":expected_new,
            "genuinely_new_automatic_four_hue_classified":counts["classified"],
            "genuinely_new_image_access_failed":counts["acquisition_failed"],
            "genuinely_new_acquired_ROI_unclassifiable":counts["acquired_unclassifiable"],
            "original_hue_counts":{key:counts.get(key,0) for key in (*HUES,"mixed_uncertain")},
            "one_photo_per_original_observation_in_new_nested_wave":True,
            "original_photo_full_prefix5_5_over_N_over_M_selection_probability_retained":True,
            "source_original_photo_manifest_SHA256":source_sha,
            "source_all_137_completed_work_unit_terminal_SHA_proofs":audit["completed_unit_audit"],
            "source_new_photo_Parquet_SHA256":sha(parquet),
            "original_image_downloads_or_ROI_reinferences_in_this_collector":0,
            "botanically_confirmed_wild_native_open_blooming_individuals":False,
            "globally_unbiased_white_flora_rate_estimated":False,
        }
        (stage/"receipt.json").write_text(
            json.dumps(receipt,sort_keys=True,indent=2)+"\n",encoding="utf-8")
        os.replace(stage,output)
    except BaseException:
        shutil.rmtree(stage,ignore_errors=True)
        raise
    print("MOIS_FCP_U100_NESTED5_ALL_GENUINE_NEW_ORIGINAL_PIXELS",json.dumps({
        "new_original_image_attempts":expected_new,
        "new_classified_four_hue_ROIs":counts["classified"],
        "failed_photo_access":counts["acquisition_failed"],
        "acquired_but_classification_unknown":counts["acquired_unclassifiable"],
        "old_original_image_attempts_reused":expected_reused,
        "original_images_redownloaded_in_collector":0,
    },sort_keys=True),flush=True)
    return receipt


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("source","workers","output"):
        p.add_argument("--"+name,type=Path,required=True)
    a=p.parse_args()
    run(a.source,a.workers,a.output)


if __name__=="__main__":
    main()
