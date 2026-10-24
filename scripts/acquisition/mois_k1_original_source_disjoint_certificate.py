#!/usr/bin/env python3
"""One-time SHA-authenticated old-versus-K1 source photo-ID exclusion certificate.

Reads existing signed historical ROI Parquets and the complete *precolour*
K1 original-photo work shards once; emits a small reusable gzip photo-ID
roster for further, cheaper genuine K1 global cumulative updates.
No network image requests or new flower-colour inference.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import os
from pathlib import Path
import shutil
import tempfile
from collections import Counter

from scripts.acquisition.mois_completed_pixel_checkpoint import require, sha, HUES
from scripts.acquisition.mois_93m_balanced_prior_ROI_reuse_bridge import (
    PRIOR_MAIN, PRIOR_OLD, PRIOR_R1, PRIOR_R2, prior_original_photo_index,
)
from scripts.acquisition.mois_fcp_k1_append_completed_original_photo_group import (
    K1_SOURCE_SHA, PARENT_SOURCE_SHA, SOURCE_ARM, SOURCE_PLANNED,
    SELECTED_SOURCE_OBSERVATIONS, REUSED_OLD_ORIGINALS, ORIGINAL_RUN,
)

SCHEMA="mois_sha_certified_929138_historical_disjoint_76753_original_K1_ID_source_v1"
OLD_COUNTS={"main":PRIOR_MAIN, "prior_group00":PRIOR_OLD,
            "prior_R1":PRIOR_R1, "prior_R2":PRIOR_R2}
HISTORICAL_COUNT=929138
HISTORICAL_HUES=sum(x[1] for x in OLD_COUNTS.values())
GROUPS=75
UNITS=300
COLUMNS=("photo_id","observation_id","inat_taxon_id","measurement_id",
         "selection_probability","original_worker_group")


def certify_source(source, old_pixels, old_obs, original_audits, output, *,
                   expected_old=HISTORICAL_COUNT, expected_hues=HISTORICAL_HUES,
                   expected_new=SOURCE_PLANNED, expected_sha=K1_SOURCE_SHA):
    """Verify every source-only photo ID is new before ANY new image analysis."""
    source,output=Path(source),Path(output)
    require(not output.exists(),"Signed K1 source photo ID ledger cannot overwrite old proof")
    require(len(old_pixels)==expected_old
            and sum(x[2] in HUES for x in old_pixels.values())==expected_hues,
            "Previously processed historical photo-ID/hue source count mismatch")
    manifest=source/"manifest.json"
    require(sha(manifest)==expected_sha,"Pinned source-only K1 manifest SHA changed")
    m=json.loads(manifest.read_text(encoding="utf-8"))
    require(
        m.get("schema")=="mois_full_source_observation_sampling_cohort_v1"
        and m.get("full_frame_photo_ids")==93209351
        and m.get("source_arm")==SOURCE_ARM
        and m.get("source_precolour_balanced_manifest_SHA256")==PARENT_SOURCE_SHA
        and m.get("new_photo_attempts_planned")==expected_new
        and m.get("all_balanced_selected_original_observations")==SELECTED_SOURCE_OBSERVATIONS
        and m.get("reused_original_photo_ids")==REUSED_OLD_ORIGINALS
        and m.get("work_units")==UNITS
        and m.get("worker_ids")==list(range(UNITS))
        and m.get("group_ids")==list(range(GROUPS))
        and m.get("units_per_group")==4
        and m.get("original_model_batch_size")==256
        and m.get("known_previous_full_photo_and_observation_IDs_source_SHA_verified") is True
        and m.get("previously_counted_photo_ID_duplicates_are_NOT_new_inferences") is True
        and m.get("original_photo_source_probabilities_unchanged_by_prior_pixel_reuse") is True
        and m.get("technical_unknown_not_other_colour") is True
        and m.get("images_downloaded_by_this_source_bridge")==0
        and m.get("new_colour_ROIs_inferred_by_this_source_bridge")==0,
        "Wrong full-frame original K1 probability source, group or history contract")
    rows={}
    obs=set()
    overlap_obs=set()
    groups=Counter()
    for unit in m["worker_ids"]:
        file=source/"shards"/f"shard-{unit:03d}.csv"
        require(file.is_file() and sha(file)==m["photo_shard_sha256"][str(unit)],
                "Original K1 precolour photo shard SHA mismatch")
        with file.open(newline="",encoding="utf-8") as stream:
            reader=csv.DictReader(stream)
            required={"photo_id","observation_id","inat_taxon_id","measurement_id",
                      "selection_probability","source_arm","shard"}
            require(required <= set(reader.fieldnames or []),
                    "K1 probability original-photo columns incomplete")
            shard_rows=list(reader)
        require(0<len(shard_rows)<=256,"Invalid signed original photo unit size")
        for x in shard_rows:
            photo=x["photo_id"]
            source_obs=x["observation_id"]
            taxon=x["inat_taxon_id"]
            require(photo.isdecimal() and taxon.isdecimal() and bool(source_obs)
                    and photo not in rows and photo not in old_pixels
                    and source_obs not in obs
                    and x["source_arm"]==SOURCE_ARM and x["shard"]==str(unit),
                    "Historical or new original K1 photo-ID/observation collision")
            try:
                prob=float(x["selection_probability"])
            except (TypeError,ValueError) as exc:
                raise ValueError("Invalid K1 original per-photo inclusion probability") from exc
            require(math.isfinite(prob) and 0<prob<=1,
                    "Invalid K1 original per-photo inclusion probability")
            group=unit//4
            rows[photo]=(photo,source_obs,taxon,x["measurement_id"],
                         x["selection_probability"],group)
            obs.add(source_obs)
            if source_obs in old_obs:
                overlap_obs.add(source_obs)
            groups[group]+=1
    require(len(rows)==expected_new==len(obs)
            and set(groups)==set(range(GROUPS))
            and all(groups[g]==min(1024,expected_new-g*1024) for g in range(GROUPS)),
            "Entire K1 selected-original photo group count not conserved")
    output.parent.mkdir(parents=True,exist_ok=True)
    stage=Path(tempfile.mkdtemp(prefix="mois-original-K1-disjoint-",dir=output.parent))
    try:
        compressed=stage/"K1_all_original_new_photo_source_IDs.csv.gz"
        with compressed.open("wb") as binary:
            with gzip.GzipFile(fileobj=binary,mode="wb",filename="",mtime=0) as gz:
                with io.TextIOWrapper(gz,encoding="utf-8",newline="") as stream:
                    writer=csv.writer(stream,lineterminator="\n")
                    writer.writerow(COLUMNS)
                    for photo in sorted(rows,key=int):
                        writer.writerow(rows[photo])
        info={
            "schema":SCHEMA,
            "original_full93m_photo_ID_catalogue_metadata_only":93209351,
            "source_K1_image_workflow_run":ORIGINAL_RUN,
            "historical_verified_unique_original_photo_IDs":len(old_pixels),
            "historical_verified_four_hue_original_photo_ROIs":expected_hues,
            "historical_original_photo_arm_SHA_receipts":original_audits,
            "frozen_K1_original_source_manifest_SHA256":expected_sha,
            "source_parent_probability_SHA256":PARENT_SOURCE_SHA,
            "source_arm":SOURCE_ARM,
            "K1_all_new_precolour_source_original_photo_IDs":len(rows),
            "K1_all_new_source_original_observations":len(obs),
            "K1_exact_historical_photo_ID_collisions":0,
            "K1_historical_original_observation_ID_overlaps_are_not_independence":
                len(overlap_obs),
            "K1_source_group_count":len(groups),
            "K1_new_source_photo_counts_by_group":
                {str(g):groups[g] for g in range(GROUPS)},
            "K1_new_source_photo_roster":compressed.name,
            "K1_new_source_photo_roster_sha256":sha(compressed),
            "source_selected_before_any_K1_colour_inference":True,
            "original_photo_joint_selection_probability_retained":True,
            "prior_original_pixels_recomputed_or_refetched":False,
            "actual_photo_bytes_downloaded_or_model_inference_in_this_certificate":0,
            "independently_verified_open_flower_native_wild_status":False,
            "global_white_flora_fraction_unweighted_estimated":False,
        }
        (stage/"K1_original_source_disjointness_certificate.json").write_text(
            json.dumps(info,indent=2,sort_keys=True)+"\n",encoding="utf-8")
        os.replace(stage,output)
    except BaseException:
        shutil.rmtree(stage,ignore_errors=True)
        raise
    print("MOIS_K1_COMPLETE_ORIGINAL_PHOTO_ID_DISJOINTNESS_AUDITED",json.dumps({
        "historic_photo_IDs":len(old_pixels),
        "historical_four_hue_ROIs":expected_hues,
        "K1_original_precolour_photo_IDs":len(rows),
        "exact_photo_ID_overlaps":0,
        "completed_and_pending_K1_source_groups":len(groups),
        "new_original_images_downloaded":0,
    },sort_keys=True),flush=True)
    return info


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("source","prior","output"):
        p.add_argument("--"+key,required=True,type=Path)
    a=p.parse_args()
    sources=[
        ("main",a.prior/"main"/"combined_actual_photo_checkpoint.parquet",
         a.prior/"main"/"checkpoint_receipt.json"),
        ("prior_group00",a.prior/"group00"/"measured_photos.parquet",
         a.prior/"group00"/"receipt.json"),
        ("prior_R1",a.prior/"r1"/"measured_photos.parquet",a.prior/"r1"/"receipt.json"),
        ("prior_R2",a.prior/"r2"/"measured_photos.parquet",a.prior/"r2"/"receipt.json"),
    ]
    old,old_obs,audits,overlap,_,_=prior_original_photo_index(sources,OLD_COUNTS)
    require(overlap["unique_old_photo_IDs"]==HISTORICAL_COUNT
            and not overlap.get("old_cross_arm_same_photo_evidence_conflicts",0),
            "Historical duplicate/contradictory original photo evidence")
    certify_source(a.source,old,old_obs,audits,a.output)


if __name__=="__main__":
    main()
