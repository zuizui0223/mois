#!/usr/bin/env python3
"""Preflight exactly 3141 still-unmeasured FCP U100 first-photo sources BEFORE K1 ends.

Data-only, outcome-blind identity/probability/sha contract:
* Audit all 20,242 one-observation-per-taxon original-photo decisions.
* Compare 15,483 old-photo reuses, 1,618 K1-reserved photo IDs and exactly
  3,141 NEW original photos against the immutable K1 76,753-photo source.
* Verify each of the 13 frozen 256-photo YOLO11/EfficientSAM work shards.
No original photo HTTP requests, model imports/inference or image labels.
Do not claim native/wild, independently visible blooming or unbiased prevalence.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import tempfile

from scripts.acquisition.mois_completed_pixel_checkpoint import sha, require
from scripts.acquisition.mois_fcp_u100_collect_new_first_photo_ROIs import (
    PREFIX_PHOTOS_SHA, SOURCE_SHA, NEW, PRIOR_REUSE, RESERVED_K1,
    check_prefix_probabilities,
)
from scripts.acquisition.mois_fcp_u100_depth1_prior_k1_reuse_queue import (
    PREFIX_ARM, FULL_TARGET, canonical_prefix_photo, ACTIONS,
)
from scripts.acquisition.mois_k1_original_source_disjoint_certificate import (
    SCHEMA as K1_CERT_SCHEMA,
)

SCHEMA = "mois_source_only_20242_FCP_U100_prefix1_3141_new_K1_disjoint_preflight_v1"
K1_CERT_SHA = "47fe5775475f7346cd6f331c928bfe20d3217152a45124481eee0c76726c692d"
K1_ROSTER_SHA = "b620e56cc8386670d183f7a7a74f707cd22de5ccdbb55c1dace2c4a9617a0386"
DISPOSITIONS_SHA = "dacf3a8df021bac2de9f6f02fa779b246e71e544be4e473296f90ad19f66e63f"
SOURCE_TAR_SHA = "48b31262c20db2501889106f8ea12c650fd875cf85ed1570c10f8bdc986747a6"


def preflight(source, prefix_csv, dispositions_csv, k1_certificate,
              k1_roster_gz, output, *,
              source_sha=SOURCE_SHA, prefix_sha=PREFIX_PHOTOS_SHA,
              dispositions_sha=DISPOSITIONS_SHA,
              k1_certificate_sha=K1_CERT_SHA, k1_roster_sha=K1_ROSTER_SHA,
              eligible=FULL_TARGET, reuse=PRIOR_REUSE, reserved=RESERVED_K1,
              new=NEW, expected_work_units=13, expected_groups=4,
              expected_K1_source=76753):
    source, output=Path(source),Path(output)
    require(not output.exists(), "Refuse to overwrite a previous signed U100 preflight")
    require(eligible==reuse+reserved+new,
            "Every eligible FCP U100 nominal taxon must have a source disposition")
    manifest=source/"manifest.json"
    require(
        sha(manifest)==source_sha and sha(prefix_csv)==prefix_sha
        and sha(dispositions_csv)==dispositions_sha
        and sha(k1_certificate)==k1_certificate_sha
        and sha(k1_roster_gz)==k1_roster_sha,
        "Frozen original 93m prefix/K1 source manifest or photo ID bytes altered")
    doc=json.loads(manifest.read_text(encoding="utf-8"))
    require(doc.get("schema")=="mois_full_source_observation_sampling_cohort_v1"
            and doc.get("full_frame_photo_ids")==93209351
            and doc.get("source_arm")==PREFIX_ARM
            and doc.get("selected_FCP_U100_eligible_source_taxa")==eligible
            and doc.get("selected_original_photo_observations")==eligible
            and doc.get("previously_measured_original_photo_ids_reused")==reuse
            and doc.get("inflight_K1_photo_IDs_reserved_without_redownload")==reserved
            and doc.get("new_photo_attempts_planned")==new
            and doc.get("work_units")==expected_work_units
            and doc.get("worker_ids")==list(range(expected_work_units))
            and doc.get("group_ids")==list(range(expected_groups))
            and doc.get("units_per_group")==4
            and doc.get("original_model_batch_size")==256
            and doc.get("technical_pilot_photo_IDs_not_refetched")==0
            and doc.get("unsupported_original_photo_format_selected_as_unknown")==0
            and doc.get("source_photo_inclusion_probability_is_per_original_photo") is True
            and doc.get("photo_already_reserved_by_K1_never_fetched_twice") is True
            and doc.get("FCP_top100_photo_prefix_inclusion_probability_corrected_to_1_over_N_M") is True
            and doc.get("historic_pixels_reused_at_exact_original_photo_ID_only") is True
            and doc.get("native_flowering_botanical_and_syndrome_status_verified") is False
            and doc.get("colour_measured") is False
            and doc.get("original_images_downloaded_by_source_bridge")==0
            and doc.get("no_original_image_ROI_inference_by_source_bridge") is True
            and doc.get("new_FCP_photo_ROI_jobs_not_launched_until_separate_resource_and_duplicate_gate") is True
            and doc.get("frozen_all_selected_FCP_photo_dispositions_SHA256")
                ==dispositions_sha
            and doc.get("source_full93m_prefix1_selected_csv_SHA256")==prefix_sha,
            "FCP U100 actual 3141-photo image-free source plan was modified")
    k1_doc=json.loads(Path(k1_certificate).read_text(encoding="utf-8"))
    require(k1_doc.get("schema")==K1_CERT_SCHEMA
            and k1_doc.get("K1_all_new_precolour_source_original_photo_IDs")
                ==expected_K1_source
            and k1_doc.get("K1_exact_historical_photo_ID_collisions")==0
            and k1_doc.get("K1_new_source_photo_roster_sha256")==k1_roster_sha
            and k1_doc.get("original_photo_joint_selection_probability_retained") is True
            and k1_doc.get("actual_photo_bytes_downloaded_or_model_inference_in_this_certificate")==0,
            "Signed 76753 K1 new-photo source exclusion certificate invalid")
    prior_source=check_prefix_probabilities(prefix_csv,prefix_sha=prefix_sha,
                                            expected_eligible=eligible)
    prefix_photo={}
    with Path(prefix_csv).open(newline="",encoding="utf-8") as stream:
        reader=csv.DictReader(stream)
        require({"photo_id","observation_id","inat_taxon_id",
                 "photo_url_large","photo_license","original_image_extension",
                 "parent_U100_joint_original_photo_inclusion_probability",
                 "depth1_joint_original_photo_inclusion_probability",
                 "progressive_source_arm"}<=set(reader.fieldnames or ()),
                "Original FCP source photo metadata missing")
        for row in reader:
            pid=row["photo_id"]
            obs,taxon,p=prior_source[pid]
            require(row["observation_id"]==obs
                    and row["inat_taxon_id"]==taxon
                    and row["progressive_source_arm"]==PREFIX_ARM,
                    "Original prefix photo/taxon/provenance mismatch")
            parent=float(row["parent_U100_joint_original_photo_inclusion_probability"])
            require(math.isfinite(parent)
                    and math.isclose(parent/p,100,rel_tol=1e-12),
                    "U100 original 100-observation parent photo probability changed")
            ext=canonical_prefix_photo(
                pid,row["photo_license"],row["original_image_extension"],
                row["photo_url_large"])
            prefix_photo[pid]=(obs,taxon,p,ext)
    require(len(prefix_photo)==eligible, "Original eligible taxon prefix lost a photo")
    k1={}
    with gzip.open(k1_roster_gz,"rt",encoding="utf-8",newline="") as stream:
        reader=csv.DictReader(stream)
        require({"photo_id","observation_id","inat_taxon_id",
                 "selection_probability","original_worker_group"}<=set(reader.fieldnames or ()),
                "Original K1 source photo ID roster missing selection metadata")
        for row in reader:
            pid=row["photo_id"]
            require(pid.isdecimal() and pid not in k1,
                    "Original K1 source contains a repeated photo ID")
            k1[pid]=(row["observation_id"],row["inat_taxon_id"])
    require(len(k1)==expected_K1_source,
            "Not all 76753 distinct precolour new K1 photo IDs were checked")

    disposition={}
    actions={a:0 for a in ACTIONS}
    seen_taxa=set()
    seen_obs=set()
    with Path(dispositions_csv).open(newline="",encoding="utf-8") as stream:
        reader=csv.DictReader(stream)
        required={"photo_id","observation_id","inat_taxon_id",
                  "depth1_joint_original_photo_inclusion_probability",
                  "source_fcp_u100_parent_joint_photo_probability",
                  "source_arm","reuse_or_new_photo_outcome",
                  "prior_original_ROI_morph","prior_photo_arm"}
        require(required<=set(reader.fieldnames or ()),
                "Frozen photo disposition identity or unknown columns missing")
        for row in reader:
            pid=row["photo_id"]
            require(pid in prefix_photo and pid not in disposition,
                    "FCP U100 original selected photo ID missing/repeated")
            obs,taxon,p,supported=prefix_photo[pid]
            action=row["reuse_or_new_photo_outcome"]
            require(row["observation_id"]==obs
                    and row["inat_taxon_id"]==taxon
                    and row["source_arm"]==PREFIX_ARM
                    and action in actions
                    and obs not in seen_obs and taxon not in seen_taxa,
                    "FCP U100 original photo/observation/taxon source disposition altered")
            try:
                inexact=float(row["depth1_joint_original_photo_inclusion_probability"])
                parent=float(row["source_fcp_u100_parent_joint_photo_probability"])
            except (ValueError,TypeError) as err:
                raise ValueError("Missing exact original FCP prefix probability") from err
            require(math.isclose(inexact,p,rel_tol=1e-12)
                    and math.isclose(parent,100*p,rel_tol=1e-12),
                    "FCP photo inclusion probability differs from signed prefix1")
            if action=="reserved_by_inflight_K1_original_photo":
                require(pid in k1 and k1[pid]==(obs,taxon),
                        "K1 reserved original image ID not in its frozen source")
            else:
                require(pid not in k1,
                        "An already active K1 original image would be downloaded twice")
            if action=="new_original_photo_ROI_queued":
                require(supported is not None and bool(prefix_photo[pid][0])
                        and not row["prior_photo_arm"]
                        and not row["prior_original_ROI_morph"],
                        "FCP proposed new original photo was already measured or inaccessible")
            disposition[pid]=(obs,taxon,p,action)
            actions[action]+=1
            seen_taxa.add(taxon)
            seen_obs.add(obs)
    require(len(disposition)==eligible
            and actions["prior_original_photo_ROI_reuse"]==reuse
            and actions["reserved_by_inflight_K1_original_photo"]==reserved
            and actions["new_original_photo_ROI_queued"]==new
            and actions["prior_independent_technical_photo_holdout"]==0
            and actions["unsupported_original_photo_format_technical_unknown"]==0
            # Signed Counter snapshots omit zero-valued categories.
            # Reconcile only nonzero actions, while separately requiring both
            # technical/unsupported counts to be exactly zero above.
            and {k:v for k,v in actions.items() if v}
                ==doc["photo_ID_technical_unknown_and_reuse_dispositions"],
            "Not all 20242 FCP-U100 selection/old/K1/new original-photo outcomes conserved")

    shards={}
    sharded=set()
    sharded_taxa=set()
    sharded_obs=set()
    for unit in range(expected_work_units):
        file=source/"shards"/f"shard-{unit:03d}.csv"
        require(file.is_file() and sha(file)==doc["photo_shard_sha256"][str(unit)],
                "Frozen original FCP U100 3141 image source shard SHA changed")
        with file.open(newline="",encoding="utf-8") as stream:
            reader=csv.DictReader(stream)
            required={"photo_id","observation_id","inat_taxon_id",
                      "selection_probability","source_arm","shard",
                      "photo_url_large","measurement_id"}
            require(required<=set(reader.fieldnames or ()),
                    "New-original-photo worker shard missing essential ID/ROI columns")
            rows=list(reader)
        require(0<len(rows)<=256
                and len(rows)==min(256,new-256*unit),
                "FCP-U100 original 256-photo work shard size incorrect")
        for row in rows:
            pid=row["photo_id"]
            require(pid in disposition and pid not in sharded
                    and row["observation_id"] not in sharded_obs
                    and row["inat_taxon_id"] not in sharded_taxa
                    and row["source_arm"]==PREFIX_ARM
                    and row["shard"]==str(unit),
                    "FCP first-photo work packets reused or altered original ID")
            obs,taxon,p,action=disposition[pid]
            try: actual_p=float(row["selection_probability"])
            except (ValueError,TypeError) as err:
                raise ValueError("Missing per-photo original worker probability") from err
            require(action=="new_original_photo_ROI_queued"
                    and row["observation_id"]==obs and row["inat_taxon_id"]==taxon
                    and math.isclose(actual_p,p,rel_tol=1e-12)
                    and bool(row["photo_url_large"]) and bool(row["measurement_id"]),
                    "An old/K1 original photo was scheduled for second image inference")
            sharded.add(pid)
            sharded_obs.add(obs)
            sharded_taxa.add(taxon)
        shards[str(unit)]={"original_new_photo_rows":len(rows),"original_CSV_sha256":sha(file)}
    expected_new_ids={pid for pid,d in disposition.items()
                      if d[3]=="new_original_photo_ROI_queued"}
    require(sharded==expected_new_ids and len(sharded)==new
            and set(shards)==set(doc["photo_shard_sha256"]),
            "Frozen FCP U100 work queue is not the 3141 unique never-measured photos")
    report={
        "schema":SCHEMA,
        "signed_original_source_manifest_sha256":source_sha,
        "signed_original_prefix_taxon_photo_SHA256":prefix_sha,
        "signed_20242_photo_dispositions_sha256":dispositions_sha,
        "signed_K1_all76753_new_original_source_roster_sha256":k1_roster_sha,
        "signed_K1_old_historical_photo_ID_disjoint_certificate_sha256":k1_certificate_sha,
        "upstream_original_photo_IDs_metadata_only":93209351,
        "eligible_nominal_FCP_U100_source_taxa":eligible,
        "source_original_observations_one_per_taxon":eligible,
        "already_measured_old_photo_IDs_to_reuse_NOT_RECOMPUTED":reuse,
        "already_running_K1_new_photo_IDs_to_reuse_NOT_RECOMPUTED":reserved,
        "never_before_measured_FCP_U100_original_photo_ROI_target_IDs":new,
        "signed_13_original_work_unit_audit":shards,
        "four_independent_groups_planned":expected_groups,
        "new_original_photos_overlap_with_K1":0,
        "old_historical_original_photo_IDs_recomputed":0,
        "technical_nonresponse_unknowns_kept_as_unknowns":True,
        "all_photo_inclusion_probabilities_individually_source_checked":True,
        "all_photo_taxtype_observation_IDs_individually_source_checked":True,
        "new_original_image_bytes_downloaded_here":0,
        "flower_colour_models_run_here":0,
        "new_original_photo_colour_results_generated_here":0,
        "true_wild_native_flowering_status_not_independently_validated":True,
        "ready_for_inference_ONLY_AFTER_K1_75_WORKERS_AND_75_SIGNED_GROUP_RELEASES":True,
    }
    output.parent.mkdir(parents=True,exist_ok=True)
    stage=Path(tempfile.mkdtemp(prefix="mois-u100-source-preflight-",dir=output.parent))
    try:
        (stage/"fcp_u100_3141_photo_only_source_preflight_receipt.json").write_text(
            json.dumps(report,sort_keys=True,indent=2)+"\n",encoding="utf-8")
        os.replace(stage,output)
    except BaseException:
        shutil.rmtree(stage,ignore_errors=True)
        raise
    print("MOIS_FCP_U100_3141_ORIGINAL_PHOTO_PRECOLOUR_PREFLIGHT",json.dumps({
        "full_candidate_metadata":93209351,
        "eligible_taxa":eligible,"already_old_RGB_photo_IDs":reuse,
        "reserved_K1_image_photo_IDs":reserved,
        "exact_new_original_photo_IDs_undownloaded":new,
        "source_verified_frozen_units":len(shards),
        "duplicate_new_K1_image_downloads":0,
        "original_image_downloads_or_colour_inferences":0,
    },sort_keys=True),flush=True)
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ("source","prefix","dispositions","k1-certificate",
                "k1-roster","output"):
        parser.add_argument("--"+key,type=Path,required=True)
    a=parser.parse_args()
    preflight(a.source,a.prefix,a.dispositions,
              a.k1_certificate,a.k1_roster,a.output)


if __name__=="__main__":
    main()
