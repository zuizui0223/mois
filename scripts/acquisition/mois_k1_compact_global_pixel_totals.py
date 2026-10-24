#!/usr/bin/env python3
"""Cheap exact worldwide original-photo counts from 929138 historical+sealed K1.

Historical 929138 genuine image-photo IDs and all 76753 PRECOLOUR K1 photo IDs
were already independently SHA-joined once.  Read that immutable source
certificate/3.95MB compressed roster plus SHA-certified completed K1 original
ROI receipts, not historical 0.9m-row Parquets or raw pixels again.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import tempfile

from scripts.acquisition.mois_completed_pixel_checkpoint import HUES, require, sha
from scripts.acquisition.mois_k1_original_source_disjoint_certificate import (
    SCHEMA as CERT_SCHEMA, HISTORICAL_COUNT, HISTORICAL_HUES, GROUPS, COLUMNS,
)
from scripts.acquisition.mois_fcp_k1_cumulative_receipts import (
    SCHEMA as CUM_SCHEMA, SOURCE_RUN, SOURCE_ARM, NEW_TARGET, NUM_GROUPS,
    REUSED_OLD_ORIGINALS, SOURCE_SHA,
)

SCHEMA = "mois_K1_global_original_photo_ID_exact_cumulative_from_frozen_source_v1"
CERT_SHA = "47fe5775475f7346cd6f331c928bfe20d3217152a45124481eee0c76726c692d"
ROSTER_SHA = "b620e56cc8386670d183f7a7a74f707cd22de5ccdbb55c1dace2c4a9617a0386"
REFERENCE_SHA = "ae25b993ea1f2f69dcfa187a2a925794210f8ce1823bb262dd3e70db4a391ffd"
REFERENCE_GROUPS = 16
REFERENCE_ORIGINALS = 16_384
REFERENCE_CLASSIFIED = 5_983


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def aggregate(certificate_path, roster_path, cumulative_path, old16_path,
              output, *,
              certificate_sha=CERT_SHA, roster_sha=ROSTER_SHA,
              reference_sha=REFERENCE_SHA,
              expected_historical=HISTORICAL_COUNT,
              expected_hues=HISTORICAL_HUES,
              expected_new_source=NEW_TARGET, expected_groups=GROUPS,
              reference_groups=REFERENCE_GROUPS):
    output=Path(output)
    require(not output.exists(),
            "Never overwrite SHA-frozen cumulative genuine original photo counts")
    require(sha(certificate_path)==certificate_sha
            and sha(roster_path)==roster_sha
            and sha(old16_path)==reference_sha,
            "Pinned independently published old/new source or reference bytes changed")
    cert=read_json(certificate_path)
    k1=read_json(cumulative_path)
    reference=read_json(old16_path)
    require(
        cert.get("schema")==CERT_SCHEMA
        and cert.get("frozen_K1_original_source_manifest_SHA256")==SOURCE_SHA
        and cert.get("historical_verified_unique_original_photo_IDs")==expected_historical
        and cert.get("historical_verified_four_hue_original_photo_ROIs")==expected_hues
        and cert.get("K1_all_new_precolour_source_original_photo_IDs")==expected_new_source
        and cert.get("K1_all_new_source_original_observations")==expected_new_source
        and cert.get("K1_exact_historical_photo_ID_collisions")==0
        and cert.get("K1_source_group_count")==expected_groups
        and cert.get("K1_new_source_photo_roster_sha256")==roster_sha
        and cert.get("source_arm")==SOURCE_ARM
        and cert.get("source_selected_before_any_K1_colour_inference") is True
        and cert.get("original_photo_joint_selection_probability_retained") is True
        and cert.get("prior_original_pixels_recomputed_or_refetched") is False
        and cert.get("actual_photo_bytes_downloaded_or_model_inference_in_this_certificate")==0
        and cert.get("independently_verified_open_flower_native_wild_status") is False
        and cert.get("global_white_flora_fraction_unweighted_estimated") is False,
        "Entire precolour K1 no-old-duplicate original source certificate invalid")
    count_by_group=Counter()
    photo_ids=set()
    obs_ids=set()
    last=-1
    with gzip.open(roster_path,"rt",encoding="utf-8",newline="") as file:
        reader=csv.DictReader(file)
        require(tuple(reader.fieldnames or ())==COLUMNS,
                "Wrong signed selected K1 original photo ID/probability roster schema")
        for item in reader:
            pid=item["photo_id"]
            obs=item["observation_id"]
            taxon=item["inat_taxon_id"]
            group_string=item["original_worker_group"]
            require(pid.isdecimal() and int(pid)>last and taxon.isdecimal()
                    and bool(obs) and bool(item["measurement_id"])
                    and group_string.isdecimal(),
                    "Selected K1 source photo IDs not original sorted unique observations")
            last=int(pid)
            group=int(group_string)
            require(group in range(expected_groups) and pid not in photo_ids
                    and obs not in obs_ids,
                    "Duplicate original K1 photo/observation or unexpected group")
            try:
                p=float(item["selection_probability"])
            except (TypeError,ValueError) as exc:
                raise ValueError("Missing original K1 source photograph probability") from exc
            require(math.isfinite(p) and 0<p<=1,
                    "K1 original source inclusion probability invalid")
            photo_ids.add(pid)
            obs_ids.add(obs)
            count_by_group[group]+=1
    require(len(photo_ids)==expected_new_source==len(obs_ids)
            and dict(count_by_group)=={int(g):n for g,n in
                cert["K1_new_source_photo_counts_by_group"].items()},
            "Complete signed K1 original source/photo group inclusion roster changed")
    require(
        k1.get("schema")==CUM_SCHEMA
        and k1.get("source_run")==SOURCE_RUN
        and k1.get("new_original_photo_target_total")==expected_new_source
        and k1.get("previously_measured_K1_photo_IDs_reused_elsewhere_NOT_NEW")==REUSED_OLD_ORIGINALS
        and k1.get("total_scheduled_original_worker_groups")==expected_groups
        and k1.get("source_precolour_probability_photo_sampling_preserved") is True
        and k1.get("within_K1_original_source_group_units_disjoint") is True
        and k1.get("new_original_images_downloaded_or_colour_inferred_in_THIS_ledger")==0
        and k1.get("this_is_NOT_a_crossarm_photo_ID_deduplicated_global_total") is True
        and k1.get("global_wild_flora_or_blooming_petal_status_not_independently_verified") is True,
        "Signed genuine completed-only K1 original-photo source receipt changed")
    items=k1.get("per_group")
    groups=k1.get("completed_original_worker_group_IDs")
    require(isinstance(items,list) and items and isinstance(groups,list)
            and sorted(groups)==groups
            and len(items)==len(groups)==len(set(groups))
            ==k1["completed_and_SHA_sealed_original_worker_groups"]
            and len(groups)<=expected_groups,
            "Pending or repeated original K1 work groups counted as measured")
    counts=Counter()
    per_group={}
    for item in items:
        group=item.get("original_worker_group")
        require(type(group) is int and group in count_by_group and group not in per_group,
                "Original source K1 completed group duplicated or outside source")
        n=item.get("new_original_photo_attempts")
        hue=item.get("four_hue_classified")
        missing=item.get("image_acquisition_failed")
        other=item.get("acquired_but_unclassifiable")
        h=item.get("original_hue_counts")
        digest=item.get("original_group_receipt_SHA256")
        require(type(n) is int and n==count_by_group[group]
                and all(type(x) is int and x>=0 for x in (hue,missing,other))
                and hue+missing+other==n
                and isinstance(h,dict) and
                    set(h)==set(HUES)|{"mixed_uncertain"}
                and sum(h[k] for k in HUES)==hue
                and h["mixed_uncertain"]==missing+other
                and isinstance(digest,str) and len(digest)==64
                and all(z in "0123456789abcdef" for z in digest)
                and item.get("distinct_conditional_joint_photo_probabilities",0)>=1,
                "K1 completed original photo terminal/hue/unknown/probability changed")
        per_group[group]=item
        counts.update(attempted=n,classified=hue,failed=missing,uncertain=other)
    require(sorted(per_group)==groups
            and counts["attempted"]==k1["completed_real_new_original_photo_attempts"]
            and counts["classified"]==k1["completed_automatic_four_hue_classifications"]
            and counts["failed"]==k1["completed_original_acquisition_failures"]
            and counts["uncertain"]==
                k1["completed_acquired_but_automatic_four_hue_unclassifiable"],
            "New original photo group evidence differs from K1 published cumulative counts")
    group_digest=hashlib.sha256("\n".join(
        f'{g}:{per_group[g]["original_group_receipt_SHA256"]}' for g in groups
    ).encode()).hexdigest()
    require(group_digest==k1["source_original_group_receipt_aggregate_SHA256"],
            "Original K1 signed group receipt set SHA differs from cumulative evidence")
    # Exact agreement with earlier expensive full 929138-row photo-ID union.
    if reference_groups:
        must=list(range(reference_groups))
        require(set(must)<=set(groups) and
                reference.get("K1_source_finished_original_worker_groups")==must
                and reference.get("historical_actual_unique_original_photo_IDs")==expected_historical
                and reference.get("historical_automatic_four_hue_photo_ROIs")==expected_hues
                and reference.get("K1_vs_old_crossarm_exact_photo_ID_overlaps")==0,
                "Independent old full-Parquet deduplicated source comparison absent")
        prefix_n=sum(per_group[g]["new_original_photo_attempts"] for g in must)
        prefix_h=sum(per_group[g]["four_hue_classified"] for g in must)
        require(prefix_n==reference["new_K1_sha_verified_original_photo_IDs_disjoint_from_history"]
                and prefix_h==reference["new_K1_frozen_four_hue_ROIs"]
                and expected_historical+prefix_n==
                    reference["exact_global_historical_plus_K1_distinct_original_photo_IDs"]
                and expected_hues+prefix_h==
                    reference["exact_global_historical_plus_K1_automatic_four_hue_photo_ROIs"],
                "Compressed fast K1 source union does not reproduce costly independent full-ID audit")
    report={
        "schema":SCHEMA,
        "full_worldwide_upstream_candidate_photo_IDs_METADATA_ONLY":93209351,
        "signed_historical_true_distinct_original_photo_IDs":expected_historical,
        "signed_historical_original_four_hue_ROIs":expected_hues,
        "source_K1_all_precolour_original_photo_IDs_disjoint_from_historical":
            expected_new_source,
        "source_original_K1_all_groups_photo_coverage":expected_groups,
        "genuine_K1_completed_and_SHA_sealed_groups":len(groups),
        "genuine_K1_completed_group_IDs":groups,
        "K1_original_new_photo_attempts":counts["attempted"],
        "K1_actual_original_four_hue_ROIs":counts["classified"],
        "K1_original_photo_acquisition_failed":counts["failed"],
        "K1_original_acquired_unclassifiable":counts["uncertain"],
        "exact_worldwide_historical_plus_completed_K1_distinct_original_photo_IDs":
            expected_historical+counts["attempted"],
        "exact_worldwide_historical_plus_completed_K1_four_hue_ROIs":
            expected_hues+counts["classified"],
        "already_old_pixel_K1_reuses_not_counted_new":REUSED_OLD_ORIGINALS,
        "zero_all_source_K1_vs_historical_photo_ID_overlap":True,
        "source_precolour_76753_ID_gzip_SHA256":roster_sha,
        "independent_initial_full_16group_worldwide_ID_union_SHA256":reference_sha,
        "frozen_original_K1_complete_source_disjoint_certificate_SHA256":certificate_sha,
        "sealed_original_K1_group_receipt_set_SHA256":group_digest,
        "counts_reproduce_independent_full_929138_ID_union":bool(reference_groups),
        "K1_photo_observation_sample_dependency_not_assumed_independent":True,
        "K1_original_model_source_provenance_and_unknowns_preserved":True,
        "additional_raw_image_downloads_or_ROI_inference_during_this_update":0,
        "historical_929138_photo_Parquets_re_downloaded_for_this_update":False,
        "independent_botanical_open_blooming_wild_native_photo_status_verified":False,
        "global_flora_white_rate_unweighted_estimated":False,
    }
    output.parent.mkdir(parents=True,exist_ok=True)
    stage=Path(tempfile.mkdtemp(prefix="mois-k1-compact-global-",dir=output.parent))
    try:
        (stage/"mois_original_photo_global_compact_cumulative_count.json").write_text(
            json.dumps(report,sort_keys=True,indent=2)+"\n",encoding="utf-8")
        os.replace(stage,output)
    except BaseException:
        shutil.rmtree(stage,ignore_errors=True)
        raise
    print("MOIS_K1_WORLDWIDE_UNIQUE_REAL_ORIGINAL_PIXELS_FROM_SOURCE_CERT",json.dumps({
        "global_original_unique_photos":
            report["exact_worldwide_historical_plus_completed_K1_distinct_original_photo_IDs"],
        "global_automatic_four_hue_ROIs":
            report["exact_worldwide_historical_plus_completed_K1_four_hue_ROIs"],
        "K1_finished_photo_worker_groups":len(groups),
        "K1_new_original_photo_attempts":counts["attempted"],
        "historical_Parquet_downloads":0,
        "new_original_image_downloads":0,
    },sort_keys=True),flush=True)
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("certificate","roster","cumulative","old16","output"):
        p.add_argument("--"+key,required=True,type=Path)
    a=p.parse_args()
    aggregate(a.certificate,a.roster,a.cumulative,a.old16,a.output)


if __name__=="__main__":
    main()
