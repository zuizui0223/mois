#!/usr/bin/env python3
"""FCP 93m U100 prefix1 -> NEW original ROI work, exact historical/K1 reuse.

The global FCP U100 sample is 20,242 original photo IDs from one independent
observation per taxon (selected source probability 1/N_taxon/M_photo). This
phase reads only source-verified Parquets from four previously measured FCP/MOIS
arms and the 80,035 already reserved, actively running balanced K1 PHOTO IDs.

Historical exact-photo identity => reuse unchanged prior ROI as an observation
with the NEW target's FCP prefix1 probability (not the old arm's propensity).
Active K1 photo identity => DEFER; do not fetch a second identical image while
the K1 wave is running. Photo ID, not observation ID, is the deduplication key:
a different photo of the same observation remains a valid random sample under
its original sampling design and is not silently replaced or omitted.

Unsupported selected image extensions and technical holdout records are kept
as signed technical unknowns, not assigned any visible flower colours.

This produces source-only frozen 256-photo original FCP ROI-v4 work units, but
does NOT launch ANY new S3 photo fetch while the K1 8-worker acquisition runs.
After K1 completion the ready shards can be separately enabled under an
explicit new-work, bounded-provider-budget workflow.
"""
from __future__ import annotations
import argparse
from collections import Counter
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import tempfile

from scripts.acquisition.mois_onepass_full93m_target_tiers import (
    require, sha256,
)
from scripts.acquisition.mois_queue_from_full93m_target_tiers import WORKER_COLUMNS
from scripts.acquisition.mois_93m_balanced_prior_ROI_reuse_bridge import (
    PRIOR_MAIN,PRIOR_OLD,PRIOR_R1,PRIOR_R2,
    prior_original_photo_index, prior_source_observation_relation,
    SOURCE_SAMPLE_SHA, SOURCE_SELECTED_SHA,
)
from scripts.acquisition.mois_fcp_u100_progressive_depth1_photo_source import (
    SCHEMA as PREFIX_SCHEMA,
    ARM as PREFIX_ARM,
    make_first_row as validate_exact_source_photo,
)

SCHEMA="mois_FCP_U100_full93m_progressive_prefix1_prior_and_active_K1_photo_reuse_v1"
OUTPUT_SOURCE="mois_full_source_observation_sampling_cohort_v1"
PREFIX_RECEIPT_SHA="2119e3e2f4c769401da5694ae8e243a409006c4c495de4c8f4de20446bd5e0c0"
PREFIX_PHOTOS_SHA="0403905f04d2def07d3f2e0caee28728a54923ef187f1fc2c63318857012b713"
FULL_TARGET=20_242
K1_TOTAL=80_035
UPSTREAM_PHOTOS=93_209_351
TECH_SHA="7038130d6a3cc45100f11a3994da7d6735c554dcb8e6bc829c27f5bef82f2639"

def technical_reservations(path, *, expected_sha=TECH_SHA):
    """Do not re-download the independently archived genuine 16-photo pilot."""
    path=Path(path)
    require(sha256(path)==expected_sha,"Original photo technical holdout source SHA changed")
    seen=set()
    observations=set()
    with path.open(newline="",encoding="utf-8") as f:
        reader=csv.DictReader(f)
        require({"photo_id","observation_id","acquired","status"}<=set(reader.fieldnames or ()),
                "Original photo technical trial needs real photo and observation IDs")
        for row in reader:
            pid=row["photo_id"]
            obs=row["observation_id"]
            require(pid.isdecimal() and bool(obs) and pid not in seen,
                    "Duplicate or malformed previously benchmarked original photo IDs")
            seen.add(pid)
            observations.add(obs)
    require(bool(seen),"No original technical photo identity reservation present")
    return seen,observations


ACTIONS=("prior_original_photo_ROI_reuse",
         "reserved_by_inflight_K1_original_photo",
         "prior_independent_technical_photo_holdout",
         "unsupported_original_photo_format_technical_unknown",
         "new_original_photo_ROI_queued")

def canonical_prefix_photo(photo, photo_license, ext, raw_url):
    # Generic check even for old/technical source: unsupported photos are not
    # treated as image-free certified petal colours.
    from scripts.acquisition.mois_onepass_full93m_target_tiers import normalized_extension
    supported=normalized_extension(ext)
    if supported:
        require(raw_url==f"https://inaturalist-open-data.s3.amazonaws.com/photos/{photo}/large.{supported}",
                "Unsupported or mismatched source photo URL")
    else:
        require(raw_url=="","Original unsupported image source URL must be missing")
    return supported


def make_work(prefix_csv, prefix_receipt,
              old_sources, k1_selected, k1_manifest, technical_csv,
              output, *,
              expected_prefix=FULL_TARGET,
              expected_prefix_receipt_sha=PREFIX_RECEIPT_SHA,
              expected_prefix_selected_sha=PREFIX_PHOTOS_SHA,
              expected_k1=K1_TOTAL,
              expected_k1_selected_sha=SOURCE_SELECTED_SHA,
              expected_k1_manifest_sha=SOURCE_SAMPLE_SHA,
              expected_old=None,
              expected_technical_sha=TECH_SHA,
              max_new=FULL_TARGET):
    prefix_csv,prefix_receipt,k1_selected,k1_manifest,output=map(
        Path,(prefix_csv,prefix_receipt,k1_selected,k1_manifest,output))
    require(not output.exists(),"Never overwrite signed FCP U100 first-photo pixel reuse queue")
    if expected_old is None:
        expected_old=dict(main=PRIOR_MAIN,prior_group00=PRIOR_OLD,
                          prior_R1=PRIOR_R1,prior_R2=PRIOR_R2)
    require(sha256(prefix_csv)==expected_prefix_selected_sha
            and sha256(prefix_receipt)==expected_prefix_receipt_sha,
            "Globally frozen full93m FCP-U100 photo prefix SHA differs")
    signed=json.loads(prefix_receipt.read_text(encoding="utf-8"))
    require(signed.get("schema")==
                 "mois_full93m_fcp_U100_prefix1_complete_20242_nominal_taxa_v1"
            and signed.get("first_original_photo_candidates")==expected_prefix
            and signed.get("full_93m_original_photo_id_opportunities")==UPSTREAM_PHOTOS
            and signed.get("selected_original_photo_candidates_SHA256")==sha256(prefix_csv)
            and signed.get("original_photo_inclusion_probability")==
                "1/N_taxon/M_original_photos_per_observation, conditional on 100-obs source eligibility",
            "No source-audited first one-photo-per-FCP-U100-eligible-taxon sample")
    require(sha256(k1_selected)==expected_k1_selected_sha
            and sha256(k1_manifest)==expected_k1_manifest_sha,
            "Actively running K1 original-photo source SHA changed")
    k1=json.loads(k1_manifest.read_text(encoding="utf-8"))
    require(k1.get("schema")==OUTPUT_SOURCE
            and k1.get("sampled_distinct_original_observations")==expected_k1
            and k1.get("sampling_schema")==
                "mois_93m_tier06_taxon_month_lat_stratified_onephoto_v1"
            and k1.get("source_sample_decided_without_photo_pixels") is True,
            "K1 photo queue is not the signed 80,035-observation source")
    k1_photos={}
    k1_obs=set()
    with k1_selected.open(newline="",encoding="utf-8") as f:
        reader=csv.DictReader(f)
        require({"photo_id","observation_id","photo_source_disposition",
                 "joint_photo_selection_probability"} <=set(reader.fieldnames or ()),
                "K1 exact photo-ID source exposure missing")
        for row in reader:
            pid=row["photo_id"]
            obs=row["observation_id"]
            require(pid.isdecimal() and bool(obs) and pid not in k1_photos,
                    "K1 previously selected original photo IDs not unique")
            k1_photos[pid]=row
            k1_obs.add(obs)
    require(len(k1_photos)==expected_k1,
            "K1 original photo source sample changed before prior reuse")

    prior,prior_observations,prior_metadata,prior_overlap,prior_duplicates,prior_conflicts = (
        prior_original_photo_index(old_sources,expected_old)
    )
    technical,technical_obs=technical_reservations(
        technical_csv,expected_sha=expected_technical_sha
    )
    require(not technical & set(prior),
            "Original photo technical holdout unexpectedly reused in prior pixel checkpoint")

    output.parent.mkdir(parents=True,exist_ok=True)
    staged=Path(tempfile.mkdtemp(prefix="mois-fcp-U100-prefix1-reuse-",dir=output.parent))
    try:
        chosen=staged/"all_FCP_U100_prefix1_photo_source_dispositions.csv"
        fields=[
            "photo_id","observation_id","inat_taxon_id","species",
            "recorded_observation_month","latitude_band",
            "depth1_joint_original_photo_inclusion_probability",
            "source_fcp_u100_parent_joint_photo_probability",
            "photo_source_disposition","source_arm",
            "prior_photo_arm","prior_original_ROI_morph",
            "prior_original_ROI_status","prior_source_observation_relation",
            "reuse_or_new_photo_outcome",
        ]
        new=[]
        counts=Counter()
        uid=set()
        seen_taxa=set()
        overlaps=Counter()
        with prefix_csv.open("r",newline="",encoding="utf-8") as src,\
             chosen.open("w",newline="",encoding="utf-8") as dst:
            reader=csv.DictReader(src)
            require({"photo_id","observation_id","inat_taxon_id",
                     "depth1_joint_original_photo_inclusion_probability",
                     "parent_U100_joint_original_photo_inclusion_probability",
                     "photos_in_original_observation",
                     "full_93m_distinct_original_observations_of_taxon",
                     "photo_source_disposition","photo_url_large"} <=set(reader.fieldnames or ()),
                    "FCP-U100 ranked prefix original ID/probability columns missing")
            writer=csv.DictWriter(dst,fieldnames=fields,lineterminator="\n")
            writer.writeheader()
            for row in reader:
                pid=row["photo_id"];obs=row["observation_id"]
                taxon=row["inat_taxon_id"]
                require(pid.isdecimal() and bool(obs) and taxon.isdecimal()
                        and pid not in uid and taxon not in seen_taxa,
                        "Selected FCP one-observation photo/taxon source not unique")
                uid.add(pid);seen_taxa.add(taxon)
                n=int(row["full_93m_distinct_original_observations_of_taxon"])
                m=int(row["photos_in_original_observation"])
                pi=float(row["depth1_joint_original_photo_inclusion_probability"])
                parentpi=float(row["parent_U100_joint_original_photo_inclusion_probability"])
                require(n>=100 and m>=1 and
                        0<pi<=1 and math.isclose(pi,1/n/m,rel_tol=1e-12) and
                        math.isclose(parentpi/pi,100,rel_tol=1e-12) and
                        row["progressive_source_arm"]==PREFIX_ARM,
                        "True FCP prefix sampling inclusion probability/provenance changed")
                normalized=canonical_prefix_photo(
                    pid,row["photo_license"],row["original_image_extension"],
                    row["photo_url_large"])
                old=prior.get(pid)
                previous_arm=previous_morph=previous_status=relationship=""
                if old is not None:
                    previous_arm,old_obs,previous_morph,previous_status,_sha=old
                    relationship=prior_source_observation_relation(old_obs,obs)
                    action="prior_original_photo_ROI_reuse"
                elif pid in k1_photos:
                    require(k1_photos[pid]["observation_id"]==obs,
                            "K1 in-flight photo ID has different signed observation UUID")
                    action="reserved_by_inflight_K1_original_photo"
                elif pid in technical:
                    action="prior_independent_technical_photo_holdout"
                elif normalized is None:
                    action="unsupported_original_photo_format_technical_unknown"
                else:
                    action="new_original_photo_ROI_queued"
                    new.append(row)
                require(action in ACTIONS,"Unexpected original photo source disposition")
                counts[action]+=1
                if obs in k1_obs and pid not in k1_photos:
                    overlaps["K1_observation_other_photo_preserved_as_separate_random_photo"]+=1
                if obs in prior_observations and pid not in prior:
                    overlaps["previous_observation_other_photo_preserved_as_separate_random_photo"]+=1
                writer.writerow({
                    "photo_id":pid,"observation_id":obs,
                    "inat_taxon_id":taxon,"species":row["species"],
                    "recorded_observation_month":row["recorded_month"],
                    "latitude_band":row["latitude_band"],
                    "depth1_joint_original_photo_inclusion_probability":row[
                        "depth1_joint_original_photo_inclusion_probability"],
                    "source_fcp_u100_parent_joint_photo_probability":row[
                        "parent_U100_joint_original_photo_inclusion_probability"],
                    "photo_source_disposition":row["photo_source_disposition"],
                    "source_arm":PREFIX_ARM,
                    "prior_photo_arm":previous_arm,
                    "prior_original_ROI_morph":previous_morph,
                    "prior_original_ROI_status":previous_status,
                    "prior_source_observation_relation":relationship,
                    "reuse_or_new_photo_outcome":action,
                })
        require(len(uid)==len(seen_taxa)==expected_prefix
                and sum(counts.values())==expected_prefix
                and len(new)==counts["new_original_photo_ROI_queued"]
                and len(new)<=max_new,
                "FCP original photo overlap and next ROI target accounting changed")
        new.sort(key=lambda row:(
            hashlib.sha256(("mois-fcp-prefix1-balanced-job-shuffle-v1|"+row["observation_id"]).encode()).hexdigest(),
            int(row["photo_id"])))
        folder=staged/"shards"
        folder.mkdir()
        hashes={}
        for part in range((len(new)+255)//256):
            file=folder/f"shard-{part:03d}.csv"
            with file.open("w",newline="",encoding="utf-8") as f:
                w=csv.DictWriter(f,fieldnames=WORKER_COLUMNS,lineterminator="\n")
                w.writeheader()
                for z in new[part*256:(part+1)*256]:
                    pid=z["photo_id"]
                    mid="MOIS-"+hashlib.sha256(
                        ("mois-original-global-RGB-v1|"+pid).encode()).hexdigest()[:24].upper()
                    w.writerow({
                        "measurement_id":mid,"photo_id":pid,
                        "inat_taxon_id":z["inat_taxon_id"],"species":z["species"],
                        "observation_id":z["observation_id"],"observer_id":z["observer_id"],
                        "observed_on":z["observed_on"],"latitude":z["latitude"],
                        "longitude":z["longitude"],"positional_accuracy":z["positional_accuracy"],
                        "photo_license":z["photo_license"],"photo_url_large":z["photo_url_large"],
                        "image_filename":mid+".jpg","source_arm":PREFIX_ARM,
                        "selection_probability":z["depth1_joint_original_photo_inclusion_probability"],
                        "source_snapshot_date":"2026-09-27","shard":str(part),
                    })
            hashes[str(part)]=sha256(file)
        n_units=len(hashes)
        report={
            "schema":OUTPUT_SOURCE,
            "fcp_prefix_source_schema":SCHEMA,
            "full_frame_photo_ids":UPSTREAM_PHOTOS,
            "source_arm":PREFIX_ARM,
            "selection_probability":None,
            "source_photo_inclusion_probability_is_per_original_photo":True,
            "selected_FCP_U100_eligible_source_taxa":expected_prefix,
            "selected_original_photo_observations":expected_prefix,
            "previously_measured_original_photo_ids_reused":counts["prior_original_photo_ROI_reuse"],
            "inflight_K1_photo_IDs_reserved_without_redownload":
                counts["reserved_by_inflight_K1_original_photo"],
            "technical_pilot_photo_IDs_not_refetched":
                counts["prior_independent_technical_photo_holdout"],
            "unsupported_original_photo_format_selected_as_unknown":
                counts["unsupported_original_photo_format_technical_unknown"],
            "new_photo_attempts_planned":len(new),
            "frozen_FCP_original_ROI_source_work_units":n_units,
            "work_units":n_units,
            "group_ids":list(range((n_units+3)//4)),
            "worker_ids":list(range(n_units)),
            "units_per_group":4,
            "original_model_batch_size":256,
            "photo_shard_sha256":hashes,
            "photo_ID_technical_unknown_and_reuse_dispositions":dict(counts),
            "other_photo_same_original_observation_kept_for_design":dict(overlaps),
            "prior_independent_real_photo_source_info":prior_metadata,
            "prior_exact_photo_crossarm_overlap_summary":prior_overlap,
            "prior_old_original_photo_duplicate_rows":len(prior_duplicates),
            "prior_old_original_photo_conflicting_evidence_rows":len(prior_conflicts),
            "source_full93m_prefix1_receipt_SHA256":sha256(prefix_receipt),
            "source_full93m_prefix1_selected_csv_SHA256":sha256(prefix_csv),
            "source_original_K1_inflight_photo_IDs_SHA256":sha256(k1_selected),
            "prior_technical_original_photo_reservations_SHA256":sha256(technical_csv),
            "frozen_all_selected_FCP_photo_dispositions_SHA256":sha256(chosen),
            "FCP_top100_photo_prefix_inclusion_probability_corrected_to_1_over_N_M":True,
            "original_sampled_observation_not_replaced_when_already_measured_in_another_arm":True,
            "historic_pixels_reused_at_exact_original_photo_ID_only":True,
            "photo_already_reserved_by_K1_never_fetched_twice":True,
            "native_flowering_botanical_and_syndrome_status_verified":False,
            "colour_measured":False,
            "original_images_downloaded_by_source_bridge":0,
            "no_original_image_ROI_inference_by_source_bridge":True,
            "not_globally_unbiased_native_flora_white_flower_rate":True,
            "new_FCP_photo_ROI_jobs_not_launched_until_separate_resource_and_duplicate_gate":True,
        }
        (staged/"manifest.json").write_text(
            json.dumps(report,indent=2,sort_keys=True)+"\n",encoding="utf-8")
        os.replace(staged,output)
        print("MOIS_FCP_U100_PREFIX1_SOURCE_REUSE_WITH_K1_INFLIGHT_SAFETY",
              json.dumps({k:report[k] for k in (
                  "selected_original_photo_observations",
                  "previously_measured_original_photo_ids_reused",
                  "inflight_K1_photo_IDs_reserved_without_redownload",
                  "technical_pilot_photo_IDs_not_refetched",
                  "unsupported_original_photo_format_selected_as_unknown",
                  "new_photo_attempts_planned")},sort_keys=True),flush=True)
        return report
    except BaseException:
        shutil.rmtree(staged,ignore_errors=True)
        raise


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for flag in ("prefix-csv","prefix-receipt","k1-selected","k1-manifest",
                 "main-parquet","main-receipt",
                 "prior-group00-parquet","prior-group00-receipt",
                 "prior-r1-parquet","prior-r1-receipt",
                 "prior-r2-parquet","prior-r2-receipt","technical","output"):
        ap.add_argument("--"+flag,type=Path,required=True)
    a=ap.parse_args()
    sources=[
        ("main",a.main_parquet,a.main_receipt),
        ("prior_group00",a.prior_group00_parquet,a.prior_group00_receipt),
        ("prior_R1",a.prior_r1_parquet,a.prior_r1_receipt),
        ("prior_R2",a.prior_r2_parquet,a.prior_r2_receipt),
    ]
    make_work(a.prefix_csv,a.prefix_receipt,sources,
              a.k1_selected,a.k1_manifest,a.technical,a.output)


if __name__=="__main__":
    main()
