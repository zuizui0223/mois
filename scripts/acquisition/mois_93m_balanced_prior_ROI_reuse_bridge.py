#!/usr/bin/env python3
"""Build original FCP ROI-v4 work ONLY for new signed balanced-sample photos.

The full 93,209,351 upstream photo-ID opportunities and the authenticated
597,070-photo tier06 observation probability sample remain immutable. A separate
source-only design already selected exactly 80,035 distinct original
observations/photos with known inclusion probabilities by nominal taxon ×
recorded month × latitude band (K1).

This bridge only removes EXACT photo IDs already measured in prior true-pixel
SHA-certified arms; it does NOT replace the selected photo by another photo
of the same observation or silently drop previously seen observations, since
that would alter the sampling probability. Prior exact-photo ROIs are saved in
a re-use sidecar; the new signed ROI queue contains only never-downloaded
original photo IDs. Previously measured image labels and the new source
sample's inclusion probabilities remain separate until a full method audit.

No photo is downloaded or inferred by this script.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import shutil
import tempfile
from collections import Counter

from scripts.acquisition.mois_onepass_full93m_target_tiers import sha256, require
from scripts.acquisition.mois_queue_from_full93m_target_tiers import WORKER_COLUMNS

SOURCE_SAMPLE_SHA="c0ba3a316a336816065590540e75034bd8c2b9b7362034163673702103c44977"
SOURCE_SELECTED_SHA="b282295807fd48125ce5646385b534651da83fa7b0102634deb2d9941359ea55"
SOURCE_EXPECTED=80_035
SOURCE_ARM="mois_global93m_tier06_taxon_month_lat_probability_onephoto_K1_v1"
PRIOR_MAIN=(911_767,410_323)
PRIOR_OLD=(4_138,1_573)
PRIOR_R1=(4_468,1_760)
PRIOR_R2=(8_765,3_259)
PRIOR_SOURCE_KEYS=("main","prior_group00","prior_R1","prior_R2")
EXPECTED_UPSTREAM=93_209_351
TECH_SOURCE_SHA="7038130d6a3cc45100f11a3994da7d6735c554dcb8e6bc829c27f5bef82f2639"
SCHEMA="mois_93m_balanced_observation_onephoto_prior_original_ROI_reuse_bridge_v1"
SOURCE_SCHEMA="mois_93m_tier06_taxon_month_lat_stratified_onephoto_v1"
HUES=("white","yellow_orange","red_pink","blue_purple")


def prior_original_photo_index(sources, expected_counts):
    """Authenticated prior source arms, checksums and photo+observation IDs."""
    import pandas as pd
    ids={}
    obs_to_arm={}
    source_info=[]
    for arm,path,receipt in sources:
        require(arm in PRIOR_SOURCE_KEYS and arm not in {x["source_arm"] for x in source_info},
                "Unknown or duplicated previously signed original-photo source arm")
        path,receipt=map(Path,(path,receipt))
        r=json.loads(receipt.read_text(encoding="utf-8"))
        n,hues=expected_counts[arm]
        if arm=="main":
            require(r.get("schema")=="mois_sealed_actual_pixel_checkpoint_v1"
                    and r.get("combined_unique_photo_attempts")==n
                    and r.get("combined_four_hue_classified")==hues
                    and r["output_sha256"][path.name]==sha256(path),
                    "Previous complete six-arm original-photo Parquet or SHA changed")
        else:
            require(r.get("schema")=="mois_global93m_observation_baseline_real_RGB_photo_results_v1"
                    and r.get("attempted")==n and r.get("classifiable")==hues
                    and r["measured_photo_parquet_sha256"]==sha256(path),
                    "Previous independent original-image ROI source or SHA changed")
        columns=["photo_id","observation_id","morph","measurement_status","image_sha256"]
        frame=pd.read_parquet(path,columns=columns).fillna("")
        require(len(frame)==n and frame.photo_id.astype(str).str.fullmatch(r"\d+").all()
                and frame.photo_id.astype(str).is_unique,
                "Corrupt previously classified original photo source ID inventory")
        count=0
        for record in frame.itertuples(index=False,name=None):
            pid,obs,morph,status,image_sha=map(str,record)
            require(pid not in ids and morph in (*HUES,"mixed_uncertain"),
                    "Overlapping original photo ID across previous source sampling arms")
            known=morph in HUES
            if status:
                require(known==(status=="classified_four_state_morph"),
                        "Prior technical unknown falsely recorded as petal hue")
            if obs and obs not in ("nan","None","null"):
                require(obs not in obs_to_arm or obs_to_arm[obs]==arm,
                        "Previous cross-arm original observation unexpectedly duplicated")
                obs_to_arm[obs]=arm
            ids[pid]=(arm,obs,morph,status,image_sha)
            count+=known
        require(count==hues,"Signed historical original photo hue sum changed")
        source_info.append({
            "source_arm":arm,
            "original_photo_attempts":n,
            "four_hue_classified":hues,
            "photo_parquet_SHA256":sha256(path),
            "receipt_SHA256":sha256(receipt),
        })
    require({x["source_arm"] for x in source_info}==set(PRIOR_SOURCE_KEYS),
            "Missing prior authentic original photo source-arm ledger")
    return ids,obs_to_arm,source_info


def build(sample_manifest, sample_selected, prior_sources, technical,
          output, *,expected_sample=SOURCE_EXPECTED,
          expected_sample_sha=SOURCE_SAMPLE_SHA,
          expected_selected_sha=SOURCE_SELECTED_SHA,
          expected_prior_counts=None,
          expected_technical_sha=TECH_SOURCE_SHA,
          max_work_photos=SOURCE_EXPECTED):
    output=Path(output)
    require(not output.exists(),"Sampling selection and prior photo reuse map is immutable")
    if expected_prior_counts is None:
        expected_prior_counts={
            "main":PRIOR_MAIN,"prior_group00":PRIOR_OLD,
            "prior_R1":PRIOR_R1,"prior_R2":PRIOR_R2}
    require(sha256(sample_manifest)==expected_sample_sha
            and sha256(sample_selected)==expected_selected_sha,
            "Original full 93m-source probability sample manifest or selected roster SHA changed")
    source=json.loads(Path(sample_manifest).read_text(encoding="utf-8"))
    require(source.get("schema")=="mois_full_source_observation_sampling_cohort_v1"
            and source.get("sampling_schema")==SOURCE_SCHEMA
            and source.get("source_arm")==SOURCE_ARM
            and source.get("all_upstream_original_photo_IDs")==EXPECTED_UPSTREAM
            and source.get("original_tier06_source_photo_IDs_accounted")==597070
            and source.get("sampled_distinct_original_observations")==expected_sample
            and source.get("sampled_per_stratum_K")==1
            and source.get("source_sample_decided_without_photo_pixels") is True
            and source.get("selection_probability_is_per_photo_in_signed_shards") is True
            and source["original_source_sample_full_table_SHA256"]==sha256(sample_selected),
            "Not the exact SHA-signed full-tier observation-balanced precolour source")
    prior_ids,prior_obs,prior_audit=prior_original_photo_index(
        prior_sources,expected_prior_counts)
    if expected_technical_sha is not None:
        require(sha256(technical)==expected_technical_sha,
                "Original tier06 independent technical photo reservation has been altered")
    with Path(technical).open(newline="",encoding="utf-8") as f:
        reader=csv.DictReader(f)
        require({"photo_id","observation_id"}<=set(reader.fieldnames or ()),
                "Original pixel technical holdout identity absent")
        technical_ids={row["photo_id"] for row in reader}
    require(all(t.isdecimal() for t in technical_ids),
            "Invalid independently source-reserved original photo IDs")
    require(not technical_ids.intersection(prior_ids),
            "Technical original photos already counted as source-certified biological images")

    output.parent.mkdir(parents=True,exist_ok=True)
    stage=Path(tempfile.mkdtemp(prefix="mois-balanced-source-ROI-new-only-",dir=output.parent))
    counts=Counter()
    seen=set();observations=set()
    kept=[];reused=[]
    try:
        with Path(sample_selected).open(newline="",encoding="utf-8") as f:
            reader=csv.DictReader(f)
            require({"photo_id","observation_id","inat_taxon_id","species",
                     "photo_source_disposition",
                     "joint_photo_selection_probability",
                     "photos_in_original_observation",
                     "observations_in_original_source_stratum"}<=set(reader.fieldnames or ()),
                    "Original stratum/photo probability not retained")
            for row in reader:
                pid,obs=row["photo_id"],row["observation_id"]
                require(pid.isdecimal() and bool(obs)
                        and pid not in seen and obs not in observations
                        and row["source_arm"]==SOURCE_ARM,
                        "Duplicate or mismatched original sample photo/observation identity")
                seen.add(pid);observations.add(obs)
                m=int(row["photos_in_original_observation"])
                n=int(row["observations_in_original_source_stratum"])
                inclusion=float(row["joint_photo_selection_probability"])
                require(m>=1 and n>=1 and 0<inclusion<=1
                        and abs(inclusion/(0.0064/n/m)-1)<1e-12,
                        "Source sample photo multipicity/stratum probability changed")
                if row["photo_source_disposition"]!="supported":
                    counts["unsupported_original_format_unknown"]+=1
                    continue
                if pid in technical_ids:
                    counts["technical_holdout_photo_not_redownloaded"]+=1
                    continue
                if pid in prior_ids:
                    arm, old_obs, old_hue, old_status, old_sha=prior_ids[pid]
                    require(not old_obs or old_obs==obs,
                            "Reused genuine photo ID has conflicting original observation UUID")
                    reused.append({
                        "photo_id":pid,"observation_id":obs,
                        "source_arm":SOURCE_ARM,
                        "balanced_joint_photo_probability":row["joint_photo_selection_probability"],
                        "prior_measured_arm":arm,
                        "prior_automatic_ROI_morph":old_hue,
                        "prior_ROI_status":old_status,
                        "prior_original_image_SHA256":old_sha,
                    })
                    counts["exact_previous_photo_IDs_reused"]+=1
                    continue
                counts["fresh_photo_observation_already_measured_other_image"]+=obs in prior_obs
                kept.append(row)
                counts["never_before_downloaded_original_photo_ids"]+=1
        require(len(seen)==expected_sample
                and counts["never_before_downloaded_original_photo_ids"]==
                    len(kept)
                and len(reused)==counts["exact_previous_photo_IDs_reused"]
                and len(kept)<=max_work_photos,
                "All sampled 93m-original-photo ID outcomes not accounted for")
        require(len(kept)>0,"No new original photo pixels left to measure in this true-probability source")
        # Same source manifest/photo grain expected by the original frozen ROI
        # evaluator. Do not change photo probability when some previously
        # classified selected photos are re-used instead of downloaded.
        shards=stage/"shards";shards.mkdir()
        shard_sha={}
        batch_size=256
        for unit,start in enumerate(range(0,len(kept),batch_size)):
            path=shards/f"shard-{unit:03d}.csv"
            with path.open("w",encoding="utf-8",newline="") as out:
                writer=csv.DictWriter(out,fieldnames=WORKER_COLUMNS,lineterminator="\n")
                writer.writeheader()
                for row in kept[start:start+batch_size]:
                    pid=row["photo_id"]
                    mid="MOIS-"+hashlib.sha256(
                        ("mois-original-global-RGB-v1|"+pid).encode()
                    ).hexdigest()[:24].upper()
                    writer.writerow({
                        "measurement_id":mid,
                        "photo_id":pid,
                        "inat_taxon_id":row["inat_taxon_id"],
                        "species":row["species"],
                        "observation_id":row["observation_id"],
                        "observer_id":row["observer_id"],
                        "observed_on":row["observed_on"],
                        "latitude":row["latitude"],
                        "longitude":row["longitude"],
                        "positional_accuracy":row["positional_accuracy"],
                        "photo_license":row["photo_license"],
                        "photo_url_large":row["photo_url_large"],
                        "image_filename":mid+".jpg",
                        "source_arm":SOURCE_ARM,
                        "selection_probability":row["joint_photo_selection_probability"],
                        "source_snapshot_date":"2026-09-27",
                        "shard":str(unit),
                    })
            shard_sha[str(unit)]=sha256(path)
        reused_file=stage/"already_measured_original_photo_reuse.csv"
        columns=("photo_id","observation_id","source_arm",
                 "balanced_joint_photo_probability","prior_measured_arm",
                 "prior_automatic_ROI_morph","prior_ROI_status",
                 "prior_original_image_SHA256")
        with reused_file.open("w",encoding="utf-8",newline="") as f:
            writer=csv.DictWriter(f,fieldnames=columns,lineterminator="\n")
            writer.writeheader();writer.writerows(reused)
        info={
            "schema":"mois_full_source_observation_sampling_cohort_v1",
            "source_sampling_schema":SCHEMA,
            "full_frame_photo_ids":EXPECTED_UPSTREAM,
            "source_tier":6,
            "source_group":"balanced_taxon_month_lat_K1",
            "source_arm":SOURCE_ARM,
            "source_precolour_balanced_manifest_SHA256":expected_sample_sha,
            "source_precolour_balanced_selected_original_photos_SHA256":expected_selected_sha,
            "sampling_nested_joint_probability_per_photo":True,
            "selection_probability":None,
            "colour_measured":False,
            "new_photo_attempts_planned":len(kept),
            "photo_shard_sha256":shard_sha,
            "work_units":len(shard_sha),
            "worker_ids":list(range(len(shard_sha))),
            "group_ids":list(range(math.ceil(len(shard_sha)/4))),
            "units_per_group":4,
            "original_model_batch_size":batch_size,
            "all_balanced_selected_original_observations":expected_sample,
            "known_previous_full_photo_and_observation_IDs_source_SHA_verified":True,
            "prior_original_pixel_arms":prior_audit,
            "reused_original_photo_ids":len(reused),
            "reused_exact_photo_results_SHA256":sha256(reused_file),
            "previously_measured_original_observations_with_different_sample_photo":
                counts["fresh_photo_observation_already_measured_other_image"],
            "technical_holdout_images_not_redownloaded":
                counts["technical_holdout_photo_not_redownloaded"],
            "unsupported_original_formats_technical_unknown":
                counts["unsupported_original_format_unknown"],
            "original_photo_source_probabilities_unchanged_by_prior_pixel_reuse":True,
            "all_original_photo_ids_in_master_probability_frame_preserved":True,
            "technical_unknown_not_other_colour":True,
            "images_downloaded_by_this_source_bridge":0,
            "new_colour_ROIs_inferred_by_this_source_bridge":0,
            "is_a_stratified_photo_probability_sample_NOT_unweighted_global_prevalence":True,
            "photo_specific_anthesis_native_wild_not_independently_reviewed":True,
        }
        require((len(kept)+len(reused)+
                 counts["unsupported_original_format_unknown"]+
                 counts["technical_holdout_photo_not_redownloaded"])==expected_sample,
                "Some full-tier observation/photo probability sampling candidates lost")
        (stage/"manifest.json").write_text(
            json.dumps(info,indent=2,sort_keys=True)+"\n",encoding="utf-8")
        os_replace=__import__("os").replace
        os_replace(stage,output)
        print("MOIS_BALANCED_SAMPLE_ONLY_NEW_ORIGINAL_ROI_FROZEN_SHARDS",
              json.dumps({k:info[k] for k in (
                  "all_balanced_selected_original_observations",
                  "new_photo_attempts_planned",
                  "reused_original_photo_ids",
                  "previously_measured_original_observations_with_different_sample_photo",
                  "technical_holdout_images_not_redownloaded",
                  "unsupported_original_formats_technical_unknown",
                  "work_units")},sort_keys=True),flush=True)
        return info
    except BaseException:
        shutil.rmtree(stage,ignore_errors=True)
        raise


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for k in ("sample-manifest","sample-selected","main-parquet",
              "main-receipt","old-parquet","old-receipt",
              "r1-parquet","r1-receipt","r2-parquet","r2-receipt",
              "technical","output"):
        ap.add_argument("--"+k,type=Path,required=True)
    a=ap.parse_args()
    sources=(
        ("main",a.main_parquet,a.main_receipt),
        ("prior_group00",a.old_parquet,a.old_receipt),
        ("prior_R1",a.r1_parquet,a.r1_receipt),
        ("prior_R2",a.r2_parquet,a.r2_receipt),
    )
    build(a.sample_manifest,a.sample_selected,sources,a.technical,a.output)


if __name__=="__main__":
    main()
