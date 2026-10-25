#!/usr/bin/env python3
"""Exactly deduplicate FCP U100 ranks2..5 against finished genuine MOIS pixels.

Never put a photograph into a NEW original image ROI worker if its photo ID
already occurs in any genuine SHA-verified historical/K1/FCP original terminal,
even when the old terminal hue is technically unknown. Preserve nested original
observation sampling, with rank2..5 incremental photo inclusion probability
4 / N_original_observations_taxon / M_photos_in_chosen_observation.

Source and per-taxon selection are blind to colour. Missing flowering/native/
wild status is never inferred and missing hues remain unknown.
"""
from __future__ import annotations

import argparse
from collections import Counter,defaultdict
import csv
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import tempfile

from scripts.acquisition.mois_completed_pixel_checkpoint import HUES, sha, require
from scripts.acquisition.mois_fcp_u100_depth5_global_source import (
    SCHEMA as GLOBAL_SCHEMA,STAGE_INCREMENTAL,GLOBAL_ELIGIBLE,
)
from scripts.acquisition.mois_fcp_u100_depth5_source import (
    ARM as NESTED_ARM,FIELDS as SOURCE_FIELDS,DEPTH,
)
from scripts.acquisition.mois_93m_balanced_prior_ROI_reuse_bridge import (
    PRIOR_MAIN,PRIOR_OLD,PRIOR_R1,PRIOR_R2,
    prior_original_photo_index,prior_source_observation_relation,
)
from scripts.acquisition.mois_onepass_full93m_target_tiers import (
    normalized_extension,
)
from scripts.acquisition.mois_queue_from_full93m_target_tiers import WORKER_COLUMNS

SCHEMA="mois_FCP_U100_rank2to5_new_original_photo_pixels_after_all_prior_ID_dedup_v1"
GLOBAL_STAGE5_RECEIPT_SHA="037e74dd236d04c44f0d84c7b9d0b3e436e72fe6ddd470cc315c1dc56b0e197f"
GLOBAL_RANK2TO5_SHA="e25ee074ca844ab0ec411991376ef9788625ed1d5f0b6355cde0be986fcffb13"
K1_PARQUET_SHA="0d3de3374a06b1e7b9ad0281e2ab7844632ed139db871488524061191c8f34c1"
K1_RECEIPT_SHA="e86463bc61f9eebd8cbdbe01269f767c140bf243f9e374027008e9f095eb9044"
FCP_PARQUET_SHA="ee79c60671dda65a3d2f4f9319af952dbe81ab805ce269fe32b3bb3448b0d6b8"
FCP_RECEIPT_SHA="5dbd4590ae96729cd9469fc8403a8ee4528d18e1fb5b928d5a108483ea29150b"
HISTORICAL_COUNTS={
    "main":PRIOR_MAIN,"prior_group00":PRIOR_OLD,
    "prior_R1":PRIOR_R1,"prior_R2":PRIOR_R2,
}


def prior_actual_images(history,k1,fcp, *,
                        k1_sha=K1_PARQUET_SHA,k1_receipt_sha=K1_RECEIPT_SHA,
                        fcp_sha=FCP_PARQUET_SHA,fcp_receipt_sha=FCP_RECEIPT_SHA):
    """Build read-only actual original-photo ID terminal reuse with strict SHAs."""
    import pandas as pd
    history=Path(history);k1=Path(k1);fcp=Path(fcp)
    sources=[
        ("main",history/"main"/"combined_actual_photo_checkpoint.parquet",
         history/"main"/"checkpoint_receipt.json"),
        ("prior_group00",history/"group00"/"measured_photos.parquet",
         history/"group00"/"receipt.json"),
        ("prior_R1",history/"r1"/"measured_photos.parquet",
         history/"r1"/"receipt.json"),
        ("prior_R2",history/"r2"/"measured_photos.parquet",
         history/"r2"/"receipt.json"),
    ]
    already,oldobs,source_audit,overlap,_,_=prior_original_photo_index(
        sources,HISTORICAL_COUNTS)
    require(len(already)==929138 and overlap["unique_old_photo_IDs"]==929138
            and overlap["signed_arm_photo_rows"]==929138
            and not overlap.get("old_cross_arm_same_photo_evidence_conflicts",0),
            "Historical 929138 photo ROI source ID proof changed")
    for arm,root,filename,receipt_sha,parquet_sha,expected,known in [
        ("K1",k1,"measured_photos.parquet",k1_receipt_sha,k1_sha,76753,27689),
        ("FCP_first1_new3141",fcp,"new_FCP_depth1_original_photo_ROIs.parquet",
         fcp_receipt_sha,fcp_sha,3141,1201),
    ]:
        parquet=root/filename
        receipt=root/"receipt.json"
        require(sha(receipt)==receipt_sha and sha(parquet)==parquet_sha,
                "Previously measured K1/FCP photo Parquet or receipt SHA mismatch")
        report=json.loads(receipt.read_text(encoding="utf-8"))
        if arm=="K1":
            require(report.get("attempted")==expected and report.get("classifiable")==known
                    and report.get("measured_photo_parquet_sha256")==parquet_sha,
                    "K1 completed original photo result changed")
        else:
            require(report.get("newly_true_image_measured_original_photo_attempts")==expected
                    and report.get("newly_true_image_FCP_four_hue_classified")==known
                    and report.get("actual_new_image_parquet_SHA256")==parquet_sha,
                    "FCP completed first original photo result changed")
        frame=pd.read_parquet(parquet,columns=[
            "photo_id","observation_id","morph","measurement_status","image_sha256"
        ]).fillna("")
        require(len(frame)==expected and frame.photo_id.astype(str).is_unique
                and int(frame.morph.isin(HUES).sum())==known,
                "Measured actual original photo IDs/terminal four-hue outcome changed")
        for item in frame.itertuples(index=False,name=None):
            pid,obs,morph,status,image_sha=map(str,item)
            require(pid.isdecimal() and pid not in already
                    and bool(obs) and morph in (*HUES,"mixed_uncertain")
                    and (morph in HUES)==(status=="classified_four_state_morph"),
                    "Old, K1 and FCP source original photo IDs or technical labels collided")
            already[pid]=(arm,obs,morph,status,image_sha)
    require(len(already)==1_009_032,"Exact already measured source photos not conserved")
    return already,{"prior_photo_attempts":len(already),
                    "source_old_historical_audits":source_audit,
                    "previous_original_photo_IDs_actually_measured_and_reused":True}


def freeze_new_rank2to5(source_gz,receipt_path,already,output, *,
                       expected_delta=STAGE_INCREMENTAL,
                       expected_taxa=GLOBAL_ELIGIBLE,
                       expected_source_sha=GLOBAL_RANK2TO5_SHA,
                       expected_receipt_sha=GLOBAL_STAGE5_RECEIPT_SHA):
    """Only truly new original image source candidates become photo work."""
    source_gz,receipt_path,output=map(Path,(source_gz,receipt_path,output))
    require(not output.exists(),"The FCP U100 source photo queue is immutable")
    require(sha(source_gz)==expected_source_sha
            and sha(receipt_path)==expected_receipt_sha,
            "Immutable global 80968 nested original photo source checksum changed")
    doc=json.loads(receipt_path.read_text(encoding="utf-8"))
    require(doc.get("schema")==GLOBAL_SCHEMA
            and doc.get("eligible_taxa_with_at_least100_original_observations")==expected_taxa
            and doc.get("global_incremental_rank2to5_candidate_photo_IDs_NOT_NEW_INFERENCES")==expected_delta
            and doc.get("global_first5_nested_source_photo_IDs")==DEPTH*expected_taxa
            and doc.get("global_first1_original_photo_IDs_unchanged")==expected_taxa
            and doc.get("ranks2to5_genuine_never_measured_photo_IDs_not_yet_deduplicated") is True
            and doc.get("source_probabilities_original_5_over_N_over_M_and_delta_4_over_N_over_M") is True
            and doc.get("new_image_downloads_and_original_ROI_inferences_performed")==0
            and doc["csv_outputs_SHA256"][
                "global_FCP_U100_incremental_rank2to5_candidates_NEED_DEDUP.csv.gz"]
                ==expected_source_sha,
            "Signed full93m global FCP U100 first5 photo ranks source contract changed")
    require(isinstance(already,dict),"Photo-ID complete source index missing")
    observed=set()
    selected_ids=set()
    taxa=defaultdict(set)
    new=[]
    reused=[]
    unsupported=[]
    counts=Counter()
    with gzip.open(source_gz,"rt",encoding="utf-8",newline="") as stream:
        r=csv.DictReader(stream)
        require(tuple(r.fieldnames or ())==SOURCE_FIELDS,
                "Photo rank2to5 original-source CSV data field order changed")
        for x in r:
            pid=x["photo_id"]
            obs=x["observation_id"]
            taxon=x["inat_taxon_id"]
            rank=x["U100_source_original_observation_rank"]
            require(pid.isdecimal() and taxon.isdecimal() and obs
                    and pid not in selected_ids and obs not in observed
                    and rank in {"2","3","4","5"}
                    and x["incremental_since_depth1"]=="true"
                    and x["source_arm"]==NESTED_ARM
                    and x["independent_accepted_species_and_native_stage_verified"]=="false",
                    "Duplicate source photograph/observation or nonincremental FCP rank")
            selected_ids.add(pid)
            observed.add(obs)
            require(rank not in taxa[taxon],
                    "FCP depth5 original observation rank repeated within taxon")
            taxa[taxon].add(rank)
            nobs=int(x["full_93m_distinct_original_observations_of_taxon"])
            m=int(x["photos_in_original_observation"])
            require(nobs>=100 and m>=1,"Not genuine FCP U100-eligible original taxon")
            nested=float(x["cumulative_prefix5_joint_original_photo_probability"])
            increment=float(x["rank2to5_incremental_joint_original_photo_probability"])
            prior1=float(x["depth1_joint_original_photo_inclusion_probability"])
            require(
                math.isclose(prior1,1/nobs/m,rel_tol=1e-12)
                and math.isclose(nested,5/nobs/m,rel_tol=1e-12)
                and math.isclose(increment,4/nobs/m,rel_tol=1e-12)
                and 0<increment<nested<=1,
                "FCP rank2..5 original photo selection probability changed")
            source_ext=normalized_extension(x["original_image_extension"])
            require(x["photo_source_disposition"]==(
                "supported" if source_ext else "unsupported_original_extension"),
                "FCP original image technical extension state changed")
            require(x["photo_url_large"]==(
                f"https://inaturalist-open-data.s3.amazonaws.com/photos/{pid}/large.{source_ext}"
                if source_ext else ""),
                "Original selected photo URL does not match signed public-photo ID")
            old=already.get(pid)
            if old is not None:
                prior_arm,prior_obs,morph,status,photo_sha=old
                relation=prior_source_observation_relation(prior_obs,obs)
                reused.append({
                    "photo_id":pid,"observation_id":obs,"inat_taxon_id":taxon,
                    "U100_original_rank":rank,
                    "U100_prefix5_joint_photo_probability":format(nested,".17g"),
                    "incremental_rank2to5_joint_photo_probability":format(increment,".17g"),
                    "prior_image_arm":prior_arm,
                    "prior_automatic_four_hue_morph":morph,
                    "prior_ROI_terminal_status":status,
                    "prior_observation_identity_relation":relation,
                    "prior_image_original_SHA256":photo_sha,
                })
                counts["reused_actual_genuine_original_photo_IDs_NOT_NEW"]+=1
            elif source_ext is None:
                unsupported.append({
                    "photo_id":pid,"observation_id":obs,"inat_taxon_id":taxon,
                    "U100_original_rank":rank,
                    "technical_status":"unsupported_original_image_extension",
                    "rank2to5_joint_photo_probability":format(increment,".17g"),
                })
                counts["unsupported_source_photo_pixels_NOT_WHITE"]+=1
            else:
                new.append(x)
                counts["never_previously_measured_new_original_photo_ROI_targets"]+=1
    require(len(selected_ids)==expected_delta
            and len(taxa)==expected_taxa
            and all(ranks=={"2","3","4","5"} for ranks in taxa.values())
            and len(reused)+len(unsupported)+len(new)==expected_delta,
            "FCP original stage5 source lacks four independent ranked observations per taxon")
    require(new, "No genuinely new source photo pixels left after full photo-ID reuse")
    # Blind source hash-based shuffle improves worker geography/taxon balance.
    new.sort(key=lambda x:(
        hashlib.sha256(("mois-FCP-depth5-photo-blind-shuffle-v1|"+
                        x["observation_id"]).encode()).hexdigest(),
        int(x["photo_id"])
    ))
    output.parent.mkdir(parents=True,exist_ok=True)
    stage=Path(tempfile.mkdtemp(prefix="mois-next-FCP-stage5-",dir=output.parent))
    try:
        units={}
        for unit,start in enumerate(range(0,len(new),256)):
            file=stage/"shards"/f"shard-{unit:03d}.csv"
            file.parent.mkdir(exist_ok=True)
            with file.open("w",encoding="utf-8",newline="") as stream:
                writer=csv.DictWriter(stream,fieldnames=WORKER_COLUMNS,
                                      lineterminator="\n")
                writer.writeheader()
                for x in new[start:start+256]:
                    pid=x["photo_id"]
                    measurement_id="MOIS-"+hashlib.sha256(
                        ("mois-original-global-RGB-v1|"+pid).encode()
                    ).hexdigest()[:24].upper()
                    writer.writerow({
                        "measurement_id":measurement_id,"photo_id":pid,
                        "inat_taxon_id":x["inat_taxon_id"],"species":x["species"],
                        "observation_id":x["observation_id"],
                        "observer_id":x["observer_id"],"observed_on":x["observed_on"],
                        "latitude":x["latitude"],"longitude":x["longitude"],
                        "positional_accuracy":x["positional_accuracy"],
                        "photo_license":x["photo_license"],
                        "photo_url_large":x["photo_url_large"],
                        "image_filename":measurement_id+".jpg","source_arm":NESTED_ARM,
                        "selection_probability":x["rank2to5_incremental_joint_original_photo_probability"],
                        "source_snapshot_date":"2026-09-27","shard":str(unit),
                    })
            units[str(unit)]=sha(file)
        def sidecar(name,headers,rows):
            file=stage/name
            with file.open("w",encoding="utf-8",newline="") as stream:
                writer=csv.DictWriter(stream,fieldnames=headers,lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)
            return sha(file)
        reuse_headers=("photo_id","observation_id","inat_taxon_id","U100_original_rank",
                       "U100_prefix5_joint_photo_probability",
                       "incremental_rank2to5_joint_photo_probability",
                       "prior_image_arm","prior_automatic_four_hue_morph",
                       "prior_ROI_terminal_status","prior_observation_identity_relation",
                       "prior_image_original_SHA256")
        reuse_sha=sidecar("previous_original_image_photo_ID_reuse.csv",
                           reuse_headers,reused)
        unsupported_headers=("photo_id","observation_id","inat_taxon_id",
                             "U100_original_rank","technical_status",
                             "rank2to5_joint_photo_probability")
        unsupported_sha=sidecar("unsupported_original_image_source_unknown.csv",
                                unsupported_headers,unsupported)
        manifest={
            "schema":SCHEMA,
            "full_global_upstream_original_photo_IDs_metadata_only":93209351,
            "global_source_nominal_taxa":42111,
            "U100_eligible_nominal_taxa":expected_taxa,
            "first5_independent_original_photo_observation_rank_count":5,
            "new_rank2to5_source_photo_candidates":expected_delta,
            "old_existing_global_real_original_photo_attempts":len(already),
            "new_original_photo_IDs_selected_never_previously_measured":len(new),
            "original_photo_IDs_with_existing_ROI_reused":len(reused),
            "original_photo_unsupported_formats_technical_unknown":len(unsupported),
            "source_U100_incremental_rank2to5_photo_probability_is_4_over_N_over_M":True,
            "source_U100_cumulative_prefix5_probability_is_5_over_N_over_M":True,
            "source_photo_sample_decided_without_photo_pixels_or_colour":True,
            "original_unclassifiable_ROIs_from_prior_stay_unknown":True,
            "original_images_downloaded_by_this_source_process":0,
            "new_four_hue_ROIs_inferred_by_this_source_process":0,
            "source_full_80k_incremental_candidate_SHA256":expected_source_sha,
            "source_full_global_nested5_receipt_SHA256":expected_receipt_sha,
            "worker_photo_shard_SHA256":units,
            "number_original_photo_work_units":len(units),
            "worker_group_ids":list(range((len(units)+3)//4)),
            "original_model_batch_size":256,
            "original_units_per_group":4,
            "source_worker_no_lowres_flower_exclusion":True,
            "historical_old_and_K1_FCP_real_pixel_counts_not_doubled":True,
            "photo_ID_reuse_sidecar_sha256":reuse_sha,
            "unsupported_image_unknown_sidecar_sha256":unsupported_sha,
            "true_wild_native_blooming_plant_individuals_not_verified":True,
            "probability_unbiased_global_white_flower_rate_not_claimed":True,
        }
        (stage/"manifest.json").write_text(
            json.dumps(manifest,sort_keys=True,indent=2)+"\n",encoding="utf-8")
        os.replace(stage,output)
    except BaseException:
        shutil.rmtree(stage,ignore_errors=True)
        raise
    print("MOIS_FCP_U100_NESTED5_GENUINE_NEW_ORIGINAL_ROI_TARGET_QUEUE",json.dumps({
        "source_incremental_candidates":expected_delta,
        "previous_real_photo_IDs":len(already),
        "reused_prior_original_photo_IDs":len(reused),
        "unsupported_source_original_photo_unknowns":len(unsupported),
        "never_measured_original_photo_target_IDs":len(new),
        "source_shards":len(units),"original_image_downloads":0,"colour_inferences":0,
    },sort_keys=True),flush=True)
    return manifest


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("global-gz","global-receipt","historical","k1","fcp","output"):
        p.add_argument("--"+name,type=Path,required=True)
    a=p.parse_args()
    old,source_receipt=prior_actual_images(a.historical,a.k1,a.fcp)
    result=freeze_new_rank2to5(a.global_gz,a.global_receipt,old,a.output)
    report_path=a.output/"source_real_1009032_photo_ID_prior_audit.json"
    report_path.write_text(json.dumps(source_receipt,sort_keys=True,indent=2)+"\n",
                           encoding="utf-8")
    print("MOIS_FCP_STAGE5_TRUE_SOURCE_PRIOR_ORIGINAL_SHA_AUDITED",
          len(old),flush=True)


if __name__=="__main__":
    main()
