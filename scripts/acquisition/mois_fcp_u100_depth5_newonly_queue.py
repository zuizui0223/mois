#!/usr/bin/env python3
"""Deduplicate FCP U100 rank2..5 SOURCE photo IDs against 1,009,032 real prior pixels.

Image work starts ONLY from authentic never-before-inferred supported original
photo IDs, not all 80,968 prefix5 incremental metadata candidates. Preserve
the complete nested 101,210-photo observation/taxon source including old colour
unknowns, first1+next4 inclusion probabilities and unsupported image formats.

No original photo downloaded, no new flower classifier run here.
"""
from __future__ import annotations

import argparse
from collections import Counter
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

from scripts.acquisition.mois_completed_pixel_checkpoint import HUES, require, sha
from scripts.acquisition.mois_fcp_u100_depth5_source import (
    ARM, SCHEMA as SOURCE_GROUP_SCHEMA, FIELDS as STAGE5_FIELDS,
)
from scripts.acquisition.mois_fcp_u100_depth5_global_source import (
    SCHEMA as GLOBAL_STAGE5_SCHEMA,
)
from scripts.acquisition.mois_93m_balanced_prior_ROI_reuse_bridge import (
    PRIOR_MAIN, PRIOR_OLD, PRIOR_R1, PRIOR_R2,
    prior_original_photo_index, prior_source_observation_relation,
)
from scripts.acquisition.mois_queue_from_full93m_target_tiers import (
    WORKER_COLUMNS, OUT_SCHEMA,
)

SCHEMA="mois_FCP_U100_20242_taxa_nest5_exact_new_original_photo_ROI_queue_v1"
SHA_STAGE5="037e74dd236d04c44f0d84c7b9d0b3e436e72fe6ddd470cc315c1dc56b0e197f"
SHA_SOURCE="a2df91d66748a3fd7aab0d03d0eedebd3ca73d659a3950a8c1ef3ad9cf6c7bd8"
SHA_DELTA="e25ee074ca844ab0ec411991376ef9788625ed1d5f0b6355cde0be986fcffb13"
SHA_K1="0d3de3374a06b1e7b9ad0281e2ab7844632ed139db871488524061191c8f34c1"
SHA_FCP="ee79c60671dda65a3d2f4f9319af952dbe81ab805ce269fe32b3bb3448b0d6b8"
HISTORIC=929_138
HISTORIC_HUE=416_915
NEW_K1=76_753
NEW_K1_HUE=27_689
NEW_FCP=3_141
NEW_FCP_HUE=1_201
TOTAL_PRIOR=1_009_032
TOTAL_HUES=445_805
N_TAXA=20_242
N_FIRST5=101_210
N_DELTA=80_968
PRIOR_COUNTS={
    "main":PRIOR_MAIN,"prior_group00":PRIOR_OLD,
    "prior_R1":PRIOR_R1,"prior_R2":PRIOR_R2,
}
DECISION_FIELDS=(
    "photo_id","inat_taxon_id","observation_id","original_rank","source_arm",
    "first5_photo_probability_5_over_N_M","incremental_rank2to5_probability_4_over_N_M",
    "already_measured_source_arm","already_measured_hue_or_unknown",
    "already_measured_status","observation_ID_relation_to_old_photo",
    "photo_decision",
)
WORK_BATCH=256
WORK_GROUP_UNITS=4


def prior_from_signed_sources(historical, k1, fcp):
    import pandas as pd

    paths=[
        ("main",historical/"main"/"combined_actual_photo_checkpoint.parquet",
         historical/"main"/"checkpoint_receipt.json"),
        ("prior_group00",historical/"group00"/"measured_photos.parquet",
         historical/"group00"/"receipt.json"),
        ("prior_R1",historical/"r1"/"measured_photos.parquet",
         historical/"r1"/"receipt.json"),
        ("prior_R2",historical/"r2"/"measured_photos.parquet",
         historical/"r2"/"receipt.json"),
    ]
    previous,_,sources,history,_,conflicts=prior_original_photo_index(paths,PRIOR_COUNTS)
    require(len(previous)==HISTORIC and not conflicts
            and history["unique_old_photo_IDs"]==HISTORIC
            and sum(x[2] in HUES for x in previous.values())==HISTORIC_HUE,
            "Historic source-SHA-verified original photo ID/hue baseline changed")
    known=dict(previous)
    for filepath,expected_sha,attempts,hues,arm in (
        (k1,SHA_K1,NEW_K1,NEW_K1_HUE,"K1_signed_finished_new_original_photo"),
        (fcp,SHA_FCP,NEW_FCP,NEW_FCP_HUE,"FCP_depth1_signed_new_original_photo"),
    ):
        require(sha(filepath)==expected_sha,
                "Original image finished K1/FCP photo Parquet SHA changed")
        df=pd.read_parquet(filepath,columns=[
            "photo_id","observation_id","morph","measurement_status"]).fillna("")
        require(len(df)==attempts and df.photo_id.astype(str).str.fullmatch(r"\d+").all()
                and df.photo_id.astype(str).is_unique
                and int(df.morph.isin(HUES).sum())==hues,
                "Previously processed original photo image count/hues changed")
        for pid,obs,morph,status in df.itertuples(index=False,name=None):
            pid,obs,morph,status=map(str,(pid,obs,morph,status))
            require(pid not in known and bool(obs)
                    and morph in (*HUES,"mixed_uncertain")
                    and (morph in HUES)==(status=="classified_four_state_morph"),
                    "Old, K1 and FCP identical original photo IDs or unknown status conflict")
            known[pid]=(arm,obs,morph,status,"")
    require(len(known)==TOTAL_PRIOR
            and sum(x[2] in HUES for x in known.values())==TOTAL_HUES,
            "Photo-ID-authenticated FCP/K1/history full one-million old source changed")
    return known,sources


def reconcile(prefix5,delta,receipt_path,known,output, *,
              expected_source_sha=SHA_SOURCE, expected_delta_sha=SHA_DELTA,
              expected_receipt_sha=SHA_STAGE5,
              expected_taxa=N_TAXA, expected_first5=N_FIRST5,
              expected_delta=N_DELTA, expected_prior=TOTAL_PRIOR):
    output=Path(output)
    require(not output.exists(),
            "Cannot replace previously decided nested5 original-image photo source")
    require(sha(prefix5)==expected_source_sha and sha(delta)==expected_delta_sha
            and sha(receipt_path)==expected_receipt_sha,
            "Signed original nested FCP5 prefix+incremental source SHA differs")
    manifest=json.loads(Path(receipt_path).read_text(encoding="utf-8"))
    require(manifest.get("schema")==GLOBAL_STAGE5_SCHEMA
            and manifest.get("global_first5_nested_source_photo_IDs")==expected_first5
            and manifest.get("global_incremental_rank2to5_candidate_photo_IDs_NOT_NEW_INFERENCES")==expected_delta
            and manifest.get("eligible_taxa_with_at_least100_original_observations")==expected_taxa
            and manifest.get("ranks2to5_genuine_never_measured_photo_IDs_not_yet_deduplicated") is True
            and manifest.get("csv_outputs_SHA256",{}).get(
                "global_FCP_U100_nested_first5_original_photo_candidates.csv.gz")
                ==expected_source_sha
            and manifest.get("csv_outputs_SHA256",{}).get(
                "global_FCP_U100_incremental_rank2to5_candidates_NEED_DEDUP.csv.gz")
                ==expected_delta_sha
            and manifest.get("new_image_downloads_and_original_ROI_inferences_performed")==0,
            "Immutable worldwide FCP-U100 depth5 source totals changed")
    require(len(known)==expected_prior,
            "Genuine prior original image pixels missing from source reuse gate")

    source_rows={}
    per_taxon={}
    ranks=Counter()
    first1=set()
    seen_obs=set()
    with gzip.open(prefix5,"rt",encoding="utf-8",newline="") as f:
        reader=csv.DictReader(f)
        require(set(("photo_id","observation_id","inat_taxon_id",
                     "U100_source_original_observation_rank",
                     "source_rank1_original_photo_id",
                     "cumulative_prefix5_joint_original_photo_probability",
                     "rank2to5_incremental_joint_original_photo_probability",
                     "photo_source_disposition","photo_url_large","source_arm"))
                <=set(reader.fieldnames or ()),
                "Frozen original FCP depth5 photo/observation inclusion source incomplete")
        for r in reader:
            pid,obs,taxon=r["photo_id"],r["observation_id"],r["inat_taxon_id"]
            rank=int(r["U100_source_original_observation_rank"])
            require(pid.isdecimal() and bool(obs) and taxon.isdecimal()
                    and rank in (1,2,3,4,5)
                    and pid not in source_rows and obs not in seen_obs
                    and r["source_arm"]==ARM,
                    "Nested5 source photo, independent observation or rank repeated")
            n=int(r["full_93m_distinct_original_observations_of_taxon"])
            m=int(r["photos_in_original_observation"])
            p=5/n/m
            require(n>=100 and m>=1
                    and math.isclose(float(r["cumulative_prefix5_joint_original_photo_probability"]),
                                     p,rel_tol=1e-12)
                    and math.isclose(float(r["rank2to5_incremental_joint_original_photo_probability"]),
                                     4/n/m,rel_tol=1e-12),
                    "FCP original photograph nested5 probability changed")
            if rank==1:
                require(pid in known,"Already measured FCP U100 first photo unexpectedly unavailable")
                first1.add(taxon)
            source_rows[pid]=r
            per_taxon.setdefault(taxon,[]).append(r)
            ranks[taxon]+=1
            seen_obs.add(obs)
    require(len(source_rows)==expected_first5 and len(ranks)==expected_taxa
            and len(first1)==expected_taxa and all(x==5 for x in ranks.values()),
            "World FCP U100 nested5 original observation-photo/taxon coverage drift")
    original_delta={}
    with gzip.open(delta,"rt",encoding="utf-8",newline="") as f:
        reader=csv.DictReader(f)
        require(set(("photo_id","observation_id","inat_taxon_id",
                     "U100_source_original_observation_rank"))<=set(reader.fieldnames or ()),
                "Original incremental source photo/observation grain missing")
        for row in reader:
            pid=row["photo_id"]
            require(pid in source_rows and pid not in original_delta
                    and source_rows[pid]==row
                    and int(row["U100_source_original_observation_rank"]) in (2,3,4,5),
                    "Wrong nested first5 original photo in rank2to5 incremental source")
            original_delta[pid]=row
    require(len(original_delta)==expected_delta and len(set(original_delta))==expected_delta,
            "Expected complete candidate source ranks2to5 for 20242 taxa")
    output.parent.mkdir(parents=True,exist_ok=True)
    staged=Path(tempfile.mkdtemp(prefix="mois-fcp-depth5-onlynew-",dir=output.parent))
    try:
        decisions=staged/"original_FCP5_all_first5_photo_reuse_and_new_dispositions.csv.gz"
        fresh=[]
        unsupported=[]
        counts=Counter()
        eligible_source_taxa=set()
        with decisions.open("wb") as raw:
            with gzip.GzipFile(fileobj=raw,mode="wb",filename="",mtime=0) as compressed:
                with io.TextIOWrapper(compressed,encoding="utf-8",newline="") as f:
                    writer=csv.DictWriter(f,fieldnames=DECISION_FIELDS,lineterminator="\n")
                    writer.writeheader()
                    for taxon in sorted(ranks,key=int):
                        group=sorted(per_taxon[taxon],
                            key=lambda row:int(row["U100_source_original_observation_rank"]))
                        # Every original observation remains in its hash-ranked
                        # source position; no replacing unknown images with bloom photos.
                        first=group[0]["photo_id"]
                        require(all(x["source_rank1_original_photo_id"]==first for x in group),
                                "Nested FCP original photo hash-ranked first1 changed")
                        eligible_source_taxa.add(taxon)
                        for row in group:
                            pid=row["photo_id"]
                            prev=known.get(pid)
                            oldarm=oldmorph=oldstatus=relation=""
                            if prev is not None:
                                oldarm,oldobs,oldmorph,oldstatus,_=prev
                                relation=prior_source_observation_relation(oldobs,row["observation_id"])
                                action="prior_original_photo_exact_ID_reuse"
                            elif row["photo_source_disposition"]!="supported":
                                require(not row["photo_url_large"],
                                        "Unsupported original image supplied an apparent URL")
                                action="original_photo_unsupported_extension_unknown"
                                unsupported.append(row)
                            else:
                                from scripts.acquisition.mois_onepass_full93m_target_tiers import normalized_extension
                                ext=normalized_extension(row["original_image_extension"])
                                require(ext is not None and row["photo_url_large"]==
                                        f"https://inaturalist-open-data.s3.amazonaws.com/photos/{pid}/large.{ext}",
                                        "Previously unmeasured original image lacks immutable S3 source")
                                action="new_original_photo_ROI_to_infer"
                                fresh.append(row)
                            counts[action]+=1
                            writer.writerow({
                                "photo_id":pid,
                                "inat_taxon_id":taxon,
                                "observation_id":row["observation_id"],
                                "original_rank":row["U100_source_original_observation_rank"],
                                "source_arm":ARM,
                                "first5_photo_probability_5_over_N_M":
                                    row["cumulative_prefix5_joint_original_photo_probability"],
                                "incremental_rank2to5_probability_4_over_N_M":
                                    row["rank2to5_incremental_joint_original_photo_probability"],
                                "already_measured_source_arm":oldarm,
                                "already_measured_hue_or_unknown":oldmorph,
                                "already_measured_status":oldstatus,
                                "observation_ID_relation_to_old_photo":relation,
                                "photo_decision":action,
                            })
        require(len(eligible_source_taxa)==expected_taxa
                and counts["prior_original_photo_exact_ID_reuse"]>=expected_taxa
                and counts["new_original_photo_ROI_to_infer"]+counts["original_photo_unsupported_extension_unknown"]+
                    counts["prior_original_photo_exact_ID_reuse"]==expected_first5,
                "FCP depth5 photo reuse/unknown/inference decision no longer conserves source")
        new=sorted(fresh,key=lambda r:(
            hashlib.sha256(("mois-fcp-u100-nested5-work-shuffle-v1|"+r["observation_id"]).encode()).digest(),
            int(r["photo_id"])))
        folder=staged/"shards";folder.mkdir()
        hashes={}
        photo_shard_unique=set()
        for part in range((len(new)+WORK_BATCH-1)//WORK_BATCH):
            file=folder/f"shard-{part:03d}.csv"
            with file.open("w",newline="",encoding="utf-8") as f:
                w=csv.DictWriter(f,fieldnames=WORKER_COLUMNS,lineterminator="\n")
                w.writeheader()
                for row in new[part*WORK_BATCH:(part+1)*WORK_BATCH]:
                    pid=row["photo_id"]
                    require(pid not in known and pid not in photo_shard_unique,
                            "Already inferred original photo ID inadvertently scheduled twice")
                    photo_shard_unique.add(pid)
                    mid="MOIS-"+hashlib.sha256(
                        ("mois-original-global-RGB-v1|"+pid).encode()).hexdigest()[:24].upper()
                    w.writerow({
                        "measurement_id":mid,"photo_id":pid,
                        "inat_taxon_id":row["inat_taxon_id"],"species":row["species"],
                        "observation_id":row["observation_id"],
                        "observer_id":row["observer_id"],"observed_on":row["observed_on"],
                        "latitude":row["latitude"],"longitude":row["longitude"],
                        "positional_accuracy":row["positional_accuracy"],
                        "photo_license":row["photo_license"],"photo_url_large":row["photo_url_large"],
                        "image_filename":mid+".jpg","source_arm":ARM,
                        "selection_probability":row["cumulative_prefix5_joint_original_photo_probability"],
                        "source_snapshot_date":"2026-09-27","shard":str(part),
                    })
            hashes[str(part)]=sha(file)
        require(len(photo_shard_unique)==len(fresh),"New FCP5 original photo work units missing IDs")
        source_manifest={
            "schema":OUT_SCHEMA,
            "original_sampling_schema":SCHEMA,
            "full_frame_photo_ids":93_209_351,
            "source_arm":ARM,
            "source_group":"FCP_U100_nested5_rank2to5_only_NOMINAL_ELIGIBLE",
            "colour_measured":False,
            "new_photo_attempts_planned":len(fresh),
            "original_photo_technical_unsupported_new_photos_not_downloaded":len(unsupported),
            "previously_measurement_outcomes_reused_without_inference":
                counts["prior_original_photo_exact_ID_reuse"],
            "original_first5_photo_source_IDs":expected_first5,
            "original_incremental_rank2to5_candidate_photo_IDs":expected_delta,
            "source_original_eligible_taxa":expected_taxa,
            "original_photo_nested5_source_receipt_SHA256":expected_receipt_sha,
            "previously_measured_historical_K1_FCP_1009032_photo_IDs_accounted":expected_prior,
            "selection_probability":None,
            "individual_original_photo_nested_prefix_inclusion_probabilities_preserved":True,
            "pending_image_work_is_only_source_photo_IDs_not_previously_inferred":True,
            "old_unknown_ROI_not_treated_as_flower_colour":True,
            "new_image_downloads_or_ROI_inferences_by_this_source_step":0,
            "source_metadata_wild_native_bloom_not_independently_verified":True,
            "photo_shard_sha256":hashes,
            "work_units":len(hashes),
            "worker_ids":list(range(len(hashes))),
            "group_ids":list(range((len(hashes)+WORK_GROUP_UNITS-1)//WORK_GROUP_UNITS)),
            "units_per_group":WORK_GROUP_UNITS,
            "original_model_batch_size":WORK_BATCH,
            "technical_unknowns_preserved":True,
            "full_old_K1_FCP_reuse_photo_count_conservation":
                counts["prior_original_photo_exact_ID_reuse"]+counts["new_original_photo_ROI_to_infer"]+
                counts["original_photo_unsupported_extension_unknown"],
            "photo_dispositions":dict(counts),
            "all_source_photo_dispositions_SHA256":sha(decisions),
        }
        (staged/"manifest.json").write_text(
            json.dumps(source_manifest,sort_keys=True,indent=2)+"\n",encoding="utf-8")
        os.replace(staged,output)
    except BaseException:
        shutil.rmtree(staged,ignore_errors=True)
        raise
    print("MOIS_REAL_NEW_FCP_U100_DEPTH5_ORIGINAL_IMAGE_SOURCE_QUEUE",json.dumps({
        "source_first5_nominal_photo_candidates":expected_first5,
        "rank2to5_new_source_candidates_before_pixel_reuse":expected_delta,
        "old_exact_photo_IDs_reused":counts["prior_original_photo_exact_ID_reuse"],
        "technical_unknown_unsupported_original_images":len(unsupported),
        "new_original_photo_IDs_to_download_and_infer":len(fresh),
        "new_original_photo_inferences_performed":0,
        "original_image_downloads":0,
    },sort_keys=True),flush=True)
    return source_manifest


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ("prefix5","delta","stage5-receipt","historical","k1","fcp","output"):
        parser.add_argument("--"+key,required=True,type=Path)
    a=parser.parse_args()
    old,_=prior_from_signed_sources(a.historical,a.k1,a.fcp)
    reconcile(a.prefix5,a.delta,a.stage5_receipt,old,a.output)


if __name__=="__main__":
    main()
