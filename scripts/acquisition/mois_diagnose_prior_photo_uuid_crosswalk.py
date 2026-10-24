#!/usr/bin/env python3
"""Diagnose exact-photo historical observation IDs against signed 93m UUIDs.

**Read only** 80,035 original-photo probability-sample rows and previously
measured source-authenticated photo/observation Parquets. Check original
photo-ID overlap and distinguish truly missing historical observation fields,
iNaturalist numeric observation IDs, UUIDs and disagreeing genuine UUIDs.

No original images downloaded, no ROI inferred, no colours changed, no source
candidate photo dropped. Raw photo/observation identifiers are not published in
the diagnostic; small example hashes and format signatures suffice to design
a canonical crosswalk with original provenance.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
import re
from pathlib import Path

from scripts.acquisition.mois_onepass_full93m_target_tiers import require, sha256

SCHEMA="mois_signed_original_photo_observation_UUID_crosswalk_diagnosis_v1"
SOURCE_SELECTED_SHA="b282295807fd48125ce5646385b534651da83fa7b0102634deb2d9941359ea55"
SOURCE_MANIFEST_SHA="c0ba3a316a336816065590540e75034bd8c2b9b7362034163673702103c44977"
SOURCE_SELECTED_OBS=80035
SOURCE_PHOTO_OPPORTUNITIES=93_209_351
SOURCE_ARMS=("main","group00","R1","R2")
REQUIRED=("photo_id","observation_id")
MISSING=frozenset(("", "nan","none","null","na","n/a"))
BASE_COUNTS={"main":(911767,410323),"group00":(4138,1573),
             "R1":(4468,1760),"R2":(8765,3259)}
UUID=re.compile(r"^[0-9a-fA-F]{8}-(?:[0-9a-fA-F]{4}-){3}[0-9a-fA-F]{12}$")
HEX32=re.compile(r"^[0-9a-fA-F]{32}$")
HEX64=re.compile(r"^[0-9a-fA-F]{64}$")


def field_type(value):
    s=str(value or "").strip()
    if s.lower() in MISSING:
        return "missing"
    if UUID.fullmatch(s):return "uuid"
    if s.isdecimal():return "numeric"
    if HEX32.fullmatch(s):return "hex32"
    if HEX64.fullmatch(s):return "hex64"
    return "other"


def hashed_example(photo_id,source_obs,prior_obs):
    return {
        "photo_hash12":hashlib.sha256(str(photo_id).encode()).hexdigest()[:12],
        "source_observation_format":field_type(source_obs),
        "prior_observation_format":field_type(prior_obs),
        "source_observation_length":len(str(source_obs)),
        "prior_observation_length":len(str(prior_obs)),
        "same_value_after_strip":str(source_obs).strip()==str(prior_obs).strip(),
    }


def probe(sample,manifest,prior_arms,output,*,expected_source=SOURCE_SELECTED_OBS,
          expected_counts=None,expected_sample_sha=SOURCE_SELECTED_SHA,
          expected_manifest_sha=SOURCE_MANIFEST_SHA):
    import pandas as pd

    sample,manifest,output=map(Path,(sample,manifest,output))
    require(not output.exists(),"Observation-ID photo overlap audit is immutable")
    require(sha256(sample)==expected_sample_sha and
            sha256(manifest)==expected_manifest_sha,
            "Input 93m signed original-photo sample changed")
    info=json.loads(manifest.read_text())
    require(info["schema"]=="mois_full_source_observation_sampling_cohort_v1"
            and info["sampled_distinct_original_observations"]==expected_source
            and info["all_upstream_original_photo_IDs"]==SOURCE_PHOTO_OPPORTUNITIES
            and info["source_sample_decided_without_photo_pixels"] is True,
            "Input is not original FCP-style 93m precolour probability sample")
    source=pd.read_csv(sample,usecols=["photo_id","observation_id","inat_taxon_id"],dtype=str,
                       keep_default_na=False)
    require(len(source)==expected_source and source.photo_id.is_unique
            and source.observation_id.is_unique
            and source.photo_id.str.fullmatch(r"\d+").all()
            and source.observation_id.str.strip().ne("").all(),
            "Input original sample photo/observation UUIDs duplicated or empty")
    expected_counts=expected_counts if expected_counts is not None else BASE_COUNTS
    require(set(prior_arms)==set(SOURCE_ARMS) and
            set(expected_counts)==set(SOURCE_ARMS),
            "Must audit all four signed previously attempted source arms")
    stats={}
    examples=[]
    source_index=source.set_index("photo_id")
    for arm in SOURCE_ARMS:
        path,receipt=map(Path,prior_arms[arm])
        old=json.loads(receipt.read_text())
        n,h=expected_counts[arm]
        if arm=="main":
            require(old.get("schema")=="mois_sealed_actual_pixel_checkpoint_v1"
                    and old.get("combined_unique_photo_attempts")==n
                    and old.get("combined_four_hue_classified")==h
                    and old["output_sha256"][path.name]==sha256(path),
                    "Old six-arm true original pixel/photo count or source SHA changed")
        else:
            require(old.get("schema")=="mois_global93m_observation_baseline_real_RGB_photo_results_v1"
                    and old.get("attempted")==n and old.get("classifiable")==h
                    and old.get("measured_photo_parquet_sha256")==sha256(path),
                    "Old signed tier06 original photo/colour source changed")
        prior=pd.read_parquet(path,columns=list(REQUIRED))
        prior["photo_id"]=prior["photo_id"].astype(str)
        prior["observation_id"]=prior["observation_id"].fillna("").astype(str)
        require(len(prior)==n and prior.photo_id.is_unique
                and prior.photo_id.str.fullmatch(r"\d+").all(),
                "Previously measured original photo IDs duplicated")
        matches=prior[prior.photo_id.isin(source_index.index)].copy()
        matches["source_obs"]=source_index.reindex(matches.photo_id).observation_id.to_numpy()
        found=Counter()
        local_examples=[]
        for pid,old_obs,source_obs in zip(matches.photo_id,
                                         matches.observation_id,matches.source_obs):
            oldkind,newkind=field_type(old_obs),field_type(source_obs)
            equal=str(old_obs).strip()==str(source_obs).strip()
            category=("same" if equal else
                      "prior_missing" if oldkind=="missing" else
                      f"different_{oldkind}_versus_{newkind}")
            found[category]+=1
            found[f"old_format_{oldkind}"]+=1
            found[f"source_format_{newkind}"]+=1
            if not equal and len(local_examples)<12:
                local_examples.append(hashed_example(pid,source_obs,old_obs))
        require(sum(v for k,v in found.items() if k=="same" or k=="prior_missing"
                    or k.startswith("different_"))==len(matches),
                "Original photo overlap categories failed count conservation")
        stats[arm]={
            "previous_original_photo_attempts":n,
            "previous_four_hue_classified":h,
            "photo_ID_overlap_with_new_probability_sample":len(matches),
            "exact_same_observation_ID":found["same"],
            "historical_missing_observation_ID":found["prior_missing"],
            "historical_actual_ID_disagreement":
                len(matches)-found["same"]-found["prior_missing"],
            "identity_format_counts":dict(sorted(found.items())),
            "prior_original_photo_Parquet_SHA256":sha256(path),
            "prior_original_receipt_SHA256":sha256(receipt),
            "hashed_type_examples":local_examples,
        }
        examples.extend({"arm":arm,**example} for example in local_examples[:3])
    total_matches=sum(x["photo_ID_overlap_with_new_probability_sample"] for x in stats.values())
    report={
        "schema":SCHEMA,
        "input_full_93m_photo_ID_opportunities":SOURCE_PHOTO_OPPORTUNITIES,
        "source_93m_original_sample_distinct_observations":expected_source,
        "source_sample_SHA256":sha256(sample),
        "source_manifest_SHA256":sha256(manifest),
        "total_per_arm_overlapping_photo_rows_including_any_cross_arm_duplicates":
            total_matches,
        "per_arm_diagnosis":stats,
        "hashed_format_examples":examples[:12],
        "no_original_photo_images_downloaded":True,
        "no_FCP_ROI_hue_inference_performed":True,
        "no_original_photo_IDs_or_early_hue_labels_mutated":True,
        "ambiguity_requires_canonical_observation_crosswalk_or_photo_only_technical_unknown":True,
        "cannot_reuse_disagreeing_genuine_observation_as_same_independent_observation":True,
    }
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(report,sort_keys=True,indent=2)+"\n")
    print("MOIS_FCP_GLOBAL93M_PRIOR_ORIGINAL_OBSERVATION_ID_DIAGNOSIS",
          json.dumps({key:{
              "overlap":value["photo_ID_overlap_with_new_probability_sample"],
              "same":value["exact_same_observation_ID"],
              "missing":value["historical_missing_observation_ID"],
              "disagreement":value["historical_actual_ID_disagreement"],
              "formats":{k:v for k,v in value["identity_format_counts"].items()
                         if k.startswith("different_")}
          } for key,value in stats.items()},sort_keys=True),flush=True)
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("sample","manifest","main-parquet","main-receipt",
                 "group00-parquet","group00-receipt","r1-parquet","r1-receipt",
                 "r2-parquet","r2-receipt","output"):
        p.add_argument("--"+name,type=Path,required=True)
    a=p.parse_args()
    probe(a.sample,a.manifest,
          {"main":(a.main_parquet,a.main_receipt),
           "group00":(a.group00_parquet,a.group00_receipt),
           "R1":(a.r1_parquet,a.r1_receipt),
           "R2":(a.r2_parquet,a.r2_receipt)},a.output)


if __name__=="__main__":
    main()
