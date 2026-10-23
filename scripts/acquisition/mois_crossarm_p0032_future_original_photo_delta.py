#!/usr/bin/env python3
"""Append only certified p0032 original ROI groups133..144 to cross-arm registry.

Prior 910542-photo cross-arm certificate combines the immutable historic
six-arm pixel frame, tier06 group00, first all16/bucket0..3, p0032 batch001
groups121..132, and wave2 groups0..1. Here ONLY newly completed original
p0032 groups133..144 are read; verify SHA/true photo terminal and reject any
photo or observation collision with *other original sampling arms*, while
allowing the same observation in different compute groups of p0032 itself.

No historical 0.9m-row Parquet rewrite, no old source image redownload.
"""
from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path
import re

from scripts.acquisition.mois_completed_pixel_checkpoint import sha, require, HUES
from scripts.acquisition.mois_append_p0032_postcheckpoint_source_worker import (
    SCHEMA as P32_SCHEMA, RUN as P32_RUN,
    SOURCE_SHA as P32_SOURCE_SHA, BEGIN_GROUP, END_GROUP,
)

SCHEMA="mois_original_photo_crossarm_910542_plus_p0032_postdelta_source_SHA_v1"
OLD_SCHEMA="mois_crossarm_r1_r2_p0032_original_photo_source_arm_incremental_v1"
OLD_TOTAL=910542
OLD_HUES=409768
ALREADY_P32_DELTA=12288


def identity(path,expected=None):
    import pandas as pd
    d=pd.read_parquet(path,columns=["photo_id","observation_id"])
    ids=d.photo_id.astype(str)
    obs=d.observation_id.fillna("").astype(str)
    require(ids.str.fullmatch(r"[0-9]+").all()
            and ids.is_unique and obs.str.strip().ne("").all(),
            "Original photo or observation source-ID inventory corrupted")
    if expected is not None:
        require(len(ids)==expected,"Prior known original-photo source arm count changed")
    return set(ids),set(obs)


def validate(prior_registry, p32delta, early, broad, wave2_groups, new_groups, output,
             *, prior_total=OLD_TOTAL, prior_hues=OLD_HUES,
             p32_delta_n=ALREADY_P32_DELTA,
             expected_early=4138, expected_broad=4468):
    import pandas as pd

    output=Path(output)
    require(not output.exists(),"Immutable source photo cross-arm registry cannot be overwritten")
    prior=json.loads(Path(prior_registry).read_text(encoding="utf-8"))
    require(prior.get("schema")==OLD_SCHEMA
            and prior.get("full_distinct_original_photo_attempts")==prior_total
            and prior.get("full_automatically_classified_four_hue_original_photo_ROIs")==prior_hues
            and prior.get("prior_separate_arms_and_new_wave2_no_duplicate_photo_or_observation") is True
            and prior.get("original_images_fetched_and_ROI_inferred_during_this_ledger")==0
            and prior.get("all_wave2_group_photo_and_source_coverage_SHA_verified") is True,
            "Prior original-image crossarm registry was not completely source/photo-certified")
    require(new_groups,"No genuinely new p0032 original photo group Parquet released")
    early_id,early_obs=identity(early,expected_early)
    broad_id,broad_obs=identity(broad,expected_broad)
    p32_id,p32_obs=identity(p32delta,p32_delta_n)
    other_id=early_id|broad_id
    other_obs=early_obs|broad_obs
    # All old tier06 original photographs must be independently photo-disjoint
    # from old p0032 delta; history was already completely verified.
    require(not other_id.intersection(p32_id)
            and not other_obs.intersection(p32_obs),
            "Signed prior tier06 and p0032 original source arms unexpectedly overlap")

    actual_wave2_groups = {int(x["group"]):x for x in prior["per_group_original_photo_source_audit"]}
    require(len(actual_wave2_groups)==len(wave2_groups)
            and set(actual_wave2_groups)==set(prior["wave2_groups_completed"]),
            "Signed prior tier06 wave2 work groups missing/repeated")
    wave2_photo=set()
    wave2_obs=set()
    for folder in wave2_groups:
        path=Path(folder)
        r=json.loads((path/"original_photo_group_receipt.json").read_text(encoding="utf-8"))
        group=r.get("original_source_warm_group")
        require(type(group) is int and group in actual_wave2_groups
                and sha(path/"original_photo_group_receipt.json") ==
                    actual_wave2_groups[group]["group_source_SHA256"]
                and sha(path/"original_photo_ROIs.parquet") ==
                    actual_wave2_groups[group]["group_photo_Parquet_SHA256"],
                "Prior signed wave2 original-photo source worker receipt changed")
        ids,obs=identity(path/"original_photo_ROIs.parquet",
                         actual_wave2_groups[group]["original_photo_attempts"])
        require(not wave2_photo.intersection(ids)
                and not other_id.intersection(ids)
                and not other_obs.intersection(obs),
                "Wave2 certified photo IDs or independent source observation identities changed")
        wave2_photo.update(ids)
        wave2_obs.update(obs)
    require(len(wave2_photo)==prior["wave2_original_photo_attempts"],
            "Prior wave2 genuine original pixel counts lost")
    other_id.update(wave2_photo)
    other_obs.update(wave2_obs)

    n=Counter();groups_seen=set();old_original_p32_ids=p32_id.copy()
    new_photo=set();new_obs=set();out=[]
    for group_dir in new_groups:
        group_dir=Path(group_dir)
        f=group_dir/"original_photo_ROIs.parquet"
        r_path=group_dir/"p0032_original_group_receipt.json"
        art=group_dir/"original_worker_artifact_zip_identity.json"
        c=group_dir/"original_photo_hue_counts.csv"
        rec=json.loads(r_path.read_text(encoding="utf-8"))
        att=json.loads(art.read_text(encoding="utf-8"))
        group=rec.get("original_frozen_source_worker_group")
        require(type(group) is int and BEGIN_GROUP<=group<END_GROUP
                and group not in groups_seen
                and rec.get("schema")==P32_SCHEMA
                and rec.get("signed_original_source_run")==P32_RUN
                and rec.get("source_manifest_sha256")==P32_SOURCE_SHA
                and rec.get("original_source_photo_opportunities")==148071
                and rec.get("frozen_original_fcp_roi_v4_attested") is True
                and rec.get("original_photo_and_terminal_SHA_certified") is True
                and rec.get("new_S3_photo_downloads_or_model_reinferences")==0
                and rec.get("other_prior_original_worker_ZIPs_downloaded")==0,
                "Unrecognized/tampered p0032 future original photo work group")
        groups_seen.add(group)
        expected_units=list(range(group*4,min((group+1)*4,579)))
        require(rec.get("original_completed_source_units")==expected_units,
                "Original 256-image p0032 source work unit roster altered")
        require(sha(f)==rec["actual_output_sha256"][f.name]
                and sha(c)==rec["actual_output_sha256"][c.name]
                and att.get("original_p0032_group")==group
                and att.get("source_run")==P32_RUN
                and att.get("original_source_SHA256")==P32_SOURCE_SHA
                and att.get("original_terminal_group_receipt_SHA256")==sha(r_path)
                and att.get("artifact_name")==f"mois-93m-p0032-RGB-group-{group}"
                and re.fullmatch(r"sha256:[0-9a-f]{64}",
                                 str(att.get("artifact_zip_sha256") or "")) is not None,
                "Original new p0032 worker ZIP/dataset/receipt SHA does not match frozen evidence")
        columns=["photo_id","observation_id","morph","measurement_status",
                 "source_arm","selection_probability"]
        d=pd.read_parquet(f,columns=columns)
        ids=d.photo_id.astype(str);obs=d.observation_id.fillna("").astype(str)
        classified=d.morph.isin(HUES)
        fail=d.measurement_status.eq("image_acquisition_failed")
        require(ids.str.fullmatch(r"[0-9]+").all()
                and ids.is_unique and obs.str.strip().ne("").all()
                and not new_photo.intersection(ids)
                and not old_original_p32_ids.intersection(ids)
                and not other_id.intersection(ids)
                and not other_obs.intersection(obs)
                and d.morph.isin((*HUES,"mixed_uncertain")).all()
                and classified.eq(d.measurement_status.eq("classified_four_state_morph")).all()
                and pd.to_numeric(d.selection_probability,errors="coerce").between(
                    0,1,inclusive="right").all(),
                "New p0032 authentic photo IDs/observations clash or hue/propensities changed")
        require(len(d)==rec.get("original_true_photo_attempts")
                and int(classified.sum())==rec.get("original_FCP_roi_four_hue_classified")
                and int(fail.sum())==rec.get("original_photo_image_acquisition_failed")
                and len(set(obs))==rec.get("original_distinct_observations")
                and int(classified.sum())+int(fail.sum())+
                    rec.get("original_acquired_unclassifiable")==len(d)
                and all(rec["hue_counts"].get(h)==int((d.morph==h).sum())
                        for h in (*HUES,"mixed_uncertain")),
                "New source-original p0032 photo count or per-hue terminal receipts changed")
        new_photo.update(ids)
        new_obs.update(obs)
        n["photos"]+=len(d)
        n["four_hue"]+=int(classified.sum())
        n["acquisition_failed"]+=int(fail.sum())
        out.append({
            "original_group":group,
            "original_photo_attempts":len(d),
            "four_hue_classified":int(classified.sum()),
            "original_worker_artifact_ZIP_identity_SHA256":sha(art),
            "original_photo_parquet_SHA256":sha(f),
            "source_terminal_receipt_SHA256":sha(r_path),
        })
    out.sort(key=lambda x:x["original_group"])
    require(len(out)==len(new_groups),"New group duplicated/replayed in original photo crossarm")
    report={
        "schema":SCHEMA,
        "previous_original_photo_crossarm_verified_SHA256":sha(prior_registry),
        "historical_authenticated_unique_photo_IDs":prior_total,
        "historical_authenticated_four_hue_classified":prior_hues,
        "new_p0032_postdelta_original_worker_groups":sorted(groups_seen),
        "new_p0032_postdelta_authentic_original_photo_attempts":n["photos"],
        "new_p0032_postdelta_four_hue_classified":n["four_hue"],
        "new_p0032_original_acquisition_failed":n["acquisition_failed"],
        "new_p0032_distinct_observation_IDs":len(new_obs),
        "historical_and_new_genuine_original_photo_attempts":prior_total+n["photos"],
        "historical_and_new_four_hue_classifications":prior_hues+n["four_hue"],
        "new_original_photo_ID_overlap_with_old_photo_source_arms":0,
        "new_original_observation_ID_overlap_with_independent_tier06_arms":0,
        "original_p0032_photo_and_source_manifest_unique_across_group_boundaries":True,
        "historical_prior_p0032_group_artifacts_immutable_SHA_snapshot_verified":True,
        "new_image_downloads_or_ROI_colour_inference_performed":False,
        "large_historical_pixel_Parquet_rewritten":False,
        "mixed_sampling_arms_not_unweighted_global_white_flower_frequency":True,
        "taxon_and_wild_flowering_stage_not_independently_verified":True,
        "new_group_source_SHA_attestations":out,
    }
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print("MOIS_GENUINE_P0032_NEW_PHOTO_GROUPS_CROSSARM_SHA_VERIFIED",
          json.dumps({k:report[k] for k in (
              "new_p0032_postdelta_original_worker_groups",
              "new_p0032_postdelta_authentic_original_photo_attempts",
              "new_p0032_postdelta_four_hue_classified",
              "historical_and_new_genuine_original_photo_attempts",
              "historical_and_new_four_hue_classifications")},sort_keys=True),flush=True)
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("prior-registry","p32delta","early","broad","output"):
        p.add_argument("--"+key,type=Path,required=True)
    p.add_argument("--wave2",type=Path,action="append",required=True)
    p.add_argument("--new-p32",type=Path,action="append",required=True)
    a=p.parse_args()
    validate(a.prior_registry,a.p32delta,a.early,a.broad,
             a.wave2,a.new_p32,a.output)


if __name__=="__main__":
    main()
