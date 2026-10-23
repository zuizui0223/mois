#!/usr/bin/env python3
"""Seal one new, post-checkpoint p0032 worker as immutable original-photo Parquet.

Historical 887600-photo main source pixels + the signed 12288-photo p0032
batch001 already cover original worker groups 0..132 (source run 39401484641).
Only as-yet unpublished original groups 133..144 may be appended here.
The original FCP YOLO11/EfficientSAM ROI-v4 method/terminal image bytes are
validated from stored original worker artifacts, WITHOUT opening S3 photos,
re-running the detector or recomputing prior biological flower-colour labels.

No independent flower-specific anthesis/wild/native status is established.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import json
import os
from pathlib import Path
import shutil
import tempfile

from scripts.acquisition.mois_completed_pixel_checkpoint import (
    WAVES, HUES, sha, require, collect_wave,
)
from scripts.analysis.mois_reassemble_original_new4096_measurements import METHOD

SCHEMA="mois_p0032_post_checkpoint_immutable_original_ROI_worker_group_v1"
RUN,EXPECTED_PHOTOS,SOURCE_SHA=WAVES["p0032"]
BEGIN_GROUP=133
END_GROUP=145


def verify_frozen_original_artifacts(earlier_snapshot, live_snapshot):
    """Require all previously processed true p0032 worker ZIPs stay unchanged."""
    prior=json.loads(Path(earlier_snapshot).read_text(encoding="utf-8"))
    live=json.loads(Path(live_snapshot).read_text(encoding="utf-8"))
    require(prior.get("source_run_id")==RUN and live.get("source_run_id")==RUN,
            "p0032 original source run changed")
    def canonical(d):
        assets={}
        for x in d["artifacts"]:
            name=x.get("name","")
            pref="mois-93m-p0032-RGB-group-"
            require(name.startswith(pref) and name[len(pref):].isdecimal()
                    and 0 <= int(name[len(pref):]) < END_GROUP
                    and isinstance(x.get("digest"),str)
                    and x["digest"].startswith("sha256:")
                    and len(x["digest"]) == 71
                    and isinstance(x.get("id"),int) and x["id"] > 0
                    and x.get("expired") is not True,
                    "Unrecognized original-photo worker ZIP SHA inventory")
            group=int(name[len(pref):])
            require(group not in assets,"Duplicate genuine original RGB group artifact")
            assets[group]={k:x[k] for k in ("id","name","digest","size_in_bytes")}
        return assets
    old,new=canonical(prior),canonical(live)
    require(len(old)==BEGIN_GROUP and set(old)==set(range(BEGIN_GROUP)),
            "Old signed batch001 source does not cover groups0..132")
    for group,item in old.items():
        require(group in new and new[group]==item,
                "Previously source-verified original p0032 ZIP changed or disappeared")
    return {
        "schema":"mois_original_p0032_future_groups_only_new_SHA_gate_v1",
        "original_p0032_run":RUN,
        "prior_completed_original_photo_work_groups":sorted(old),
        "new_completed_original_photo_work_groups":sorted(set(new)-set(old)),
        "prior_source_worker_SHA_records_unchanged":True,
        "new_original_photo_downloads_or_roi_reinferences":0,
    }


def seal(source, workers, group, output, *, expected_sha=SOURCE_SHA,
         expected_attempts=EXPECTED_PHOTOS, expected_run=RUN):
    import pandas as pd

    source,workers,output=map(Path,(source,workers,output))
    require(not output.exists(),"Refuse to overwrite an immutable signed original-photo worker")
    require(type(group) is int and BEGIN_GROUP <= group < END_GROUP,
            "This pipeline is exclusively for the still unmeasured p0032 workers133..144")
    mf=source/"manifest.json"
    require(sha(mf)==expected_sha,"Original p0032 source/observation CSV manifest SHA mismatch")
    manifest=json.loads(mf.read_text(encoding="utf-8"))
    require(manifest.get("schema")=="mois_full_source_observation_sampling_cohort_v1"
            and manifest.get("full_frame_photo_ids")==93209351
            and manifest.get("new_photo_attempts_planned")==expected_attempts
            and manifest.get("group_ids")==list(range(END_GROUP))
            and manifest.get("worker_ids")==list(range(manifest["work_units"])),
            "Unrecognized original 148071-photo source cohort")
    expected_units=list(range(group*4,min((group+1)*4,manifest["work_units"])))
    require(expected_units,"Missing source photo units for p0032 worker group")
    rows,audit=collect_wave(source,workers,METHOD)
    require(audit["completed_units"]==expected_units
            and audit["attempted"]==len(rows)>0
            and audit["planned_photo_attempts"]==expected_attempts
            and set(audit["missing_units"])==
                set(manifest["worker_ids"])-set(expected_units),
            "Completed genuine source photo work units differ from frozen p0032 source")
    ids=[r["photo_id"] for r in rows]
    obs={r["observation_id"] for r in rows}
    require(len(ids)==len(set(ids))
            and all(s.isdecimal() for s in ids)
            and all(bool(x) for x in obs),
            "Original photo or observation identity missing or repeated inside work group")
    hues=Counter(r["morph"] for r in rows)
    require(all(r["morph"] in (*HUES,"mixed_uncertain") for r in rows)
            and sum(hues[h] for h in HUES)==audit["classifiable"]
            and (audit["attempted"]==
                 audit["acquisition_failed"]+audit["classifiable"]+
                 audit["unclassifiable_acquired"]),
            "Original photo pixel class/unknown or download outcome altered")
    dest=output.parent
    dest.mkdir(parents=True,exist_ok=True)
    temp=Path(tempfile.mkdtemp(prefix="mois-p0032-photo-delta-",dir=dest))
    try:
        table=pd.DataFrame(rows)
        table=table.sort_values("photo_id",key=lambda s:s.map(int)).reset_index(drop=True)
        photo=temp/"original_photo_ROIs.parquet"
        table.to_parquet(photo,index=False,compression="zstd")
        summary=temp/"original_photo_hue_counts.csv"
        with summary.open("w",newline="",encoding="utf-8") as f:
            w=csv.writer(f,lineterminator="\n")
            w.writerow(("morph","true_original_photo_attempts"))
            for c in (*HUES,"mixed_uncertain"):w.writerow((c,hues[c]))
        receipt={
            "schema":SCHEMA,
            "signed_original_source_run":expected_run,
            "source_manifest_sha256":expected_sha,
            "original_source_photo_opportunities":expected_attempts,
            "original_frozen_source_worker_group":group,
            "original_completed_source_units":expected_units,
            "original_terminal_source_units_SHA_attestations":audit["completed_unit_audit"],
            "original_true_photo_attempts":len(rows),
            "original_FCP_roi_four_hue_classified":audit["classifiable"],
            "original_photo_image_acquisition_failed":audit["acquisition_failed"],
            "original_acquired_unclassifiable":audit["unclassifiable_acquired"],
            "hue_counts":{c:hues[c] for c in (*HUES,"mixed_uncertain")},
            "original_distinct_observations":len(obs),
            "frozen_original_fcp_roi_v4_attested":True,
            "original_photo_and_terminal_SHA_certified":True,
            "frozen_original_p0032_precolour_source_probability_unchanged":True,
            "other_prior_original_worker_ZIPs_downloaded":0,
            "new_S3_photo_downloads_or_model_reinferences":0,
            "standalone_p0032_group_not_yet_fully_global_photo_deduplicated":True,
            "globally_biological_flower_hue_rate_unestimated":True,
            "actual_output_sha256":{
                photo.name:sha(photo),
                summary.name:sha(summary),
            },
        }
        (temp/"p0032_original_group_receipt.json").write_text(
            json.dumps(receipt,indent=2,sort_keys=True)+"\n")
        os.replace(temp,output)
        print("MOIS_P0032_GENUINE_ORIGINAL_NEW_APPEND_ONLY_WORKER",
              json.dumps({k:receipt[k] for k in (
                  "original_frozen_source_worker_group",
                  "original_true_photo_attempts",
                  "original_FCP_roi_four_hue_classified",
                  "original_photo_image_acquisition_failed")},sort_keys=True),flush=True)
        return receipt
    except BaseException:
        shutil.rmtree(temp,ignore_errors=True)
        raise


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("source","workers","output"):
        p.add_argument("--"+k,type=Path,required=True)
    p.add_argument("--group",type=int,required=True)
    args=p.parse_args()
    seal(args.source,args.workers,args.group,args.output)


if __name__=="__main__":
    main()
