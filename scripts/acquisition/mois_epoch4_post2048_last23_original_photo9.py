#!/usr/bin/env python3
"""Terminal epoch4 partial 23-taxon cohort after signed offset2048.

Never re-query the previous 17 disjoint 128-taxon cohorts. This is a 23-taxon
source *only* until the same-run original-RGB worker, independent audit and
formal append complete. All historic photo IDs / observation UUIDs are checked;
all-history image raw SHA identity is NOT established.
"""
import argparse, csv, hashlib, json, re, time
from pathlib import Path
import pandas as pd
from scripts.acquisition import mois_epoch4_post1792_rebased_photo9 as prev

m=prev.m
LAST_DELTA="old/epoch4_offset2048/new_real_original_nine_photo9_delta.parquet"
LAST_DELTA_SHA="255e1813f808c392b03fbefec11269dae9618ff4194d26d70176f47756b9cd7c"
LAST_ATLAS_SHA="20c4692a821b2cfbaf41e1f32bf0817e718165870fe1e1ee5e60a461b94ae8ee"
m.SCHEMA="mois_incremental_20261116_epoch4_last23_post2048_true_original_photo9_v1"
m.SOURCE_ARM="MOIS_EPOCH4_LAST23_FULL_ID_UUID_RECENT_RAW_DISJOINT"
m.ATLAS_SHA=LAST_ATLAS_SHA
m.OLD_PHOTOS=1231345
m.OLD_NUMERIC_PHOTOS=467566
m.OLD_NUMERIC_TAXA=39912
m.OLD_CLEAR_TAXA=39173
m.MAX_SOURCE_TAXA=23
m.DELTA_FILES=m.DELTA_FILES+(LAST_DELTA,)
m.DELTA_N=m.DELTA_N+(1,)

def history_last23():
    for p,expected in (
        (prev.EPOCH3_DELTA_PATH,prev.EPOCH3_DELTA_SHA),
        (prev.PREV_DELTA_PATH,prev.PREV_DELTA_SHA),
        (prev.NEW_DELTA_PATH,prev.NEW_DELTA_SHA),
        (LAST_DELTA,LAST_DELTA_SHA)):
        m.require(Path(p).is_file() and m.sha256(p)==expected,
                  "historical signed independent audit delta missing/altered: "+p)
    ids,obs,raw=prev._full_history()
    m.require(len(ids)==1231345 and len(raw)==101 and
        "763811321" in ids and
        "183a6f88-4367-417b-a0b3-404d0a0ccd0b" in obs and
        "fb25d663cd95983e6e44078efe8a29504a2f3ad1d1371a9bedb550913ebffe56" in raw,
        "offset2048 signed original photo, independent observation, raw SHA absent")
    return ids,obs,raw
m.history=history_last23

def select_last23(atlas,n=23,offset=2176):
    z=atlas.loc[atlas.after_strict9.eq(0)].copy()
    m.require(len(atlas)==42111 and len(z)==2199 and
      atlas.inat_taxon_id.is_unique and z.inat_taxon_id.str.fullmatch(r"\d+").all() and
      n==23 and offset==2176 and offset+n==len(z),
      "last 23 frozen zero-ROI taxon identities or partial offset changed")
    z["priority"]=z.inat_taxon_id.map(lambda s:hashlib.sha256(
      ("MOIS-20261115-EPOCH4-HASH-DISJOINT-2199|"+s).encode()).hexdigest())
    ordered=z.sort_values("priority",kind="stable").inat_taxon_id.tolist()
    m.require(len(ordered)==len(set(ordered))==2199,
              "duplicate frozen taxon in full 2199 selection")
    previous={x for i in range(0,2176,128) for x in ordered[i:i+128]}
    selected=ordered[2176:]
    m.require(len(previous)==2176 and len(selected)==23 and
              not set(selected)&previous, "last23 collide with previous 2176 taxon IDs")
    return selected
m.select_taxa=select_last23

def source_last23(atlas,out,pace=0.8,quota=23,offset=2176):
    out=Path(out)
    m.require(not out.exists() and 0.6<=pace<=4 and quota==23 and offset==2176,
              "last23 one-time partial source invalid or attempted duplicate")
    frame=m.read_atlas(atlas)
    ids,obs,raw=history_last23()
    taxa=select_last23(frame,quota,offset)
    selected=[];stat=[]
    for i,tid in enumerate(taxa):
        st={"inat_taxon_id":tid,"status":"source_queried","accepted_unused_photo":0}
        try:
            payload=m.query(tid)
            recent=[x for x in payload["results"] if m.valid_new_observation(x)]
            valid,_=m.source_rows(recent,tid)
            valid=[p for p in valid if p["photo_id"] not in ids and
                   p["observation_id"] not in obs]
            valid.sort(key=lambda p:hashlib.sha256((
                "MOIS-20261114-SECOND-INDEPENDENT-ORIGINAL|"+p["photo_id"]+"|"+
                p["observation_id"]).encode()).hexdigest())
            if valid:
                selected.append(valid[0]);st["accepted_unused_photo"]=1
            st["API_recent_records"]=len(recent)
        except Exception as exc:
            st["status"]="API_ERROR_REVIEW";st["error"]=str(exc)[:200]
        stat.append(st)
        time.sleep(pace)
    m.require(len(stat)==23 and not any(x["status"]=="API_ERROR_REVIEW" for x in stat),
              "partial source API incomplete, prohibit worker and append")
    p=pd.DataFrame(selected)
    if not p.empty:
        m.require(p.photo_id.is_unique and p.observation_id.is_unique and
            p.inat_taxon_id.is_unique and
            not set(p.photo_id)&ids and not set(p.observation_id)&obs and
            p.photo_license.isin(m.CC).all() and
            p.photo_url_large.str.startswith(
                "https://inaturalist-open-data.s3.amazonaws.com/photos/").all() and
            set(p.inat_taxon_id).issubset(set(taxa)),
            "source candidate reused, unlicensed, wrong host, or non-independent")
    out.mkdir(parents=True)
    pd.DataFrame(stat).to_csv(out/"source_taxon_status.csv",index=False)
    p.to_csv(out/"source_candidates.csv",index=False)
    report={"schema":m.SCHEMA,"stage":"last23_source_metadata_only",
        "old_original_photo_attempt_ID_count":m.OLD_PHOTOS,
        "old_original_photo_obs_ID_overlap":0,
        "old_all_image_raw_SHA_disjointness_proven":False,
        "recent_nonempty_raw_image_SHA_count":len(raw),
        "searched_zero_numeric_taxa":23,"hash_disjoint_taxon_offset":2176,
        "new_independent_taxa":len(p),
        "new_independent_original_photo_candidates":len(p),
        "source_API_error_count":0,"created_utc_lower_bound":m.MIN_CREATED,
        "original_pixels_fetched":0,"strict9_numeric_photos":0,
        "source_status_SHA256":m.sha256(out/"source_taxon_status.csv"),
        "source_candidates_SHA256":m.sha256(out/"source_candidates.csv")}
    if not p.empty:
        p=p.copy()
        p["measurement_id"]="MOIS-NEW1411-"+p.photo_id
        p["image_filename"]=p.measurement_id+".jpg";p["shard"]="0"
        (out/"shards").mkdir()
        shard=out/"shards/shard-000.csv"
        p.to_csv(shard,index=False,lineterminator="\n")
        manifest={"schema":m.WORKER_SOURCE_SCHEMA,"source_arm":m.SOURCE_ARM,
          "colour_measured":False,"work_units":1,"worker_ids":[0],
          "photo_shard_sha256":{"0":m.sha256(shard)},
          "new_photo_attempts_planned":len(p),
          "source_previous_attempts":m.OLD_PHOTOS,
          "source_original_42111_photo9_atlas_SHA256":m.ATLAS_SHA,
          "source_historical_original_photo_and_observation_IDs_checked":True,
          "source_nominal_taxa":42111,"raw_image_pixels_downloaded":False,
          "numeric_nine_fraction_ROIs_measured":False,
          "previously_measured_photo_IDs_rerun":0}
        (out/"manifest.json").write_text(json.dumps(manifest,sort_keys=True,indent=2)+"\n")
        report["source_shard_SHA256"]=m.sha256(shard)
        report["source_manifest_SHA256"]=m.sha256(out/"manifest.json")
    (out/"receipt.json").write_text(json.dumps(report,sort_keys=True,indent=2)+"\n")
    print("MOIS_EPOCH4_LAST23_SOURCE_ONLY",json.dumps(report,sort_keys=True),flush=True)
m.source=source_last23

def main():
    parser=argparse.ArgumentParser()
    stages=parser.add_subparsers(dest="stage",required=True)
    for stage,flags in (
       ("source",("atlas","output")),("audit",("source","worker","output")),
       ("append",("atlas","source","worker","audit","output"))):
        sub=stages.add_parser(stage)
        for flag in flags:sub.add_argument("--"+flag,required=True)
        if stage=="source":sub.add_argument("--offset",type=int,default=2176)
    args=parser.parse_args()
    if args.stage=="source":source_last23(args.atlas,args.output,offset=args.offset)
    elif args.stage=="audit":m.audit(args.source,args.worker,args.output)
    else:m.append(args.atlas,args.source,args.worker,args.audit,args.output)

if __name__=="__main__":main()
