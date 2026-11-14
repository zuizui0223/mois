#!/usr/bin/env python3
"""One-shot post-source-2026-11-14 incremental original-photo ROI9 pipeline.

This scans only genuinely *newly created* iNaturalist observations since the
last 1909-taxon source scan ended, for 128 previously zero-numeric taxa.
All 1,231,339 previous photo IDs and observation UUIDs must be disjoint.
Source-only never counts as pixel measurement; byte SHA is checked only
against the 95 recent measured photos, not the full older 1.23m images.
"""
from __future__ import annotations
import argparse, hashlib, json, re, time
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen
from urllib.error import HTTPError, URLError

import numpy as np
import pandas as pd
from scripts.analysis.mois_species_first_palette9 import PALETTE, class_labels, sha256
from scripts.analysis.mois_species_first_wave3_incremental import add_to_catalogue
from scripts.acquisition.mois_zero9_recent_inat_new_created_photo_source import (
    strict_new_photo, source_rows, CC)
from scripts.acquisition.mois_full_index_real_original_photo_measure import (
    SCHEMA as WORKER_SOURCE_SCHEMA)

SCHEMA="mois_incremental_20261114_new_created_cc_original_photo9_v1"
SOURCE_ARM="MOIS_POSTSOURCE_20261114_NEW_CREATED_INDEPENDENT_ORIGINALS"
ATLAS_SHA="2929b68758f000906ac303b9edd26ddf60bcdc317ddf7f4185e6a36c6bb65380"
EXPECTED_ATLAS=42111
OLD_PHOTOS=1231339
OLD_NUMERIC_PHOTOS=467565
OLD_NUMERIC_TAXA=39911
OLD_CLEAR_TAXA=39172
MIN_CREATED="2026-11-14T05:35:00+00:00"
CREATED_D1="2026-11-14"
HUES=("white","yellow_orange","red_pink","blue_purple")
MAX_SOURCE_TAXA=128
LIMIT=40
API="https://api.inaturalist.org/v1/observations"
DELTA_FILES=("old/new11/mois_new11_iNat_original_photo9_audited_delta.parquet",
             "old/new64/mois_postsnap1909_authentic_original_photo9_delta.parquet",
             "old/new18/mois_postsnap1909_remaining18_original_photo9_delta.parquet",
             "old/offset1152/new_real_original_nine_photo9_delta.parquet",
             "old/offset1536/new_real_original_nine_photo9_delta.parquet")
DELTA_N=(11,64,18,1,1)
HISTORY_PATH="old/history/historical_1231244_original_photo_observation_taxon_ID_registry.csv.gz"

def require(x, message):
    if not bool(x): raise ValueError("MOIS_20261114_INCREMENTAL: "+message)

def read_atlas(path):
    require(sha256(path)==ATLAS_SHA,"latest after offset1536 signed atlas bytes changed")
    a=pd.read_csv(path,compression="infer",dtype={"inat_taxon_id":str})
    require(len(a)==EXPECTED_ATLAS and a.inat_taxon_id.is_unique,"42,111 nominal taxa not unique")
    require(int(a.after_strict9.gt(0).sum())==OLD_NUMERIC_TAXA,"previous numeric taxon tally drift")
    require(int(a.dominant_hue9.fillna("").ne("").sum())==OLD_CLEAR_TAXA,"previous clear hue drift")
    return a

def history():
    old=pd.read_csv(HISTORY_PATH,compression="gzip",
                    usecols=["photo_id","observation_id"],dtype=str,keep_default_na=False)
    require(len(old)==1231244 and old.photo_id.is_unique,"full historical photo-ID census absent")
    allframes=[old]
    recent_sha=[]
    for path,n in zip(DELTA_FILES,DELTA_N):
        if path=="old/offset1152/new_real_original_nine_photo9_delta.parquet":
            require(sha256(path)=="9d11bbfd5bd46f82b41b887f16e0ca07abaeccbd246bb9b804673abdf27483b5",
                    "offset1152 independent audit delta bytes do not match signed SHA")
        d=pd.read_parquet(path)
        if path=="old/offset1536/new_real_original_nine_photo9_delta.parquet":
            require(sha256(path)=="2273b25f50054e38fd5e28d8621609d97c532d907717da641d9bd095a3b97657",
                    "offset1536 independent audit delta bytes do not match signed SHA")
        require(len(d)==n,"missing signed recent image delta "+path)
        require({"photo_id","observation_id","image_sha256"}.issubset(d.columns),
                "recent image ID/SHA information missing "+path)
        d=d.fillna({"photo_id":"","observation_id":"","image_sha256":""})
        d[["photo_id","observation_id"]]=d[["photo_id","observation_id"]].astype(str)
        require(d.photo_id.is_unique and d.observation_id.is_unique,
                "recent original photo or independent observation repeats "+path)
        allframes.append(d[["photo_id","observation_id"]])
        recent_sha.extend(s for s in d.image_sha256.astype(str) if re.fullmatch("[0-9a-f]{64}",s))
    p=pd.concat(allframes,ignore_index=True)
    require(len(p)==OLD_PHOTOS and p.photo_id.is_unique,
            "1,231,339 previous photo IDs not complete")
    # Historical observations may contain several photos, but each recent
    # cohort must use a different observation from all preceding cohorts.
    old_obs=set(old.observation_id)-{""}
    for d in allframes[1:]:
        new_obs=set(d.observation_id)-{""}
        require(not new_obs.intersection(old_obs),"previous recently added observation overlaps")
        old_obs.update(new_obs)
    require(len(recent_sha)==len(set(recent_sha)),
            "nonempty recent 95 original byte SHA collision")
    return set(p.photo_id),old_obs,set(recent_sha)

def query(taxon):
    query={"taxon_id":str(taxon),"created_d1":CREATED_D1,
           "photos":"true","photo_licensed":"true","per_page":str(LIMIT),
           "page":"1","order_by":"created_at","order":"desc","captive":"false"}
    req=Request(API+"?"+urlencode(query),
                headers={"Accept":"application/json",
                         "User-Agent":"MOIS-unique-original-photo-ecology/2.0"})
    for trial in range(3):
        try:
            with urlopen(req,timeout=25) as r: payload=json.loads(r.read(2_000_000))
            require(isinstance(payload,dict) and isinstance(payload.get("results"),list),
                    "unexpected API response")
            return payload
        except (HTTPError,URLError,TimeoutError,ValueError) as e:
            if isinstance(e,HTTPError) and e.code not in (429,500,502,503,504):
                raise
            if trial==2: raise
            time.sleep(2**trial)

def valid_new_observation(x):
    try:
        t=datetime.fromisoformat(str(x.get("created_at","")).replace("Z","+00:00"))
        return t.tzinfo is not None and t.astimezone(timezone.utc)>=datetime.fromisoformat(MIN_CREATED)
    except (TypeError,ValueError,AttributeError):
        return False

def select_taxa(a,n=MAX_SOURCE_TAXA,offset=0):
    z=a.loc[a.after_strict9.eq(0)].copy()
    require(len(z)==2200 and z.inat_taxon_id.str.fullmatch(r"\d+").all(),
            "after offset1536 exactly 2,200 taxon IDs must remain without numeric ROI")
    z["priority"]=z.inat_taxon_id.map(lambda s:hashlib.sha256(
        ("MOIS-20261114-INDEPENDENT-CREATED-POST-LAST-SOURCE|"+s).encode()).hexdigest())
    return z.sort_values("priority",kind="stable").iloc[offset:offset+n].inat_taxon_id.tolist()

def source(atlas,out,pace=0.8,quota=MAX_SOURCE_TAXA,offset=0):
    out=Path(out)
    require(not out.exists() and 0.6<=pace<=4 and quota==MAX_SOURCE_TAXA,
            "source rerun or unbounded API request forbidden")
    require(isinstance(offset,int) and offset>=0 and offset%MAX_SOURCE_TAXA==0
            and offset+quota<=2200,"invalid hash-disjoint source offset")
    atlas=read_atlas(atlas)
    oldids,oldobs,rawhash=history()
    taxa=select_taxa(atlas,quota,offset)
    selected=[];stat=[]
    for index,tid in enumerate(taxa):
        st={"inat_taxon_id":tid,"status":"source_queried","accepted_unused_photo":0}
        try:
            payload=query(tid)
            recent=[x for x in payload["results"] if valid_new_observation(x)]
            valid,_=source_rows(recent,tid)
            # Only candidate *new* observations, not existing old source IDs.
            valid=[p for p in valid
                   if p["photo_id"] not in oldids and p["observation_id"] not in oldobs]
            valid.sort(key=lambda p:hashlib.sha256(
                ("MOIS-20261114-SECOND-INDEPENDENT-ORIGINAL|"+p["photo_id"]+
                 "|"+p["observation_id"]).encode()).hexdigest())
            if valid:
                selected.append(valid[0])
                st["accepted_unused_photo"]=1
            st["API_recent_records"]=len(recent)
        except Exception as e:
            st["status"]="API_ERROR_REVIEW"
            st["error"]=str(e)[:200]
        stat.append(st)
        if (index+1)%16==0: print("MOIS_INCREMENTAL_NEW_OBSERVATION_METADATA",index+1,flush=True)
        time.sleep(pace)
    require(not any(x["status"]=="API_ERROR_REVIEW" for x in stat),
            "incomplete source API cohort; do not freeze or pixel-measure")
    p=pd.DataFrame(selected)
    require(len(p)<=quota,"candidate quota changed")
    if not p.empty:
        require(p.photo_id.is_unique and p.observation_id.is_unique and
                p.inat_taxon_id.is_unique,"new source candidate not independent")
        require(not set(p.photo_id)&oldids and not set(p.observation_id)&oldobs,
                "new source repeats entire historical photo or observation")
        require(p.photo_license.isin(CC).all() and
                p.photo_url_large.str.startswith(
                 "https://inaturalist-open-data.s3.amazonaws.com/photos/").all(),
                "original image host or photo CC licence missing")
        require(set(p.inat_taxon_id).issubset(set(taxa)),
                "new photo belongs to taxon not sampled")
    out.mkdir(parents=True)
    pd.DataFrame(stat).to_csv(out/"source_taxon_status.csv",index=False)
    p.to_csv(out/"source_candidates.csv",index=False)
    report={"schema":SCHEMA,"stage":"post_last_source_API_verified_original_candidates",
            "old_original_photo_attempt_ID_count":OLD_PHOTOS,
            "old_original_photo_obs_ID_overlap":0,
            "old_all_image_raw_SHA_disjointness_proven":False,
            "recent_95_nonempty_raw_image_SHA_count":len(rawhash),
            "searched_zero_numeric_taxa":len(taxa),
            "hash_disjoint_taxon_offset":offset,
            "new_independent_taxa":len(p),
            "new_independent_original_photo_candidates":len(p),
            "source_API_error_count":0,
            "created_utc_lower_bound":MIN_CREATED,
            "original_pixels_fetched":0,"strict9_numeric_photos":0,
            "source_status_SHA256":sha256(out/"source_taxon_status.csv"),
            "source_candidates_SHA256":sha256(out/"source_candidates.csv")}
    if not p.empty:
        p=p.copy()
        p["measurement_id"]="MOIS-NEW1411-"+p.photo_id
        p["image_filename"]=p.measurement_id+".jpg"
        p["shard"]="0"
        (out/"shards").mkdir()
        shard=out/"shards/shard-000.csv"
        p.to_csv(shard,index=False,lineterminator="\n")
        manifest={"schema":WORKER_SOURCE_SCHEMA,"source_arm":SOURCE_ARM,
                  "colour_measured":False,"work_units":1,"worker_ids":[0],
                  "photo_shard_sha256":{"0":sha256(shard)},
                  "new_photo_attempts_planned":len(p),
                  "source_previous_attempts":OLD_PHOTOS,
                  "source_original_42111_photo9_atlas_SHA256":ATLAS_SHA,
                  "source_historical_original_photo_and_observation_IDs_checked":True,
                  "source_nominal_taxa":42111,"raw_image_pixels_downloaded":False,
                  "numeric_nine_fraction_ROIs_measured":False,
                  "previously_measured_photo_IDs_rerun":0}
        (out/"manifest.json").write_text(json.dumps(manifest,sort_keys=True,indent=2)+"\n")
        report["source_shard_SHA256"]=sha256(shard)
        report["source_manifest_SHA256"]=sha256(out/"manifest.json")
    (out/"receipt.json").write_text(json.dumps(report,sort_keys=True,indent=2)+"\n")
    print("MOIS_INCREMENTAL_SOURCE_ONLY",json.dumps(report,sort_keys=True),flush=True)

def audited_frame(src,worker):
    src=Path(src);worker=Path(worker)
    m=json.loads((src/"manifest.json").read_text())
    r=json.loads((worker/"shard_000/receipt.json").read_text())
    shard=src/"shards/shard-000.csv"
    t=worker/"shard_000/photos.csv"
    require(m["schema"]==WORKER_SOURCE_SCHEMA and m["source_arm"]==SOURCE_ARM and
            m["source_previous_attempts"]==OLD_PHOTOS and
            0<m["new_photo_attempts_planned"]<=MAX_SOURCE_TAXA and
            sha256(shard)==m["photo_shard_sha256"]["0"],
            "original source manifest/source roster changed")
    require(r["schema"]=="mois_full_original93m_real_RGB_photo_unit_v1"
            and r["photo_colours_measured"] is True
            and r["attempted"]==m["new_photo_attempts_planned"]
            and r["source_sha256"]==sha256(shard)
            and r["photos_sha256"]==sha256(t),
            "genuine FCP pixel receipt and terminal source SHA failed")
    original=pd.read_csv(shard,dtype=str,keep_default_na=False)
    measured=pd.read_csv(t,dtype=str,keep_default_na=False)
    require(len(original)==len(measured)==m["new_photo_attempts_planned"] and
            original.photo_id.is_unique and measured.photo_id.is_unique and
            original.observation_id.is_unique and original.inat_taxon_id.is_unique and
            set(original.photo_id)==set(measured.photo_id),
            "source terminal photo ID or row count overlap")
    for c in ("inat_taxon_id","observation_id","photo_url_large","photo_license"):
        require(original.set_index("photo_id").sort_index()[c].equals(
                measured.set_index("photo_id").sort_index()[c]),
                "worker changed original image identity, observation or licence")
    values=measured[["flower_fraction_"+q for q in PALETTE]].apply(
        pd.to_numeric,errors="coerce").to_numpy(float)
    numeric=(measured.morph.isin(HUES).to_numpy(bool) &
             np.isfinite(values).all(axis=1) &
             (values>=0).all(axis=1)&(values<=1.00001).all(axis=1) &
             (np.abs(np.nan_to_num(values).sum(axis=1)-1)<=.001))
    require(int(numeric.sum())==r["classifiable"],"worker classifiable strict9 changed")
    labels=np.full(len(measured),"",dtype=object)
    if numeric.any(): labels[numeric]=class_labels(values[numeric])
    res=measured[["inat_taxon_id","photo_id","observation_id","image_sha256"]].copy()
    res["strict9"]=numeric;res["dominant_label9"]=labels
    for i,col in enumerate(PALETTE):
        res["strict9_"+col]=np.where(numeric,values[:,i],np.nan)
    require(res.loc[~res.strict9,["strict9_"+q for q in PALETTE]].isna().all().all(),
            "technical unknown nine numeric fractions must remain missing")
    require(res.loc[res.strict9,"image_sha256"].astype(str).str.fullmatch(
        r"[0-9a-f]{64}").all(),"successful original byte image SHA missing")
    return res,{"worker_receipt_sha256":sha256(worker/"shard_000/receipt.json"),
                "worker_terminal_sha256":sha256(t),"source_shard_sha256":sha256(shard)}

def audit(source_root,worker_root,out):
    out=Path(out)
    require(not out.exists(),"refuse independent audit overwrite")
    d,sha=audited_frame(source_root,worker_root)
    ids,obs,oldraw=history()
    require(not set(d.photo_id)&ids and not set(d.observation_id)&obs,
            "new actual worker reuses an old photo ID or observation")
    newraw=[s for s in d.image_sha256.fillna("").astype(str) if re.fullmatch(
        r"[0-9a-f]{64}",s)]
    require(len(newraw)==len(set(newraw)) and not set(newraw)&oldraw,
            "new byte-identical image collides with previous 95 originals")
    out.mkdir()
    delta=out/"new_real_original_nine_photo9_delta.parquet"
    d.to_parquet(delta,index=False)
    receipt={"schema":SCHEMA,"stage":"genuine_pixels_independently_audited",
             "new_authentic_original_attempts":len(d),
             "new_strict9_numeric_photos":int(d.strict9.sum()),
             "new_numeric_first_taxa":int(d.loc[d.strict9,"inat_taxon_id"].nunique()),
             "new_clear_dominant_hue_taxa":int(d.dominant_label9.ne("").sum()),
             "new_technical_unknowns":int((~d.strict9).sum()),
             "all_historical_photo_and_observation_ID_collisions":0,
             "raw_byte_SHA_collisions_with_last95":0,
             "raw_byte_SHA_collisions_with_old1231244_certified":False,
             "technical_unknown_nine_fractions_all_missing":True,
             "source_and_worker_sha256":sha,
             "verified_delta_sha256":sha256(delta)}
    (out/"receipt.json").write_text(json.dumps(receipt,sort_keys=True,indent=2)+"\n")
    print("MOIS_INCREMENTAL_GENUINE_NINE_AUDIT",json.dumps(receipt,sort_keys=True),flush=True)

def append(atlas_path,source_root,worker_root,audit_root,out):
    out=Path(out)
    require(not out.exists(),"no duplicate atlas append")
    atlas=read_atlas(atlas_path)
    d,sha=audited_frame(source_root,worker_root)
    signed=json.loads((Path(audit_root)/"receipt.json").read_text())
    delta=Path(audit_root)/"new_real_original_nine_photo9_delta.parquet"
    require(signed["schema"]==SCHEMA and signed["stage"]=="genuine_pixels_independently_audited"
            and signed["verified_delta_sha256"]==sha256(delta)
            and signed["new_authentic_original_attempts"]==len(d),
            "independent genuine-pixel audit missing")
    actual=pd.read_parquet(delta)
    require(set(actual.photo_id.astype(str))==set(d.photo_id.astype(str)) and
            int(actual.strict9.sum())==int(d.strict9.sum()),
            "audit delta not identical to pixel source")
    require(atlas.set_index("inat_taxon_id").loc[d.inat_taxon_id,
              "after_strict9"].eq(0).all(),
            "cannot append a photo to taxon already numerically covered")
    now=add_to_catalogue(atlas,d,expected_photos=len(d))
    require(len(now)==EXPECTED_ATLAS and int(now.after_strict9.gt(0).sum())==
            OLD_NUMERIC_TAXA+int(d.strict9.sum()) and
            int(now.dominant_hue9.fillna("").ne("").sum())==
            OLD_CLEAR_TAXA+int(d.dominant_label9.ne("").sum()),
            "numerical or clear-hue coverage changed unexpectedly")
    out.mkdir()
    target=out/"mois_42111_post20261114_new_independent_original_photo9.csv.gz"
    now.to_csv(target,index=False,compression="gzip")
    receipt={"schema":SCHEMA,"stage":"signed_new_original_photo9_appended",
             "new_original_photo_attempts":len(d),
             "new_strict9_photo_ROI9":int(d.strict9.sum()),
             "new_first_numeric_taxa":int(d.strict9.sum()),
             "new_clear_sampled_photo_hue_taxa":int(d.dominant_label9.ne("").sum()),
             "new_technical_unknown_photos":int((~d.strict9).sum()),
             "cumulative_photo_attempts":OLD_PHOTOS+len(d),
             "cumulative_genuine_numeric_photos":OLD_NUMERIC_PHOTOS+int(d.strict9.sum()),
             "cumulative_numeric_taxa":OLD_NUMERIC_TAXA+int(d.strict9.sum()),
             "cumulative_clear_sampled_photo_taxa":OLD_CLEAR_TAXA+
                int(d.dominant_label9.ne("").sum()),
             "zero_numeric_remaining":EXPECTED_ATLAS-OLD_NUMERIC_TAXA-int(d.strict9.sum()),
             "old1231244_raw_SHA_fully_compared":False,
             "derived_atlas_sha256":sha256(target),
             "audit_receipt_sha256":sha256(Path(audit_root)/"receipt.json")}
    (out/"receipt.json").write_text(json.dumps(receipt,sort_keys=True,indent=2)+"\n")
    print("MOIS_INCREMENTAL_SIGNED_ROI9_APPEND",json.dumps(receipt,sort_keys=True),flush=True)

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    sub=parser.add_subparsers(dest="stage",required=True)
    for stage,flags in [
        ("source",["atlas","output"]),
        ("audit",["source","worker","output"]),
        ("append",["atlas","source","worker","audit","output"])]:
        q=sub.add_parser(stage)
        for f in flags:q.add_argument("--"+f,required=True)
        if stage=="source":q.add_argument("--offset",type=int,default=0)
    args=parser.parse_args()
    if args.stage=="source": source(args.atlas,args.output,offset=args.offset)
    elif args.stage=="audit":audit(args.source,args.worker,args.output)
    else:append(args.atlas,args.source,args.worker,args.audit,args.output)
