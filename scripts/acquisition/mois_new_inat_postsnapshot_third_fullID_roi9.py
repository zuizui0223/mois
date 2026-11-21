#!/usr/bin/env python3
"""Third hash-disjoint post-snapshot new iNaturalist photo wave.

Requires both earlier signed source/audit pairs and all 1,231,244 historic
photo IDs. Never downloads old pixels, never auto-appends the old atlas.
"""
from __future__ import annotations
import argparse, json, shutil
from pathlib import Path
from datetime import datetime, timezone
from urllib.parse import urlsplit
import pandas as pd
from scripts.acquisition.mois_new_inat_postsnapshot_verified_original_roi9 import (
    check, write_json, UUID, ARM, EXPECTED_ATTEMPTED, WORKER_SCHEMA,
    audit as base_audit,
)
from scripts.acquisition.mois_zero9_recent_inat_new_created_photo_source import (
    ATLAS_SHA, source_sha, SCHEMA as API_SCHEMA, CC,
)
from scripts.acquisition.mois_full_index_real_original_photo_measure import sha as sha256

OFFSET = 128
HIST_NAME = "historical_1231244_original_photo_observation_taxon_ID_registry.csv.gz"
AUD_NAME = "mois_new_inat_post20261028_actual_original_nine_roi.csv.gz"
def load_prior(src, audited, expected_count, offset=None):
    src, audited = Path(src), Path(audited)
    m = json.loads((src/"manifest.json").read_text())
    g = json.loads((src/"gate_receipt.json").read_text())
    a = json.loads((audited/"receipt.json").read_text())
    reg=src/HIST_NAME
    check(reg.is_file() and m.get("schema")==WORKER_SCHEMA
          and m.get("source_arm")==ARM
          and m.get("source_previous_attempts")==EXPECTED_ATTEMPTED
          and m.get("source_prior_no_reused_original_photo_or_observation") is True
          and m.get("source_previous_ID_registry_SHA256")==sha256(reg)
          and g.get("historical_verified_photo_ID_count")==EXPECTED_ATTEMPTED
          and g.get("historical_full_registry_sha256")==sha256(reg)
          and g.get("source_manifest_sha256")==sha256(src/"manifest.json")
          and m.get("new_photo_attempts_planned")==expected_count
          and a.get("source_photo_manifest_sha256")==sha256(src/"manifest.json")
          and a.get("historical_full_ID_registry_sha256")==sha256(reg)
          and a.get("original_licensed_photo_attempts")==expected_count
          and a.get("photo_observation_taxon_collision_count")==0
          and a.get("actual_image_downloaded")==expected_count
          and a.get("actual_genuine_strict9_numeric_success")+
              a.get("actual_technical_unknowns")==expected_count
          and a.get("official_iNat_42111_original_frame_auto_appended") is False
          and a.get("audited_derivative_sha256")==sha256(audited/AUD_NAME),
          "Previous source/independent original pixel audit not signed complete")
    if offset is not None:
        check(m.get("source_fixed_hash_offset")==offset
              and g.get("new_CC_photo_opportunities")==expected_count,
              "Prior disjoint wave offset/count wrong")
    shard=src/"shards"/"shard-000.csv"
    check(sha256(shard)==m["photo_shard_sha256"]["0"],
          "Prior frozen original photo shard SHA mismatch")
    ids=pd.read_csv(shard,dtype=str,keep_default_na=False)
    check(len(ids)==expected_count and ids.photo_id.is_unique
          and ids.observation_id.is_unique and ids.inat_taxon_id.is_unique,
          "Prior independent photo/observation/taxon missing")
    raw=a.get("new_photo_image_SHA256",{})
    check(set(raw)==set(ids.photo_id)
          and len([v for v in raw.values() if v])==expected_count
          and len(set(raw.values()))==expected_count
          and all(len(v)==64 and int(v,16)>=0 for v in raw.values()),
          "Signed prior actual raw image SHA incomplete or collided")
    return m,g,a,reg,ids,raw

def gate(meta,receipt,atlas,first_source,first_audit,second_source,second_audit,output):
    meta,receipt,atlas,output=map(Path,(meta,receipt,atlas,output))
    check(not output.exists(),"Never overwrite/re-run same signed source")
    r=json.loads(receipt.read_text())
    check(r.get("schema")==API_SCHEMA and r.get("sampled_zero_ROI9_source_taxa")==64
          and r.get("source_deterministic_hash_offset")==OFFSET
          and r.get("original_source_zero_ROI9_taxa")==2229
          and r.get("source_taxa_API_errors")==0
          and r.get("old_photo_ID_1_231_244_ledger_disjointness_verified") is False
          and r.get("original_photo_pixels_downloaded")==0
          and r.get("genuine_nine_fraction_ROI_measured")==0,
          "Unsigned or unsuitable offset128 source receipt")
    check(source_sha(meta)==r["output_sha256"][meta.name]
          and source_sha(atlas)==ATLAS_SHA,
          "External metadata or frozen old atlas bytes changed")
    census=pd.read_csv(atlas,dtype={"inat_taxon_id":str})
    check(len(census)==42111 and census.inat_taxon_id.is_unique
          and int(census.after_strict9.eq(0).sum())==2229,
          "Old atlas species frame incompatible")
    first=load_prior(first_source,first_audit,2)
    second=load_prior(second_source,second_audit,1,offset=64)
    fm,fg,fa,fr,fi,fb=first
    sm,sg,sa,sr,si,sb=second
    check(sha256(fr)==sha256(sr)
          and sm.get("source_first_new_photo_IDs")==fi.photo_id.tolist()
          and sm.get("source_first_new_observation_UUIDs")==fi.observation_id.tolist()
          and sm.get("source_first_new_image_SHA256")==fb
          and sa.get("additional_first_new_photo_IDs_excluded")==2
          and sa.get("additional_first_new_raw_image_SHA_collisions")==0,
          "Second-source lineage does not match first signed source/audit")
    check(not(set(fi.photo_id)&set(si.photo_id))
          and not(set(fi.observation_id)&set(si.observation_id))
          and not(set(fi.inat_taxon_id)&set(si.inat_taxon_id))
          and not(set(fb.values())&set(sb.values())),
          "Duplicate ID/taxon/raw pixel SHA among prior two waves")
    history=pd.read_csv(fr,dtype=str,keep_default_na=False,
                        usecols=["photo_id","observation_id","inat_taxon_id"])
    check(len(history)==EXPECTED_ATTEMPTED and history.photo_id.is_unique
          and history.photo_id.ne("").all(),
          "Full historical 1,231,244 original photo IDs not intact")
    candidate=pd.read_csv(meta,dtype=str,keep_default_na=False)
    check(len(candidate)==r["source_taxa_with_new_post_snapshot_API_photo_opportunities"]
          and len(candidate)==r["new_independent_original_iNaturalist_photo_IDs_seen"]
          and 0<=len(candidate)<=64,
          "Offset128 API candidate row/receipt mismatch")
    check(candidate.photo_id.is_unique and candidate.observation_id.is_unique
          and candidate.inat_taxon_id.is_unique,
          "Offset128 photo/observation/taxon duplicate within source")
    check(candidate.inat_taxon_id.isin(
          census.loc[census.after_strict9.eq(0),"inat_taxon_id"]).all(),
          "Candidate taxon already has old numerical ROI9")
    prev_ids=pd.concat([fi,si],ignore_index=True)
    check(not(set(candidate.photo_id)&(set(history.photo_id)|set(prev_ids.photo_id)))
          and not(set(candidate.observation_id)&(
              set(history.observation_id)|set(prev_ids.observation_id)))
          and not(set(candidate.inat_taxon_id)&set(prev_ids.inat_taxon_id)),
          "Candidate reused any previously attempted ID or post-snapshot taxon")
    for t in candidate.itertuples(index=False):
        u=urlsplit(t.photo_url_large)
        check(t.photo_id.isdigit() and UUID.fullmatch(t.observation_id) is not None
              and t.photo_license in CC
              and u.scheme=="https"
              and u.hostname=="inaturalist-open-data.s3.amazonaws.com"
              and u.path.startswith("/photos/"+t.photo_id+"/")
              and u.path.lower().endswith((".jpg",".jpeg",".png")),
              "Invalid CC/taxon/URL/observation provenance")
        dt=datetime.fromisoformat(t.source_observation_created_at.replace("Z","+00:00"))
        check(dt.tzinfo is not None and
              dt.astimezone(timezone.utc)>=datetime(2026,10,28,tzinfo=timezone.utc),
              "Old observation created before snapshot cutoff")
    output.mkdir(parents=True)
    shutil.copyfile(fr,output/HIST_NAME)
    (output/"shards").mkdir()
    if len(candidate):
        shard=output/"shards"/"shard-000.csv"
        current=candidate.copy()
        current["measurement_id"]="MOIS-POST20261028-OFFSET128-"+current.photo_id
        current["image_filename"]=current.measurement_id+".jpg"
        current["shard"]="0"
        current.to_csv(shard,index=False,lineterminator="\n")
        hashes={"0":sha256(shard)}
    else:hashes={}
    prev_raw={**fb,**sb}
    m={"schema":WORKER_SCHEMA,"source_arm":ARM,"colour_measured":False,
       "source_previous_attempts":EXPECTED_ATTEMPTED,
       "source_previous_ID_registry_SHA256":sha256(output/HIST_NAME),
       "source_prior_no_reused_original_photo_or_observation":True,
       "source_prior_post_snapshot_photo_IDs":prev_ids.photo_id.tolist(),
       "source_prior_post_snapshot_observation_UUIDs":prev_ids.observation_id.tolist(),
       "source_prior_post_snapshot_original_image_SHA256":prev_raw,
       "source_first_audit_receipt_sha256":sha256(Path(first_audit)/"receipt.json"),
       "source_second_audit_receipt_sha256":sha256(Path(second_audit)/"receipt.json"),
       "source_post_snapshot_API_receipt_SHA256":source_sha(receipt),
       "source_post_snapshot_API_metadata_SHA256":source_sha(meta),
       "source_fixed_hash_offset":OFFSET,
       "new_photo_attempts_planned":len(candidate),
       "work_units":int(len(candidate)>0),
       "worker_ids":[0] if len(candidate) else [],
       "photo_shard_sha256":hashes,
       "raw_image_pixels_downloaded":False,"new_nineROI_measured":False,
       "official_iNat_42111_original_frame_auto_appended":False,
       "species_wide_fixed_colour_proven":False}
    write_json(output/"manifest.json",m)
    g={"schema":"mois_postsnapshot_offset128_hist_plus2plus1_signed_gate_v1",
       "historical_verified_photo_ID_count":len(history),
       "previous_new_photo_source_count":len(prev_ids),
       "prior_raw_original_SHA_count":len(prev_raw),
       "new_CC_photo_opportunities":len(candidate),
       "new_photo_IDs_historical_collision_count":0,
       "new_observation_UUID_historical_collision_count":0,
       "source_manifest_sha256":sha256(output/"manifest.json"),
       "historical_full_registry_sha256":sha256(output/HIST_NAME),
       "source_only":True,"strict9_measured":0,"new_photo_pixels_downloaded":0}
    write_json(output/"gate_receipt.json",g)
    print("MOIS_OFFSET128_SIGNED_FULL_ID_PLUS_THREE_NEW_RAW_SHA_GATE",
          json.dumps(g,sort_keys=True),flush=True)

def audit_third(source,worker,output):
    source,worker,output=map(Path,(source,worker,output))
    m=json.loads((source/"manifest.json").read_text())
    g=json.loads((source/"gate_receipt.json").read_text())
    check(g.get("schema")=="mois_postsnapshot_offset128_hist_plus2plus1_signed_gate_v1"
          and m.get("source_fixed_hash_offset")==OFFSET
          and g.get("source_manifest_sha256")==sha256(source/"manifest.json")
          and g.get("historical_full_registry_sha256")==sha256(source/HIST_NAME)
          and len(m["source_prior_post_snapshot_photo_IDs"])==3
          and len(m["source_prior_post_snapshot_original_image_SHA256"])==3,
          "Incomplete signed prior wave provenance")
    now=pd.read_csv(source/"shards"/"shard-000.csv",dtype=str,keep_default_na=False)
    check(not(set(now.photo_id)&set(m["source_prior_post_snapshot_photo_IDs"]))
          and not(set(now.observation_id)&set(m["source_prior_post_snapshot_observation_UUIDs"])),
          "Third wave repeats first/second source IDs")
    b=pd.read_csv(worker/"shard_000"/"photos.csv",dtype=str,keep_default_na=False)
    sha=set(v for v in b.image_sha256 if v)
    check(not(sha&set(m["source_prior_post_snapshot_original_image_SHA256"].values())),
          "Original image byte SHA matches prior new photo (duplicate upload)")
    base_audit(source,worker,output)
    a=json.loads((output/"receipt.json").read_text())
    a.update({"additional_prior_new_photo_IDs_excluded":3,
              "additional_prior_new_observation_UUID_collisions":0,
              "additional_prior_new_raw_image_SHA_collisions":0,
              "full_historical_raw_image_SHA_collision_check_possible":False,
              "official_iNat_42111_original_frame_auto_appended":False})
    write_json(output/"receipt.json",a)
    print("MOIS_OFFSET128_REAL_PIXEL_INDEPENDENT_AUDIT",json.dumps({
        k:a[k] for k in ("original_licensed_photo_attempts",
                        "actual_image_downloaded","actual_genuine_strict9_numeric_success",
                        "actual_clear_photo_nine_palette_labels",
                        "actual_technical_unknowns")},sort_keys=True),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    subs=p.add_subparsers(dest="command",required=True)
    g=subs.add_parser("gate")
    for k in ("meta","receipt","atlas","first-source","first-audit",
              "second-source","second-audit","output"):
        g.add_argument("--"+k,required=True)
    a=subs.add_parser("audit")
    for k in ("source","worker","output"):a.add_argument("--"+k,required=True)
    v=p.parse_args()
    if v.command=="gate":
        gate(v.meta,v.receipt,v.atlas,v.first_source,v.first_audit,
             v.second_source,v.second_audit,v.output)
    else:audit_third(v.source,v.worker,v.output)
