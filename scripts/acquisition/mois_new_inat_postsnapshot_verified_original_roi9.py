#!/usr/bin/env python3
"""Source-hash pinned post-snapshot iNaturalist originals: full historic ID gate and independent pixel ROI9 audit.

Only created-after-2026-10-27 licensed iNat S3 photos; never re-download
the original 1,231,244 attempted IDs. Original pixel bytes are not distributed.
"""
from __future__ import annotations
import argparse,hashlib,json,re
from pathlib import Path
from urllib.parse import urlsplit
from datetime import datetime,timezone
import numpy as np
import pandas as pd
from scripts.acquisition.mois_full_index_real_original_photo_measure import sha as sha256
from scripts.analysis.mois_species_first_palette9 import PALETTE,class_labels
from scripts.acquisition.mois_zero9_recent_inat_new_created_photo_source import (
    ATLAS_SHA,source_sha,SCHEMA as API_SCHEMA,CC,CUTOFF,
)
from scripts.acquisition.mois_full_index_real_original_photo_measure import SCHEMA as WORKER_SCHEMA

ARM="MOIS_NEW_INAT_CREATED_20261028_ZERO9_DISJOINT"
EXPECTED_ATTEMPTED=1231244
PIX=["flower_fraction_"+x for x in PALETTE]
HUES={"white","yellow_orange","red_pink","blue_purple"}
UUID=re.compile(r"[0-9a-fA-F]{8}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{12}")
def check(ok,why):
    if not ok:raise ValueError(why)
def write_json(path,data):
    Path(path).write_text(json.dumps(data,sort_keys=True,indent=2)+"\n")

def gate(meta,receipt,atlas,history_paths,output):
    output=Path(output);check(not output.exists(),"Refusing to overwrite a source manifest")
    r=json.loads(Path(receipt).read_text())
    check(r.get("schema")==API_SCHEMA
          and r.get("sampled_zero_ROI9_source_taxa")==64
          and r.get("source_deterministic_hash_offset")==0
          and r.get("original_source_zero_ROI9_taxa")==2229
          and r.get("original_photo_pixels_downloaded")==0
          and r.get("genuine_nine_fraction_ROI_measured")==0
          and r.get("old_photo_ID_1_231_244_ledger_disjointness_verified") is False,
          "Unsigned or previously classified new iNaturalist source")
    check(source_sha(meta)==r["output_sha256"][Path(meta).name],
          "External source metadata SHA mismatch")
    check(source_sha(atlas)==ATLAS_SHA,"Frozen after56 species atlas SHA mismatch")
    census=pd.read_csv(atlas,dtype={"inat_taxon_id":str})
    check(len(census)==42111 and census.inat_taxon_id.is_unique
          and int(census.after_strict9.eq(0).sum())==2229,
          "42,111 taxon old original photo atlas not intact")
    new=pd.read_csv(meta,dtype=str,keep_default_na=False)
    check(len(new)==r["source_taxa_with_new_post_snapshot_API_photo_opportunities"]
          and len(new)==r["new_independent_original_iNaturalist_photo_IDs_seen"]
          and 0<=len(new)<=64,
          "Candidate count mismatches post-2026-10-28 source receipt")
    if len(new):
        check(new.photo_id.is_unique and new.observation_id.is_unique
              and new.inat_taxon_id.is_unique,
              "photo UUID/taxon repeated in the new source")
        check(new.inat_taxon_id.isin(
            census.loc[census.after_strict9.eq(0),"inat_taxon_id"]).all(),
            "New photo for existing measured taxon is forbidden")
        for t in new.itertuples():
            url=urlsplit(t.photo_url_large)
            check(t.photo_id.isdigit() and UUID.fullmatch(t.observation_id) is not None
                  and t.photo_license in CC and
                  url.scheme=="https"
                  and url.hostname=="inaturalist-open-data.s3.amazonaws.com"
                  and url.path.startswith("/photos/"+t.photo_id+"/")
                  and url.path.lower().endswith((".jpg",".jpeg",".png")),
                  "Invalid independent photo ID/CC license/UUID/source S3 URL")
            when=datetime.fromisoformat(t.source_observation_created_at.replace("Z","+00:00"))
            check(when.tzinfo is not None
                  and when.astimezone(timezone.utc)>=datetime(2026,10,28,tzinfo=timezone.utc),
                  "Old observation outside post-freeze created_at cutoff")
    paths=[Path(v.strip()) for v in Path(history_paths).read_text().splitlines() if v.strip()]
    check(len(paths)==57 and len(set(map(str,paths)))==57,
          "Complete historical original photo registry should contain baseline+56 waves")
    photo_ids=set();obs_ids=set();n=0;input_signatures={}
    frames=[]
    for p in paths:
        check(p.is_file(),"Missing old signed ID ledger "+str(p))
        d=pd.read_parquet(p,columns=["photo_id","observation_id","inat_taxon_id"])
        d=d.fillna("").astype(str)
        check(d.photo_id.ne("").all() and d.photo_id.is_unique
              and not set(d.photo_id).intersection(photo_ids),
              "Historical photo ID collided, not a valid signed union")
        photo_ids.update(d.photo_id)
        obs_ids.update(d.observation_id[d.observation_id.ne("")])
        n+=len(d)
        input_signatures[str(p)]=sha256(p)
        frames.append(d)
    check(n==EXPECTED_ATTEMPTED and len(photo_ids)==EXPECTED_ATTEMPTED,
          "Historical registry differs from 1,231,244 actually attempted original photos")
    new_photo_overlaps=set(new.photo_id).intersection(photo_ids)
    new_obs_overlaps=set(new.observation_id).intersection(obs_ids)
    check(not new_photo_overlaps and not new_obs_overlaps,
          "Newly offered photo/observation is in the actual historical photo ID registry")
    historical=pd.concat(frames,ignore_index=True)
    output.mkdir(parents=True)
    registry=output/"historical_1231244_original_photo_observation_taxon_ID_registry.csv.gz"
    historical.to_csv(registry,index=False,compression="gzip")
    sh=output/"shards";sh.mkdir()
    if len(new):
        new=new.copy()
        new["measurement_id"]="MOIS-POST20261028-"+new.photo_id
        new["image_filename"]=new.measurement_id+".jpg"
        new["shard"]="0"
        p=sh/"shard-000.csv"
        new.to_csv(p,index=False,lineterminator="\n")
        shard_sha={"0":sha256(p)}
    else:shard_sha={}
    manifest={
        "schema":WORKER_SCHEMA,"source_arm":ARM,"colour_measured":False,
        "source_previous_attempts":EXPECTED_ATTEMPTED,
        "source_previous_ID_registry_SHA256":sha256(registry),
        "source_previous_complete_56wave_atlas_SHA256":ATLAS_SHA,
        "source_post_snapshot_API_metadata_SHA256":source_sha(meta),
        "source_post_snapshot_API_receipt_SHA256":source_sha(receipt),
        "source_created_d1":CUTOFF,
        "source_strict_license_taxon_UUID_gate":True,
        "source_prior_no_reused_original_photo_or_observation":True,
        "new_photo_attempts_planned":len(new),
        "work_units":int(len(new)>0),"worker_ids":[0] if len(new)>0 else [],
        "photo_shard_sha256":shard_sha,
        "raw_image_pixels_downloaded":False,"new_nineROI_measured":False,
        "species_wide_fixed_colour_proven":False,
    }
    write_json(output/"manifest.json",manifest)
    result={"schema":"mois_post20261028_new_iNat_frozen_photo_ID_disjointness_gate_v1",
            "source_API_taxa_queried":64,"new_CC_photo_opportunities":len(new),
            "new_photo_IDs_historical_collision_count":len(new_photo_overlaps),
            "new_observation_UUID_historical_collision_count":len(new_obs_overlaps),
            "historical_verified_photo_ID_count":len(photo_ids),
            "historical_source_wave_deltas_checked":len(paths),
            "historical_full_registry_sha256":sha256(registry),
            "historical_inputs_sha256":input_signatures,
            "source_manifest_sha256":sha256(output/"manifest.json"),
            "frozen_source_shard_sha256":shard_sha,
            "new_photo_pixels_downloaded":0,"strict9_measured":0,
            "source_only":True}
    write_json(output/"gate_receipt.json",result)
    print("MOIS_POST_SNAPSHOT_INAT_FULL_1231244_ID_DISJOINTNESS_GATE",
          json.dumps({k:result[k] for k in (
              "new_CC_photo_opportunities","historical_verified_photo_ID_count",
              "new_photo_IDs_historical_collision_count",
              "new_observation_UUID_historical_collision_count")},sort_keys=True),flush=True)
    return result

def audit(source,worker,output):
    source=Path(source);worker=Path(worker);output=Path(output)
    check(not output.exists(),"Refuse to overwrite completed authentic ROI audit")
    m=json.loads((source/"manifest.json").read_text())
    g=json.loads((source/"gate_receipt.json").read_text())
    check(m.get("schema")==WORKER_SCHEMA and m.get("source_arm")==ARM
          and m.get("source_prior_no_reused_original_photo_or_observation") is True
          and m.get("source_previous_attempts")==EXPECTED_ATTEMPTED
          and g.get("historical_verified_photo_ID_count")==EXPECTED_ATTEMPTED
          and g.get("historical_full_registry_sha256")==
            sha256(source/"historical_1231244_original_photo_observation_taxon_ID_registry.csv.gz")
          and g.get("source_manifest_sha256")==sha256(source/"manifest.json")
          and 1<=m.get("new_photo_attempts_planned",0)<=64
          and m.get("work_units")==1,"Missing full ID-gated original photo source")
    a_path=source/"shards"/"shard-000.csv"
    check(sha256(a_path)==m["photo_shard_sha256"]["0"],"Original frozen source shard modified")
    receipt=worker/"shard_000"/"receipt.json"
    terminal=worker/"shard_000"/"photos.csv"
    check(receipt.is_file() and terminal.is_file(),"Missing genuine original image pixel worker")
    r=json.loads(receipt.read_text())
    check(r["schema"]=="mois_full_original93m_real_RGB_photo_unit_v1"
          and r["photo_colours_measured"] is True
          and r["source_sha256"]==sha256(a_path)
          and r["photos_sha256"]==sha256(terminal)
          and r["attempted"]==m["new_photo_attempts_planned"],
          "Worker photo receipt SHA or number of processed pixels mismatched")
    a=pd.read_csv(a_path,dtype=str,keep_default_na=False)
    b=pd.read_csv(terminal,dtype=str,keep_default_na=False)
    check(len(a)==len(b)==m["new_photo_attempts_planned"]
          and a.photo_id.is_unique and a.observation_id.is_unique and a.inat_taxon_id.is_unique
          and b.photo_id.is_unique
          and set(a.photo_id)==set(b.photo_id),
          "Missing/duplicated true original photo worker outcome")
    left=a.set_index("photo_id").sort_index()
    right=b.set_index("photo_id").sort_index()
    for col in ("inat_taxon_id","observation_id","photo_license","photo_url_large"):
        check(left[col].equals(right[col]),"Frozen source provenance changed: "+col)
    check(set(PIX).issubset(b.columns),"Not all nine ROI-v4 measured colour dimensions present")
    rgb=b[PIX].apply(pd.to_numeric,errors="coerce").to_numpy(float)
    strict=b.morph.isin(HUES).to_numpy(bool)&np.isfinite(rgb).all(axis=1)&(rgb>=0).all(axis=1)&(rgb<=1.00001).all(axis=1)&(np.abs(np.nan_to_num(rgb).sum(axis=1)-1)<=.001)
    labels=np.full(len(b),"",dtype=object)
    labels[strict]=class_labels(rgb[strict])
    check(int(b.measurement_status.eq("image_acquisition_failed").sum())==
          r["acquisition_failed"],"Image acquisition failures not counted as technical unknown")
    acquired=~b.measurement_status.eq("image_acquisition_failed")
    check(b.loc[acquired,"image_sha256"].ne("").all(),
          "Downloaded original raw pixels missing SHA256")
    goodsha=b.loc[acquired,"image_sha256"]
    check(goodsha.is_unique,"Byte-identical new images found between two distinct photo URLs")
    registry=pd.read_csv(source/"historical_1231244_original_photo_observation_taxon_ID_registry.csv.gz",
                         dtype=str,keep_default_na=False,
                         usecols=["photo_id","observation_id"])
    check(not set(a.photo_id).intersection(set(registry.photo_id))
          and not set(a.observation_id).intersection(set(registry.observation_id)),
          "Historical photo ID or observation UUID was reused")
    audited=b.copy()
    audited["strict9_original_photo_numeric"]=strict
    audited["dominant_hue9_label"]=labels
    audited.loc[~strict,PIX]=""
    check(audited.loc[~strict,PIX].eq("").all().all(),
          "Unresolved technical image gained invented 9-component hue")
    check(np.allclose(rgb[strict].sum(axis=1),1,atol=.001),
          "strict9 authentic numerical fractions do not sum to 1")
    output.mkdir()
    path=output/"mois_new_inat_post20261028_actual_original_nine_roi.csv.gz"
    audited.to_csv(path,index=False,compression="gzip")
    report={"schema":"mois_post20261028_fullID_gated_actual_RGB9_independent_audit_v1",
            "original_licensed_photo_attempts":len(b),
            "actual_image_downloaded":int(acquired.sum()),
            "actual_genuine_strict9_numeric_success":int(strict.sum()),
            "actual_clear_photo_nine_palette_labels":int((labels!="").sum()),
            "new_first_strict9_nominal_taxa":int(strict.sum()),
            "actual_technical_unknowns":int((~strict).sum()),
            "photo_observation_taxon_collision_count":0,
            "source_photo_manifest_sha256":sha256(source/"manifest.json"),
            "historical_full_ID_registry_sha256":sha256(source/"historical_1231244_original_photo_observation_taxon_ID_registry.csv.gz"),
            "original_photo_worker_receipt_sha256":sha256(receipt),
            "original_photo_worker_terminal_sha256":sha256(terminal),
            "audited_derivative_sha256":sha256(path),
            "new_photo_image_SHA256":dict(zip(b.photo_id,b.image_sha256)),
            "full_historical_raw_image_SHA_collision_check_possible":False,
            "GBIF_col_xr_new_accepted_revalidated":False,
            "official_iNat_42111_original_frame_auto_appended":False,
            "species_wide_fixed_flower_colour_proven":False,
            "GBIF_occurrence_SQL_authenticated_downloaded":False}
    write_json(output/"receipt.json",report)
    print("MOIS_POSTSNAPSHOT_ACTUAL_NEW_INAT_ORIGINAL_FCP_ROI9_AUDIT",
          json.dumps({k:report[k] for k in (
              "original_licensed_photo_attempts","actual_image_downloaded",
              "actual_genuine_strict9_numeric_success",
              "actual_clear_photo_nine_palette_labels","actual_technical_unknowns")},sort_keys=True),flush=True)
    return report

def main():
    p=argparse.ArgumentParser()
    sub=p.add_subparsers(dest="command",required=True)
    a=sub.add_parser("gate")
    for key in ("meta","receipt","atlas","history-paths","output"):
        a.add_argument("--"+key,required=True)
    b=sub.add_parser("audit")
    for key in ("source","worker","output"):
        b.add_argument("--"+key,required=True)
    v=p.parse_args()
    if v.command=="gate":gate(v.meta,v.receipt,v.atlas,v.history_paths,v.output)
    else:audit(v.source,v.worker,v.output)
if __name__=="__main__":main()
