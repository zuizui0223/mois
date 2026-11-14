#!/usr/bin/env python3
"""Fresh POST-2026-10-28 iNat originals from 1909 previously unqueried MOIS taxa.

Three stages: source (post-freeze authentic photo metadata + all 1,231,255
prior original IDs), audit (actual pinned FCP YOLO11/EfficientSAM nine fractions
only), append (42,111 taxon atlas only after signed full photo+observation QC).
No historical image is re-downloaded. Source-only candidate counts are NOT
measured photos. Fresh original image SHA overlap against all historical raw
image bytes is NOT proven by metadata IDs alone.
"""
from __future__ import annotations
import argparse,hashlib,json
from datetime import datetime,timezone
from pathlib import Path
from urllib.parse import urlsplit

import numpy as np
import pandas as pd

from scripts.analysis.mois_species_first_palette9 import PALETTE,COLUMNS,class_labels,sha256
from scripts.analysis.mois_species_first_wave3_incremental import add_to_catalogue
from scripts.acquisition.mois_zero9_recent_inat_new_created_photo_source import CC
from scripts.acquisition.mois_full_index_real_original_photo_measure import SCHEMA as WORKER_SOURCE_SCHEMA

SCHEMA="mois_remaining1909_new_postsnapshot_original_photo_ROI9_pipeline_v1"
SOURCE_SCHEMA="mois_post20261028_new_inat_other1909_zero_numeric_taxa_source_union_v1"
SOURCE_ARM="MOIS_NEW_INAT_POST20261028_REMAINING1909_BREADTH"
ATLAS_SHA="8d9a683fffd73c3a67fe63552e1490d4f629b688bb1eba0d1b4b1bafd143c22f"
PRIOR_ATTEMPTS=1231255
PRIOR_NUMERIC=467539
PRIOR_SUPPORTED=39885
PRIOR_HUE=39147
BASE_OLD_PHOTOS=1231244
MAX_PHOTOS=64
HUES=("white","yellow_orange","red_pink","blue_purple")


def require(value,reason):
    if not bool(value):raise ValueError("MOIS_FRESH1909: "+reason)


def load_prior(atlas_path,history_path,new11_delta_path):
    require(sha256(atlas_path)==ATLAS_SHA,
            "Signed 42111 species atlas after eleven true new original photos drifted")
    atlas=pd.read_csv(atlas_path,compression="infer",dtype={"inat_taxon_id":str})
    require(len(atlas)==42111 and atlas.inat_taxon_id.is_unique
            and int(atlas.after_strict9.gt(0).sum())==PRIOR_SUPPORTED
            and int(atlas.dominant_hue9.fillna("").ne("").sum())==PRIOR_HUE,
            "Authentic iNaturalist original photo9 taxon coverage changed")
    prior=pd.read_csv(history_path,compression="gzip",
                      usecols=["photo_id","observation_id"],dtype=str,keep_default_na=False)
    delta=pd.read_parquet(new11_delta_path,columns=["photo_id","observation_id"])
    delta=delta.fillna("").astype(str)
    require(len(prior)==BASE_OLD_PHOTOS and len(delta)==11 and
            prior.photo_id.is_unique and delta.photo_id.is_unique,
            "Original 1231244 + 11 photo ID old source census absent")
    both=pd.concat([prior,delta],ignore_index=True)
    require(len(both)==PRIOR_ATTEMPTS and both.photo_id.is_unique and
            both.observation_id.loc[both.observation_id.ne("")].is_unique,
            "Historical original photo or independent observation UUID duplicated")
    return atlas,both


def pick_source(photo_source,prior_atlas,prior,*,quota=MAX_PHOTOS):
    photo_source=Path(photo_source)
    require(1<=quota<=MAX_PHOTOS,"Bound actual original new photo cohort to <=64")
    m=json.loads((photo_source/"receipt.json").read_text())
    require(m.get("schema")==SOURCE_SCHEMA
            and m.get("new_source_taxa_API_queried")==1909
            and m.get("original_source_zero_numeric_at_pre20261028")==2229
            and m.get("actual_new_photo_pixels_downloaded")==0
            and m.get("new_strict9_numerical_ROI")==0,
            "Post-snapshot source not authenticated as 1909 new metadata taxa")
    p=photo_source/"MOIS_post_snapshot_new1909_unmeasured_taxa_CC_photo_ID_opportunities.csv"
    require(sha256(p)==m["results_SHA256"][p.name],
            "New source original photo-ID metadata SHA changed")
    candidates=pd.read_csv(p,dtype=str,keep_default_na=False)
    require(candidates.photo_id.is_unique and candidates.observation_id.is_unique
            and candidates.inat_taxon_id.is_unique and
            candidates.photo_id.str.fullmatch(r"\d+").all(),
            "New source metadata photo, observation or taxon repeated")
    remain=set(prior_atlas.loc[prior_atlas.after_strict9.eq(0),"inat_taxon_id"])
    priorphotos=set(prior.photo_id.astype(str))
    priorobs=set(prior.observation_id.astype(str))-{""}
    frames=[]
    for r in candidates.itertuples(index=False):
        if r.inat_taxon_id not in remain or r.photo_id in priorphotos or r.observation_id in priorobs:
            continue
        u=urlsplit(r.photo_url_large)
        if r.photo_license not in CC or u.scheme!="https" or u.hostname!="inaturalist-open-data.s3.amazonaws.com":
            continue
        if not u.path.startswith("/photos/"+str(r.photo_id)+"/"):continue
        try:
            created=datetime.fromisoformat(str(r.source_observation_created_at).replace("Z","+00:00"))
            if created.tzinfo is None or created.astimezone(timezone.utc)<datetime(2026,10,28,tzinfo=timezone.utc):
                continue
        except (ValueError,TypeError):continue
        frames.append(r._asdict())
    if not frames:return pd.DataFrame(),int(len(candidates))
    allowed=pd.DataFrame(frames)
    allowed["priority"]=allowed.apply(
        lambda r:hashlib.sha256((
            "MOIS-POSTSNAP-REAL-1909-NEW-FLOWER-ORIGINAL|"+r.photo_id+"|"+
            r.observation_id).encode()).hexdigest(),axis=1)
    allowed=allowed.sort_values("priority",kind="stable").head(quota).drop(columns="priority")
    require(allowed.inat_taxon_id.is_unique and allowed.photo_id.is_unique
            and allowed.observation_id.is_unique,
            "Some current original photo source candidate repeated")
    return allowed.reset_index(drop=True),int(len(candidates))


def source(photo_source,atlas_path,historical_path,new11_path,out,*,quota=MAX_PHOTOS):
    out=Path(out);require(not out.exists(),"Never overwrite frozen new original photo source")
    atlas,previous=load_prior(atlas_path,historical_path,new11_path)
    selected,seen=pick_source(photo_source,atlas,previous,quota=quota)
    out.mkdir(parents=True)
    if selected.empty:
        report={"schema":SCHEMA,"stage":"no_unused_photo_source",
                "input_candidate_photo_IDs":seen,"independent_new_eligible_photos":0,
                "all_1231255_original_photo_ID_and_observation_UUID_checked":True,
                "new_raw_photo_pixels_downloaded":0,"new_numeric9_ROI":0}
        (out/"receipt.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
        print("MOIS_POSTSNAP_1909_NO_NEW_DISJOINT_ORIGINAL_PHOTOS",json.dumps(report),flush=True)
        return report
    selected=selected.copy()
    selected["measurement_id"]="MOIS-FRESH1909-"+selected.photo_id
    selected["image_filename"]=selected.measurement_id+".jpg"
    selected["shard"]="0"
    (out/"shards").mkdir()
    path=out/"shards"/"shard-000.csv"
    selected.to_csv(path,index=False,lineterminator="\n")
    manifest={
        "schema":WORKER_SOURCE_SCHEMA,
        "source_arm":SOURCE_ARM,
        "colour_measured":False,
        "work_units":1,"worker_ids":[0],
        "photo_shard_sha256":{"0":sha256(path)},
        "new_photo_attempts_planned":len(selected),
        "source_previous_attempts":PRIOR_ATTEMPTS,
        "source_original_42111_photo9_atlas_SHA256":ATLAS_SHA,
        "source_historical_original_photo_and_observation_IDs_checked":True,
        "source_nominal_taxa":42111,
        "source_original_full_photo_candidates":93209351,
        "source_new_api_metadata_candidates":seen,
        "source_strict_CC_new_observation_created_date_verified":True,
        "raw_image_pixels_downloaded":False,
        "numeric_nine_fraction_ROIs_measured":False,
        "species_wide_photo_colour_fixed_proven":False,
    }
    (out/"manifest.json").write_text(json.dumps(manifest,indent=2,sort_keys=True)+"\n")
    report={"schema":SCHEMA,"stage":"new_real_photo_ID_source_frozen",
            "input_candidate_photo_IDs":seen,
            "independent_new_eligible_photos":len(selected),
            "source_manifest_SHA256":sha256(out/"manifest.json"),
            "source_photo_shard_SHA256":sha256(path),
            "source_old_photo_ID_overlap":0,"source_old_observation_UUID_overlap":0,
            "new_raw_photo_pixels_downloaded":0,"new_numeric9_ROI":0}
    (out/"receipt.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print("MOIS_POSTSNAP_NEW1909_SIGNED_FRESH_PHOTO_SOURCE",json.dumps(report,sort_keys=True),flush=True)
    return report


def genuine_frame(source_root,worker_root):
    src=Path(source_root);w=Path(worker_root)
    m=json.loads((src/"manifest.json").read_text())
    require(m.get("schema")==WORKER_SOURCE_SCHEMA and m.get("source_arm")==SOURCE_ARM
            and m.get("source_previous_attempts")==PRIOR_ATTEMPTS
            and 0<m.get("new_photo_attempts_planned",0)<=MAX_PHOTOS
            and m.get("source_historical_original_photo_and_observation_IDs_checked") is True,
            "Real 1909 next photo source manifold not signed")
    source_path=src/"shards"/"shard-000.csv"
    measured=w/"shard_000"/"photos.csv"
    receipt=w/"shard_000"/"receipt.json"
    require(sha256(source_path)==m["photo_shard_sha256"]["0"],
            "New source photo URL/ID roster SHA changed")
    r=json.loads(receipt.read_text())
    require(r.get("schema")=="mois_full_original93m_real_RGB_photo_unit_v1"
            and r.get("photo_colours_measured") is True
            and r.get("attempted")==m["new_photo_attempts_planned"]
            and r.get("source_sha256")==sha256(source_path)
            and r.get("photos_sha256")==sha256(measured),
            "Real FCP original colour terminal SHA or source identity changed")
    original=pd.read_csv(source_path,dtype=str,keep_default_na=False)
    terminal=pd.read_csv(measured,dtype=str,keep_default_na=False)
    require(len(original)==len(terminal)==m["new_photo_attempts_planned"]
            and original.photo_id.is_unique and terminal.photo_id.is_unique
            and original.observation_id.is_unique and original.inat_taxon_id.is_unique
            and set(original.photo_id)==set(terminal.photo_id),
            "New original photo worker missing or duplicating one sourced image")
    for c in ("inat_taxon_id","observation_id","photo_url_large","photo_license"):
        require(original.set_index("photo_id").sort_index()[c].equals(
                terminal.set_index("photo_id").sort_index()[c]),
                "Real original photo source taxon/obs/URL/license changed")
    columns=["flower_fraction_"+p for p in PALETTE]
    vals=terminal[columns].apply(pd.to_numeric,errors="coerce").to_numpy(float)
    good=(terminal.morph.isin(HUES).to_numpy(bool)&
          np.isfinite(vals).all(axis=1)&(vals>=0).all(axis=1)&
          (vals<=1.00001).all(axis=1)&
          (abs(np.nan_to_num(vals).sum(axis=1)-1)<=.001))
    require(int(good.sum())==r.get("classifiable",0),
            "Original numerically classifiable 9-component ROI mismatch against terminal")
    labels=np.full(len(terminal),"",dtype=object)
    if good.any():labels[good]=class_labels(vals[good])
    require(terminal.loc[good,"image_sha256"].ne("").all(),
            "Real source image SHA missing for successful nine-palette ROI")
    frame=terminal[["inat_taxon_id","photo_id","observation_id","image_sha256"]].copy()
    frame["strict9"]=good;frame["dominant_label9"]=labels
    for i,c in enumerate(PALETTE):
        frame["strict9_"+c]=np.where(good,vals[:,i],np.nan)
    require(frame.photo_id.is_unique and frame.observation_id.is_unique
            and frame.inat_taxon_id.is_unique,
            "Duplicate observation or taxon among post-snapshot true photo terminals")
    return frame,{"raw_worker_receipt_sha256":sha256(receipt),
                  "worker_photo_terminal_sha256":sha256(measured),
                  "source_photo_SHA256":sha256(source_path)}


def audit(source_root,worker_root,out):
    out=Path(out);require(not out.exists(),"Refuse overwrite true photo ROI audit")
    frame,proof=genuine_frame(source_root,worker_root)
    out.mkdir(parents=True)
    csv=out/"mois_new1909_genuine_original_photo9_and_unknowns.csv.gz"
    frame.to_csv(csv,index=False,compression="gzip")
    result={"schema":SCHEMA,"stage":"actual_new_independent_FCP_nineROI_audited",
        "new_authentic_original_photo_attempts":len(frame),
        "new_original_photo_numeric9_success":int(frame.strict9.sum()),
        "new_clear_dominant_hue9":int(frame.dominant_label9.ne("").sum()),
        "new_technical_unknowns":int((~frame.strict9).sum()),
        "new_source_unique_taxa":frame.inat_taxon_id.nunique(),
        "image_SHA_collisions_with_all_1231255_prior_originals_certified":False,
        "source_and_worker_photo_ID_UUID_collisions":0,
        "photo_pixels_persisted":False,
        "GBIF_occurrence_hue_measured":False,
        "original_source_and_worker_SHA256":proof,
        "derived_numeric9_csv_SHA256":sha256(csv)}
    (out/"receipt.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print("MOIS_REAL_NEW1909_SOURCE_FCP_PHOTO_NINE_ROI_AUDIT",json.dumps(result,sort_keys=True),flush=True)
    return result


def append(atlas_path,old_history,new11_delta,source_root,worker_root,audit_root,out):
    out=Path(out);require(not out.exists(),"No silent repeated original photo append")
    atlas,history=load_prior(atlas_path,old_history,new11_delta)
    frame,proof=genuine_frame(source_root,worker_root)
    r=json.loads((Path(audit_root)/"receipt.json").read_text())
    csv=Path(audit_root)/"mois_new1909_genuine_original_photo9_and_unknowns.csv.gz"
    require(r["schema"]==SCHEMA and
            r["stage"]=="actual_new_independent_FCP_nineROI_audited"
            and r["derived_numeric9_csv_SHA256"]==sha256(csv),
            "Independent signed real original RGB9 audit missing")
    actual=pd.read_csv(csv,compression="gzip",dtype={"inat_taxon_id":str,"photo_id":str,"observation_id":str})
    require(set(actual.photo_id)==set(frame.photo_id) and
            int(actual.strict9.sum())==int(frame.strict9.sum()) and
            actual.dominant_label9.fillna("").tolist()==frame.dominant_label9.tolist(),
            "Actual audited nine-palette numeric colours drifted from original worker")
    require(not set(frame.photo_id).intersection(set(history.photo_id))
            and not set(frame.observation_id).intersection(set(history.observation_id)),
            "New iNaturalist photo ID or observation UUID was already acquired")
    require(atlas.set_index("inat_taxon_id").loc[frame.inat_taxon_id,
               "after_strict9"].eq(0).all(),
            "New post-snapshot source previously had measured ROI9")
    new=add_to_catalogue(atlas,frame,expected_photos=len(frame))
    success=int(frame.strict9.sum())
    clear=int(frame.dominant_label9.ne("").sum())
    require(int(new.after_strict9.gt(0).sum())==PRIOR_SUPPORTED+success and
            int(new.dominant_hue9.fillna("").ne("").sum())==PRIOR_HUE+clear,
            "Original photo9 numeric/trait taxon atlas append cannot be audited")
    out.mkdir(parents=True)
    p1=out/"mois_42111_taxa_postsnap_new1909_original_photo9.csv.gz"
    p2=out/"mois_postsnap1909_authentic_original_photo9_delta.parquet"
    new.to_csv(p1,index=False,compression="gzip")
    frame.to_parquet(p2,index=False)
    report={"schema":SCHEMA,"stage":"full_42111_taxon_atlas_new1909_true_photo9_appended",
       "new_original_photo_attempts":len(frame),
       "new_strict9_numeric_original_photos":success,
       "new_clear_photo_hue_taxa":clear,
       "cumulative_original_photo_attempts":PRIOR_ATTEMPTS+len(frame),
       "cumulative_original_nine_ROI_photos":PRIOR_NUMERIC+success,
       "nominal_taxa_with_genuine_nineROI":PRIOR_SUPPORTED+success,
       "nominal_taxa_with_clear_sampled_photo_hue":PRIOR_HUE+clear,
       "nominal_taxa_without_numeric_nineROI":42111-PRIOR_SUPPORTED-success,
       "prior_original_photo_and_observation_UUID_collisions":0,
       "full_historical_image_content_SHA_exclusion_proven":False,
       "source_original_93m_frame_unchanged":True,
       "GBIF_new_accepted_names_not_yet_reconciled":True,
       "signed_audit_receipt_SHA256":sha256(Path(audit_root)/"receipt.json"),
       "output_SHA256":{str(p.name):sha256(p) for p in (p1,p2)}}
    (out/"receipt.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print("MOIS_REAL_NEW1909_ADDITIONAL_ORIGINAL_NINE_PALETTE_ATLAS_APPEND",
          json.dumps(report,sort_keys=True),flush=True)
    return report


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    sub=p.add_subparsers(dest="mode",required=True)
    a=sub.add_parser("source")
    for k in ("candidates","atlas","historical","new11","output"):
        a.add_argument("--"+k,required=True)
    a.add_argument("--quota",type=int,default=MAX_PHOTOS)
    b=sub.add_parser("audit")
    for k in ("source","workers","output"):b.add_argument("--"+k,required=True)
    c=sub.add_parser("append")
    for k in ("atlas","historical","new11","source","workers","audit","output"):
        c.add_argument("--"+k,required=True)
    args=p.parse_args()
    if args.mode=="source":source(args.candidates,args.atlas,args.historical,args.new11,args.output,quota=args.quota)
    elif args.mode=="audit":audit(args.source,args.workers,args.output)
    else:append(args.atlas,args.historical,args.new11,args.source,args.workers,args.audit,args.output)
