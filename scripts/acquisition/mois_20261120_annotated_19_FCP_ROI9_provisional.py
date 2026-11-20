#!/usr/bin/env python3
"""Source-pinned 19 visual perianth candidates -> original FCP ROI-v4 technical
measurements. Not botanically verified colour, not fixed species trait.

Preselection was a 440px contact-sheet visual screen of 46 distinct genuinely
flower-annotated original iNaturalist CC photos. Original full source SHA and
per-photo raw SHA are required. All outputs are provisional pixel diagnostics,
never permitted into the global nine-hue accepted-species atlas automatically.
"""
from __future__ import annotations
import argparse,hashlib,io,json
from pathlib import Path
from zipfile import ZipFile
import numpy as np
import pandas as pd

from scripts.acquisition.mois_full_index_real_original_photo_measure import (
    SCHEMA as FCP_SOURCE_SCHEMA, sha as sha_file
)
from scripts.analysis.mois_species_first_palette9 import PALETTE,class_labels

SOURCE_SHA="c2ae31d6d492ba08b07457ff68dfc39c76870723974383219700d7ec7ed5395e"
SOURCE_CSV="flowering_annotation_original_photo_SHA_UUID_sources.csv"
TRIAGE_SHA="5f580ace30f8b359e3b9b7efdc22bdd90f57fb810a09d4d517aca7459399e572"
N_SOURCE=46
N_CANDIDATE=19

def check(ok,msg):
    if not bool(ok): raise ValueError("MOIS_FLOWER_FCP_CANDIDATE19: "+msg)

def sha(raw):
    return hashlib.sha256(raw).hexdigest()

def sources(archive,triage_path):
    check(sha_file(archive)==SOURCE_SHA,"True 46-annotated-photo ZIP digest changed")
    check(sha_file(triage_path)==TRIAGE_SHA,"Independent visual triage source digest changed")
    with ZipFile(archive) as z:
        receipt=json.loads(z.read("receipt.json"))
        check(receipt["new_annotated_original_photo_SHA_checked"]==N_SOURCE,
              "Previous actual original photo count changed")
        for name,digest in receipt["outputs_SHA256"].items():
            check(sha(z.read(name))==digest,"Source member SHA invalid "+name)
        for name,digest in receipt["contact_SHA256"].items():
            check(sha(z.read(name))==digest,"Source contact SHA invalid")
        recorded=pd.read_csv(io.BytesIO(z.read(SOURCE_CSV)),
                             dtype={"photo_id":str,"observation_id":str})
    triage=pd.read_csv(triage_path,dtype={
        "photo_id":str,"observation_uuid":str,"original_SHA256":str
    })
    check(len(recorded)==len(triage)==N_SOURCE and
          recorded.photo_id.is_unique and triage.photo_id.is_unique,
          "46 real original image IDs did not match per-image review")
    d=recorded.merge(triage,on="photo_id",how="inner",validate="one_to_one")
    check(len(d)==N_SOURCE and
          d.observation_id.eq(d.observation_uuid).all() and
          d.original_image_sha256.eq(d.original_SHA256).all() and
          d.candidate_species.eq(d.species).all() and
          d.original_url.eq(d.original_photo_url).all(),
          "Original per-image SHA/observation/taxon provenance inconsistent")
    check(d.iNaturalist_flowers_annotation_12_13_in_response.eq(True).all(),
          "Source 12|13 flowering annotation is not verified for every photo")
    count=d.thumbnail_visual_triage.value_counts().to_dict()
    check(count=={"clear_candidate":19,"uncertain":8,"not_measurable":19},
          "Independent source visual triage unexpectedly changed")
    keep=d.loc[d.thumbnail_visual_triage.eq("clear_candidate")].copy()
    check(len(keep)==N_CANDIDATE and keep.candidate_species.nunique()==9,
          "19 candidate visibly floral source photos or nine species changed")
    check(keep.fullres_botanical_petal_verification.eq("NOT_DONE").all()
          and keep.nine_hue_measurement_eligible.eq(False).all()
          and keep.manual_organ_mask_verified.eq(False).all(),
          "A provisional 440px visual candidate was misrepresented as botanical truth")
    return keep.reset_index(drop=True),receipt

def prep(args):
    target=Path(args.out)
    check(not target.exists(),"Cannot overwrite signed FCP source manifest")
    d,prior=sources(args.source,args.triage)
    out=d[["photo_id","photo_url_large","photo_license","observation_id",
           "inat_taxon_id","candidate_species","original_image_sha256",
           "thumbnail_visual_triage"]].copy()
    out["measurement_id"]=d.photo_id.map(lambda x:
        "MOIS-"+sha(("FLOWERING-CANDIDATE-19|"+x).encode()).upper()[:24])
    out["image_filename"]=out.measurement_id+".jpg"
    out["shard"]="0"
    check(out.photo_id.is_unique and out.observation_id.is_unique
          and out.measurement_id.is_unique
          and out.photo_license.isin(
              ["cc0","cc-by","cc-by-sa","cc-by-nc","cc-by-nc-sa"]).all(),
          "Invalid accepted original image source/license/UUID")
    source_cols=["measurement_id","image_filename","photo_id",
       "observation_id","inat_taxon_id","candidate_species",
       "photo_url_large","photo_license","shard",
       "original_image_sha256","thumbnail_visual_triage"]
    out=out[source_cols]
    shardfolder=target/"shards"
    shardfolder.mkdir(parents=True)
    shard=shardfolder/"shard-000.csv"
    out.to_csv(shard,index=False,lineterminator="\n")
    manifest={
      "schema":FCP_SOURCE_SCHEMA,"source_arm":"MOIS_20261120_FLOWERS_ANNOTATED_THUMBNAIL_SCREEN_ONLY",
      "colour_measured":False,"work_units":1,"worker_ids":[0],
      "photo_shard_sha256":{"0":sha_file(shard)},
      "new_photo_attempts_planned":len(out),
      "source_previous_actual_CC_original48_ZIP_SHA256":SOURCE_SHA,
      "source_visual_thumbnail_triage_SHA256":TRIAGE_SHA,
      "source_flowers_annotation_12_13_all":True,
      "source_original_full_bytes_SHA_audited":True,
      "src_440px_candidate_not_verified_petal_organ":True,
      "botanically_confirmed_flowers":0,
      "global_photo9_atlas_append_allowed":False,
      "source_no_floral_hue_used_to_choose_individual_photo":True,
      "new_pixel_ROIs_measured":False
    }
    (target/"manifest.json").write_text(json.dumps(manifest,indent=2,sort_keys=True)+"\n")
    audit={"n_prior_source_rows":N_SOURCE,
       "n_thumbnail_clearly_floral_candidates":N_CANDIDATE,
       "n_source_species":out.candidate_species.nunique(),
       "provisional_source_no_true_perianth_validation":True,
       "source_sha256":SOURCE_SHA,"triage_sha256":TRIAGE_SHA,
       "shard_sha256":sha_file(shard),
       "manifest_sha256":sha_file(target/"manifest.json")}
    (target/"source_receipt.json").write_text(json.dumps(audit,indent=2)+"\n")
    print("MOIS_FCP19_FROZEN_SOURCE_SUCCESS",json.dumps(audit,sort_keys=True),flush=True)

def audit(args):
    root=Path(args.out)
    check(not root.exists(),"Refuse overwriting actual pixel ROI diagnostic")
    d,rec=sources(args.source,args.triage)
    src=Path(args.manifest)
    manifest=json.loads((src/"manifest.json").read_text())
    shard=src/"shards/shard-000.csv"
    check(manifest["schema"]==FCP_SOURCE_SCHEMA
          and manifest["src_440px_candidate_not_verified_petal_organ"]
          and sha_file(shard)==manifest["photo_shard_sha256"]["0"],
          "Actual FCP source shard not signed or falsely botanically confirmed")
    expected=pd.read_csv(shard,dtype=str,keep_default_na=False)
    folder=Path(args.worker)/"shard_000"
    receipt=json.loads((folder/"receipt.json").read_text())
    table=folder/"photos.csv"
    check(receipt["schema"]=="mois_full_original93m_real_RGB_photo_unit_v1"
          and receipt["photo_colours_measured"]
          and receipt["attempted"]==N_CANDIDATE
          and receipt["source_sha256"]==sha_file(shard)
          and receipt["photos_sha256"]==sha_file(table),
          "Actual FCP original pixel worker source receipt did not verify")
    measured=pd.read_csv(table,dtype=str,keep_default_na=False)
    check(len(expected)==len(measured)==N_CANDIDATE
          and expected.measurement_id.is_unique
          and measured.measurement_id.is_unique
          and set(measured.measurement_id)==set(expected.measurement_id),
          "FCP terminal source measurement IDs do not match all 19 source candidates")
    out=expected.merge(measured,on="measurement_id",how="inner",validate="one_to_one")
    for k in ["photo_id","observation_id","inat_taxon_id","photo_license"]:
        orig=out[k+"_x"] if k+"_x" in out else out[k]
        reported=out[k+"_y"] if k+"_y" in out else out[k]
        check(orig.astype(str).equals(reported.astype(str)),
              "FCP worker changed original image, photo license or observation identity: "+k)
    fractions=["flower_fraction_"+k for k in PALETTE]
    check(set(fractions).issubset(out.columns),
          "Measured nine-palette original pixel ROI fractions missing")
    numeric=out[fractions].apply(pd.to_numeric,errors="coerce").to_numpy(float)
    four=out.morph.isin(["white","yellow_orange","red_pink","blue_purple"]).to_numpy(bool)
    good=(four & np.isfinite(numeric).all(axis=1)
          & (numeric>=0).all(axis=1)&(numeric<=1.001).all(axis=1)
          & (abs(np.nansum(numeric,axis=1)-1)<.001))
    check(int(good.sum())==receipt["classifiable"],
          "Strict measured technical pixel denominator drift")
    label=np.full(len(out),"",dtype=object)
    if good.any():label[good]=class_labels(numeric[good])
    out["numeric9_ROI_fraction_valid"]=good
    out["provisional_photo9_dominant_label"]=label
    out["manual_botanical_petal_organ_verified"]=False
    out["species_fixed_petal_colour_inferred"]=False
    out["analysis_sample_inclusion_authorized"]=False
    out["photo_FCP_pixel_ROI_not_necessarily_true_petal"]=True
    root.mkdir(parents=True)
    file=root/"source19_unverified_corolla_FCP_true_RGB9_diagnostics.csv"
    out.to_csv(file,index=False,lineterminator="\n")
    report={"schema":"MOIS_FCP19_SOURCE_BACKED_BOTANICALLY_UNVERIFIED_DIAGNOSTIC",
       "source_true_original_flowers_annotated_ZIP_SHA256":SOURCE_SHA,
       "source_thumbnail_screen_SHA256":TRIAGE_SHA,
       "FCP_method_from_worker":receipt["method"],
       "source_photos":N_CANDIDATE,
       "n_original_worker_downloaded":int(receipt["downloaded"]),
       "n_original_worker_classifiable_four_morph":int(receipt["classifiable"]),
       "n_9_component_numeric_ROI":int(good.sum()),
       "n_provisional_clear_9hue_labels":int((label!="").sum()),
       "n_acquisition_failed":int(receipt["acquisition_failed"]),
       "n_not_classifiable":int(receipt["unclassifiable_acquired"]),
       "all_botanical_petal_verification_still_false":True,
       "no_global_photo9_atlas_append":True,
       "all_original_prior_photo_ID_history_after_base2026_not_certified":True,
       "FCP_output_csv_sha256":sha_file(file)}
    (root/"receipt.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print("MOIS_FCP19_REAL_RGB9_UNVERIFIED_PIXEL_DIAGNOSTIC",
          json.dumps({k:v for k,v in report.items() if k!="FCP_method_from_worker"},
                     sort_keys=True),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    cmd=p.add_subparsers(dest="cmd",required=True)
    first=cmd.add_parser("source")
    first.add_argument("--source",required=True)
    first.add_argument("--triage",required=True)
    first.add_argument("--out",required=True)
    second=cmd.add_parser("audit")
    for key in ("source","triage","manifest","worker","out"):
        second.add_argument("--"+key,required=True)
    args=p.parse_args()
    if args.cmd=="source":prep(args)
    else:audit(args)
