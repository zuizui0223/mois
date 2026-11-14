#!/usr/bin/env python3
"""Freeze GBIF HUMAN_OBSERVATION CC still-image URLs for a separate FCP ROI9 pilot.

This source can NEVER append blindly to the frozen 93,209,351 iNaturalist
photo ID population. Field basis is NOT proof the photo actually shows petals.
Pixels are obtained only by a subsequent pinned FCP ROI-v4 worker after the
independent media and botanical source gates have passed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

from scripts.acquisition.mois_zero9_gbif_other_cc_media_pilot import (
    BASELINE, safe_image, source_sha, require
)
from scripts.acquisition.mois_external_gbif_CC_original_photo_roi_v4 import (
    COHORT, validate_https_public_host
)

SCHEMA="mois_CCGBIF_field_human_original_CC_source_for_FCP_ROI9_v1"
FIELD_SCHEMA="mois_GBIF_COLXR_289taxa_field_original_CC_media_union_v1"
FIELD_FILE="mois_CCGBIF_COLXR_289taxa_human_observation_original_CC_photo_URLs.csv.gz"


def select(field_union,host_allowlist,*,quota=32,exclude_URLs=()):
    root=Path(field_union)
    info=json.loads((root/"receipt.json").read_text())
    require(info.get("schema")==FIELD_SCHEMA
            and info.get("source_accepted_species_without_iNat_numeric_ROI")==289
            and info.get("new_image_pixels_downloaded")==0
            and info.get("new_nine_fraction_ROI9_measurements")==0,
            "GBIF field-only media union must be a signed zero-photo source")
    p=root/FIELD_FILE
    require(p.is_file() and source_sha(p)==info["outputs_sha256"]["field"],
            "Field-only exact COL XR media-source URL hash missing or changed")
    require(1<=quota<=96 and host_allowlist,"Choose an explicitly approved non-iNat HTTPS host")
    permitted={s.lower() for s in host_allowlist}
    blocked=set(exclude_URLs)
    d=pd.read_csv(p,compression="gzip",dtype=str,keep_default_na=False)
    require(d.image_url.is_unique and d.gbif_basisOfRecord.eq("HUMAN_OBSERVATION").all(),
            "Source field-photo identity or independent URL uniqueness lost")
    exact=d.GBIF_media_attached_to_correct_species.str.lower().eq("true")
    require(exact.all(),"Source GBIF occurrence/media to accepted species link drifted")
    chosen=[]
    for row in d.to_dict("records"):
        if row["image_url"] in blocked:continue
        if not safe_image(row["image_url"],row["image_license"]):continue
        try: validate_https_public_host(row["image_url"],permitted)
        except ValueError:continue
        if not row["inat_taxon_id"] or not row["gbif_occurrence_key"] or not row["gbif_dataset_key"]:
            continue
        row["fixed_priority"]=hashlib.sha256((
           "MOIS-GBIF-CC-FIELD-ORIGINAL-NINEROI|"+
           row["inat_taxon_id"]+"|"+row["gbif_occurrence_key"]+"|"+
           row["image_url"]).encode()).hexdigest()
        chosen.append(row)
    if not chosen:
        return pd.DataFrame()
    selected=pd.DataFrame(chosen)
    selected=selected.sort_values("fixed_priority",kind="stable").drop_duplicates(
        "inat_taxon_id").head(quota).reset_index(drop=True)
    require(selected.inat_taxon_id.is_unique and selected.image_url.is_unique and
            selected.gbif_occurrence_key.is_unique,
            "A field photo source taxon/image/GBIF occurrence is repeated")
    selected["measurement_id"]=selected.image_url.map(
        lambda x:"MOIS-GBIF-CC-HUMAN-ROI9-"+hashlib.sha256(x.encode()).hexdigest()[:20])
    selected["image_filename"]=selected.measurement_id+".jpg"
    selected["photo_url_large"]=selected.image_url
    selected["photo_license"]=selected.image_license
    selected["observation_id"]="GBIF-occurrence:"+selected.gbif_occurrence_key
    selected["photo_id"]=selected.image_url.map(
        lambda x:"GBIF-CC-human:"+hashlib.sha256(x.encode()).hexdigest()[:22])
    selected["shard"]="0"
    selected=selected.drop(columns=["fixed_priority"])
    return selected


def run(field_union,output,hosts,*,quota=32,excluded_source=None):
    output=Path(output);require(not output.exists(),"Do not overwrite immutable approved field photo source")
    ex=set()
    if excluded_source is not None:
        p=Path(excluded_source)
        prior=pd.read_csv(p,dtype=str,keep_default_na=False)
        require(prior.image_url.is_unique,"Source external previously attempted photo URLs are not unique")
        ex=set(prior.image_url)
    rows=select(field_union,hosts,quota=quota,exclude_URLs=ex)
    output.mkdir(parents=True)
    if rows.empty:
        state={"schema":SCHEMA,"eligible_direct_HTTPS_original_field_photos":0,
               "GBIF_field_source_SHA256":source_sha(Path(field_union)/"receipt.json"),
               "old_GBIF_photo_URLs_excluded":len(ex),
               "photo_pixels_downloaded":0,
               "new_9fraction_ROI":0,
               "new_photo_run_authorized":False}
        (output/"no_eligible_field_photos.json").write_text(json.dumps(state,indent=2)+"\n")
        print("MOIS_CCGBIF_FIELD_NO_GENUINE_NEW_HOST_APPROVED_IMAGES",
              json.dumps(state,sort_keys=True),flush=True)
        return state
    (output/"shards").mkdir()
    source=output/"shards"/"shard-000.csv"
    rows.to_csv(source,index=False,lineterminator="\n")
    manifest={
        "schema":COHORT,
        "source_arm":"MOIS_CCGBIF_EXPLICIT_HUMAN_OBSERVATION_ONLY",
        "source_field_photo_provenance":SCHEMA,
        "source_GBIF_field_origin_receipt_SHA256":source_sha(Path(field_union)/"receipt.json"),
        "source_real_original_photo_urls_SHA256":source_sha(source),
        "source_historical_original_MOIS_photo_taxon_atlas_SHA256":BASELINE,
        "source_taxon_count":len(rows),
        "original_external_image_count_planned":len(rows),
        "source_manual_approved_public_image_hosts":sorted(set(x.lower() for x in hosts)),
        "source_unique_GBIF_occurrence_per_species":True,
        "source_non_HUMAN_OBSERVATION_photos_allowed":False,
        "prior_external_photo_urls_excluded":len(ex),
        "photo_colour_measured":False,
        "botanical_flower_visibility_checked":False,
        "original_image_sha_overlap_with_MOIS_1231244_prior_pixels_checked":False,
        "downstream_appending_to_iNaturalist_original_frame_allowed":False,
        "raw_image_pixels_persisted":False,
    }
    (output/"manifest.json").write_text(json.dumps(manifest,indent=2,sort_keys=True)+"\n")
    print("MOIS_CCGBIF_FIELD_ONLY_SHA_FROZEN_FCP_ORIGINAL_SOURCE",json.dumps({
       "unique_field_original_photo_taxa":len(rows),
       "permitted_public_HTTPS_hosts":manifest["source_manual_approved_public_image_hosts"],
       "previous_external_URLs_excluded":len(ex),
       "new_numerical_ROI9_measured":0},sort_keys=True),flush=True)
    return manifest


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--field-union",required=True)
    p.add_argument("--output",required=True)
    p.add_argument("--allow-host",action="append",required=True)
    p.add_argument("--quota",type=int,default=32)
    p.add_argument("--exclude-source")
    a=p.parse_args()
    run(a.field_union,a.output,a.allow_host,quota=a.quota,excluded_source=a.exclude_source)
