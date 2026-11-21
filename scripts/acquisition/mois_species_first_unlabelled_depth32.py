#!/usr/bin/env python3
"""Refresh original-photo fallback depth from signed 93.2m metadata for unlabelled MOIS taxa.

The 1,181,146-photo baseline remains immutable. We scan actual original-photo
metadata once per 16 disjoint source groups, selecting up to 32 independently
observed CC-licensed original photos *per currently unlabelled taxon*.
This is source-only, with no image-colour inference. The eventual acquisition
must independently exclude all recently attempted photo/observation IDs too.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from scripts.acquisition.mois_species_first_zero9_source import (
    PHOTO_SHA, group_source, output_rows, check
)
from scripts.acquisition.mois_onepass_full93m_target_tiers import sha256

SCHEMA="mois_species_first_93m_current_unlabelled_depth32_original_photo_source_v1"
ATLAS_TWELFTH_SHA="48a6b76750da2e19227bc3b21bf9b1120216fc79bb7ff6c28c1bac5fb7650989"
FRAME=42111
BASE_PHOTOS=1181146

def choose_unlabelled(atlas, group):
    """Include never-strict9 and numeric-but-mixed, never already clear hues."""
    check(0<=group<16 and len(atlas)==FRAME and atlas.inat_taxon_id.is_unique,
          "Latest known 42111 original-photo taxon atlas incomplete")
    needed={"inat_taxon_id","hue_status","dominant_hue9"}
    check(needed <= set(atlas.columns), "Missing botanical taxon/hue status")
    c=atlas.copy()
    c.inat_taxon_id=c.inat_taxon_id.fillna("").astype(str)
    c.hue_status=c.hue_status.fillna("").astype(str)
    c.dominant_hue9=c.dominant_hue9.fillna("").astype(str)
    unknown=c.dominant_hue9.eq("") & c.hue_status.isin(
        ("no_strict9","mixed_image_unresolved"))
    mask=c.inat_taxon_id.map(lambda t:t.isdecimal() and int(t)%128//8==group)
    out=set(c.loc[unknown&mask,"inat_taxon_id"])
    check(all(int(k)%128//8==group for k in out),
          "A source taxon escaped its signed disjoint original group")
    return out

def run(registry, group_receipt, archives, base_photos, current_atlas,
        output,*,group=0,limit=32,verify_sha=True):
    output=Path(output)
    check(not output.exists() and 0<=group<16 and 11<=limit<=48,
          "Depth refresh must be bigger than 10; output must not exist")
    if verify_sha:
        check(sha256(base_photos)==PHOTO_SHA,
              "Signed original 1181146 photo ID ledger was modified")
        check(sha256(current_atlas)==ATLAS_TWELFTH_SHA,
              "Signed post-twelfth 42111 source photo-hue atlas changed")
    prior=pd.read_parquet(base_photos,columns=["photo_id","observation_id"])
    check(len(prior)==BASE_PHOTOS and prior.photo_id.is_unique,
          "Base original photo attempts are not distinct")
    ids=set(prior.photo_id.fillna("").astype(str))
    obs=set(prior.observation_id.fillna("").astype(str))-{""}
    atlas=pd.read_csv(current_atlas,compression="infer",dtype={"inat_taxon_id":str})
    wanted=choose_unlabelled(atlas,group)
    selections,stats=group_source(registry,group_receipt,archives,group,
                                  wanted,ids,obs,limit=limit)
    queue=output_rows(selections)
    if not queue.empty:
        check(queue.photo_id.is_unique and queue.observation_id.is_unique and
              queue.groupby("inat_taxon_id").size().max()<=limit and
              queue.inat_taxon_id.isin(wanted).all(),
              "Source candidates are repeated or reintroduce a hue-labelled taxon")
        check(not set(queue.photo_id)&ids and not set(queue.observation_id)&obs,
              "Baseline original photo/observation was reselected")
    first=queue.loc[queue.source_observation_rank==1].copy()
    check(first.inat_taxon_id.is_unique and first.photo_id.is_unique and
          first.observation_id.is_unique,
          "Disjoint first photo per unlabelled source species selection failed")
    output.mkdir(parents=True)
    prefix=f"mois_unlabelled_depth{limit}_group{group:02d}"
    allfile=output/(prefix+"_original_photo_options.csv.gz")
    firstfile=output/(prefix+"_first_photos.csv")
    queue.to_csv(allfile,index=False,compression="gzip")
    first.to_csv(firstfile,index=False)
    still_no=sorted(wanted-set(selections),key=int)
    (output/"unlabelled_taxa_without_new_original_photos.json").write_text(
        json.dumps(still_no)+"\n")
    note={
        "schema":SCHEMA,
        "input_original_photo_metadata_count":93209351,
        "source_total_nominal_taxa":FRAME,
        "source_group":group,
        "source_target_unlabelled_taxa":len(wanted),
        "target_taxa_with_one_or_more_eligible_original_photo":len(selections),
        "target_taxa_with_no_eligible_original_photo":len(still_no),
        "candidate_original_photos_all_depth":len(queue),
        "source_max_observations_per_taxon":limit,
        "frozen_base_original_photo_IDs_excluded":BASE_PHOTOS,
        "recent_incremental_photo_IDs_not_excluded_in_source":True,
        "recheck_all_wave_photo_and_observation_IDs_before_pixels":True,
        "source_signed_full93m_original_group_receipt_SHA256":sha256(group_receipt),
        "source_signed_full93m_metadata_registry_SHA256":sha256(registry),
        "target_current_wave12_ninehue_atlas_SHA256":sha256(current_atlas),
        "source_group_photo_rows_scanned":int(stats["source_rows"]),
        "new_original_photo_downloads":0,
        "new_nine_palette_ROIs":0,
        "GBIF_occurrences_downloaded":0,
        "outputs_SHA256":{f.name:sha256(f) for f in (allfile,firstfile)},
    }
    (output/"receipt.json").write_text(json.dumps(note,indent=2,sort_keys=True)+"\n")
    print("MOIS_FULL93M_DEPTH32_UNLABELLED_SOURCE",json.dumps(note,sort_keys=True),flush=True)
    return note

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for arg in ("registry","group-receipt","archives","base-photos",
                "current-atlas","output"):
        p.add_argument("--"+arg,required=True)
    p.add_argument("--group",type=int,required=True)
    p.add_argument("--limit",type=int,default=32)
    a=p.parse_args()
    run(a.registry,a.group_receipt,a.archives,a.base_photos,a.current_atlas,
        a.output,group=a.group,limit=a.limit)

if __name__=="__main__":
    main()
