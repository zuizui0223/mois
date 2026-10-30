#!/usr/bin/env python3
"""Frozen ROI9 photo/observation identity for MOIS automatic bias work.

No human annotation, no image pixels, no botanical truth inference.
"""
from __future__ import annotations
import hashlib
from pathlib import Path

import pandas as pd

SEED="MOIS_20261029_PREOUTCOME_OBSERVATION_FRAME_BIAS_AUDIT_V1"
ARMS=(
    "historical_K1_FCP_source_heterogeneous",
    "bias_refill_2048", "p2_first_4096", "bulk128_131072",
)
EXPECTED_PHOTOS=1181146
EXPECTED_SUCCESS=451643
COLS=(
    "photo_id","observation_id","inat_taxon_id","source_design",
    "observed_on","latitude","longitude","strict9",
)

def require(ok,message):
    if not bool(ok):
        raise ValueError("MOIS_AUTO_OBSERVATION_FRAME: "+message)

def sha_file(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda:f.read(4*1024*1024),b""):
            h.update(block)
    return h.hexdigest()

def canonical_rank(observation_id,photo_id):
    return hashlib.sha256(
        (SEED+"|canonical|"+str(observation_id)+"|"+str(photo_id)).encode()
    ).hexdigest()

def coarse_zone(zone):
    if zone=="N_tropical":return "N_tropical"
    if zone=="S_tropical":return "S_tropical"
    if zone=="N_temperate":return "N_extratropical"
    if zone=="S_temperate":return "S_extratropical"
    return "polar_or_unknown"

def observation_frame(photo):
    """One hash-minimized original photo per observation (outcome-blind).

    Retain original source-row indices so the numerical hue can be joined
    back *only* after selection. This is not a sample of global wild plants.
    """
    require(set(COLS).issubset(photo.columns),"Frozen source columns missing")
    d=photo.loc[:,COLS].copy()
    for c in ("photo_id","observation_id","inat_taxon_id","source_design"):
        d[c]=d[c].fillna("").astype(str).str.strip()
    require(d.photo_id.is_unique and d.photo_id.str.fullmatch(r"[0-9]+").all(),
            "Original photo IDs must be unique numeric strings")
    require(d.source_design.isin(ARMS).all(),"Unrecognized source selection arm")
    require(d.strict9.notna().all() and
            pd.api.types.is_bool_dtype(d.strict9.dtype),
            "Missing/invalid original strict9 flag")
    without_obs=int(d.observation_id.eq("").sum())
    available=d.loc[d.observation_id.ne("")].copy()
    require(len(available)>0,"No authentic observation UUIDs")
    available["_canonical_rank"]=[
        canonical_rank(o,p) for o,p in
        zip(available.observation_id,available.photo_id)
    ]
    available=available.sort_values(
        ["observation_id","_canonical_rank","photo_id"],kind="mergesort"
    )
    chosen=available.drop_duplicates("observation_id",keep="first").copy()
    require(chosen.observation_id.is_unique,"Original observation dedup failed")
    return chosen,{
        "source_original_photo_attempts":int(len(d)),
        "photos_without_observation_uuid":without_obs,
        "source_photos_with_uuid":int(len(available)),
        "distinct_represented_observations":int(len(chosen)),
        "additional_photos_of_already_represented_observations":
            int(len(available)-len(chosen)),
        "canonical_success_count":int(chosen.strict9.sum()),
        "canonical_unknown_count":int((~chosen.strict9).sum()),
        "canonical_photo_chosen_outcome_blind":True,
    }
