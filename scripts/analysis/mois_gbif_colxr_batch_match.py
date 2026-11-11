#!/usr/bin/env python3
"""GBIF v2 COL XR batch name matching from an authenticated MOIS 9-hue atlas.

Official GBIF matching-ws POST /v2/species/match accepts <=1000 NameUsageQuery
objects and preserves input order. Explicit calibration against 138 prior
GET name matches precedes all batch work. Conservative exact/accepted/species
requires botanical kingdom, canonical agreement; no fuzzy/synonym promotion.
"""
from __future__ import annotations
import argparse
import gzip
import hashlib
import json
import re
import time
from pathlib import Path
from urllib.error import HTTPError,URLError
from urllib.request import Request,urlopen

import pandas as pd
from scripts.analysis.mois_gbif_colxr_nine_palette_name_pilot import (
    COL_XR, API, PALETTE, NAME_RULE, inspect_response, require
)

SCHEMA="mois_gbif_colxr_v2_batch_exact_species_20261111_v1"
CALIBRATION_RUN=42931715566


def source_sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda:f.read(1<<20),b""):h.update(block)
    return h.hexdigest()


def match_batch(names,*,tries=3):
    require(0<len(names)<=1000,"GBIF batch request exceeds official documented maximum")
    payload=json.dumps([{"scientificName":x} for x in names],
                       ensure_ascii=False,separators=(",",":")).encode("utf-8")
    # Checklist defaults to COL XR on current GBIF v2 match service.
    # The explicit 138-species calibration rejects any default dataset drift.
    req=Request(API+"?checklistKey="+COL_XR,data=payload,method="POST",headers={
        "Accept":"application/json","Content-Type":"application/json",
        "User-Agent":"MOIS-floral-colour-global-species-name-research/1.0"})
    last=None
    for n in range(tries):
        try:
            with urlopen(req,timeout=90) as f:
                result=json.loads(f.read(35000000))
            require(isinstance(result,list) and len(result)==len(names),
                    "GBIF POST matching batch response must be an ordered complete array")
            return result
        except HTTPError as e:
            last=e
            if e.code not in (429,500,502,503,504):break
        except (URLError,TimeoutError,ValueError) as e:
            last=e
        if n+1<tries:time.sleep(2**n)
    raise RuntimeError("GBIF POST v2 batch refused/error: "+str(last))


def parse_row(row, response):
    info=inspect_response(response,row.species)
    classification=response.get("classification") or []
    kingdoms=[c.get("name") for c in classification
              if isinstance(c,dict) and c.get("rank")=="KINGDOM"]
    if info["status"]=="colxr_exact_accepted_species" and "Plantae" not in kingdoms:
        info["status"]="kingdom_context_review"
    return {"inat_taxon_id":row.inat_taxon_id,
            "species":row.species,
            "source_hue9":row.dominant_hue9,
            "hue_status":row.hue_status,
            "COL_XR_checklistKey":COL_XR,
            "GBIF_v2_batch_query":"NameUsageQuery.scientificName",
            "GBIF_kingdom":"|".join(str(x) for x in kingdoms),
            **info}


def verify_calibration(source,*,minimum=130):
    known=pd.read_csv(source,dtype=str,keep_default_na=False)
    good=known.loc[known.status.eq("colxr_exact_accepted_species")].copy()
    require(len(good)==138 and good.key.is_unique and good.inat_taxon_id.is_unique,
            "Pinned 138 known COL XR key table is invalid")
    inputs=good.source_species.tolist()
    response=match_batch(inputs)
    matches=[]
    mis=[]
    for z,r in zip(good.itertuples(index=False),response):
        value=inspect_response(r,z.source_species)
        accepted=(value["status"]=="colxr_exact_accepted_species"
                  and value["key"]==z.key)
        matches.append(accepted)
        if not accepted:mis.append(dict(name=z.source_species,key_get=z.key,
            key_post=value["key"],status=value["status"]))
    # No allowance for mis-identification: a changed COL release is a halt.
    require(sum(matches)>=minimum and not mis,
            "GBIF batch COL XR taxonomy does not reproduce 138 independent known GET accepted species keys: "+repr(mis[:5]))
    return {"schema":"mois_gbif_colxr_batch_calibration_v1",
            "known_get_species":len(good),
            "post_identical_exact_species_keys":sum(matches),
            "known_get_sha256":source_sha(source),
            "dataset":COL_XR,
            "mismatch_count":len(mis)}


def select_nominal_species(source):
    catalogue=pd.read_csv(source,compression="infer",dtype=str,keep_default_na=False)
    require({"inat_taxon_id","species","dominant_hue9","hue_status"}<=set(catalogue.columns)
            and catalogue.inat_taxon_id.is_unique,"Original full MOIS species catalogue invalid")
    requested=catalogue.loc[catalogue.dominant_hue9.isin(PALETTE) &
                             catalogue.species.str.match(NAME_RULE)].copy()
    requested=requested.sort_values("inat_taxon_id",kind="stable").reset_index(drop=True)
    require(len(requested)>=20000,"Expected source-wide species with nine-colour labels absent")
    return requested


def batch_atlas(source,calibration,output,*,chunk_size=250,pause=1.0):
    output=Path(output)
    require(not output.exists() and 1<=chunk_size<=250 and 1<=pause<=10,
            "Unsafe GBIF batch size/pacing or overwrite of prior run")
    require(Path(calibration).is_file(),"Prior official GBIF v2 calibration source required")
    rows=select_nominal_species(source)
    output.mkdir(parents=True)
    shards=output/"shards";shards.mkdir()
    n_chunks=(len(rows)+chunk_size-1)//chunk_size
    counts={}
    for i in range(n_chunks):
        slice_=rows.iloc[i*chunk_size:(i+1)*chunk_size]
        names=slice_.species.tolist()
        try:
            results=match_batch(names)
            report=[parse_row(row,r) for row,r in
                    zip(slice_.itertuples(index=False),results)]
        except Exception as e:
            # Do not treat unavailable API matches as biological absence.
            report=[{"inat_taxon_id":r.inat_taxon_id,"species":r.species,
                     "source_hue9":r.dominant_hue9,"hue_status":r.hue_status,
                     "COL_XR_checklistKey":COL_XR,
                     "status":"api_error_review","key":"",
                     "error_note":str(e)[:300]}
                    for r in slice_.itertuples(index=False)]
        path=shards/f"gbif_colxr_batch_{i:03d}.csv.gz"
        pd.DataFrame(report).to_csv(path,index=False,compression="gzip")
        counts[path.name]=source_sha(path)
        print("MOIS_GBIF_COLXR_MATCH_BATCH",i+1,n_chunks,
              "rows",len(report),"accepted",sum(r["status"]=="colxr_exact_accepted_species" for r in report),
              flush=True)
        time.sleep(pause)
    frames=[pd.read_csv(shards/f"gbif_colxr_batch_{i:03d}.csv.gz",dtype=str,keep_default_na=False)
            for i in range(n_chunks)]
    full=pd.concat(frames,ignore_index=True)
    require(len(full)==len(rows) and full.inat_taxon_id.is_unique and
            set(full.inat_taxon_id)==set(rows.inat_taxon_id),
            "GBIF batch matched catalogue lost or repeated source plant taxa")
    good=full.loc[full.status.eq("colxr_exact_accepted_species")].copy()
    collision_keys=set(good.loc[good.key.duplicated(keep=False),"key"])
    if collision_keys:
        full.loc[full.key.isin(collision_keys),"status"]="duplicate_accepted_species_key_review"
    # Unresolved source taxa remain in full output with explicit status.
    full.to_csv(output/"mois_28312_named_hue_taxa_gbif_colxr_match.csv.gz",
                index=False,compression="gzip")
    allowed=full.loc[full.status.eq("colxr_exact_accepted_species")].copy()
    allowed.to_csv(output/"mois_gbif_colxr_exact_accepted_species_hue_roster.csv.gz",
                   index=False,compression="gzip")
    receipt={"schema":SCHEMA,
        "verified_138_colxr_get_calibration_reference_sha256":source_sha(calibration),
        "source_nine_colour_catalogue_sha256":source_sha(source),
        "actual_requested_nominal_taxa":len(rows),
        "successful_exact_accepted_colxr_species":len(allowed),
        "unresolved_or_conflicting_taxa":len(rows)-len(allowed),
        "palette_count_in_source":int(rows.dominant_hue9.nunique()),
        "match_status_counts":{str(k):int(v) for k,v in full.status.value_counts().items()},
        "shard_sha256":counts,
        "taxon_key_scope":"COL_XR_not_legacy_integer_GBIF",
        "GBIF_occurrence_records_downloaded":0,
        "GBIF_download_DOI_created":False,
        "species_range_photo_hue_fixed_proven":False}
    (output/"receipt.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print("MOIS_GBIF_COLXR_FULL9_PALETTE_NAME_MATCH",json.dumps({
        k:receipt[k] for k in ("actual_requested_nominal_taxa",
                              "successful_exact_accepted_colxr_species",
                              "unresolved_or_conflicting_taxa","match_status_counts")},sort_keys=True),flush=True)
    return receipt


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--mode",choices=["calibrate","full"],required=True)
    p.add_argument("--source",required=True);p.add_argument("--output",required=True)
    p.add_argument("--calibration",default=None)
    p.add_argument("--chunk-size",type=int,default=250)
    p.add_argument("--pause",type=float,default=1.0)
    args=p.parse_args()
    if args.mode=="calibrate":
        r=verify_calibration(args.source)
        path=Path(args.output)
        require(not path.exists(),"Refuse overwrite")
        path.parent.mkdir(parents=True,exist_ok=True)
        path.write_text(json.dumps(r,indent=2,sort_keys=True)+"\n")
        print("MOIS_GBIF_COLXR_POST_V2_BATCH_CALIBRATED",
              json.dumps(r,sort_keys=True),flush=True)
    else:
        require(args.calibration and Path(args.calibration).is_file(),
                "Full match batch requires successful independent GET-key calibration")
        cal=json.loads(Path(args.calibration).read_text())
        require(cal.get("known_get_species")==138 and
                cal.get("post_identical_exact_species_keys")==138 and
                cal.get("mismatch_count")==0,
                "GBIF POST does not preserve known COL XR accepted taxon-key identities")
        batch_atlas(args.source,args.calibration,args.output,
                    chunk_size=args.chunk_size,pause=args.pause)


if __name__=="__main__":
    main()
