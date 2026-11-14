#!/usr/bin/env python3
"""Merge independently acquired GBIF COL XR CC-media metadata batches for MOIS.

Only a source-opportunity index: no flowers/petals are verified, no pixels
retrieved, no colour inferred. Review host and image URL reuse before acquisition.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
from urllib.parse import urlparse
import pandas as pd

from scripts.acquisition.mois_zero9_gbif_other_cc_media_pilot import (
    SCHEMA as PILOT_SCHEMA, BASELINE, safe_image, source_sha, require
)

SCHEMA="mois_zero9_304taxa_noninat_CC_GBIF_source_union_v1"
CANDIDATE="mois_zero9_noninat_CC_GBIF_image_metadata_candidates.csv.gz"
STATUS="mois_zero9_noninat_CC_GBIF_name_and_media_status.csv"


def union(roots,out):
    out=Path(out)
    require(not out.exists(),"Do not overwrite earlier provenance-certified CC media roster")
    require(len(roots)==5,"Require prior 48 plus 4 completely separate 64-taxon batches")
    allphotos=[];allnames=[];sha={}
    offsets=[0,48,112,176,240]
    sizes=[48,64,64,64,64]
    for ix,root in enumerate(roots):
        root=Path(root)
        rpath=root/"receipt.json"
        receipt=json.loads(rpath.read_text())
        require(receipt.get("schema")==PILOT_SCHEMA
                and receipt.get("source_signed_original_atlas_sha256")==BASELINE
                and receipt.get("original_photo_pixels_downloaded")==0
                and receipt.get("nine_fraction_flower_colour_measurements")==0
                and receipt.get("selected_zero9_binomial_taxa_sample")==sizes[ix]
                and receipt.get("sample_source_fixed_hash_offset",0)==offsets[ix],
                "A GBIF source sampling batch changed taxon counts / frozen source identity")
        pp=root/CANDIDATE;ss=root/STATUS
        for file in (pp,ss):
            require(file.is_file() and
                    source_sha(file)==receipt["outputs_sha256"][file.name],
                    "GBIF CC-media source results not byte-SHA authenticated")
        photos=pd.read_csv(pp,compression="gzip",dtype=str,keep_default_na=False)
        taxa=pd.read_csv(ss,dtype=str,keep_default_na=False)
        require(len(taxa)==sizes[ix] and taxa.inat_taxon_id.is_unique,
                "Sample source taxa were repeated or omitted")
        require(len(photos)==int(receipt["noninat_CC_candidate_media_metadata"]),
                "GBIF original media candidates changed since signed terminal")
        if not photos.empty:
            require(photos.inat_taxon_id.isin(taxa.inat_taxon_id).all() and
                    photos.gbif_occurrence_key.ne("").all(),
                    "GBIF media may be assigned only to source-audited taxa/occurrences")
            require(photos.apply(
                lambda row:safe_image(row.image_url,row.image_license),axis=1).all(),
                "External image URL or CC license changed after prior screening")
        sha[str(offsets[ix])]={"receipt":source_sha(rpath),
                               "media":source_sha(pp),"taxa":source_sha(ss)}
        taxa["pilot_source_offset"]=offsets[ix]
        allphotos.append(photos)
        allnames.append(taxa)
    taxa=pd.concat(allnames,ignore_index=True)
    photos=pd.concat(allphotos,ignore_index=True)
    require(len(taxa)==304 and taxa.inat_taxon_id.is_unique,
            "Pilot input sample taxon hash windows overlapped")
    distinct=photos.drop_duplicates(subset=["image_url"],keep=False).copy()
    reused=photos.loc[photos.duplicated("image_url",keep=False)].copy()
    # Same external image may be present against several GBIF occurrences or
    # species. Do not let that URL be downloaded for any taxon automatically.
    if not distinct.empty:
        distinct["media_host"]=distinct.image_url.map(
            lambda x:(urlparse(x).hostname or "").lower())
        require(distinct.image_url.is_unique and distinct.media_host.ne("").all(),
                "Unsafe URL or duplicate external source photo retained")
        require(distinct.inat_taxon_id.isin(taxa.inat_taxon_id).all(),
                "Exact GBIF taxon association missing from immutable trial roster")
    hosts=(distinct.groupby("media_host",dropna=False)
          .agg(unique_candidate_urls=("image_url","nunique"),
               source_taxa=("inat_taxon_id","nunique"))
          .sort_values("unique_candidate_urls",ascending=False).reset_index()
          if not distinct.empty else
          pd.DataFrame(columns=["media_host","unique_candidate_urls","source_taxa"]))
    out.mkdir(parents=True)
    sources=out/"mois_CC_noninat_GBIF_unique_image_metadata_304taxa.csv.gz"
    allstatus=out/"mois_CC_noninat_GBIF_304taxa_exact_name_status.csv"
    dup=out/"mois_CC_noninat_GBIF_duplicate_image_URLs_for_manual_review.csv"
    perhost=out/"mois_CC_noninat_GBIF_public_image_host_inventory.csv"
    distinct.to_csv(sources,index=False,compression="gzip")
    taxa.to_csv(allstatus,index=False)
    reused.to_csv(dup,index=False)
    hosts.to_csv(perhost,index=False)
    receipt={
        "schema":SCHEMA,
        "source_original_42111_taxon_photo_atlas_sha256":BASELINE,
        "source_zero9_nominal_taxa_total":2229,
        "stratified_hash_draw_taxa":len(taxa),
        "stratified_source_hash_offsets":offsets,
        "official_COL_XR_exact_accepted_taxa_count":int(taxa.name_status.eq("accepted_exact").sum()),
        "taxa_with_any_CC_noninat_media_metadata":int(
            photos.inat_taxon_id.nunique()) if not photos.empty else 0,
        "distinct_taxa_with_unique_external_photo_url":int(
            distinct.inat_taxon_id.nunique()) if not distinct.empty else 0,
        "raw_candidate_url_rows":len(photos),
        "eligible_unique_external_URL_metadata_count":len(distinct),
        "duplicate_URL_rows_held_for_review":len(reused),
        "image_hosts_seen":len(hosts),
        "source_artifacts_sha256":sha,
        "photo_pixels_downloaded":0,
        "floral_organ_visible_independently_verified":False,
        "external_media_raw_image_sha_overlap_with_prior_1231244_not_verified":True,
        "new_strict9_numerical_photo_measurements":0,
        "new_colour_sampled_taxa":0,
        "GBIF_occurrence_geographic_colour_measured":False,
        "outputs_sha256":{p.name:source_sha(p) for p in
                          (sources,allstatus,dup,perhost)}
    }
    (out/"receipt.json").write_text(json.dumps(receipt,sort_keys=True,indent=2)+"\n")
    print("MOIS_ZERO9_GBIF_COLXR_304TAXA_CC_MEDIA_UNION",json.dumps({
       k:receipt[k] for k in ("stratified_hash_draw_taxa",
          "official_COL_XR_exact_accepted_taxa_count",
          "distinct_taxa_with_unique_external_photo_url",
          "eligible_unique_external_URL_metadata_count",
          "duplicate_URL_rows_held_for_review","image_hosts_seen")},sort_keys=True),flush=True)
    return receipt


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root",action="append",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    union(a.root,a.output)
