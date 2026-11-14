#!/usr/bin/env python3
"""Find GBIF CC field-observation photos for MOIS species with no original ROI9.

Unlike specimen-first GBIF searches, explicitly ask the official GBIF
occurrence API for HUMAN_OBSERVATION + StillImage at an exact COL XR
accepted species. Validate returned basisOfRecord, species name, license
and non-iNaturalist media URL. Retain source uncertainty; NO image bytes
downloaded and NO hue measured in this source step.
"""
from __future__ import annotations
import argparse,hashlib,json,time
from pathlib import Path
from urllib.parse import urlparse

import pandas as pd

from scripts.acquisition.mois_zero9_gbif_other_cc_media_pilot import (
    BASELINE, COL_XR, request_json, safe_image, source_sha, require
)
from scripts.acquisition.mois_zero9_gbif_other_cc_media_union import SCHEMA as UNION_SCHEMA

SCHEMA="mois_noninat_GBIF_COLXR_HUMAN_OBSERVATION_CC_289taxa_source_v1"
PALETTE_SOURCE="mois_CC_noninat_GBIF_304taxa_exact_name_status.csv"
FIELD_SOURCE="https://api.gbif.org/v1/occurrence/search"


def signed_accepted(union_root,group):
    root=Path(union_root)
    r=json.loads((root/"receipt.json").read_text())
    require(r.get("schema")==UNION_SCHEMA
            and r.get("source_original_42111_taxon_photo_atlas_sha256")==BASELINE
            and r.get("stratified_hash_draw_taxa")==304
            and r.get("official_COL_XR_exact_accepted_taxa_count")==289
            and r.get("photo_pixels_downloaded")==0,
            "Original signed 304-taxon zero9 GBIF source unavailable")
    p=root/PALETTE_SOURCE
    require(source_sha(p)==r["outputs_sha256"][p.name],
            "Original accepted COL XR source taxonomy roster SHA drift")
    d=pd.read_csv(p,dtype=str,keep_default_na=False)
    require(len(d)==304 and d.inat_taxon_id.is_unique,"Full 304 source-taxonomic identity list missing")
    accepted=d.loc[d.name_status.eq("accepted_exact") &
                   d.accepted_colxr_key.ne("")].copy()
    require(len(accepted)==289 and accepted.accepted_colxr_key.is_unique,
            "289 exact COL XR accepted species source identities must be unique")
    require(group in range(4),"Only four disjoint index strata accepted")
    chosen=accepted.loc[accepted.inat_taxon_id.map(
        lambda x:int(hashlib.sha256(
            ("MOIS-GBIF-FIELD-VS-HERBARIUM-289|"+x).encode()).hexdigest(),16)%4==group)]
    return chosen.sort_values("inat_taxon_id",kind="stable")


def field_media(hit,*,inat_taxon_id,binomial,accepted_key):
    """One actual HUMAN_OBSERVATION record; never infer flower from its basis."""
    species=" ".join(str(hit.get("species") or "").split()[:2]).casefold()
    if species!=str(binomial).casefold() or str(hit.get("basisOfRecord"))!="HUMAN_OBSERVATION":
        return [],1
    occurrence=str(hit.get("key") or "")
    dataset=str(hit.get("datasetKey") or "")
    if not occurrence or not dataset:
        return [],1
    records=[];review=0
    for photo in hit.get("media",[]) or []:
        if not isinstance(photo,dict):continue
        if str(photo.get("type") or photo.get("mediaType") or "").lower() not in (
            "stillimage","image"):continue
        url=str(photo.get("identifier") or "")
        license_=str(photo.get("license") or "")
        if not safe_image(url,license_):
            review+=1;continue
        records.append({
            "inat_taxon_id":str(inat_taxon_id),
            "source_binomial":str(binomial),
            "gbif_COL_XR_accepted_species_key":str(accepted_key),
            "gbif_occurrence_key":occurrence,
            "gbif_dataset_key":dataset,
            "image_url":url,
            "image_license":license_,
            "media_host":(urlparse(url).hostname or "").lower(),
            "gbif_basisOfRecord":"HUMAN_OBSERVATION",
            "GBIF_media_attached_to_correct_species":True,
            "floral_organ_visible_verified":False,
            "source_photo_URL_sha256":hashlib.sha256(url.encode()).hexdigest(),
            "pixel_bytes_downloaded":False,
            "numeric9_flower_ROI_measured":False,
            "field_observation_local_hue_verified":False,
        })
    return records,review


def run(union,out,*,group,pace=.3,limit=60,max_pages=2):
    out=Path(out)
    require(not out.exists() and .2<=pace<=5 and 10<=limit<=100 and
            max_pages in (1,2),"Bound taxon-level public GBIF field photo API query")
    taxa=signed_accepted(union,group)
    require(len(taxa)>0,"This exact COL XR source stratum has no taxa")
    found=[];statuses=[]
    for n,row in enumerate(taxa.itertuples(index=False)):
        info={"inat_taxon_id":str(row.inat_taxon_id),
              "source_binomial":str(row.source_binomial),
              "accepted_colxr_key":str(row.accepted_colxr_key),
              "source_stratum":group,
              "HUMAN_OBSERVATION_qualifying_media":0,
              "search_status":"unqueried"}
        acquired=[];rejected=0
        count=0
        try:
            for page in range(max_pages):
                results=request_json(FIELD_SOURCE,{
                    "checklistKey":COL_XR,
                    "taxonKey":row.accepted_colxr_key,
                    "basisOfRecord":"HUMAN_OBSERVATION",
                    "mediaType":"StillImage",
                    "occurrenceStatus":"PRESENT",
                    "limit":limit,"offset":limit*page,
                })
                count=int(results.get("count",0))
                hits=results.get("results",[])
                for hit in hits:
                    media,bad=field_media(hit,inat_taxon_id=row.inat_taxon_id,
                                         binomial=row.source_binomial,
                                         accepted_key=row.accepted_colxr_key)
                    acquired.extend(media);rejected+=bad
                if acquired or results.get("endOfRecords") or len(hits)<limit:
                    break
                time.sleep(pace)
            by_url={}
            for v in acquired:
                by_url.setdefault(v["image_url"],v)
            acquired=list(by_url.values())
            info.update(search_status="retrieved",
                        field_search_GBIF_occurrences_total=count,
                        source_pagination_not_exhaustive=count>(page+1)*limit,
                        HUMAN_OBSERVATION_qualifying_media=len(acquired),
                        source_media_or_taxon_rejected=rejected)
            found.extend(acquired)
        except Exception as exc:
            info.update(search_status="api_error_review",error=str(exc)[:180],
                        HUMAN_OBSERVATION_qualifying_media=0)
        statuses.append(info)
        if (n+1)%20==0:
            print("MOIS_GBIF_COLXR_FIELD_ONLY_SOURCE_PROGRESS",group,n+1,len(taxa),flush=True)
        time.sleep(pace)
    out.mkdir(parents=True)
    main=out/f"MOIS_CCGBIF_HUMAN_OBSERVATION_source_photo_group{group}.csv.gz"
    status=out/f"MOIS_CCGBIF_HUMAN_OBSERVATION_search_status_group{group}.csv"
    cols=["inat_taxon_id","source_binomial","gbif_COL_XR_accepted_species_key",
          "gbif_occurrence_key","gbif_dataset_key","image_url","image_license",
          "media_host","gbif_basisOfRecord","GBIF_media_attached_to_correct_species",
          "floral_organ_visible_verified","source_photo_URL_sha256",
          "pixel_bytes_downloaded","numeric9_flower_ROI_measured",
          "field_observation_local_hue_verified"]
    pd.DataFrame(found,columns=cols).to_csv(main,index=False,compression="gzip")
    pd.DataFrame(statuses).to_csv(status,index=False)
    report={
       "schema":SCHEMA,
       "source_old_iNaturalist_photo_atlas_SHA256":BASELINE,
       "source_304_taxa_union_receipt_SHA256":source_sha(Path(union)/"receipt.json"),
       "stratum":group,"source_accepted_species_queried":len(taxa),
       "successful_GBIF_field_API_taxa":sum(x["search_status"]=="retrieved" for x in statuses),
       "api_error_unknown_taxa":sum(x["search_status"]=="api_error_review" for x in statuses),
       "source_taxa_with_human_observation_CC_still_media":len({x["inat_taxon_id"] for x in found}),
       "licensed_noninat_field_still_image_URL_metadata":len(found),
       "image_pixels_downloaded":0,
       "numeric9_ROI_measured":0,
       "human_observation_basis_is_not_flower_visibility_proof":True,
       "GBIF_occurrence_local_hue_verified":False,
       "outputs_sha256":{x.name:source_sha(x) for x in (main,status)}
    }
    (out/"receipt.json").write_text(json.dumps(report,sort_keys=True,indent=2)+"\n")
    print("MOIS_REAL_GBIF_COLXR_FIELD_ONLY_CC_IMAGE_SOURCE",
          json.dumps(report,sort_keys=True),flush=True)
    return report


def combine(roots,out):
    out=Path(out)
    require(not out.exists() and len(roots)==4,"Four authenticated taxon strata mandatory")
    sources=[];statuses=[];receipt_sha={}
    for i,root in enumerate(roots):
        root=Path(root)
        m=json.loads((root/"receipt.json").read_text())
        img=root/f"MOIS_CCGBIF_HUMAN_OBSERVATION_source_photo_group{i}.csv.gz"
        stat=root/f"MOIS_CCGBIF_HUMAN_OBSERVATION_search_status_group{i}.csv"
        require(m.get("schema")==SCHEMA and m.get("stratum")==i
                and source_sha(img)==m["outputs_sha256"][img.name]
                and source_sha(stat)==m["outputs_sha256"][stat.name],
                "Missing or modified GBIF field photo source signed group")
        sources.append(pd.read_csv(img,compression="gzip",dtype=str,keep_default_na=False))
        statuses.append(pd.read_csv(stat,dtype=str,keep_default_na=False))
        receipt_sha[str(i)]=source_sha(root/"receipt.json")
    taxa=pd.concat(statuses,ignore_index=True)
    imgs=pd.concat(sources,ignore_index=True)
    require(len(taxa)==289 and taxa.inat_taxon_id.is_unique,
            "289 accepted source photo taxa missing/overlap across groups")
    dup=imgs.loc[imgs.duplicated("image_url",keep=False)].copy()
    unique=imgs.loc[~imgs.duplicated("image_url",keep=False)].copy()
    require(unique.image_url.is_unique and
            unique.inat_taxon_id.isin(taxa.inat_taxon_id).all() and
            unique.gbif_basisOfRecord.eq("HUMAN_OBSERVATION").all(),
            "GBIF field still original-photo source contains repeated URL or specimen")
    hosts=unique.groupby("media_host").agg(
        URL_count=("image_url","nunique"),
        distinct_taxa=("inat_taxon_id","nunique")).sort_values(
            "distinct_taxa",ascending=False).reset_index()
    out.mkdir()
    files={
       "field":out/"mois_CCGBIF_COLXR_289taxa_human_observation_original_CC_photo_URLs.csv.gz",
       "taxa":out/"mois_CCGBIF_COLXR_289taxa_human_observation_source_status.csv",
       "hosts":out/"mois_CCGBIF_field_CC_media_host_counts.csv",
       "quarantined":out/"mois_CCGBIF_field_duplicate_URL_review.csv",
    }
    unique.to_csv(files["field"],index=False,compression="gzip")
    taxa.to_csv(files["taxa"],index=False)
    hosts.to_csv(files["hosts"],index=False)
    dup.to_csv(files["quarantined"],index=False)
    report={
       "schema":"mois_GBIF_COLXR_289taxa_field_original_CC_media_union_v1",
       "source_accepted_species_without_iNat_numeric_ROI":len(taxa),
       "source_taxa_with_actual_human_observation_CC_still_image":int(
           unique.inat_taxon_id.nunique()),
       "noninat_public_CC_human_observation_media_URLs_unique":len(unique),
       "duplicate_URL_rows_held_for_review":len(dup),
       "noninat_human_observation_photo_hosts":len(hosts),
       "taxa_with_metadata_API_errors":int(taxa.search_status.eq("api_error_review").sum()),
       "image_content_verified_flower_only":False,
       "new_image_pixels_downloaded":0,
       "new_nine_fraction_ROI9_measurements":0,
       "GBIF_spatially_observed_flower_hues":False,
       "signed_group_receipts_sha256":receipt_sha,
       "source_external_image_sha_duplicate_vs_old_iNaturalist_unverified":True,
       "outputs_sha256":{k:source_sha(v) for k,v in files.items()},
    }
    (out/"receipt.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print("MOIS_REAL_GBIF_289TAXA_FIELD_CC_IMAGE_SOURCE_UNION",
          json.dumps(report,sort_keys=True),flush=True)
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    sub=p.add_subparsers(dest="mode",required=True)
    a=sub.add_parser("group")
    a.add_argument("--union",required=True)
    a.add_argument("--group",required=True,type=int)
    a.add_argument("--pace",type=float,default=.3)
    a.add_argument("--output",required=True)
    b=sub.add_parser("combine")
    b.add_argument("--root",required=True,action="append")
    b.add_argument("--output",required=True)
    args=p.parse_args()
    if args.mode=="group":
        run(args.union,args.output,group=args.group,pace=args.pace)
    else:
        combine(args.root,args.output)


if __name__=="__main__":
    main()
