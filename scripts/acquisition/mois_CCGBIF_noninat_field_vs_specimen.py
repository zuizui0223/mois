#!/usr/bin/env python3
"""Classify GBIF original occurrence-media identity before external flower ROI9.

This audits basisOfRecord and exact URL attribution on already-source-signed
non-iNaturalist CC photos; does NOT retrieve image bytes or infer petal colours.
An herbarium specimen image is not treated as a field flower observation.
"""
from __future__ import annotations
import argparse,hashlib,json,time
from pathlib import Path
from urllib.request import Request,urlopen
from urllib.error import HTTPError,URLError
from urllib.parse import quote
import pandas as pd
from scripts.acquisition.mois_zero9_gbif_other_cc_media_pilot import (
    BASELINE, safe_image, source_sha, require
)
from scripts.acquisition.mois_zero9_gbif_other_cc_media_union import SCHEMA as UNION

SCHEMA="mois_CCGBIF_noninat_occurrence_basis_two_photos_per_taxon_source_v1"
FIELD={"HUMAN_OBSERVATION","MACHINE_OBSERVATION","OBSERVATION","LIVING_SPECIMEN"}
SPECIMEN={"PRESERVED_SPECIMEN","FOSSIL_SPECIMEN","MATERIAL_SAMPLE","MATERIAL_CITATION"}


def signed_candidates(union_dir,group=0):
    union_dir=Path(union_dir)
    m=json.loads((union_dir/"receipt.json").read_text())
    require(m.get("schema")==UNION and
            m.get("source_original_42111_taxon_photo_atlas_sha256")==BASELINE and
            m.get("eligible_unique_external_URL_metadata_count")==3365 and
            m.get("distinct_taxa_with_unique_external_photo_url")==217 and
            m.get("photo_pixels_downloaded")==0,
            "Only the signed real 304-taxon GBIF COL XR media ledger is permitted")
    path=union_dir/"mois_CC_noninat_GBIF_unique_image_metadata_304taxa.csv.gz"
    require(source_sha(path)==m["outputs_sha256"][path.name],
            "Original GBIF CC image media source SHA drifted")
    a=pd.read_csv(path,compression="gzip",dtype=str,keep_default_na=False)
    require(a.inat_taxon_id.nunique()==217 and a.image_url.is_unique and
            a.apply(lambda row:safe_image(row.image_url,row.image_license),axis=1).all(),
            "GBIF non-iNaturalist provider or CC media source invalid")
    require(group in range(4),"Four disjoint source strata mandatory")
    a=a.loc[a.inat_taxon_id.map(lambda x:
        int(hashlib.sha256(("MOIS-GBIF-CC-BASIS-SHARD|"+x).encode()).hexdigest(),16)%4==group)]
    a=a.copy()
    a["frozen_draw"]=a.apply(lambda r:hashlib.sha256((
        "MOIS-GBIF-CC-SPECIES-OBS-SAMPLE-V1|"+r.inat_taxon_id+
        "|"+r.gbif_occurrence_key+"|"+r.image_url).encode()).hexdigest(),axis=1)
    a=a.sort_values("frozen_draw",kind="stable")
    a=a.drop_duplicates(subset=["inat_taxon_id","gbif_occurrence_key"],keep="first")
    a=a.groupby("inat_taxon_id",sort=False).head(2)
    require(a.gbif_occurrence_key.is_unique,
            "GBIF observation or media URL conflicts across source taxa")
    return a.reset_index(drop=True)


def get_occurrence(gbif_id):
    require(str(gbif_id).isdigit(),"Require numeric GBIF source occurrence ID")
    url="https://api.gbif.org/v1/occurrence/"+quote(str(gbif_id),safe="")
    req=Request(url,headers={
      "User-Agent":"MOIS-CCGBIF-field-vs-specimen-flower-photo-provenance/1.0",
      "Accept":"application/json"})
    last=None
    for attempt in range(3):
        try:
            with urlopen(req,timeout=24) as response:
                return json.loads(response.read(700000))
        except HTTPError as e:
            last=e
            if e.code not in (408,429,500,502,503,504):break
        except (URLError,TimeoutError,ValueError) as e:
            last=e
        if attempt<2:time.sleep((attempt+1)*2)
    raise RuntimeError("GBIF source occurrence metadata request failed: "+str(last)[:140])


def inspect_row(row,hit):
    species=" ".join(str(hit.get("species") or "").split()[:2]).casefold()
    expected=" ".join(str(row.source_binomial).split()[:2]).casefold()
    key=str(hit.get("key") or "")
    dataset=str(hit.get("datasetKey") or "")
    images=[str(m.get("identifier") or "") for m in hit.get("media",[]) or []
            if isinstance(m,dict)]
    correct=(key==str(row.gbif_occurrence_key) and
             dataset==str(row.gbif_dataset_key) and
             species==expected and row.image_url in images)
    basis=str(hit.get("basisOfRecord") or "UNKNOWN").upper()
    if not correct: category="occurrence_url_or_taxonomy_mismatch_review"
    elif basis in FIELD:category="potential_field_observation_not_flower_verified"
    elif basis in SPECIMEN:category="preserved_or_material_specimen_not_field_flower"
    else:category="unknown_basis_review"
    return {
       "inat_taxon_id":str(row.inat_taxon_id),
       "source_binomial":str(row.source_binomial),
       "gbif_occurrence_key":str(row.gbif_occurrence_key),
       "gbif_dataset_key":str(row.gbif_dataset_key),
       "gbif_accepted_species_key":str(row.gbif_COL_XR_accepted_species_key),
       "image_url":str(row.image_url),
       "image_license":str(row.image_license),
       "original_media_host":str(row.media_host),
       "gbif_basisOfRecord":basis,
       "gbif_eventDate":str(hit.get("eventDate") or ""),
       "gbif_occurrence_media_URL_present":row.image_url in images,
       "same_source_species_and_dataset_proven":correct,
       "field_vs_specimen_gate":category,
       "flower_visible_verified":False,
       "image_pixels_downloaded":False,
       "source_original_MOIS_atlas_merged":False,
    }


def run(union_dir,out,*,group=0,pace=0.25):
    out=Path(out)
    require(not out.exists() and .2<=pace<=5,"Do not overwrite or disable polite API pacing")
    subset=signed_candidates(union_dir,group=group)
    results=[]
    for i,row in enumerate(subset.itertuples(index=False)):
        try:
            hit=get_occurrence(row.gbif_occurrence_key)
            res=inspect_row(row,hit)
        except Exception as exc:
            res={"inat_taxon_id":str(row.inat_taxon_id),
                 "source_binomial":str(row.source_binomial),
                 "gbif_occurrence_key":str(row.gbif_occurrence_key),
                 "image_url":str(row.image_url),
                 "field_vs_specimen_gate":"api_error_no_field_assumption",
                 "reason":str(exc)[:150],
                 "same_source_species_and_dataset_proven":False,
                 "flower_visible_verified":False,
                 "image_pixels_downloaded":False}
        results.append(res)
        if (i+1)%25==0:
            print("MOIS_CCGBIF_FIELD_BASIS_GENUINE_METADATA",group,i+1,len(subset),flush=True)
        time.sleep(pace)
    out.mkdir(parents=True)
    p=out/f"mois_GBIF_CC_media_occurrence_basis_group{group}.csv"
    pd.DataFrame(results).to_csv(p,index=False)
    count=pd.Series([r["field_vs_specimen_gate"] for r in results]).value_counts()
    note={
        "schema":SCHEMA,
        "source_304_taxon_metadata_receipt_SHA256":source_sha(Path(union_dir)/"receipt.json"),
        "stratum_index":group,
        "GBIF_occurrence_metadata_requests_completed":len(results),
        "distinct_zero_numeric_taxa_sampled":int(subset.inat_taxon_id.nunique()),
        "GBIF_basis_record_type_counts":{str(k):int(v) for k,v in count.items()},
        "potential_field_photo_taxa":len(set(r["inat_taxon_id"] for r in results
            if r["field_vs_specimen_gate"]=="potential_field_observation_not_flower_verified")),
        "trusted_source_taxon_occurrence_and_media_URL_joins":sum(
            r.get("same_source_species_and_dataset_proven",False) for r in results),
        "raw_original_images_downloaded":0,
        "flower_pixels_verified":0,
        "new_nine_palette_ROI_measured":0,
        "source_iNaturalist_photo_atlas_unchanged":True,
        "result_SHA256":source_sha(p),
    }
    (out/"receipt.json").write_text(json.dumps(note,sort_keys=True,indent=2)+"\n")
    print("MOIS_CCGBIF_REAL_OCCURRENCE_BASIS_GATE",json.dumps(note,sort_keys=True),flush=True)
    return note


def combine(group_roots,out):
    require(len(group_roots)==4,"Full source strata must be audited separately")
    out=Path(out);require(not out.exists(),"Do not overwrite GBIF field-vs-herbarium source audit")
    frames=[];receipts={}
    for i,path in enumerate(group_roots):
        root=Path(path);m=json.loads((root/"receipt.json").read_text())
        p=root/f"mois_GBIF_CC_media_occurrence_basis_group{i}.csv"
        require(m["schema"]==SCHEMA and m["stratum_index"]==i
                and source_sha(p)==m["result_SHA256"],
                "Missing or altered GBIF basisOfRecord source audit")
        frames.append(pd.read_csv(p,dtype=str,keep_default_na=False))
        receipts[str(i)]=source_sha(root/"receipt.json")
    all_=pd.concat(frames,ignore_index=True)
    require(all_.gbif_occurrence_key.is_unique,"Duplicate real GBIF occurrence ID after source join")
    field=all_.loc[all_.field_vs_specimen_gate.eq(
        "potential_field_observation_not_flower_verified")].copy()
    require(field.inat_taxon_id.nunique()<=217,"Outside 217-sampled exact COL XR plant taxa")
    out.mkdir()
    f=out/"mois_GBIF_CC_potential_field_observation_photos.csv"
    full=out/"mois_GBIF_CC_full_basis_304taxon_metadata_audit.csv.gz"
    field.to_csv(f,index=False)
    all_.to_csv(full,index=False,compression="gzip")
    cats=all_.field_vs_specimen_gate.value_counts()
    result={
        "schema":"mois_noninat_CC_GBIF_field_photo_basis_audit_union_v1",
        "source_304_taxa_with_CC_images":217,
        "actual_GBIF_occurrence_metadata_calls":len(all_),
        "source_taxa_examined":int(all_.inat_taxon_id.nunique()),
        "field_image_candidates_not_botanically_verified":len(field),
        "field_photo_candidate_source_taxa":int(field.inat_taxon_id.nunique()),
        "preserved_specimen_record_count":int(cats.get("preserved_or_material_specimen_not_field_flower",0)),
        "other_or_bad_provenance_count":int(len(all_)-len(field)-int(cats.get("preserved_or_material_specimen_not_field_flower",0))),
        "GBIF_field_basis_is_not_proof_of_flower_visibility":True,
        "original_images_downloaded":0,
        "numeric9_ROI_measured":0,
        "source_48_to_304_original_atlas_merged":False,
        "outputs_SHA256":{p.name:source_sha(p) for p in (f,full)},
        "signed_strata_receipt_SHA256":receipts,
    }
    (out/"receipt.json").write_text(json.dumps(result,sort_keys=True,indent=2)+"\n")
    print("MOIS_CCGBIF_ALL304_FIELD_VS_HERBARIUM_REAL_SOURCE",
          json.dumps(result,sort_keys=True),flush=True)
    return result


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    s=p.add_subparsers(dest="mode",required=True)
    group=s.add_parser("group")
    group.add_argument("--union",required=True)
    group.add_argument("--output",required=True)
    group.add_argument("--group",type=int,required=True)
    group.add_argument("--pace",type=float,default=.25)
    all_=s.add_parser("combine")
    all_.add_argument("--root",action="append",required=True)
    all_.add_argument("--output",required=True)
    a=p.parse_args()
    if a.mode=="group":run(a.union,a.output,group=a.group,pace=a.pace)
    else:combine(a.root,a.output)
