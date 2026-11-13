#!/usr/bin/env python3
"""No-image-download GBIF media pilot for unmeasured MOIS flower-colour taxa.

For a frozen hashed sample of taxa that have NO original numeric ROI9 in MOIS,
use exact accepted GBIF COL XR names and look for CC-licensed still-image media
from NON-iNaturalist providers. This is ONLY metadata opportunity discovery:
no pixels, flower hue, petal truth, or GBIF presence-local hue is inferred.
"""
from __future__ import annotations
import argparse, hashlib, json, re, time
from pathlib import Path
from urllib.parse import urlencode,urlparse
from urllib.request import Request,urlopen
from urllib.error import HTTPError,URLError
import pandas as pd
from scripts.analysis.mois_gbif_colxr_nine_palette_name_pilot import COL_XR,NAME_RULE

SCHEMA="mois_zero9_exact_colxr_gbif_other_provider_cc_image_metadata_pilot_v1"
MATCH="https://api.gbif.org/v2/species/match"
OCC="https://api.gbif.org/v1/occurrence/search"
VALID_LICENSE=("creativecommons.org/licenses/by/","creativecommons.org/licenses/by-nc/",
               "creativecommons.org/publicdomain/zero/","creativecommons.org/publicdomain/mark/")
BASELINE="0b2f4ac9b7598c89f4dc84235ed010f3e3fe91a3ab318b8fd39b66fbc0c466e4"


def require(test,text):
    if not bool(test):raise ValueError("MOIS_GBIF_MEDIA_RESCUE: "+text)


def source_sha(path):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for p in iter(lambda:f.read(1<<20),b""):h.update(p)
    return h.hexdigest()


def source_taxa(atlas,sample=48,offset=0):
    require(1<=sample<=256 and 0<=offset<=2500,"Bound GBIF API pilot to at most 256 taxa and a valid fixed source offset")
    require(len(atlas)==42111 and atlas.inat_taxon_id.is_unique,"Unrecognized immutable MOIS original-source taxa census")
    for name in ("inat_taxon_id","species","dominant_hue9","hue_status"):
        require(name in atlas,"Frozen nominal taxon and flower status missing: "+name)
    a=atlas.copy()
    a.inat_taxon_id=a.inat_taxon_id.fillna("").astype(str)
    a.species=a.species.fillna("").astype(str).str.strip()
    a.dominant_hue9=a.dominant_hue9.fillna("").astype(str).str.strip()
    a.hue_status=a.hue_status.fillna("").astype(str).str.strip()
    a["after_strict9"]=pd.to_numeric(a.after_strict9,errors="raise")
    a=a.loc[a.after_strict9.eq(0)&a.dominant_hue9.eq("")&
            a.species.str.match(NAME_RULE)].copy()
    a["priority"]=a.inat_taxon_id.map(lambda x:hashlib.sha256(
        ("MOIS-ZERO-ROI-OTHER-CC-PROVIDER-20261114|"+x).encode()).hexdigest())
    return a.sort_values("priority",kind="stable").iloc[offset:offset+sample].reset_index(drop=True)


def request_json(url,params):
    target=url+"?"+urlencode(params)
    req=Request(target,headers={"Accept":"application/json",
                "User-Agent":"MOIS-ecology-GBIF-CC-image-metadata-pilot/1.0"})
    last=""
    for attempt in range(3):
        try:
            with urlopen(req,timeout=28) as conn:
                return json.loads(conn.read(3500000))
        except (HTTPError,URLError,TimeoutError,ValueError) as exc:
            last=str(exc)
            if isinstance(exc,HTTPError) and exc.code not in (429,500,502,503,504):
                break
            if attempt<2:time.sleep(1+attempt*2)
    raise RuntimeError("Official GBIF public API failed: "+last[:160])


def is_exact_accepted_name(data,source_name):
    usage=data.get("usage") if isinstance(data,dict) else None
    if not isinstance(usage,dict):return "", "no_match"
    diag=data.get("diagnostics") or {}
    match_type=(diag.get("matchType") if isinstance(diag,dict) else None) or data.get("matchType")
    canonical=str(usage.get("canonicalName") or "")
    rank=str(usage.get("rank") or "")
    status=str(usage.get("status") or "")
    taxon_key=str(usage.get("key") or "")
    source=" ".join(source_name.split()[:2]).casefold()
    exact=" ".join(canonical.split()[:2]).casefold()
    ok=(match_type=="EXACT" and rank=="SPECIES" and status=="ACCEPTED" and
        len(source_name.split())==2 and source==exact and bool(taxon_key))
    return (taxon_key, "accepted_exact") if ok else ("","not_exact_accepted_species")


def safe_image(url,license_text):
    link=str(url or "").strip()
    lic=str(license_text or "").lower().replace("https://","").replace("http://","")
    if not link.startswith(("https://","http://")):return False
    host=(urlparse(link).hostname or "").lower()
    if not host or host.endswith("inaturalist.org") or "inaturalist" in host:
        return False
    if "inaturalist-open-data" in host:return False
    return any(lic.startswith(x) for x in VALID_LICENSE)


def candidate_media(rows,key,source_name,inat_taxon_id):
    """Accept distinct externally licensed STILL images only; not actual flowers yet."""
    out=[];review=[]
    for hit in rows:
        # Occurrences can be subspecific: accepted usage must still identify
        # the same nominated accepted species in the selected COL XR checklist.
        species=str(hit.get("species") or hit.get("scientificName") or "")
        if species and " ".join(species.split()[:2]).casefold()!=source_name.casefold():
            review.append("occurrence_name_conflict");continue
        for media in hit.get("media",[]) or []:
            if not isinstance(media,dict):continue
            if str(media.get("type") or media.get("mediaType") or "").lower() not in (
                    "stillimage","image"):continue
            url=media.get("identifier") or media.get("references") or ""
            if not safe_image(url,media.get("license")):
                review.append("unsupported_source_or_media_license");continue
            out.append({
                "inat_taxon_id":inat_taxon_id,
                "source_binomial":source_name,
                "gbif_COL_XR_accepted_species_key":str(key),
                "gbif_occurrence_key":str(hit.get("key") or ""),
                "gbif_dataset_key":str(hit.get("datasetKey") or ""),
                "image_url":str(url),
                "image_license":str(media.get("license")),
                "image_creator":str(media.get("creator") or ""),
                "image_publisher":str(media.get("publisher") or ""),
                "image_type":str(media.get("type") or media.get("mediaType")),
                "flower_visible_verified":False,
                "colour_measured":False,
                "source_taxon_id_originally_iNaturalist":True,
                "gbif_occurrence_local_flower_hue_measured":False,
            })
    return out,review


def run(atlas_path,output,*,sample=48,offset=0,page_limit=30,pages=2,pace=0.25):
    atlas_path=Path(atlas_path);output=Path(output)
    require(not output.exists(),"Do not overwrite existing authenticated source-only pilot")
    require(source_sha(atlas_path)==BASELINE,"Actual after47 42,111-taxon original-photo atlas SHA drifted")
    require(5<=page_limit<=100 and 1<=pages<=4 and 0.2<=pace<=4,
            "Bound GBIF public API load and make short pilot reproducible")
    atlas=pd.read_csv(atlas_path,compression="infer",dtype={"inat_taxon_id":str})
    taxa=source_taxa(atlas,sample=sample,offset=offset)
    require(len(taxa)>=8,"Not enough exact binomial taxa lacking original numeric ROI")
    done=[];media=[]
    for i,row in enumerate(taxa.itertuples(index=False)):
        taxon={"inat_taxon_id":row.inat_taxon_id,"source_binomial":row.species,
               "sample_hash":row.priority,"new_original_photo_pixels_downloaded":0}
        try:
            match=request_json(MATCH,{"scientificName":row.species,"checklistKey":COL_XR})
            key,state=is_exact_accepted_name(match,row.species)
            taxon["name_status"]=state
            taxon["accepted_colxr_key"]=key
            taxon["eligible_media_count"]=0
            if state=="accepted_exact":
                rejects=[]
                occurrences=[]
                for page in range(pages):
                    res=request_json(OCC,{
                        "checklistKey":COL_XR,"taxonKey":key,
                        "mediaType":"StillImage","occurrenceStatus":"PRESENT",
                        "hasCoordinate":"true","limit":page_limit,
                        "offset":page*page_limit})
                    occurrences.extend(res.get("results",[]))
                    taxon["GBIF_query_total_with_stillimage"]=int(res.get("count",0))
                    if res.get("endOfRecords") or not res.get("results"):break
                    time.sleep(pace)
                options, rejects=candidate_media(occurrences,key,row.species,row.inat_taxon_id)
                dedup={p["image_url"]:p for p in options if p["gbif_occurrence_key"]}
                taxon["eligible_media_count"]=len(dedup)
                taxon["GBIF_returned_occurrences_checked"]=len(occurrences)
                taxon["rejected_media_or_taxonomy"]=len(rejects)
                taxon["source_incomplete_after_page_limit"]=(
                    int(taxon.get("GBIF_query_total_with_stillimage",0))>len(occurrences))
                media.extend(dedup.values())
        except Exception as e:
            taxon["name_status"]="request_error_review"
            taxon["error"]=str(e)[:220]
            taxon["eligible_media_count"]=0
        done.append(taxon)
        if (i+1)%12==0:
            print("MOIS_GBIF_OTHER_PROVIDER_SOURCE_CHECK",i+1,len(taxa),flush=True)
        time.sleep(pace)
    output.mkdir(parents=True)
    source_file=output/"mois_zero9_noninat_CC_GBIF_image_metadata_candidates.csv.gz"
    status_file=output/"mois_zero9_noninat_CC_GBIF_name_and_media_status.csv"
    pd.DataFrame(media,columns=[
        "inat_taxon_id","source_binomial","gbif_COL_XR_accepted_species_key",
        "gbif_occurrence_key","gbif_dataset_key","image_url","image_license",
        "image_creator","image_publisher","image_type","flower_visible_verified",
        "colour_measured","source_taxon_id_originally_iNaturalist",
        "gbif_occurrence_local_flower_hue_measured"]).to_csv(source_file,index=False,compression="gzip")
    pd.DataFrame(done).to_csv(status_file,index=False)
    report={
        "schema":SCHEMA,"source_signed_original_atlas_sha256":source_sha(atlas_path),
        "original_source_zero_numeric_taxa_total":int(atlas.after_strict9.eq(0).sum()),
        "selected_zero9_binomial_taxa_sample":len(taxa),
        "sample_source_fixed_hash_offset":offset,
        "official_colxr_exact_accepted_species_count":sum(d.get("name_status")=="accepted_exact" for d in done),
        "selected_taxa_with_noninat_CC_stillimages":sum(d.get("eligible_media_count",0)>0 for d in done),
        "noninat_CC_candidate_media_metadata":len(media),
        "photo_candidate_urls_not_verified_to_show_flowers":True,
        "frozen_93m_ID_overlap_not_audited_for_other_source":True,
        "original_photo_pixels_downloaded":0,
        "nine_fraction_flower_colour_measurements":0,
        "new_species_photo_hue_labels":0,
        "GBIF_occurrence_local_hue_measured":False,
        "GBIF_global_SQL_downloads_submitted":0,
        "source_incomplete_when_pagination_capped":True,
        "outputs_sha256":{f.name:source_sha(f) for f in (source_file,status_file)}
    }
    (output/"receipt.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print("MOIS_GBIF_OTHER_PROVIDER_CC_PHOTO_SOURCE_PILOT",json.dumps(report,sort_keys=True),flush=True)
    return report


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--atlas",required=True)
    p.add_argument("--output",required=True)
    p.add_argument("--sample",type=int,default=48)
    p.add_argument("--offset",type=int,default=0)
    p.add_argument("--page-limit",type=int,default=30)
    p.add_argument("--pages",type=int,default=2)
    p.add_argument("--pace",type=float,default=0.25)
    a=p.parse_args()
    run(a.atlas,a.output,sample=a.sample,offset=a.offset,page_limit=a.page_limit,pages=a.pages,pace=a.pace)
