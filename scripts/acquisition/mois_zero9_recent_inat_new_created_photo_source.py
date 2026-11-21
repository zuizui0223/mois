#!/usr/bin/env python3
"""New iNaturalist source opportunity pilot AFTER frozen 93.2m source snapshot.

Query recent public observations created after 2026-10-27 for exact formerly
zero-original-ROI source taxon IDs, with independently licensed photo URLs.

IMPORTANT: This is not re-downloading any of the old 93m photos. We seek photo
IDs from a newer iNaturalist ingest period that were not accessible in the
previous frozen source. Before real image acquisition, compare all photo IDs
AND observation UUIDs with every one of the 1,231,244 old attempts, and reject
any overlap. No photo pixels or nine palette fractions are downloaded here.
"""
from __future__ import annotations
import argparse, hashlib, json, re, time
from datetime import datetime, timezone
from pathlib import Path
from urllib.error import HTTPError,URLError
from urllib.parse import urlencode,urlsplit
from urllib.request import Request,urlopen

import pandas as pd
from scripts.acquisition.mois_zero9_gbif_other_cc_media_pilot import (
    source_sha, require
)

SCHEMA="mois_new_inaturalist_created_after_20261027_zero9_photo_source_pilot_v1"
ATLAS_SHA="c504e2549f0659b9e0849a0178f689ee4caa70c7cd1ccd9ca003836c447d706b"
API="https://api.inaturalist.org/v1/observations"
CUTOFF="2026-10-28"
PHOTO_URI=re.compile(r"^/photos/(\d+)/(?:square|small|thumb|medium|large|original)\.(jpe?g|png)$")
CC={"cc-by","cc-by-nc","cc-by-sa","cc-by-nc-sa","cc0"}
SOURCE_42111=42111
DATE_CUTOFF=datetime(2026,10,28,tzinfo=timezone.utc)


def choose_zero(atlas,sample=64,offset=0):
    require(len(atlas)==SOURCE_42111 and atlas.inat_taxon_id.is_unique,
            "Full source 42,111 iNaturalist taxa must be unique")
    require(1<=sample<=128 and 0<=offset<=2229,
            "Bound API sampling so source is reproducible and rate-limited")
    x=atlas.copy()
    x["inat_taxon_id"]=x.inat_taxon_id.fillna("").astype(str)
    x["after_strict9"]=pd.to_numeric(x.after_strict9,errors="raise")
    x=x.loc[x.after_strict9.eq(0) &
            x.inat_taxon_id.str.fullmatch(r"\d+")].copy()
    require(len(x)==2229,"Genuine baseline 2229 zero-numeric-ROI source taxa changed")
    x["fixed_hash"]=x.inat_taxon_id.map(lambda tid:hashlib.sha256((
        "MOIS-NEW-INAT-POST-SNAPSHOT-20261028|"+tid).encode()).hexdigest())
    result=x.sort_values("fixed_hash",kind="stable").iloc[offset:offset+sample]
    require(result.inat_taxon_id.is_unique,"Taxonomic source ID duplicates in pinned sample")
    return result.reset_index(drop=True)


def created_after_freeze(value):
    try:
        iso=str(value or "").replace("Z","+00:00")
        parsed=datetime.fromisoformat(iso)
        require(parsed.tzinfo is not None,"Public source observation creation date timezone missing")
        return parsed.astimezone(timezone.utc)>=DATE_CUTOFF
    except (ValueError,TypeError):
        return False


def strict_new_photo(media,expected_taxon):
    """Only photo-license-approved iNat open-data source with real ID and extension."""
    if not isinstance(media,dict):return None
    photoid=media.get("id")
    if not str(photoid or "").isdigit():return None
    license_=str(media.get("license_code") or media.get("license") or "").strip().lower()
    if license_ not in CC:return None
    photo=str(media.get("url") or "")
    uri=urlsplit(photo)
    if uri.scheme!="https" or (uri.hostname or "").lower()!="inaturalist-open-data.s3.amazonaws.com":
        return None
    match=PHOTO_URI.fullmatch(uri.path)
    if match is None or match.group(1)!=str(photoid):
        return None
    ext=match.group(2)
    return {
      "photo_id":str(photoid),
      "photo_url_large":f"https://inaturalist-open-data.s3.amazonaws.com/photos/{photoid}/large.{ext}",
      "photo_license":license_,
      "photo_medium_API_URL":photo,
      "inat_taxon_id":str(expected_taxon),
    }


def source_rows(results,tid):
    """Verify exact iNaturalist taxon ID, not descendant taxon / genus / alias."""
    out=[];stats={"new_observation_records":0,"exact_taxon_records":0,
                  "license_eligible_new_photo":0,"discarded_approximate_taxon":0}
    for obs in results:
        if not isinstance(obs,dict):continue
        if not created_after_freeze(obs.get("created_at")):continue
        stats["new_observation_records"]+=1
        actual=obs.get("taxon") or {}
        if not isinstance(actual,dict) or str(actual.get("id"))!=str(tid):
            stats["discarded_approximate_taxon"]+=1;continue
        uuid=str(obs.get("uuid") or "")
        if not re.fullmatch(r"[0-9a-fA-F-]{36}",uuid) or not str(obs.get("id") or "").isdigit():
            continue
        if obs.get("captive") is True:
            continue
        stats["exact_taxon_records"]+=1
        for media in obs.get("photos",[]) or []:
            verified=strict_new_photo(media,tid)
            if verified is None:continue
            verified.update({
               "observation_id":uuid,
               "observation_numeric_id":str(obs["id"]),
               "source_observation_created_at":str(obs["created_at"]),
               "observed_on":str(obs.get("observed_on") or ""),
               "source_api":"https://api.inaturalist.org/v1/observations",
               "source_post_snapshot_date":CUTOFF,
               "source_frozen93m_photo_ID_excluded_by_time_assumption_only":False,
               "previous_1_231_244_photo_IDs_verified_disjoint":False,
               "new_raw_original_image_pixels_downloaded":False,
               "strict9_original_RGB9_measured":False,
               "species_wide_colour_fixed_proven":False,
            })
            out.append(verified);stats["license_eligible_new_photo"]+=1
    return out,stats


def request_recent(taxon_id,*,limit=20,retry=3):
    params={"taxon_id":str(taxon_id),"created_d1":CUTOFF,
            "photos":"true","photo_licensed":"true",
            "per_page":str(limit),"page":"1",
            "order_by":"created_at","order":"desc",
            "captive":"false"}
    url=API+"?"+urlencode(params)
    req=Request(url,headers={"Accept":"application/json",
        "User-Agent":"MOIS-photo-ecology-source-zeroROI-research/1.0"})
    last=None
    for trial in range(retry):
        try:
            with urlopen(req,timeout=22) as res:
                data=json.loads(res.read(2_000_000))
            require(isinstance(data,dict) and isinstance(data.get("results"),list),
                    "iNat new endpoint returned unexpected response schema")
            return data
        except (HTTPError,URLError,TimeoutError,ValueError,json.JSONDecodeError) as e:
            last=e
            if isinstance(e,HTTPError) and e.code not in (429,500,502,503,504):
                break
            if trial<retry-1:time.sleep(2**trial)
    raise RuntimeError("iNaturalist public new API failed: "+str(last)[:170])


def run(atlas_path,out,*,sample=64,offset=0,pace=.7,limit=20):
    atlas_path=Path(atlas_path);out=Path(out)
    require(not out.exists() and 0.6<=pace<=5 and 5<=limit<=40,
            "Do not overwrite original source metadata or overrun API")
    require(source_sha(atlas_path)==ATLAS_SHA,
            "True post-wave56 original photo atlas SHA mismatch")
    atlas=pd.read_csv(atlas_path,compression="infer",dtype={"inat_taxon_id":str})
    taxa=choose_zero(atlas,sample=sample,offset=offset)
    require(len(taxa)>0,"Source taxonomic draw offset outside observed 2229 zero-ROI taxa")
    stats=[];data=[]
    for i,row in enumerate(taxa.itertuples(index=False)):
        d={"inat_taxon_id":row.inat_taxon_id,"source_photo_stratum_offset":offset}
        try:
            response=request_recent(row.inat_taxon_id,limit=limit)
            options,counts=source_rows(response["results"],row.inat_taxon_id)
            d.update(counts)
            d["status"]="source_API_queried"
            d["API_reported_recent_total"]=int(response.get("total_results") or 0)
            # One new photo per independently observed source taxon in this wave.
            if options:
                options.sort(key=lambda p:hashlib.sha256((
                    "MOIS-POST20261028-SOURCE-FIXED|"+
                    p["photo_id"]+"|"+p["observation_id"]).encode()).hexdigest())
                data.append(options[0])
                d["new_api_eligible_independent_photo_opportunity"]=1
            else:d["new_api_eligible_independent_photo_opportunity"]=0
        except Exception as error:
            d["status"]="source_API_error_review"
            d["error"]=str(error)[:190]
            d["new_api_eligible_independent_photo_opportunity"]=0
        stats.append(d)
        if (i+1)%16==0:
            print("MOIS_NEW_POSTFCP_INAT_SOURCE_PILOT",i+1,len(taxa),flush=True)
        time.sleep(pace)
    out.mkdir(parents=True)
    photos=out/"new_post20261027_inat_zeroROI_photo_candidates.csv"
    statuses=out/"new_post20261027_inat_zeroROI_taxa_source_status.csv"
    frame=pd.DataFrame(data,columns=[
        "inat_taxon_id","photo_id","photo_url_large","photo_license",
        "photo_medium_API_URL","observation_id","observation_numeric_id",
        "source_observation_created_at","observed_on","source_api",
        "source_post_snapshot_date","source_frozen93m_photo_ID_excluded_by_time_assumption_only",
        "previous_1_231_244_photo_IDs_verified_disjoint",
        "new_raw_original_image_pixels_downloaded","strict9_original_RGB9_measured",
        "species_wide_colour_fixed_proven"
    ])
    require(frame.photo_id.is_unique and frame.observation_id.is_unique and
            frame.inat_taxon_id.is_unique,
            "Multiple new original photos or observations repeated across source taxon IDs")
    frame.to_csv(photos,index=False)
    pd.DataFrame(stats).to_csv(statuses,index=False)
    report={
        "schema":SCHEMA,
        "previous_full_original_photo_metadata_source_count":93209351,
        "original_source_nominal_taxa":SOURCE_42111,
        "original_source_zero_ROI9_taxa":2229,
        "post_snapshot_min_created_date":CUTOFF,
        "sampled_zero_ROI9_source_taxa":len(taxa),
        "source_deterministic_hash_offset":offset,
        "source_taxa_with_new_post_snapshot_API_photo_opportunities":len(frame),
        "new_independent_original_iNaturalist_photo_IDs_seen":len(frame),
        "source_taxa_API_errors":sum(d["status"]=="source_API_error_review" for d in stats),
        "old_photo_ID_1_231_244_ledger_disjointness_verified":False,
        "original_photo_pixels_downloaded":0,
        "genuine_nine_fraction_ROI_measured":0,
        "new_species_nine_hue_colours_proven":0,
        "source_IUCN_captive_false_filter":True,
        "output_sha256":{f.name:source_sha(f) for f in (photos,statuses)}
    }
    (out/"receipt.json").write_text(json.dumps(report,sort_keys=True,indent=2)+"\n")
    print("MOIS_NEW_INAT_POST_SNAPSHOT_64TAXA_ZERO9_SOURCE_PILOT",
          json.dumps(report,sort_keys=True),flush=True)
    return report


if __name__=="__main__":
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument("--atlas",required=True)
    a.add_argument("--output",required=True)
    a.add_argument("--sample",type=int,default=64)
    a.add_argument("--offset",type=int,default=0)
    a.add_argument("--pace",type=float,default=.7)
    a.add_argument("--limit",type=int,default=20)
    x=a.parse_args()
    run(x.atlas,x.output,sample=x.sample,offset=x.offset,pace=x.pace,limit=x.limit)
