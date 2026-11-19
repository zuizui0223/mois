#!/usr/bin/env python3
"""MOIS source-backed targeted original iNaturalist photo evidence wave.

Reuses frozen 199-genus queue from MOIS paleo PR #206, the original COLXR
37,104 species-to-iNaturalist ID crosswalk, and original 1,231,337
photo/observation history. Collects CC-licensed true photo bytes only from
exact species iNaturalist observations and exports credited contact sheets.
No automatic flower stage or botanical petal/organ identity is inferred.
No image RGB label can be promoted to confirmed flower color.
"""
from __future__ import annotations
import argparse
import hashlib
import io
import json
import re
import time
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path
from urllib.error import HTTPError,URLError
from urllib.parse import urlencode
from urllib.request import Request,urlopen
from zipfile import ZipFile

import numpy as np
import pandas as pd
from PIL import Image,ImageDraw,ImageFont,ImageOps

from scripts.acquisition.mois_zero9_recent_inat_new_created_photo_source import (
    strict_new_photo, CC,
)
from scripts.acquisition.mois_20261114_new_created_disjoint_original_photo9 import history

SCHEMA="mois_20261119_targeted_16_native_genus_flower_original_proof_v1"
SOURCE_QUEUE_BLOB="efc85eb9cc9fa842d630d69ab1962fa7a5029ab1"
SOURCE_ENV_ZIP_SHA="d8b97528188d9f8d63cb802d5b3f6e8f7f251937fd30af42aeae47e43ebad1fc"
SOURCE_ENV_MEMBER="BIO19_weather_NDVI_ESACCIv7_AGB_37104_accepted_species.parquet"
SOURCE_ID_COUNT=1231337
API="https://api.inaturalist.org/v1/observations"
TARGET_SPECIES_MAX=16
TARGET_OBS=3
DEFAULT_PAGE_LIMIT=100
MAX_IMAGE_BYTES=16_000_000
IMAGE_SIDE=440

def require(ok,message):
    if not bool(ok):
        raise ValueError("MOIS_TARGETED_PETAL_PHOTO_PROOF: "+message)

def sha_bytes(raw):
    return hashlib.sha256(raw).hexdigest()

def sha_file(path):
    return sha_bytes(Path(path).read_bytes())

def read_source_queue(path):
    path=Path(path)
    require(path.exists(),"Frozen targeted genus CSV was not restored")
    queue=pd.read_csv(path,dtype=str,keep_default_na=False)
    wanted={"rank","batch","genus","candidate_species","candidate_source_specieskey",
            "source_photo9_label","missing_confirmed_class","source_hue_status",
            "not_confirmed_petal_colour"}
    require(wanted.issubset(queue),"Genus target source fields missing")
    require(len(queue)==199 and queue["genus"].is_unique and
            queue["candidate_source_specieskey"].is_unique and
            queue.iloc[:75].batch.eq("NEXT_75_GENERA").all() and
            queue.iloc[:75].source_hue_status.eq("provisional_single").all(),
            "Frozen prioritized 199 genera/75 next targets changed")
    require(queue["not_confirmed_petal_colour"].eq("True").all(),
            "Input falsely states botanically verified petal colour")
    return queue

def source_species_map(path,queue):
    import pyarrow.parquet as pq
    require(sha_file(path)==SOURCE_ENV_ZIP_SHA,
            "37,104 accepted photo9 species ZIP SHA changed")
    with ZipFile(path) as z:
        receipt=json.loads(z.read("receipt.json"))
        raw=z.read(SOURCE_ENV_MEMBER)
    require(sha_bytes(raw)==receipt["output_SHA256"][SOURCE_ENV_MEMBER] and
            receipt["source_global_accepted_species"]==37104,
            "Source species crosswalk and internal SHA differ")
    data=pq.read_table(io.BytesIO(raw),columns=[
        "key","inat_taxon_id","species","source_hue9"]).to_pandas()
    data[["key","inat_taxon_id","species","source_hue9"]]=data[
        ["key","inat_taxon_id","species","source_hue9"]].astype(str)
    require(len(data)==37104 and data.key.is_unique and
            data.inat_taxon_id.str.fullmatch(r"\d+").all(),
            "Source exact iNaturalist taxon ID map incomplete")
    chosen=queue.iloc[:TARGET_SPECIES_MAX].merge(
        data,left_on="candidate_source_specieskey",right_on="key",
        how="left",validate="one_to_one")
    require(chosen.inat_taxon_id.notna().all() and
            chosen.candidate_species.eq(chosen.species).all() and
            chosen.source_photo9_label.eq(chosen.source_hue9).all(),
            "Crosswalk exact species identity or photo hue does not match source")
    chosen["inat_taxon_id"]=chosen.inat_taxon_id.astype(str)
    return chosen

def request_observations(taxon_id):
    params={"taxon_id":str(taxon_id),"photos":"true",
            "photo_licensed":"true","quality_grade":"research",
            "per_page":DEFAULT_PAGE_LIMIT,"page":1,
            "order_by":"observed_on","order":"desc","captive":"false"}
    url=API+"?"+urlencode(params)
    req=Request(url,headers={"Accept":"application/json",
                "User-Agent":"MOIS-academic-actual-photo-ID-petal-QA/1.0"})
    last=None
    for i in range(4):
        try:
            with urlopen(req,timeout=30) as res:
                require(res.status==200,"iNaturalist HTTP status not 200")
                raw=res.read(5_000_001)
            require(len(raw)<=5_000_000,"Overlarge iNaturalist API result")
            data=json.loads(raw)
            require(isinstance(data,dict) and isinstance(data.get("results"),list),
                    "Invalid iNaturalist result shape")
            return url,data
        except (HTTPError,URLError,TimeoutError,ValueError,json.JSONDecodeError) as exc:
            last=str(exc)
            if i<3:time.sleep(min(8,2**i))
    raise ValueError("iNaturalist source API failure: "+str(last))

def candidates(data,tid,prior_ids,prior_obs):
    options=[]
    stats=Counter()
    for obs in data["results"]:
        stats["seen_observations"]+=1
        taxon=obs.get("taxon") or {}
        if not isinstance(taxon,dict) or str(taxon.get("id"))!=str(tid):
            stats["not_exact_taxon"]+=1;continue
        if obs.get("captive") is True:
            stats["captive"]+=1;continue
        uuid=str(obs.get("uuid") or "")
        if not re.fullmatch(r"[0-9a-fA-F-]{36}",uuid):
            stats["invalid_observation_uuid"]+=1;continue
        if uuid in prior_obs:
            stats["used_prior_observation"]+=1;continue
        photo=None
        for candidate in obs.get("photos") or []:
            f=strict_new_photo(candidate,tid)
            if f is None:
                stats["no_eligible_open_CC_photo"]+=1;continue
            if f["photo_id"] in prior_ids:
                stats["used_prior_photo"]+=1;continue
            photo=f;break
        if photo is None:continue
        photo.update({
          "observation_id":uuid,
          "observation_numeric_id":str(obs.get("id","")),
          "inat_taxon_id":str(tid),
          "observed_on":str(obs.get("observed_on") or ""),
          "observer_login":str((obs.get("user") or {}).get("login") or ""),
          "observation_license":str(obs.get("license_code") or ""),
          "plant_phenology_annotations_json":json.dumps(obs.get("annotations") or []),
          "source_taxon_exact_id_verified":True,
          "source_research_grade":str(obs.get("quality_grade") or ""),
        })
        options.append(photo)
        stats["licensed_independent_options"]+=1
    options.sort(key=lambda x:(x["observed_on"],x["observation_id"]),reverse=True)
    chosen=[]
    seen_obs=set(); seen_photos=set()
    for photo in options:
        if photo["observation_id"] in seen_obs or photo["photo_id"] in seen_photos:continue
        seen_obs.add(photo["observation_id"]);seen_photos.add(photo["photo_id"])
        chosen.append(photo)
        if len(chosen)>=TARGET_OBS:break
    return chosen,dict(stats)

def fetch_original(photo):
    url=photo["photo_url_large"].replace("/large.","/original.")
    require(url.startswith("https://inaturalist-open-data.s3.amazonaws.com/photos/") and
            url.rsplit("/",1)[-1].startswith(photo["photo_id"]+"."),
            "Original photo URL source host/ID mismatch")
    req=Request(url,headers={"User-Agent":"MOIS-open-licensed-botanical-source-QA/1.0"})
    error=None
    for i in range(3):
        try:
            with urlopen(req,timeout=32) as res:
                require(res.status==200,"Original image HTTP non-200")
                raw=res.read(MAX_IMAGE_BYTES+1)
            require(0<len(raw)<=MAX_IMAGE_BYTES,"Original photo too large or empty")
            with Image.open(io.BytesIO(raw)) as im:
                im.load()
                require(im.width>=250 and im.height>=250,
                        "Original photo resolution too low")
                image=ImageOps.exif_transpose(im).convert("RGB")
            return raw,image,url
        except Exception as e:
            error=str(e)
            if i<2:time.sleep(i+1)
    raise ValueError("Original image unavailable: "+str(error))

def thumbnail(original):
    im=ImageOps.contain(original,(IMAGE_SIDE,IMAGE_SIDE))
    tile=Image.new("RGB",(IMAGE_SIDE,IMAGE_SIDE),(250,250,250))
    tile.paste(im,((IMAGE_SIDE-im.width)//2,(IMAGE_SIDE-im.height)//2))
    return tile

def contact_sheet(name,photos,target_dir):
    W=3*IMAGE_SIDE
    tile_h=IMAGE_SIDE+100
    sheet=Image.new("RGB",(W,tile_h),(255,255,255))
    draw=ImageDraw.Draw(sheet)
    for j,item in enumerate(photos):
        ix=j*IMAGE_SIDE
        sheet.paste(item["thumbnail"],(ix,0))
        text=f"{name}\nPhoto {item['photo_id']}\n{item['photo_license']} | {item['observer_login'][:28]}"
        draw.multiline_text((ix+8,IMAGE_SIDE+6),text,fill=(10,10,10),spacing=4)
    out=target_dir/"contact.jpg"
    sheet.save(out,format="JPEG",quality=82,optimize=True)
    return out

def run(args):
    root=Path(args.out)
    require(not root.exists(),"Refuse overwriting photo-identity evidence")
    queue=read_source_queue(args.targets)
    target=source_species_map(args.source_env,queue)
    old_ids,old_obs,prior_raw=history()
    require(len(old_ids)==SOURCE_ID_COUNT and len(prior_raw)>0,
            "Prior photo registry not independently recovered")
    stats=[];metadata=[];pending=[]
    root.mkdir(parents=True)
    for idx,row in enumerate(target.itertuples(index=False)):
        tid=row.inat_taxon_id
        entry={"rank":str(row.rank),"candidate_species":row.candidate_species,
               "genus":row.genus,"inat_taxon_id":tid,
               "original_photo9_proxy_label":row.source_photo9_label,
               "status":"SOURCE_SCREEN_PENDING",
               "n_images_original_bytes_verified":0,
               "n_independent_observation_photos":0,
               "already_held_colour_unverified":True}
        try:
            api_url,res=request_observations(tid)
            eligible,details=candidates(res,tid,old_ids,old_obs)
            entry["iNat_total_results_reported"]=int(res.get("total_results") or 0)
            entry.update(details)
            entry["API_url"]=api_url
            out_dir=root/("species_"+str(row.rank).zfill(3)+"_"+tid)
            out_dir.mkdir(exist_ok=True)
            images=[]
            for photo in eligible:
                try:
                    raw,img,url=fetch_original(photo)
                    digest=sha_bytes(raw)
                    if digest in prior_raw:
                        entry["raw_SHA_collision_recent93"]=True
                        continue
                    photo["original_image_sha256"]=digest
                    photo["original_image_bytes"]=len(raw)
                    photo["original_image_dimensions"]=[img.width,img.height]
                    photo["original_url"]=url
                    photo["candidate_species"]=row.candidate_species
                    photo["rank"]=int(row.rank)
                    photo["true_botanical_petal_ROI_verified"]=False
                    photo["flowering_observation_independently_verified"]=False
                    photo["raw_original_pixels_persisted"]=False
                    photo["thumbnail"]=thumbnail(img)
                    images.append(photo)
                    metadata.append({k:v for k,v in photo.items() if k!="thumbnail"})
                    if len(images)>=TARGET_OBS:break
                except Exception as exc:
                    pending.append({"rank":int(row.rank),"inat_taxon_id":tid,
                                    "photo_id":photo["photo_id"],"error":str(exc)[:170]})
            if images:
                fname=contact_sheet(row.candidate_species,images,out_dir)
                entry["contact_sheet"]=str(fname.relative_to(root))
                entry["contact_sha256"]=sha_file(fname)
                entry["n_images_original_bytes_verified"]=len(images)
                entry["n_independent_observation_photos"]=len({
                   x["observation_id"] for x in images})
                entry["status"]="REAL_CC_ORIGINALS_CONTACT_READY_REQUIRES_BOTANICAL_REVIEW"
            else:
                entry["status"]="NO_NEW_ELIGIBLE_PHOTO_REQUIRES_FURTHER_SOURCE"
        except Exception as exc:
            entry["status"]="SOURCE_OR_IDENTITY_GATE_ERROR"
            entry["error"]=str(exc)[:260]
        stats.append(entry)
        print("MOIS_TARGETED_OPEN_CC_ORIGINAL_QUERY",
              json.dumps({"rank":entry["rank"],"species":entry["candidate_species"],
                 "status":entry["status"],
                 "n_original_images":entry["n_images_original_bytes_verified"]}),
              flush=True)
        time.sleep(args.pace)
    statfile=root/"source_taxon_status.csv"
    metfile=root/"actual_CC_original_photo_UUID_SHA_metadata.csv"
    pd.DataFrame(stats).to_csv(statfile,index=False)
    pd.DataFrame(metadata).to_csv(metfile,index=False)
    (root/"photo_failure_review.json").write_text(json.dumps(pending,indent=2)+"\n")
    total=int(sum(x["n_images_original_bytes_verified"] for x in stats))
    observations=[m["observation_id"] for m in metadata]
    ids=[m["photo_id"] for m in metadata]
    require(len(ids)==len(set(ids)) and len(observations)==len(set(observations)),
            "New batch itself repeats photo or observation ID")
    digests=[m["original_image_sha256"] for m in metadata]
    require(len(digests)==len(set(digests)),"New batch repeats original image bytes")
    missing=[x["candidate_species"] for x in stats
             if x["n_images_original_bytes_verified"]<TARGET_OBS]
    contact=root/"contact_sheet_index.html"
    links=[]
    for x in stats:
        if x.get("contact_sheet"):
            links.append(f"<p>{x['candidate_species']} — {x['n_images_original_bytes_verified']} images<br><img width='950' src='{x['contact_sheet']}'></p>")
    # Static locally inspectable image index. No automated petal/anthesis label.
    contact.write_text("<!doctype html><html><meta charset='utf-8'><title>MOIS 16 Species Photo QA</title><body><h2>Unconfirmed botanical petals — real CC photos</h2>"+
            "\n".join(links)+"</body></html>")
    result={"schema":SCHEMA,"source_queue_github_blob_sha1":SOURCE_QUEUE_BLOB,
      "source_original_ENV_ZIP_SHA256":SOURCE_ENV_ZIP_SHA,
      "source_previously_attempted_photo_registry_count":len(old_ids),
      "old_observation_UUIDs_checked":len(old_obs),
      "prior_93_original_byte_digests_excluded":len(prior_raw),
      "photo_batch_rank_first":1,"n_exact_species_targeted":len(target),
      "requested_independent_observations_each":TARGET_OBS,
      "n_authentic_CC_original_image_SHA_verified":total,
      "n_species_with_at_least_one_original_photo":int(sum(
          x["n_images_original_bytes_verified"]>0 for x in stats)),
      "n_species_with_3_independent_photos":int(sum(
          x["n_independent_observation_photos"]>=3 for x in stats)),
      "incomplete_species":missing,
      "source_api_gate_errors":int(sum(x["status"]=="SOURCE_OR_IDENTITY_GATE_ERROR" for x in stats)),
      "image_fetch_errors":len(pending),
      "new_images_not_guaranteed_disjoint_from_post_1231337_subsequent_2026_waves":True,
      "full_93million_photo_IDs_not_excluded":True,
      "botanical_flowering_organ_verified":False,
      "new_RGB_9_palette_quantified":False,
      "flower_fixed_colour_or_LGM_adaptation_established":False,
      "raw_original_bytes_transient_thumbnails_only_persisted":True,
      "contact_sheets_for_human_petal_organ_review":True,
      "outputs_sha256":{
         x.name:sha_file(x) for x in [statfile,metfile,contact,root/"photo_failure_review.json"]},
      "contact_sha256":{x["contact_sheet"]:x["contact_sha256"] for x in stats if x.get("contact_sheet")}}
    (root/"receipt.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print("MOIS_TARGETED_16_PHOTO_PROOF_FINISHED",
          json.dumps(result,sort_keys=True),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--targets",required=True)
    p.add_argument("--source-env",required=True)
    p.add_argument("--out",required=True)
    p.add_argument("--pace",type=float,default=0.9)
    args=p.parse_args()
    require(.6<=args.pace<=2.0,"Provider API request pacing must remain bounded")
    run(args)
