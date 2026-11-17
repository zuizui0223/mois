#!/usr/bin/env python3
"""Blind visual *pilot* for frozen MOIS QA180; NOT a botanical truth annotation.

Select 18 photos solely by pre-frozen blind species IDs (001,011,...171),
fetch current public attributed iNaturalist large derivatives, and save two
3x3 contact sheets for manual review. No original fullsize images, no species
photo-colour predictions, no atlas updates or implied expert certification.
Photos and previews MUST NOT be committed to Git.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,time,zipfile
from pathlib import Path
from urllib.parse import urlsplit,urlunsplit
import requests
from PIL import Image, ImageDraw, ImageFont, ImageOps

def req(ok,msg):
    if not ok:raise ValueError("MOIS_QA180_BLIND_PILOT: "+str(msg))
def sha(b):return hashlib.sha256(b).hexdigest()

def draw_font(size):
    for font in ("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
                 "/usr/share/fonts/truetype/liberation2/LiberationSans-Regular.ttf"):
        try:return ImageFont.truetype(font,size)
        except OSError:pass
    return ImageFont.load_default()

def prepare(blind,metadata):
    with open(blind,newline="",encoding="utf-8-sig") as f:rows=list(csv.DictReader(f))
    req(len(rows)==538 and len({x["blind_photo_review_id"] for x in rows})==538,
        "frozen original photo sheet 538 changed")
    selected=[]
    for number in range(1,181,10):
        ident=f"MOIS-QA-{number:03d}-P01"
        row=[x for x in rows if x["blind_photo_review_id"]==ident]
        req(len(row)==1,"no exact blind photo "+ident)
        selected.append(row[0])
    req(len(selected)==18 and len({x["original_photo_page"] for x in selected})==18,
        "pilot sample not disjoint")

    with zipfile.ZipFile(metadata) as z:
        with z.open("original_photo_current_iNaturalist_API_status.csv") as f:
            values={x["blind_photo_review_id"]:x for x in csv.DictReader(
                io.TextIOWrapper(f,encoding="utf-8"))}
    req(len(values)==538,"public API QA meta not full 538")
    for x in selected:
        row=values[x["blind_photo_review_id"]]
        req(x["original_photo_page"].rsplit("/",1)[-1]==row["original_photo_id"],
            "photo identity mismatch")
        req(row["same_original_photo_in_current_observation"].lower()=="true",
            "historical photo absent from original canonical observation")
        req(row["current_photo_license_code"] not in ("","nan"),
            "photo license metadata missing")
        x["canonical_observation_id"]=row["current_numeric_observation_id"]
        x["expected_license"]=row["current_photo_license_code"]
    return selected

def get_attributed_photos(selected,session):
    ids=list(dict.fromkeys(x["canonical_observation_id"] for x in selected))
    api="https://api.inaturalist.org/v1/observations/"+",".join(ids)
    r=session.get(api,headers={"User-Agent":"MOIS-538-blind-petal-pilot-source-identity/1.0"},
                  timeout=(8,30))
    r.raise_for_status()
    j=r.json()
    req(len(j.get("results",[]))==len(ids),"wrong iNaturalist API batch size")
    lookup={str(x["id"]):x for x in j["results"]}
    req(set(lookup)==set(ids),"iNaturalist observation IDs drift")
    for item in selected:
        obs=lookup[item["canonical_observation_id"]]
        pid=item["original_photo_page"].rsplit("/",1)[-1]
        photos=[p for p in obs.get("photos",[]) if str(p.get("id"))==pid]
        req(len(photos)==1,"original photo not in current observation")
        photo=photos[0]
        src=photo.get("url") or ""
        u=urlsplit(src)
        req(u.scheme=="https" and u.hostname=="inaturalist-open-data.s3.amazonaws.com"
            and f"/photos/{pid}/" in u.path and
            u.path.rsplit(".",1)[-1] in ("jpg","jpeg","JPG") and not u.query,
            "untrusted source photo derivative URL")
        req(photo.get("license_code")==item["expected_license"],
            "photo license mismatch original metadata")
        item["source_large_url"]=urlunsplit((u.scheme,u.netloc,u.path.replace("/square.","/large."),"",""))
        req(item["source_large_url"]!=src,"source image is not square derivative URL")
        item["attribution"]=(photo.get("attribution") or "").strip()
        req(item["attribution"],"source attribution missing")
    return selected

def make_contact_sheet(items,index,path):
    page=Image.new("RGB",(3*510,3*435),(247,249,250))
    draw=ImageDraw.Draw(page)
    f1=draw_font(21);f2=draw_font(16)
    for k,item in enumerate(items):
        col=k%3;row=k//3
        x0=col*510;y0=row*435
        tile=item["thumb"].copy()
        w,h=tile.size
        page.paste(tile,(x0+(510-w)//2,y0+15+(352-h)//2))
        draw.text((x0+12,y0+371),item["blind_photo_review_id"],font=f1,fill=(22,47,63))
        draw.text((x0+12,y0+397),"photo "+item["photo_id"]+" • "+item["expected_license"],
                  font=f2,fill=(75,95,107))
        draw.rectangle((x0+2,y0+2,x0+508,y0+433),outline=(193,206,214),width=2)
    page.save(path,format="JPEG",quality=91,optimize=True)

def run(blind,metadata,outdir):
    out=Path(outdir);req(not out.exists(),"no overwrites")
    items=prepare(blind,metadata)
    with requests.Session() as s:
        get_attributed_photos(items,s)
        out.mkdir(parents=True)
        (out/"previews").mkdir()
        receipts=[]
        for k,item in enumerate(items):
            r=s.get(item["source_large_url"],timeout=(7,28),
                    headers={"User-Agent":"MOIS-photo-licensed-noncommercial-blind-review/1.0"},
                    stream=True)
            r.raise_for_status()
            data=r.content
            req(len(data)<20000000 and len(data)>1500,
                "unexpected copyrighted preview size or content")
            im=Image.open(io.BytesIO(data))
            im=ImageOps.exif_transpose(im).convert("RGB")
            w,h=im.size
            req(w>100 and h>100,"bad public source dimensions")
            thumb=im.copy();thumb.thumbnail((480,350))
            p=out/"previews"/(item["blind_photo_review_id"]+".jpg")
            thumb.save(p,format="JPEG",quality=90,optimize=True)
            item["thumb"]=thumb
            receipts.append({
              "blind_id":item["blind_photo_review_id"],"photo_id":item["original_photo_page"].rsplit("/",1)[-1],
              "current_numeric_observation_id":item["canonical_observation_id"],
              "license":item["expected_license"],"attribution":item["attribution"],
              "source_photo_page":item["original_photo_page"],
              "source_derivative_url":item["source_large_url"],"image_size":[w,h],
              "source_large_sha256":sha(data),"local_scaled_preview_SHA256":sha(p.read_bytes()),
              "petal_flower_verification":"NOT_PERFORMED"
            })
            item["photo_id"]=receipts[-1]["photo_id"]
            print("MOIS_QA180_BLIND_VISUAL_ORIGINAL_PHOTO",k+1,item["blind_photo_review_id"],w,h,flush=True)
            time.sleep(.35)
    make_contact_sheet(items[:9],0,out/"blind_contact_sheet_1.jpg")
    make_contact_sheet(items[9:],1,out/"blind_contact_sheet_2.jpg")
    (out/"SOURCE_CREDITS_AND_SCOPE.txt").write_text(
       "18 iNaturalist public photographs used as NONCOMMERCIAL temporary research preview.\n"
       "Author attribution, photo pages and current CC licenses are preserved per photo in\n"
       "source_license_manifest.json. Do not distribute outside licence conditions.\n"
       "No original fullsize photo, floral AI labels, or expert botanical truth produced.\n"
       "Source preview images are not committed to the Git repository.\n",encoding="utf-8")
    (out/"source_license_manifest.json").write_text(
        json.dumps(receipts,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    report={
      "source_blind_sheet_n":538,"pilot_n_photos":18,"pilot_n_unique_species":18,
      "all_preview_source_licenses_and_attributions_matched":True,
      "photo_selection_without_known_model_colour":True,"blinded_species":True,
      "actual_petal_colour_manual_review_completed":0,
      "photo_source_URL_based_on_current_iNaturalist_api":True,
      "pilot_selects_each_10th_blind_species_first_photo":True,
      "image_SHA256":{p.name:sha(p.read_bytes()) for p in out.glob("blind_contact_sheet_*.jpg")}
    }
    (out/"receipt.json").write_text(json.dumps(report,indent=2)+"\n")
    print("MOIS_QA180_CONTACT_SHEETS_DONE",json.dumps(report),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    for x in ("blind","meta","out"):p.add_argument("--"+x,required=True)
    a=p.parse_args()
    run(a.blind,a.meta,a.out)
