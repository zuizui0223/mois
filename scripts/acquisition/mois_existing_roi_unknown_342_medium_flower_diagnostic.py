#!/usr/bin/env python3
"""Bounded detector-only diagnosis of 342 ALREADY measured, strictROI9-unknown photos.

The single supplied SHA-locked 512-photo review table is independent of
photographic outcomes. Exactly 170 were previously screened at 500px and
are NEVER fetched by this worker; 342 previously measured image IDs need
their first low-res screen. Read medium.jpeg/.jpg/.png at most once per URL
plus bounded retries, without saving or redistributing photo RGB bytes.
Frozen FCP YOLO11 flower BOX (NOT corolla segmentation) records box confidence,
relative apparent area, and detection failure. No binary flower/not-flower
ground truth or white/nonwhite label is inferred; all prior strict ROI9 unknowns
remain unknown, and screen negatives remain eligible for botanical review.
"""
from __future__ import annotations

import argparse,csv,hashlib,io,json,math,time
from collections import Counter
from pathlib import Path
from urllib.error import HTTPError,URLError
from urllib.request import Request,urlopen

from PIL import Image,ImageOps
import numpy as np

SCHEMA="mois_source_old512_originalROI9_unknown_342_medium_RGB_flowerbox_diagnostics_v1"
MANIFEST_SCHEMA="mois_512_prior_strictROI9_unknown_existing_medium500_reuse_and_bounded_screen_v1"
MODEL_SHA="f1aaeec4664fe2c178e5cf2bc1f508977bef3e4aa7b40613026cb8ae3de789d5"
EXPECTED_TOTAL=342
GROUPS=8
BASE="https://inaturalist-open-data.s3.amazonaws.com/photos"
CONF=0.02
IOU=0.7
MAX_DET=50
DIAGNOSTIC_CLASSES=(
 "detectable_large_flower_box_ROI9_failed",
 "detectable_medium_flower_box_ROI9_failed",
 "detectable_tiny_flower_box_ROI9_failed",
 "detector_negative_requires_botanical_review",
 "fetch_or_decode_error_without_botanical_inference",
 "detector_model_error_without_botanical_inference",
)
FIELDS=(
 "photo_id","observation_id","inat_taxon_id","source_design",
 "latitude_zone","local_season","postoutcome_source_gap",
 "previous_original_photo_attempts","source_photo_opportunities",
 "screen_status","technical_diagnostic_class",
 "500px_pixel_decoded","actual_remote_medium_source_uri",
 "actual_medium_photo_sha256","source_image_width","source_image_height",
 "flower_box_count_002","largest_box_confidence","largest_box_area_fraction",
 "max_flower_box_confidence","n_boxes_area_ge_0015",
 "n_boxes_area_ge_006","image_fetch_attempted_variants",
 "error_type","strict_original_roi9_remains_unknown"
)

def require(ok,msg):
    if not bool(ok):raise ValueError("MOIS_EXISTING_342_RGB_TECH: "+msg)

def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda:f.read(4<<20),b""):h.update(b)
    return h.hexdigest()

class RequestPacer:
    def __init__(self,requests_per_second):
        require(0<requests_per_second<=1,"Per-shard S3 request rate exceeds safe bound")
        self.interval=1.0/requests_per_second
        self.last=None
    def throttle(self):
        now=time.monotonic()
        if self.last is not None:
            pause=self.last+self.interval-now
            if pause>0:time.sleep(pause)
        self.last=time.monotonic()

def _fetch_bytes(url,pacer):
    pacer.throttle()
    request=Request(url,headers={"User-Agent":"MOIS-ROI9-technical-review/1.0"})
    with urlopen(request,timeout=16) as response:
        typ=response.headers.get("Content-Type","").lower()
        if typ and not (typ.startswith("image/") or
                        typ.startswith("application/octet-stream")):
            raise ValueError("remote_source_not_image_content_type")
        data=response.read(3_000_001)
        if len(data)>3_000_000:
            raise ValueError("medium_image_larger_than_3MB")
        if not data:
            raise ValueError("empty_medium_image_response")
        return data

def photo_medium(photo_id,pacer,*,fetch=None):
    """Resolve original iNaturalist S3 extensions without guessing phenotype."""
    require(str(photo_id).isdecimal(),"Invalid signed numeric photo ID")
    fetching=_fetch_bytes if fetch is None else fetch
    counts=0
    last_error="none"
    for extension in ("jpeg","jpg","png"):
        url=f"{BASE}/{photo_id}/medium.{extension}"
        for trial in range(2):
            counts+=1
            try:
                data=fetching(url,pacer)
                photo=Image.open(io.BytesIO(data))
                photo=ImageOps.exif_transpose(photo).convert("RGB")
                photo.thumbnail((500,500),Image.Resampling.BILINEAR)
                require(min(photo.size)>0,"Image cannot have zero pixel dimension")
                return photo,url,hashlib.sha256(data).hexdigest(),counts,""
            except HTTPError as e:
                last_error=f"HTTP_{e.code}"
                if e.code in (400,401,403,404):
                    break
                if e.code not in (429,500,502,503,504):
                    break
                if trial==0:time.sleep(1.0)
            except (URLError,TimeoutError,ConnectionError,OSError,ValueError) as e:
                last_error=type(e).__name__+":"+str(e)[:80]
                if trial==0:time.sleep(1.0)
    return None,"","",counts,last_error

def classify_box_predictions(photo,prediction):
    """Numerical *image box* geometry; cannot assert open corolla or colour."""
    if getattr(prediction,"boxes",None) is None or len(prediction.boxes)==0:
        return {"screen_status":"negative",
                "technical_diagnostic_class":
                    "detector_negative_requires_botanical_review",
                "flower_box_count_002":0,
                "largest_box_confidence":0.0,"largest_box_area_fraction":0.0,
                "max_flower_box_confidence":0.0,
                "n_boxes_area_ge_0015":0,"n_boxes_area_ge_006":0}
    coords=prediction.boxes.xyxy.cpu().numpy().astype(float)
    scores=prediction.boxes.conf.cpu().numpy().astype(float)
    require(coords.ndim==2 and coords.shape[1]==4 and
            len(scores)==len(coords) and len(coords)<=MAX_DET,
            "Frozen YOLO flower detector produced invalid box geometry")
    w,h=photo.size
    a=np.clip(coords[:,2]-coords[:,0],0,w)*np.clip(coords[:,3]-coords[:,1],0,h)/(w*h)
    require(np.isfinite(a).all() and np.isfinite(scores).all() and
            (a>=0).all() and (a<=1+1e-5).all() and
            (scores>=0).all() and (scores<=1).all(),
            "Detector supplied nonphysical confidence or spatial box area")
    j=int(np.argmax(a));largest=float(a[j])
    tag=("detectable_large_flower_box_ROI9_failed" if largest>=.06 else
         "detectable_medium_flower_box_ROI9_failed" if largest>=.015 else
         "detectable_tiny_flower_box_ROI9_failed")
    return {"screen_status":"candidate","technical_diagnostic_class":tag,
            "flower_box_count_002":len(scores),
            "largest_box_confidence":float(scores[j]),
            "largest_box_area_fraction":largest,
            "max_flower_box_confidence":float(np.max(scores)),
            "n_boxes_area_ge_0015":int((a>=.015).sum()),
            "n_boxes_area_ge_006":int((a>=.06).sum())}

def run(source,weights,shard,output,*,detector=None,fetch=None,
        requests_per_second=.65,verify_model=True):
    source,weights,output=map(Path,(source,weights,output))
    require(not output.exists() and 0<=shard<GROUPS,
            "Output exists or invalid bounded source diagnostic shard")
    m=json.loads((source/"manifest.json").read_text())
    require(m.get("schema")==MANIFEST_SCHEMA and
            m.get("source_review_photo_rows")==512 and
            m.get("existing_500px_photo_screen_results_reused")==170 and
            m.get("existing_unknown_photos_needing_small_bounded_500px_recheck")==342 and
            m.get("new_photo_IDs_selected_for_world_census")==0,
            "Not original SHA-signed pre-existing photo technical review cohort")
    filename=source/"workers"/f"shard_{shard:02d}.csv"
    require(sha(filename)==m["worker_source_SHA256"][str(shard)],
            "Previously selected source original photo worker SHA changed")
    with filename.open(encoding="utf-8",newline="") as f:
        records=list(csv.DictReader(f))
    require(len(records)==m["eight_worker_photo_counts"][shard]
            and len(records) in (42,43)
            and len(set(row["photo_id"] for row in records))==len(records)
            and len(set(row["observation_id"] for row in records))==len(records),
            "Previously acquired original diagnostic photo/observation IDs inconsistent")
    if verify_model:
        require(weights.is_file() and sha(weights)==MODEL_SHA,
                "Frozen source FCP YOLO11 flower detector numeric model SHA mismatch")
    if detector is None:
        from ultralytics import YOLO
        detector=YOLO(str(weights))
    pace=RequestPacer(requests_per_second)
    rows=[]
    images=[]
    eligible=[]
    def infer_pending():
        if not images:return
        try:
            results=detector.predict(
                source=images,imgsz=512,conf=CONF,iou=IOU,
                max_det=MAX_DET,augment=False,device="cpu",verbose=False,
                save=False,save_txt=False,save_crop=False)
            require(len(results)==len(images),"YOLO lost earlier image detections")
            for index,pred in zip(eligible,results):
                d=classify_box_predictions(images[eligible.index(index)],pred)
                rows[index].update(d)
        except Exception as err:
            for index in eligible:
                rows[index].update(
                    screen_status="error",
                    technical_diagnostic_class=
                        "detector_model_error_without_botanical_inference",
                    error_type=type(err).__name__+":"+str(err)[:130])
        images.clear();eligible.clear()
    for original in records:
        photo,url,digest,attempts,problem=photo_medium(
            original["photo_id"],pace,fetch=fetch)
        row={k:"" for k in FIELDS}
        for k in ("photo_id","observation_id","inat_taxon_id",
                  "source_design","latitude_zone","local_season",
                  "postoutcome_source_gap","previous_original_photo_attempts",
                  "source_photo_opportunities"):
            row[k]=original[k]
        row["strict_original_roi9_remains_unknown"]=True
        row["image_fetch_attempted_variants"]=attempts
        if photo is None:
            row.update(screen_status="error",
              technical_diagnostic_class=
                "fetch_or_decode_error_without_botanical_inference",
              error_type=problem,**{"500px_pixel_decoded":False})
        else:
            row.update(**{
                "500px_pixel_decoded":True,
                "actual_remote_medium_source_uri":url,
                "actual_medium_photo_sha256":digest,
                "source_image_width":photo.width,
                "source_image_height":photo.height
            })
        rows.append(row)
        if photo is not None:
            images.append(photo)
            eligible.append(len(rows)-1)
        if len(images)>=12:infer_pending()
    infer_pending()
    require(len(rows)==len(records) and len(set(r["photo_id"] for r in rows))==
            len(records), "One actual photo status per photo ID is required")
    require(all(r["screen_status"] in ("candidate","negative","error")
                and r["technical_diagnostic_class"] in DIAGNOSTIC_CLASSES
                and r["strict_original_roi9_remains_unknown"] is True for r in rows),
            "Unknown source photo was dropped or given invented biological classification")
    output.mkdir(parents=True)
    outfile=output/f"shard_{shard:02d}_original_medium_diagnostic.csv"
    with outfile.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=FIELDS,lineterminator="\n")
        w.writeheader();w.writerows(rows)
    tally=Counter(r["screen_status"] for r in rows)
    cases=Counter(r["technical_diagnostic_class"] for r in rows)
    receipt={
        "schema":SCHEMA,"worker":shard,
        "original_old_photo_ids_attempted_for_first_medium_diagnostic":len(rows),
        "actual_new_source_photo_IDs_added_to_world_census":0,
        "old_source_ROI9_unknown_statuses_preserved":len(rows),
        "decoded_medium500":sum(r["500px_pixel_decoded"] is True for r in rows),
        "statuses":dict(tally),"technical_box_classes":dict(cases),
        "frozen_detector_confidence_threshold":CONF,
        "model_SHA256":MODEL_SHA,
        "original_source_manifest_SHA256":sha(source/"manifest.json"),
        "individual_original_source_worker_SHA256":sha(filename),
        "measurement_CSV_SHA256":sha(outfile),
        "no_raw_image_bytes_archived":True,
        "no_EfficientSAM_or_ROI9_inference":True,
        "YOLO_box_not_independent_botanical_flower_organ_truth":True,
        "prior_ROI9_unknowns_never_imputed_as_white_or_pigmented":True
    }
    (output/f"shard_{shard:02d}_source_receipt.json").write_text(
        json.dumps(receipt,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("MOIS_OLD_ROI9_UNKNOWN_FLOWER_BOX_DIAG_SHARD",json.dumps({
        "part":shard,"old_photo_attempts":len(rows),
        "medium_decoded":receipt["decoded_medium500"],
        "statuses":dict(tally),"source_prior_photo_RGB_reused":False},
        sort_keys=True),flush=True)
    return receipt

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("source","weights","output"):
        p.add_argument("--"+k,type=Path,required=True)
    p.add_argument("--shard",type=int,required=True)
    p.add_argument("--rps",type=float,default=.65)
    x=p.parse_args()
    run(x.source,x.weights,x.shard,x.output,requests_per_second=x.rps)

if __name__=="__main__":
    main()
