#!/usr/bin/env python3
"""Diagnostic-only reconstruction of five frozen FCP ROI-v4 masks for MOIS.

Does not change the model, colour thresholds, source colour labels or photo
catalogue. Downloads five CC-licensed provider LARGE derivatives *only* with
exact previously recorded photo-ID URL and byte SHA gates. Saves monochrome
mask-only PNGs, box geometry and per-photo receipts; never source RGB pixels.
"""
from __future__ import annotations
import argparse,hashlib,io,json,re,sys
from pathlib import Path
from urllib.parse import urlsplit
from zipfile import ZipFile
import numpy as np
import pandas as pd
from PIL import Image,ImageOps

TARGETS = {
    "164665576":"Aristolochia macroura",
    "111641760":"Aristolochia macrocarpa",
    "255188296":"Artabotrys zeylanicus",
    "446129531":"Bignonia sciuripabulum",
    "487988985":"Libidibia parvifolia",
}
PALETTE12=("white","yellow","orange","red","pink","magenta",
           "purple","blue","bronze","green","brown","black")
SOURCE_ZIP_SHA="b35e963d208a82da580dc6d093b95b000530ce9f2204791102c3a70e49e8608f"
WORKER_ZIP_SHA="8d49aeb4b9362f76a3b2b290850b6e7b0a42847052c2e4dc5e9bb7e25d79bce4"
SHARD_SHA="6c035941f6b7b77e0369ec8bbb62f6f1abf8046ad0187634a284d9aabfb9ce4b"
TERMINAL_SHA="6b3f7ea4981af8ec9546eafad4e5421765d8e674e4b64c4218194bddbed85e8a"
RECEIPT_SHA="3f59132aa9ef408b4b32bb6ecfa76ce94f21c120351990016545691cc6882684"
FROZEN_FCP_COMMIT="9fae6ccdf684a46026f72ba12e98de2c5c54bf2a"

def require(cond,msg):
    if not bool(cond):
        raise ValueError("MOIS_MASK_POSITION_QA5: "+msg)

def sha256(data):
    if isinstance(data,(str,Path)):
        data=Path(data).read_bytes()
    return hashlib.sha256(data).hexdigest()

def mask_geometry(mask):
    mask=np.asarray(mask,dtype=bool)
    require(mask.ndim==2,"mask dimension mismatch")
    yy,xx=np.nonzero(mask)
    require(len(xx)>0,"empty mask cannot be assigned flower position")
    h,w=mask.shape
    return {
      "mask_pixels":int(len(xx)),
      "bbox_xyxy":[int(xx.min()),int(yy.min()),int(xx.max()+1),int(yy.max()+1)],
      "center_xy":[float(xx.mean()),float(yy.mean())],
      "center_fraction":[float(xx.mean()/w),float(yy.mean()/h)],
      "fraction_of_frame":float(len(xx)/(w*h)),
      "image_width":int(w),"image_height":int(h)
    }

def require_exact_source_url(url,pid):
    p=urlsplit(str(url))
    require(p.scheme=="https" and p.netloc=="inaturalist-open-data.s3.amazonaws.com"
        and p.query=="" and p.fragment==""
        and re.fullmatch(r"/photos/"+re.escape(str(pid))+r"/large\.(?:jpg|jpeg|png)",p.path),
        "photo URL host, ID, resolution or suffix changed")
    return str(url)

def verify_replayed_counts(new_counts,old_terminal_row):
    names=set(new_counts)
    require(names==set(PALETTE12),"12 palette anchors are not complete")
    for c in PALETTE12:
        require(int(new_counts[c])==int(old_terminal_row["palette_count_"+c]),
                "original 12 palette counts differ: "+c)

def read_signed(source_zip,worker_zip):
    require(sha256(source_zip)==SOURCE_ZIP_SHA,"frozen 38 photo source artifact changed")
    require(sha256(worker_zip)==WORKER_ZIP_SHA,"frozen 38 result artifact changed")
    with ZipFile(source_zip) as z:
        s=z.read("shards/shard-000.csv")
        manifest=json.loads(z.read("manifest.json"))
    with ZipFile(worker_zip) as z:
        terminal=z.read("shard_000/sealed/terminal_p000_b0000000.csv")
        receipt=z.read("shard_000/receipt.json")
    require(sha256(s)==SHARD_SHA and sha256(terminal)==TERMINAL_SHA
            and sha256(receipt)==RECEIPT_SHA,"signed source/terminal SHA changed")
    require(manifest["photo_shard_sha256"]["0"]==SHARD_SHA,"38-photo manifest drift")
    worker_receipt=json.loads(receipt)
    require(worker_receipt["attempted"]==38
            and worker_receipt["source_sha256"]==SHARD_SHA
            and worker_receipt["terminal_sha256"]==TERMINAL_SHA
            and worker_receipt["method"]["source_commit"]==FROZEN_FCP_COMMIT,
            "original genuine ROI-v4 worker receipt disagreeing with source")
    source=pd.read_csv(io.BytesIO(s),dtype=str,keep_default_na=False)
    prior=pd.read_csv(io.BytesIO(terminal),dtype=str,keep_default_na=False)
    require(len(source)==len(prior)==38 and source.photo_id.is_unique
            and source.inat_taxon_id.is_unique and prior.measurement_id.is_unique,
            "nonunique 38 prior source or terminal")
    prior=prior.set_index("measurement_id")
    subset=source[source.photo_id.isin(TARGETS)].copy().sort_values("photo_id")
    require(len(subset)==len(TARGETS) and set(subset.photo_id)==set(TARGETS),
            "exact five frozen mismatch photos absent")
    for row in subset.itertuples(index=False):
        require(str(row.species)==TARGETS[row.photo_id],"photo/species source identity changed")
        require(row.photo_license=="cc-by-nc","five photo licences changed: "+row.photo_id)
        require_exact_source_url(row.photo_url_large,row.photo_id)
        require(row.measurement_id in prior.index,"missing matching measurement")
        require(prior.at[row.measurement_id,"measurement_status"]=="classified_four_state_morph",
                "a target was never strict numerical positive")
    return subset,prior,worker_receipt

def audit(*,source_zip,worker_zip,fcp_checkout,sam_dir,output):
    out=Path(output)
    require(not out.exists(),"cannot overwrite earlier frozen mask QA")
    source,prior,receipt=read_signed(source_zip,worker_zip)
    fcp=Path(fcp_checkout)
    # The shipped v4 frozen model files are hash-audited by FCP bootstrap.
    sys.path.insert(0,str(fcp))
    from fcp_pipeline.flower_roi_v4_runtime import FrozenFlowerColourEstimator,file_sha256
    from fcp_pipeline.photo_first_measurement import nearest_palette_counts
    from scripts.acquisition.run_fcp_fresh_image_partition import image_loader
    contract=json.loads((fcp/"docs/supporting/jbi_atlas_roi_estimator_contract_v4.json").read_text())
    execution=json.loads((fcp/"docs/supporting/random_photo_first_measurement_execution_v1.json").read_text())
    model=fcp/"data/atlas/qualification/roi_v4_training/jrc_yolo11n_last_v4.pt"
    require(file_sha256(model)==receipt["method"]["detector_sha256"],"detector SHA changed")
    estimator=FrozenFlowerColourEstimator(
        model,Path(sam_dir),contract,
        torch_threads=int(execution["runtime"]["torch_threads_per_worker"]))
    out.mkdir(parents=True)
    summaries=[]
    for row in source.itertuples(index=False):
        t=prior.loc[row.measurement_id]
        image_bytes=image_loader(row.photo_url_large)
        require(sha256(image_bytes)==t["image_sha256"],
                "source provider photo bytes changed, cannot compare old FCP result: "+row.photo_id)
        with Image.open(io.BytesIO(image_bytes)) as img:
            measured=estimator.measure(img)
            rgb=np.asarray(ImageOps.exif_transpose(img).convert("RGB"),dtype=np.uint8)
        mask=np.asarray(measured["flower_mask"],dtype=bool)
        require(mask.shape==rgb.shape[:2],"estimator masked another image orientation")
        geometry=mask_geometry(mask)
        require(geometry["mask_pixels"]==int(t["mask_pixels"])
                ==int(t["flower_effective_pixels"]),
                "replayed mask size changed: "+row.photo_id)
        require(measured["automated_colour_state_status"]=="automated_colour_state_admitted",
                "ROI-v4 original status no longer admitted")
        for k in ("horizontal_flip_mask_iou","horizontal_flip_colour_delta_e"):
            require(abs(float(measured[k])-float(t[k]))<1.0e-5,
                    "replayed quality metric differs for "+row.photo_id+":"+k)
        counts=nearest_palette_counts(rgb[mask])
        verify_replayed_counts(counts,t)
        mask_name="roi_mask_only_"+row.photo_id+".png"
        # Mask-only PNG does NOT include raw licensed source RGB. Authored QA derivative,
        # not a part of frozen numerical measurement or botanical expert confirmation.
        Image.fromarray(mask.astype(np.uint8)*255,mode="L").save(out/mask_name)
        boxes=measured["boxes"]
        result={
          "photo_id":str(row.photo_id),"species":str(row.species),
          "inat_taxon_id":str(row.inat_taxon_id),
          "original_observation_uuid":str(row.observation_id),
          "source_photo_page":"https://www.inaturalist.org/photos/"+str(row.photo_id),
          "source_url_large":str(row.photo_url_large),
          "photo_license":str(row.photo_license),
          "attribution":str(row.attribution),
          "source_pixel_sha256":str(t["image_sha256"]),
          "mask_png":mask_name,"mask_png_sha256":sha256(out/mask_name),
          "mask_geometry":geometry,
          "detector_boxes_xyxy_original_image":boxes,
          "detector_predictions":measured["predictions"],
          "retained_instances":int(measured["retained_instances"]),
          "flip_mask_iou":float(measured["horizontal_flip_mask_iou"]),
          "flip_colour_delta_e":float(measured["horizontal_flip_colour_delta_e"]),
          "true12_mask_palette_counts":{k:int(v) for k,v in counts.items()},
          "old_worker_mask_pixels":int(t["mask_pixels"]),
          "old_worker_measurement_status":str(t["measurement_status"]),
          "old_worker_morph":str(t["morph"]),
          "replay_matches_all_12_palette_counts":True,
          "replay_matches_mask_size_and_flip_metrics":True,
          "botanical_petal_overlap_directly_validated":False,
          "main_atlas_promoted":False,
        }
        summaries.append(result)
        print("MOIS_MASK_REPLAY_EQUIVALENT",row.photo_id,
              geometry["mask_pixels"],geometry["bbox_xyxy"],flush=True)
    pd.DataFrame([{
        "photo_id":r["photo_id"],"species":r["species"],
        "source_photo_sha256":r["source_pixel_sha256"],
        "mask_png":r["mask_png"],"mask_png_sha256":r["mask_png_sha256"],
        "mask_pixels":r["mask_geometry"]["mask_pixels"],
        "center_x":r["mask_geometry"]["center_xy"][0],
        "center_y":r["mask_geometry"]["center_xy"][1],
        "bbox_xyxy":json.dumps(r["mask_geometry"]["bbox_xyxy"]),
        "detector_box_count":len(r["detector_predictions"]),
        "detector_boxes_xyxy":json.dumps(r["detector_boxes_xyxy_original_image"]),
        "nuisance_fraction":sum(r["true12_mask_palette_counts"][k]
                                  for k in ("green","brown","black"))/
                            r["mask_geometry"]["mask_pixels"],
        "photo_license":r["photo_license"],"attribution":r["attribution"],
        "source_photo_page":r["source_photo_page"],
        "position_direct_plant_biology_review":"not_yet",
        "confidence_in_botanical_colour":"unverified",
    } for r in summaries]).to_csv(out/"five_photo_mask_geometry.csv",index=False)
    (out/"five_masks_provenance.json").write_text(
        json.dumps(summaries,sort_keys=True,indent=2,ensure_ascii=False)+"\n")
    output_sha={p.name:sha256(p) for p in out.iterdir()}
    audit_receipt={
      "schema":"mois_visual5_frozen_v4_mask_capture_qa_v1",
      "source_original_source_zip_sha256":SOURCE_ZIP_SHA,
      "source_original_worker_zip_sha256":WORKER_ZIP_SHA,
      "source_original_terminal_sha256":TERMINAL_SHA,
      "source_frozen_shard_sha256":SHARD_SHA,
      "source_FCP_method_commit":FROZEN_FCP_COMMIT,
      "detector_sha256":receipt["method"]["detector_sha256"],
      "encoder_sha256":receipt["method"]["efficient_sam_encoder_sha256"],
      "decoder_sha256":receipt["method"]["efficient_sam_decoder_sha256"],
      "replayed_photo_count":len(summaries),
      "mask_byte_exact_palette_replays":len(summaries),
      "mask_size_and_flip_replays":len(summaries),
      "new_independent_original_photo_ids":0,
      "new_nine_palette_numeric_measurements":0,
      "photographic_colour_labels_changed":0,
      "main_catalogue_changes":0,
      "raw_image_bytes_in_artifact":False,
      "mask_only_QA_derivatives":True,
      "botanical_petal_truth_verification_completed":False,
      "output_sha256":output_sha
    }
    (out/"receipt.json").write_text(json.dumps(audit_receipt,indent=2,sort_keys=True)+"\n")
    print("MOIS_FROZEN_VISUAL5_MASK_CAPTURE_REPLAY_SUCCESS",
          json.dumps(audit_receipt,sort_keys=True),flush=True)
    return audit_receipt

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("source-zip","worker-zip","fcp-checkout","sam-dir","output"):
        p.add_argument("--"+k,required=True)
    a=p.parse_args()
    audit(source_zip=a.source_zip,worker_zip=a.worker_zip,
          fcp_checkout=a.fcp_checkout,sam_dir=a.sam_dir,output=a.output)
