#!/usr/bin/env python3
"""FCP-frozen authentic original image ROI measurements for global 93m source.

Each model-bootstrap worker processes four or fewer 256-photo units with
one bounded paced prefetch queue. Prior completed unit receipts must be
reused, never downloaded again. No original copyrighted image bytes retained.
Unknown/photo-fetch failures are NOT visible-colour data.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from contextlib import nullcontext
from urllib.parse import urlsplit
import re
from pathlib import Path
import pandas as pd

HUES=("white","yellow_orange","red_pink","blue_purple")
SCHEMA="mois_full_source_observation_sampling_cohort_v1"
CANONICAL_ORIGINAL_URL=re.compile(r"^/photos/[0-9]+/large\.([a-z0-9]{1,10})$", re.I)
ORIGINAL_FROZEN_IMAGE_EXTENSIONS=frozenset(("jpg","jpeg","png"))


def original_source_loader_supported(url):
    """One pure URL preflight: never let non-FCP formats stop an RGB batch.

    All selected PHOTO IDs stay in the source queue. An unsupported URI is
    deterministically classified as a technical original-image acquisition
    failure by the wrapper below and no network request is made. The old
    signed method's three allowed extensions are not silently expanded.
    """
    u=urlsplit(str(url))
    if (u.scheme!="https" or u.netloc!="inaturalist-open-data.s3.amazonaws.com"
            or u.query or u.fragment):
        return False
    m=CANONICAL_ORIGINAL_URL.fullmatch(u.path)
    return bool(m and m.group(1).lower() in ORIGINAL_FROZEN_IMAGE_EXTENSIONS)


def original_source_image_loader(url,original_loader,permanent_error):
    if not original_source_loader_supported(url):
        raise permanent_error("Original indexed source image extension is not supported by frozen ROI-v4 photo loader; exact original source photo retained as technical unknown")
    return original_loader(url)

def require(ok,msg):
    if not ok:raise ValueError(msg)
def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for x in iter(lambda:f.read(1<<20),b""):h.update(x)
    return h.hexdigest()
def plan(source):
    return json.loads((Path(source)/"manifest.json").read_text())
def input_rows(source,index):
    info=plan(source)
    require(info["schema"]==SCHEMA and not info["colour_measured"]
            and index in info["worker_ids"],"Unknown, measured or altered source original photo work unit")
    f=Path(source)/"shards"/f"shard-{index:03d}.csv"
    require(sha(f)==info["photo_shard_sha256"][str(index)],
            "Original photo ID roster SHA changed")
    rows=pd.read_csv(f,dtype=str,keep_default_na=False)
    require(0<len(rows)<=256 and rows.photo_id.is_unique and
            rows.measurement_id.is_unique and rows.shard.eq(str(index)).all(),
            "Original signed source photos not unique and complete")
    return rows,info

class PhotoBatch:
    def __init__(self,rows,part):
        self.rows=rows
        self.part=part
    def get_batch(self,part,batch,size):
        require(part==self.part and batch==0 and size==len(self.rows),
                "Original FCP ROI work unit number and size changed")
        a=self.rows[["measurement_id","image_filename","photo_license"]].copy()
        b=self.rows[["measurement_id","image_filename","photo_license","photo_url_large"]].copy()
        return a,b,1

def measure_group(source,group,runtime,sam,out,*,runner=None,loader=None):
    source=Path(source);out=Path(out)
    info=plan(source)
    require(info["schema"]==SCHEMA,"Source receipt is not complete full-ID sampling design")
    units=list(range(group*4,min((group+1)*4,info["work_units"])))
    require(bool(units),"Requested a group outside the original source")
    done=[];pending=[]
    for index in units:
        rows,signed=input_rows(source,index)
        folder=out/f"shard_{index:03d}"
        receipt=folder/"receipt.json"
        file=folder/"photos.csv"
        if receipt.exists():
            d=json.loads(receipt.read_text())
            require(d["source_sha256"]==signed["photo_shard_sha256"][str(index)]
                    and d["photos_sha256"]==sha(file)
                    and d["attempted"]==len(rows)
                    and d["shard"]==index,
                    "Previously measured photo shard is corrupted, never silently re-infer")
            done.append(d)
        else:
            pending.append((index,rows,signed,folder))
    if not pending:return done
    production=runner is None
    if production:
        from scripts.acquisition.run_fcp_fresh_image_partition import (
            run_batch,image_loader,PermanentOriginalPhotoError)
        from scripts.acquisition.mois_prefetch import PrefetchLoader
        from scripts.analysis.mois_reassemble_original_new4096_measurements import METHOD
        runner=run_batch
        urls=[u for _,rows,_,_ in pending for u in rows.photo_url_large.tolist()]
        def loader(url):
            return original_source_image_loader(
                url,image_loader,PermanentOriginalPhotoError)
        context=PrefetchLoader(urls,loader,out/"prefetch",
            permanent_errors=(PermanentOriginalPhotoError,))
    else: context=nullcontext(loader)
    with context as prefetch:
        for index,rows,signed,folder in pending:
            folder.mkdir(parents=True,exist_ok=True)
            sealed=folder/"sealed"
            kwargs={
                "source":Path(runtime)/"scripts/analysis/measure_random_photo_first_blind_partition.py",
                "detector":Path(runtime)/"data/atlas/qualification/roi_v4_training/jrc_yolo11n_last_v4.pt",
                "sam":Path(sam),"indexed_source":PhotoBatch(rows,index)
            }
            if prefetch is not None:kwargs["loader"]=prefetch
            result=runner(Path("."),index,0,len(rows),sealed,**kwargs)
            require(result["status"]=="COMPLETE_WORK_UNIT"
                    and result["n_source_photo_ids"]==len(rows)
                    and result["n_terminal_photo_ids"]==len(rows)
                    and result["number_batches_for_this_partition"]==1
                    and not result["image_bytes_persisted"]
                    and not result["reused_historical_measurement_csv"],
                    "Genuine full original RGB/ROI source terminal did not complete")
            if production:
                require(result["method"]==METHOD,
                        "Original fixed YOLO/EfficientSAM model was changed")
            terminal=sealed/f"terminal_p{index:03d}_b0000000.csv"
            require(sha(terminal)==result["terminal_csv_sha256"],
                    "Source-original pixel colour outcome terminal SHA mismatch")
            measured=pd.read_csv(terminal,dtype=str,keep_default_na=False)
            require(len(measured)==len(rows) and measured.measurement_id.is_unique
                    and set(measured.measurement_id)==set(rows.measurement_id),
                    "Actual original image ROI output missing source IDs")
            cols=["measurement_id","morph","measurement_status","image_sha256"]
            cols += [c for c in measured.columns if c.startswith("flower_fraction_")]
            frame=rows.merge(measured[cols],on="measurement_id",
                             how="inner",validate="one_to_one",sort=False)
            known=frame.morph.isin(HUES)
            failed=frame.measurement_status.eq("image_acquisition_failed")
            require(frame.morph.isin((*HUES,"mixed_uncertain")).all()
                    and known.eq(frame.measurement_status.eq("classified_four_state_morph")).all()
                    and int(known.sum())==result["n_classifiable"]
                    and int(failed.sum())==result["n_acquisition_failed"],
                    "Image acquisition failure cannot be called a coloured flower")
            frame["classified_four_hue"]=known.map({True:"true",False:"false"})
            frame.to_csv(folder/"photos.csv",index=False,lineterminator="\n")
            summary={
                "schema":"mois_full_original93m_real_RGB_photo_unit_v1",
                "shard":index,"attempted":len(rows),
                "downloaded":int((~failed).sum()),
                "classifiable":int(known.sum()),
                "acquisition_failed":int(failed.sum()),
                "unclassifiable_acquired":int((~known & ~failed).sum()),
                "source_sha256":signed["photo_shard_sha256"][str(index)],
                "method":result["method"],
                "terminal_sha256":sha(terminal),
                "photos_sha256":sha(folder/"photos.csv"),
                "photo_colours_measured":True,
                "raw_licensed_image_bytes_persisted":False,
            }
            (folder/"receipt.json").write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n")
            done.append(summary)
    return done

def collect(source,workers,out):
    source=Path(source);workers=Path(workers);out=Path(out)
    require(not out.exists() or not any(out.iterdir()),
            "Refuse to replace signed RGB source result")
    info=plan(source)
    require(info["schema"]==SCHEMA,"Invalid complete-source cohort")
    frames=[];receipts=[]
    for index in info["worker_ids"]:
        names=list(workers.rglob(f"shard_{index:03d}/receipt.json"))
        require(len(names)==1,f"Missing/duplicate real original photo ROI work unit {index}")
        r=json.loads(names[0].read_text());p=names[0].parent/"photos.csv"
        require(r["shard"]==index and
                r["source_sha256"]==info["photo_shard_sha256"][str(index)]
                and r["photos_sha256"]==sha(p),
                "Original measured photo identity or result SHA changed")
        frame=pd.read_csv(p,dtype=str,keep_default_na=False)
        require(len(frame)==r["attempted"] and frame.photo_id.is_unique
                and int(frame.morph.isin(HUES).sum())==r["classifiable"],
                "Original ROI photo classifications mismatch")
        frames.append(frame);receipts.append(r)
    require(len({json.dumps(r["method"],sort_keys=True) for r in receipts})==1,
            "Mixed original hue classifiers across photo workers")
    combined=pd.concat(frames,ignore_index=True)
    require(len(combined)==info["new_photo_attempts_planned"]
            and combined.photo_id.is_unique,
            "New global source photo IDs duplicated or omitted")
    out.mkdir(parents=True,exist_ok=True)
    combined.to_parquet(out/"measured_photos.parquet",
                        compression="zstd",index=False)
    counters={key:sum(r[key] for r in receipts) for key in
              ("attempted","downloaded","classifiable",
               "acquisition_failed","unclassifiable_acquired")}
    receipt={
        "schema":"mois_global93m_observation_baseline_real_RGB_photo_results_v1",
        **counters,
        "full_original_upstream_photo_candidates":info["full_frame_photo_ids"],
        "observation_level_selection_probability":info["selection_probability"],
        "source_plan_sha256":sha(source/"manifest.json"),
        "measured_photo_parquet_sha256":sha(out/"measured_photos.parquet"),
        "matching_prior_measured_photos_reused":info["previously_measured_photo_ids_skipped"],
        "frozen_original_ROI_method":receipts[0]["method"],
        "technical_unknowns_not_assumed_nonwhite":True,
        "ecological_analysis_run":False,
        "raw_licensed_photo_pixels_redistributed":False,
    }
    (out/"receipt.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print("MOIS_FULL_SOURCE_OBSERVATION_BASELINE_TRUE_RGB_COMPLETE",
          json.dumps(counters,sort_keys=True),flush=True)
    return receipt

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    sub=p.add_subparsers(dest="command",required=True)
    a=sub.add_parser("work")
    for k in ("source","runtime","sam","output"):a.add_argument("--"+k,type=Path,required=True)
    a.add_argument("--group",type=int,required=True)
    a=sub.add_parser("collect")
    for k in ("source","workers","output"):a.add_argument("--"+k,type=Path,required=True)
    x=p.parse_args()
    if x.command=="work": measure_group(x.source,x.group,x.runtime,x.sam,x.output)
    else: collect(x.source,x.workers,x.output)
