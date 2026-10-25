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
from pathlib import Path
import pandas as pd

HUES=("white","yellow_orange","red_pink","blue_purple")
SCHEMA="mois_full_source_observation_sampling_cohort_v1"

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
        context=PrefetchLoader(urls,image_loader,out/"prefetch",
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


def source_prior_dedup_receipt(info):
    """Reconcile old photo-ID-only and new photo/observation dedup metadata.

    No historical photo is reread or counted as newly classified. The routed
    tier06 source froze a union of previous-photo / previous-observation
    exclusions and independent technical benchmark exclusions before pixels.
    Only *that* known signed input layout is permitted as a compatibility path.
    """
    # The 2026-10-24 full K1 75-group model run completed all 76,753
    # ORIGINAL photo terminals, but the old collector only recognized the
    # unrelated tier06 'bucket0' source layout.  This explicit, SHA-bound
    # compatibility path preserves genuine pre-colour K1 prior photo reuse
    # and observation/technical-holdout design without synthesizing counts.
    if info.get("source_group")=="balanced_taxon_month_lat_K1":
        from scripts.acquisition.mois_fcp_k1_append_completed_original_photo_group import (
            PARENT_SOURCE_SHA, SOURCE_ARM, SOURCE_PLANNED,
            SELECTED_SOURCE_OBSERVATIONS, REUSED_OLD_ORIGINALS,
        )
        previous=info.get("prior_original_crossarm_photo_overlap")
        technical=info.get("technical_holdout_images_not_redownloaded")
        unsupported=info.get("unsupported_original_formats_technical_unknown")
        require(
            info.get("schema")==SCHEMA
            and info.get("source_tier")==6
            and info.get("source_arm")==SOURCE_ARM
            and info.get("full_frame_photo_ids")==93_209_351
            and info.get("source_precolour_balanced_manifest_SHA256")==PARENT_SOURCE_SHA
            and info.get("new_photo_attempts_planned")==SOURCE_PLANNED
            and info.get("all_balanced_selected_original_observations")==
                SELECTED_SOURCE_OBSERVATIONS
            and info.get("reused_original_photo_ids")==REUSED_OLD_ORIGINALS
            and type(technical) is int and technical>=0
            and type(unsupported) is int and unsupported>=0
            and SOURCE_PLANNED+REUSED_OLD_ORIGINALS+technical+unsupported==
                SELECTED_SOURCE_OBSERVATIONS
            and isinstance(previous,dict)
            and previous.get("unique_old_photo_IDs")==929_138
            and previous.get("signed_arm_photo_rows")==929_138
            and isinstance(info.get("prior_original_pixel_arms"),list)
            and len(info["prior_original_pixel_arms"])==4
            and info.get("selection_probability") is None
            and info.get("sampling_nested_joint_probability_per_photo") is True
            and info.get("original_photo_source_probabilities_unchanged_by_prior_pixel_reuse") is True
            and info.get("known_previous_full_photo_and_observation_IDs_source_SHA_verified") is True
            and info.get("previously_counted_photo_ID_duplicates_are_NOT_new_inferences") is True
            and info.get("all_original_photo_ids_in_master_probability_frame_preserved") is True
            and info.get("technical_unknown_not_other_colour") is True
            and info.get("colour_measured") is False
            and info.get("images_downloaded_by_this_source_bridge")==0
            and info.get("new_colour_ROIs_inferred_by_this_source_bridge")==0
            and info.get("photo_specific_anthesis_native_wild_not_independently_reviewed") is True,
            "K1 source photo-ID-only vs reused-observation joint-selection denominator changed")
        return {
            "matching_prior_measured_photos_reused":REUSED_OLD_ORIGINALS,
            "source_deduplication_reporting_grain":
                "exact_original_photo_ID_reuse; source_observation_sampling_unchanged",
            "independent_original_observation_photo_targets":
                SELECTED_SOURCE_OBSERVATIONS,
            "new_original_photo_IDs_measured":SOURCE_PLANNED,
            "prior_original_photo_IDs_reused_without_new_inference":
                REUSED_OLD_ORIGINALS,
            "technical_photo_ID_holdouts_not_redownloaded":technical,
            "unsupported_original_image_formats_excluded":unsupported,
            "historical_unique_source_photo_IDs_unchanged":929_138,
            "K1_original_joint_photo_selection_probability_per_individual_photo":True,
            "historical_source_photo_ID_disjointness_preverified":True,
            "botanical_flower_native_status_unreviewed":True,
            "source_precolour_photo_dispositions_reconciled":True,
        }
    existing=info.get("previously_measured_photo_ids_skipped")
    if existing is not None:
        require(type(existing) is int and existing>=0,
                "Invalid original photo-ID skip count in source manifest")
        return {
            "matching_prior_measured_photos_reused":existing,
            "source_deduplication_reporting_grain":"prior_original_photo_id",
        }
    require(info.get("schema")==SCHEMA and info.get("source_tier")==6
            and info.get("source_group")==0
            and info.get("compute_bucket_count")==512
            and info.get("admitted_bucket_count")==64
            and info.get("not_probability_complete_global_source") is True
            and info.get("previous_source_photo_IDs_not_redownloaded") is True,
            "Unknown source metadata layout: cannot fabricate prior photo dedup counts")
    exclusions=info.get("prior_and_technical_excluded_source_photo_IDs")
    require(isinstance(exclusions,dict)
            and set(exclusions)<={
                "prior_actual_photo_or_observation",
                "prior_technical_photo_or_observation",
                "unsupported_loader_format"}
            and all(type(x) is int and x>=0 for x in exclusions.values()),
            "Invalid original source photo/observation exclusion inventory")
    prior=exclusions.get("prior_actual_photo_or_observation",0)
    technical=exclusions.get("prior_technical_photo_or_observation",0)
    unsupported=exclusions.get("unsupported_loader_format",0)
    require(unsupported==info.get("selected_new_photos_unsupported_by_frozen_loader")
            and info.get("selected_bucket_photo_IDs_before_dedup")==
                    info.get("new_photo_attempts_planned")+sum(exclusions.values()),
            "Routed source precolour prior/technical/unsupported count mismatch")
    return {
        # Original field retained for backward-compatible consumers; these
        # candidate rows can be observation-level exclusions as well as exact
        # previous-photo matches. Preserve the precise grain below.
        "matching_prior_measured_photos_reused":prior,
        "source_deduplication_reporting_grain":"prior_original_photo_or_observation",
        "prior_photo_or_observation_source_rows_excluded":prior,
        "independent_technical_photo_or_observation_rows_excluded":technical,
        "unsupported_original_image_formats_excluded":unsupported,
        "source_precolour_photo_dispositions_reconciled":True,
    }

def collect(source,workers,out):
    source=Path(source);workers=Path(workers);out=Path(out)
    require(not out.exists() or not any(out.iterdir()),
            "Refuse to replace signed RGB source result")
    info=plan(source)
    require(info["schema"]==SCHEMA,"Invalid complete-source cohort")
    skip_summary=source_prior_dedup_receipt(info)
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
        **skip_summary,
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
