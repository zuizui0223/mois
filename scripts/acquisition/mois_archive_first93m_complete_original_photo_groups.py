#!/usr/bin/env python3
"""Preserve first four ALREADY COMPLETED original-pixel groups in the full-93m wave.

This source-audited 4096-photo acquisition checkpoint is a proper SUBSET of the
18,462-new-photo original observation-probability p=.0002 wave. It is not a
new random draw, independent replication or the completed full cohort. Keep
the original source metadata and technical missingness unchanged.

No original S3 images or image pixels are fetched, stored, or redistributed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

HUES=("white","yellow_orange","red_pink","blue_purple")
FIRST_INDEX_SHA="22f7b6ddd9770e5ef51a093e547e99a25cf704b98de7c4e102d6b1bd99d6ebd9"
SOURCE_PHOTOS=93209351
DRAW_ATTEMPTED=18533
DRAW_NEW=18462
PREVIOUS_REUSED=71
SHARD_START=0
SHARD_STOP=16
EXPECTED_FROZEN_FCP={
    "detector_sha256":"f1aaeec4664fe2c178e5cf2bc1f508977bef3e4aa7b40613026cb8ae3de789d5",
    "efficient_sam_decoder_sha256":"a62f8fa5ea080447c0689418d69e58f1e83e0b7adf9c142e2bd9bcc8045c0b11",
    "efficient_sam_encoder_sha256":"84ed466ffcc5c1f8d08409bc34a23bb364ab2c15e402cb12d4335a42be0e0951",
    "measurement_git_blob":"b1d7091db44682559bbe4b48a9f3b8db648501a6",
    "source_commit":"9fae6ccdf684a46026f72ba12e98de2c5c54bf2a",
    "location_blind":True,
    "no_old_measurement_csv_inputs":True,
    "persisted_image_pixels":False,
}
OUT_PHOTOS="mois93m_first4096_completed_original_pixel_photo_source.parquet"
OUT_SOURCE="mois93m_first4096_photo_acquisition_checkpoint_receipt.json"
OUT_COUNTS="mois93m_first4096_acquired_original_hue_technical_counts.csv"


def require(ok,msg):
    if not ok:
        raise ValueError(msg)


def digest(file):
    h=hashlib.sha256()
    with Path(file).open("rb") as input:
        for chunk in iter(lambda:input.read(1<<20),b""):
            h.update(chunk)
    return h.hexdigest()


def frame(path):
    return pd.read_csv(path,dtype=str,keep_default_na=False)


def validate(queue, groups, *, first=SHARD_START,stop=SHARD_STOP):
    queue=Path(queue)
    groups=Path(groups)
    mfile=queue/"manifest.json"
    require(mfile.is_file(),"Original pre-hue p=.0002 source manifest missing")
    m=json.loads(mfile.read_text())
    require(m["schema"]=="mois_full_source_observation_sampling_cohort_v1" and
            m["full_frame_photo_ids"]==SOURCE_PHOTOS and
            m["full_source_manifest_sha256"]==FIRST_INDEX_SHA and
            m["selection_probability"]==0.0002 and
            m["selection_unit"]=="observation_uuid" and
            m["sampled_before_any_hue_measurement"] is True and
            m["source_original_images_measured"] is False and
            m["colour_measured"] is False and
            m["sampled_photo_ids_before_prior_dedup"]==DRAW_ATTEMPTED and
            m["new_photo_attempts_planned"]==DRAW_NEW and
            m["previously_measured_photo_ids_skipped"]==PREVIOUS_REUSED and
            m["selected_new_photos_unsupported_by_frozen_loader"]==0 and
            m["work_units"]==73 and m["units_per_group"]==4 and
            len(m["photo_shard_sha256"])==73 and
            m["no_per_taxon_or_observation_photo_count_cap"] is True,
            "Not the complete source-indexed first p0002 blinded original photograph roster")
    require(first>=0 and stop<=m["work_units"] and first<stop and
            first%4==0 and stop%4==0,
            "Only completed whole 4-photo-unit groups may be sealed")
    batch=[]
    unit_receipts=[]
    for n in range(first,stop):
        source_file=queue/"shards"/f"shard-{n:03d}.csv"
        raw=groups/f"group{n//4:02d}"/f"shard_{n:03d}"
        outcome_file=raw/"photos.csv"
        receipt_file=raw/"receipt.json"
        terminal_file=raw/"sealed"/f"terminal_p{n:03d}_b0000000.csv"
        require(source_file.is_file() and outcome_file.is_file() and
                receipt_file.is_file() and terminal_file.is_file(),
                "Completed actual photo/image work unit source or terminal missing: "+str(n))
        rec=json.loads(receipt_file.read_text())
        require(digest(source_file)==m["photo_shard_sha256"][str(n)]
                and rec["schema"]=="mois_full_original93m_real_RGB_photo_unit_v1"
                and rec["shard"]==n and rec["photo_colours_measured"] is True
                and rec["raw_licensed_image_bytes_persisted"] is False
                and rec["source_sha256"]==m["photo_shard_sha256"][str(n)]
                and rec["photos_sha256"]==digest(outcome_file)
                and rec["terminal_sha256"]==digest(terminal_file),
                "Original photos/signed model terminal or source shard SHA changed: "+str(n))
        require(rec["method"]==EXPECTED_FROZEN_FCP,
                "Original YOLO11/EfficientSAM fixed model source digest or threshold drifted")
        selected=frame(source_file)
        photographed=frame(outcome_file)
        require(0<len(photographed)<=256 and
                len(selected)==len(photographed)==rec["attempted"] and
                photographed.photo_id.is_unique and
                photographed.measurement_id.is_unique and
                photographed.photo_id.tolist()==selected.photo_id.tolist() and
                photographed.measurement_id.tolist()==selected.measurement_id.tolist(),
                "Authenticated source-selected original photo ID is missing/changed: "+str(n))
        for col in selected.columns:
            require(col in photographed and
                    photographed[col].tolist()==selected[col].tolist(),
                    "Source observation/taxon/date/photograph metadata drifted: "+col)
        classified=photographed.morph.isin(HUES)
        failed=photographed.measurement_status.eq("image_acquisition_failed")
        require(classified.eq(photographed.measurement_status.eq(
                    "classified_four_state_morph")).all()
                and photographed.classified_four_hue.str.lower().eq(
                    "true").eq(classified).all()
                and int(classified.sum())==rec["classifiable"]
                and int(failed.sum())==rec["acquisition_failed"]
                and int((~classified & ~failed).sum())==rec["unclassifiable_acquired"]
                and int((~failed).sum())==rec["downloaded"],
                "Actual four-colour photo ROI/acquisition failures cannot be imputed: "+str(n))
        unit_receipts.append({"shard":n,
                              "source_shard_SHA256":digest(source_file),
                              "sealed_result_SHA256":digest(outcome_file),
                              "terminal_sha256":digest(terminal_file),
                              "attempted":len(photographed),
                              "classifiable":int(classified.sum()),
                              "acquisition_failed":int(failed.sum())})
        batch.append(photographed)
    united=pd.concat(batch,ignore_index=True)
    require(united.photo_id.is_unique and united.measurement_id.is_unique,
            "Duplicate original photo IDs between genuine workers")
    require(united.source_arm.eq("global93m_probability_sample_observations").all()
            and united.selection_probability.eq("0.0002").all(),
            "Original pre-colour probability cohort or recruiting source arm was rewritten")
    outcome=united.morph.isin(HUES)
    failed=united.measurement_status.eq("image_acquisition_failed")
    require(len(united)==(stop-first)*256,
            "Not all completed first p0002 original photographed image attempts present")
    stats={
        "schema":"mois_full93m_first4096_true_RGB_partial_checkpoint_v1",
        "record_status":"AUTHENTIC_COMPLETED_FIRST_4_OF_19_GROUPS_PARTIAL_NOT_FULL_18462",
        "original_full_upstream_photo_candidates":SOURCE_PHOTOS,
        "full_source_immutable_manifest_SHA256":FIRST_INDEX_SHA,
        "sampled_source_photos_total_original":DRAW_ATTEMPTED,
        "new_source_original_photo_attempts_planned":DRAW_NEW,
        "prior_pixel_photo_ids_skipped_at_draw":PREVIOUS_REUSED,
        "source_plan_SHA256":digest(mfile),
        "preserved_first_group_count":(stop-first)//4,
        "complete_roi_work_units_in_this_partial":stop-first,
        "first_photo_unit_included":first,
        "last_photo_unit_included_inclusive":stop-1,
        "genuine_original_photo_attempts_in_partial":len(united),
        "photo_IDs_unique":True,
        "source_original_observation_ids":int(
            united.observation_id[united.observation_id.ne("")].nunique()),
        "source_nominal_taxon_ids":int(united.inat_taxon_id.nunique()),
        "downloaded_original_images":int((~failed).sum()),
        "classified_four_hues":int(outcome.sum()),
        "technical_unknown_after_download":int((~outcome&~failed).sum()),
        "verified_original_photo_acquisition_failures":int(failed.sum()),
        "hue_counts_classified":{
            k:int(united.morph.eq(k).sum()) for k in HUES},
        "photo_equal_white_fraction_among_classified_only":(
            float(united.loc[outcome,"morph"].eq("white").mean()) if outcome.any() else None),
        "original_source_scientific_model_fingerprint":EXPECTED_FROZEN_FCP,
        "original_unit_receipts":unit_receipts,
        "partial_photo_outcomes_only_not_an_unbiased_floral_census":True,
        "completion_of_remaining_14366_new_photos_NOT_claimed":True,
        "H3_native_island_or_any_abiotic_ecological_model_performed":False,
        "open_flower_petals_independently_botanically_confirmed":False,
        "raw_licensed_image_pixels_redistributed":False,
    }
    return united,stats


def run(queue,groups,output, *, first=0,stop=16):
    output=Path(output)
    require(not output.exists() or not any(output.iterdir()),
            "Do not overwrite an immutable original-photo acquisition partial result")
    united,stats=validate(queue,groups,first=first,stop=stop)
    output.mkdir(parents=True,exist_ok=True)
    united.to_parquet(output/OUT_PHOTOS,compression="zstd",index=False)
    rows=[]
    for k in HUES:
        rows.append({"status":"classifiable_four_visible_hue","morph":k,
                     "photo_attempts":int(united.morph.eq(k).sum())})
    rows += [
        {"status":"technical_unclassifiable_acquired","morph":"",
         "photo_attempts":stats["technical_unknown_after_download"]},
        {"status":"original_image_acquisition_failed","morph":"",
         "photo_attempts":stats["verified_original_photo_acquisition_failures"]}
    ]
    pd.DataFrame(rows).to_csv(output/OUT_COUNTS,index=False)
    stats["output_SHA256"]={
        OUT_PHOTOS:digest(output/OUT_PHOTOS),
        OUT_COUNTS:digest(output/OUT_COUNTS),
    }
    (output/OUT_SOURCE).write_text(json.dumps(stats,indent=2,sort_keys=True,
                                             allow_nan=False)+"\n")
    print("MOIS_FIRST93M_TRUE_RGB_FIRST4096_SOURCE_SIGNED_AND_DURABLE",
          json.dumps({k:stats[k] for k in (
              "genuine_original_photo_attempts_in_partial",
              "classified_four_hues","source_nominal_taxon_ids",
              "source_original_observation_ids",
              "downloaded_original_images",
              "verified_original_photo_acquisition_failures",
              "technical_unknown_after_download","hue_counts_classified")},
              sort_keys=True),flush=True)
    return stats


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source",type=Path,required=True)
    p.add_argument("--groups",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--first",type=int,default=0)
    p.add_argument("--stop",type=int,default=16)
    x=p.parse_args()
    run(x.source,x.groups,x.output,first=x.first,stop=x.stop)


if __name__=="__main__":
    main()
