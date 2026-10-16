#!/usr/bin/env python3
"""Compare true original photo ROI scores with LIVE observation phenology metadata.

NEVER use an observation-level iNaturalist annotation as independent proof that
the pictured original ROI shows an open petal. Mixed same-ROI flower organs and
failed image classification remain explicit photo states. Climate values, dates,
seasons and live observation annotations NEVER enter blind original ROI worker.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import zipfile
from collections import Counter
from pathlib import Path

ORIGINAL="postblind/first_original_source_matched_season_64_photo_roi_QC.csv"
ORIGINAL_RECEIPT="postblind/first_original_source_matched_season_64_roi_receipt.json"
LIVE_PREFIX="mois_2025_live_plant_phenology_supplement/"
LIVE_CSV=LIVE_PREFIX+"original_2025_64_live_iNaturalist_observation_plant_phenology_NOT_PHOTO.csv"
LIVE_RECEIPT=LIVE_PREFIX+"2025_Asclepias_live_public_observation_phenology_receipt.json"
FCP_PIXEL_SHA="1b16cdabd2936567dee87b8344019f40f315600c6510dfc1f25a03835024b6a0"
COLS=("original_photo_season","live_iNat_observation_phenology_NOT_photo",
      "automatic_original_ROI_label_NOT_botanical_anthesis","n_original_photo_observations",
      "all_photo_open_petal_tepal_unreviewed")
LABELS={"white","yellow_orange","red_pink","blue_purple"}


def sha(raw:bytes)->str:
    return hashlib.sha256(raw).hexdigest()


def csv_read(raw:bytes)->list[dict]:
    with io.StringIO(raw.decode("utf-8")) as f:
        reader=csv.DictReader(f)
        if not reader.fieldnames:
            raise ValueError("SHA-bound FCP photo or MOIS annotation table has no header")
        return list(reader)


def roi_label(row:dict)->str:
    good=row["measurement_status"]=="classified_four_state_morph"
    morph=row["morph"]
    mixed=row["same_roi_white_yellow_organ_mixed"]
    if mixed not in ("True","False"):
        raise ValueError("Original same-ROI mixed organ flag changed")
    if good and morph not in LABELS:
        raise ValueError("Original four-state ROI morph not valid")
    if not good:
        if mixed=="True":
            raise ValueError("Failed original ROI cannot assert mixed classified flower organs")
        return "ROI_failed_or_unclassifiable_UNKNOWN_flower_stage"
    if mixed=="True":
        return "same_ROI_white_yellow_floral_organs_mixed_not_morph"
    return ("white" if morph=="white" else "chromatic")+"_automatic_ROI_UNREVIEWED"


def annotation_label(source:dict)->str:
    raw=source["plant_phenology_observation_annotation_status_NOT_PHOTO"]
    known={
        "current_observation_annotation_flowering_NOT_photo",
        "current_observation_annotation_fruiting_NOT_photo",
        "current_observation_annotation_flower_budding_NOT_photo",
        "current_observation_annotation_no_evidence_of_flowering_NOT_photo",
        "no_current_observation_phenology_field_unknown",
        "no_current_plant_phenology_annotation_unknown",
        "conflicting_current_observation_phenology_NOT_photo",
        "unrecognized_current_observation_phenology_NOT_photo",
        "no_current_public_observation_unknown",
        "live_observation_date_changed_since_source_annotation_UNUSABLE",
        "original_source_photo_not_associated_with_live_observation_unknown",
    }
    if raw not in known:
        raise ValueError("Current public iNat plant-phenology status unexpectedly recoded")
    if source["photo_open_flower_petal_tepal_independently_verified"]!="False":
        raise ValueError("Live observation annotation was falsely called botanical photo proof")
    return raw


def cross(fcp_archive:Path,annotation_archive:Path,out:Path)->dict:
    if out.exists() and any(out.iterdir()):
        raise ValueError("Never overwrite archived current observational annotation source")
    with zipfile.ZipFile(fcp_archive) as z:
        photo_bytes=z.read(ORIGINAL)
        photo_receipt=json.loads(z.read(ORIGINAL_RECEIPT))
    with zipfile.ZipFile(annotation_archive) as z:
        annotation_bytes=z.read(LIVE_CSV)
        annotation_receipt=json.loads(z.read(LIVE_RECEIPT))
    if (sha(photo_bytes)!=FCP_PIXEL_SHA or
        photo_receipt.get("status")!="FIRST_ORIGINAL_SEASONAL_PIXELS_MEASURED_BOTANICAL_ANTHESIS_NOT_VALIDATED" or
        photo_receipt.get("postblind_join_sha256")!=FCP_PIXEL_SHA or
        photo_receipt.get("original_photos_measured")!=64 or
        annotation_receipt.get("status")!="CURRENT_OBSERVATION_PHENOLOGY_METADATA_ONLY_NO_PHOTO_ANTHESIS_CONFIRMED" or
        annotation_receipt.get("frozen_original_photo_receipt_sha256")!=FCP_PIXEL_SHA or
        annotation_receipt.get("original_photo_observation_live_annotation_table_sha256")!=sha(annotation_bytes) or
        annotation_receipt.get("n_original_source_photos")!=64 or
        annotation_receipt.get("n_live_v1_observation_records")!=64 or
        annotation_receipt.get("n_live_v2_uuid_resolved")!=64 or
        annotation_receipt.get("actual_photo_open_flower_petal_tepal_verified_count")!=0 or
        annotation_receipt.get("observation_annotation_does_not_certify_specific_photo_open_flower") is not True or
        annotation_receipt.get("source_original_photos_not_changed_or_replaced") is not True):
        raise ValueError("FCP 2025 64 original ROI or live iNat annotation SHA gate failed")
    photos=csv_read(photo_bytes)
    live=csv_read(annotation_bytes)
    if len(photos)!=64 or len(live)!=64:
        raise ValueError("Source original photos/live UUIDs incomplete")
    pm={p["photo_id"]:p for p in photos}
    lm={a["photo_id"]:a for a in live}
    if len(pm)!=64 or len(lm)!=64 or pm.keys()!=lm.keys():
        raise ValueError("No source-identical 1:1 original photograph to live annotation match")
    rows=Counter()
    seasons=Counter()
    annotated=Counter()
    roi=Counter()
    photo_anthesis_count=0
    records_match=0
    for pid,p in pm.items():
        a=lm[pid]
        if (p["observation_id"]!=a["original_observation_uuid"] or
            p["observed_on"]!=a["original_observation_day"] or
            p["local_season"]!=a["original_photo_season"]):
            raise ValueError("Live annotation wrongly joined across original photo UUID/month")
        if (a["live_observed_on_matches_original"]=="True" and
            a["original_photo_still_listed_with_live_observation"]=="True"):
            records_match+=1
        season=p["local_season"]
        if season not in ("summer","autumn"):
            raise ValueError("Not the exact original 2025 summer/autumn source")
        status=annotation_label(a)
        score=roi_label(p)
        rows[(season,status,score)]+=1
        seasons[season]+=1
        roi[(season,score)]+=1
        if status.startswith("current_observation_annotation_"):
            annotated[(season,status)]+=1
    if seasons!=Counter({"summer":32,"autumn":32}):
        raise ValueError("Original photo count was not balanced by source observers/season")
    out.mkdir(parents=True,exist_ok=True)
    table=out/"2025_original_photo_ROI_x_current_observation_phenology_contingency_NOT_OPEN_FLOWER.csv"
    with table.open("x",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=COLS,lineterminator="\n")
        w.writeheader()
        for (season,status,score),n in sorted(rows.items()):
            w.writerow(dict(zip(COLS,(season,status,score,n,True))))
    fruit_white=rows[("autumn","current_observation_annotation_fruiting_NOT_photo","white_automatic_ROI_UNREVIEWED")]
    flower_white=rows[("summer","current_observation_annotation_flowering_NOT_photo","white_automatic_ROI_UNREVIEWED")]
    noflower_observations=sum(n for (season,status),n in annotated.items()
                              if status=="current_observation_annotation_no_evidence_of_flowering_NOT_photo")
    report={
        "schema":"mois_original_photo_roi_by_live_observation_phenology_2025_v1",
        "status":"REAL_AUTOMATED_PHOTO_ROI_VS_LIVE_OBSERVATION_ANNOTATION_BIAS_NOT_BOTANICAL_BLOOM",
        "source_original_FCP_photo_terminal_sha256":sha(photo_bytes),
        "source_current_MOIS_annotation_csv_sha256":sha(annotation_bytes),
        "original_photo_observations":len(pm),
        "live_original_photo_and_observation_day_both_match":records_match,
        "source_photo_observations_by_season":dict(seasons),
        "iNat_term_code_12_controlled_value_21":"no_evidence_of_flowering_observation_NOT_photo",
        "observations_with_current_no_evidence_annotation":noflower_observations,
        "autumn_fruiting_annotated_original_observation_with_white_automated_ROI_UNREVIEWED":fruit_white,
        "summer_flowering_annotated_original_observation_with_white_automated_ROI_UNREVIEWED":flower_white,
        "flower_petal_tepal_BOTANICALLY_CONFIRMED_photo_count":photo_anthesis_count,
        "organ_and_live_observation_annotations_are_not_independent_photo_bloom_proof":True,
        "live_annotations_not_frozen_in_original_2026_09_27_source":True,
        "no_causal_weather_effect_or_adaptation_estimated":True,
        "no_open_flower_colour_white_frequency_computed":True,
        "contingency_csv_sha256":sha(table.read_bytes()),
        "automatic_roi_by_season_not_confirmed_open_flowers":{
            season:{score:n for (s,score),n in roi.items() if s==season}
            for season in ("summer","autumn")},
        "phenology_by_season_current_observation_only":{
            season:{label:n for (s,label),n in annotated.items() if s==season}
            for season in ("summer","autumn")
        },
    }
    (out/"signed_2025_original_roi_observation_phenology_QC_receipt.json").write_text(
        json.dumps(report,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--original-FCP-2025-photo-ROI-zip",type=Path,required=True)
    p.add_argument("--current-MOIS-2025-observation-annotations-zip",type=Path,required=True)
    p.add_argument("--outdir",type=Path,required=True)
    a=p.parse_args()
    print(json.dumps(cross(a.original_FCP_2025_photo_ROI_zip,
                           a.current_MOIS_2025_observation_annotations_zip,
                           a.outdir),sort_keys=True))


if __name__=="__main__":
    main()
