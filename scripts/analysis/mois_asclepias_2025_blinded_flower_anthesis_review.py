#!/usr/bin/env python3
"""Blinded two-botanist open-flower validation for the REAL 2025 Asclepias seasons.

Provide image URL and source photo identity to reviewers WITHOUT observation
month/season, climate, or automated ROI colour. After independent reviews,
rejoin original source/date and ROI statuses; never replace a rejected photo.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,zipfile
from collections import defaultdict
from pathlib import Path
from scripts.analysis.mois_validate_open_flower_photo_review import REVIEW_COLUMNS,decide

INPUT="postblind/first_original_source_matched_season_64_photo_roi_QC.csv"
RECEIPT="postblind/first_original_source_matched_season_64_roi_receipt.json"
EXPECTED="1b16cdabd2936567dee87b8344019f40f315600c6510dfc1f25a03835024b6a0"
RESULT_COLS=("measurement_id","photo_id","species","observed_year_month","local_season",
             "observer_id","botanical_review_status","roi_original_automated_status",
             "original_photo_roi_morph_is_not_botanical_truth",
             "roi_same_white_yellow_organs_mixed","botanical_flower_colour_eligible",
             "botanical_flower_colour_coded","denominator_count_independent_observation")
MIN_CONFIRMED_PER_SEASON=20


def sha(data:bytes):return hashlib.sha256(data).hexdigest()


def load(source:Path)->tuple[list[dict],dict]:
    with zipfile.ZipFile(source) as z:
        raw=z.read(INPUT)
        receipt=json.loads(z.read(RECEIPT))
        freeze=json.loads(z.read("source/frozen/original_season_photo_freeze_receipt.json"))
    if (receipt.get("status")!="FIRST_ORIGINAL_SEASONAL_PIXELS_MEASURED_BOTANICAL_ANTHESIS_NOT_VALIDATED"
        or receipt.get("postblind_join_sha256")!=sha(raw)
        or receipt.get("postblind_join_sha256")!=EXPECTED
        or receipt.get("original_photos_measured")!=64
        or freeze.get("observer_ids_unique_across_seasons") is not True
        or freeze.get("source_photo_queue_sha256") is not None
        or freeze.get("file_sha256",{}).get("original_source_season_frozen_128_photos.csv")!=receipt.get("source_photo_queue_sha256")):
        raise ValueError("Not intact frozen first64 FCP 2025 Asclepias pixel/source receipt")
    rows=list(csv.DictReader(io.StringIO(raw.decode("utf-8"))))
    if len(rows)!=64 or len({r["photo_id"] for r in rows})!=64:
        raise ValueError("FCP original 64 photo-ID uniqueness failed")
    if len({r["observer_id"] for r in rows})!=64 or {r["species"] for r in rows}!={"Asclepias syriaca"}:
        raise ValueError("Original source observers or taxon drift")
    if {r["local_season"] for r in rows}!={"summer","autumn"}:
        raise ValueError("Not two original photographed seasons")
    return rows,receipt


def blank(source:Path,dest:Path)->dict:
    if dest.exists():raise ValueError("Cannot rewrite botanist original-season review template")
    rows,receipt=load(source)
    with dest.open("x",newline="",encoding="utf-8") as h:
        w=csv.DictWriter(h,fieldnames=REVIEW_COLUMNS,lineterminator="\n");w.writeheader()
        for r in sorted(rows,key=lambda x:x["photo_id"]):
            w.writerow({
                "measurement_id":r["measurement_id"],"photo_id":r["photo_id"],
                "species":r["species"],"observation_url":"https://www.inaturalist.org/observations/"+r["observation_id"],
                "photo_url_large":r["photo_url_large"],"photo_license":r["photo_license"],
            })
    return {"status":"TRUE_ORIGINAL_PHOTOS_64_READY_FOR_TWO_INDEPENDENT_BOTANISTS_NOT_REVIEWED",
            "n_original_photos":64,"n_assigned_flower_phases":0,
            "photo_form_original_roi_unreviewed_sha256":sha(dest.read_bytes()),
            "original_2025_photo_source_sha256":receipt["postblind_join_sha256"],
            "reviewers_not_shown_photo_morph_season_or_climate":True}


def audit(source:Path,form:Path,outdir:Path)->dict:
    if outdir.exists() and any(outdir.iterdir()):
        raise ValueError("Cannot overwrite a botanist source review adjudication")
    original,receipt=load(source)
    known={r["measurement_id"]:r for r in original}
    observed=set()
    results=[]
    group=defaultdict(lambda:defaultdict(int))
    with form.open(newline="",encoding="utf-8") as f:
        rdr=csv.DictReader(f)
        if tuple(rdr.fieldnames or ())!=REVIEW_COLUMNS:
            raise ValueError("Human review saw forbidden date/season/colour information")
        for row in rdr:
            mid=row["measurement_id"]
            if mid not in known or mid in observed:raise ValueError("Missing/duplicate independent source photo ID")
            observed.add(mid);r=known[mid]
            if (row["photo_id"]!=r["photo_id"] or row["species"]!=r["species"] or
                row["photo_url_large"]!=r["photo_url_large"] or row["photo_license"]!=r["photo_license"]
                or row["observation_url"]!="https://www.inaturalist.org/observations/"+r["observation_id"]):
                raise ValueError("Photo identities or original image URL changed")
            a=row["r1_reviewer_id"].strip();b=row["r2_reviewer_id"].strip()
            st1=row["r1_open_flower_stage"].strip();st2=row["r2_open_flower_stage"].strip()
            t1=row["r1_roi_tissue"].strip();t2=row["r2_roi_tissue"].strip()
            if (any((st1,t1)) and not a or any((st2,t2)) and not b or a and b and a==b):
                raise ValueError("Independent botanical reviewer identities required")
            status=decide(st1,t1,st2,t2)
            c=row["adjudicator_id"].strip()
            ast=row["adjudicated_stage"].strip()
            at=row["adjudicated_roi_tissue"].strip()
            if any((c,ast,at)):
                if (not all((c,ast,at)) or c in (a,b) or status!="needs_independent_adjudication"):
                    raise ValueError("Third independent adjudicator required for conflicting full reviews")
                status=decide(ast,at,ast,at)
            good=r["measurement_status"]=="classified_four_state_morph"
            mix=r["same_roi_white_yellow_organ_mixed"]=="True"
            primary=status=="petal_tepal_open_flower_reviewed" and good and not mix
            value=r["morph"] if primary else ""
            if status=="petal_tepal_open_flower_reviewed" and not primary:
                status="reviewed_open_flower_but_unclassifiable_or_organ_mixed"
            group[r["local_season"]][status]+=1
            group[r["local_season"]]["n_original"]+=1
            group[r["local_season"]]["n_confirmed_primary_photo_colour"]+=int(primary)
            results.append(dict(zip(RESULT_COLS,(
              mid,r["photo_id"],r["species"],r["observed_year_month"],r["local_season"],
              r["observer_id"],status,r["measurement_status"],True,mix,primary,
              value,int(primary)))))
    if observed!=set(known):raise ValueError("All 64 original photos must remain in review output")
    outdir.mkdir(parents=True,exist_ok=True)
    photo=outdir/"2025_asclepias_original_photo_botanical_stage_verified_or_unknown.csv"
    with photo.open("x",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=RESULT_COLS,lineterminator="\n")
        w.writeheader();w.writerows(sorted(results,key=lambda r:int(r["photo_id"])))
    season_summary={name:dict(group[name]) for name in ("summer","autumn")}
    can_estimate=all(season_summary[s].get("n_confirmed_primary_photo_colour",0)>=MIN_CONFIRMED_PER_SEASON
                     for s in ("summer","autumn"))
    # Even after enough confirmed flowers, a model is not run here; it needs
    # explicit missingness/observer and source-climate confounding controls.
    data={"schema":"mois_summer_autumn_2025_photo_anthesis_two_independent_review_v1",
          "status":("BOTANICAL_OBSERVATION_REVIEW_PENDING" if any(
            group[s].get(c,0)>0 for s in ("summer","autumn")
            for c in ("unreviewed","needs_independent_adjudication",
                      "unresolved_flower_stage_or_roi")) else
            "BOTANICAL_OBSERVATION_REVIEW_QC_COMPLETE_NO_ASSOCIATION"),
          "source_original_2025_postblind_sha256":receipt["postblind_join_sha256"],
          "source_review_form_sha256":sha(form.read_bytes()),
          "photo_review_terminal_sha256":sha(photo.read_bytes()),
          "by_season":season_summary,
          "minimum_botanically_confirmed_primary_flower_photos_each":MIN_CONFIRMED_PER_SEASON,
          "enough_confirmed_each_season_to_consider_descriptive_rate":can_estimate,
          "no_seasonal_white_environment_regression_performed":True,
          "reviewers_not_shown_climate_date_or_automated_colour":True,
          "unreviewed_plant_never_negative_flowering_evidence":True}
    (outdir/"2025_asclepias_open_flower_review_receipt.json").write_text(
        json.dumps(data,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    return data


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fcp-season64-artifact",type=Path,required=True)
    p.add_argument("--outdir",type=Path,required=True)
    p.add_argument("--filled-independent-review-sheet",type=Path)
    a=p.parse_args()
    if a.filled_independent_review_sheet:
        print(json.dumps(audit(a.fcp_season64_artifact,a.filled_independent_review_sheet,a.outdir),sort_keys=True))
    else:
        a.outdir.mkdir(parents=True,exist_ok=True)
        sheet=a.outdir/"UNREVIEWED_2025_asclepias_two_botanist_photo_review.csv"
        report=blank(a.fcp_season64_artifact,sheet)
        (a.outdir/"review_form_source_receipt.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
        print(json.dumps(report,sort_keys=True))


if __name__=="__main__":main()
