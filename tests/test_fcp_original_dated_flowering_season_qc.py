"""Synthetic contract checks for exact-date post-ROI flowering-season tables.

These tests NEVER access a real photograph and NEVER claim a biological result.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from scripts.analysis.fcp_original_dated_flowering_season_qc import (
    analyze, season, sha, source_day,
)


def put_csv(path: Path, cols: list[str], rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        out = csv.DictWriter(f, fieldnames=cols, lineterminator="\n")
        out.writeheader()
        out.writerows(rows)


def fixtures(tmp_path: Path):
    source = tmp_path / "original-monthly-source.csv"
    terminal = tmp_path / "post-blind-original-roi.csv"
    source_receipt = tmp_path / "source-receipt.json"
    roi_receipt = tmp_path / "new-roi-receipt.json"
    src_cols = ["inat_taxon_id","species","photo_id","observation_id",
                "observer_id","latitude","longitude","observed_on",
                "source_snapshot_date"]
    entries = [
        # 1+2 are photographs of ONE source observation, not two flowers.
        (1,101,"Northern species","o1","u1",35.0,135.0,"2026-01-05","white",1,0),
        (2,101,"Northern species","o1","u1",35.0,135.0,"2026-01-05","white",1,0),
        (3,101,"Northern species","o2","u2",35.0,135.0,"2026-01-05","yellow_orange",0,1),
        # One original observation producing incompatible photo labels:
        (4,101,"Northern species","o3","u3",35.0,135.0,"2026-01-05","white",1,0),
        (5,101,"Northern species","o3","u3",35.0,135.0,"2026-01-05","yellow_orange",0,1),
        # White ray and yellow disc in SAME ROI: organ mixing, not morph evidence.
        (6,101,"Northern species","o4","u4",35.0,135.0,"2026-01-05","white",.2,.2),
        # Unclassifiable does NOT indicate a non-flowering observation.
        (7,101,"Northern species","o5","u5",35.0,135.0,"2026-01-05","mixed_uncertain",0,0),
        (8,104,"Southern species","o6","u6",-34.,20.,"2026-07-10","white",1,0),
        (9,104,"Southern species","o7","u7",-34.,20.,"2026-07-10","red_pink",0,0),
        # Tropics should use calendar month, not a guessed four-season system.
        (10,101,"Northern species","o8","u8",5.,135.,"2026-01-05","white",1,0),
        # Original observation date unknown / later than dated snapshot.
        (11,101,"Northern species","o9","u9",35.,135.,"2026-01","white",1,0),
        (12,101,"Northern species","o10","u10",35.,135.,"2026-10-05","white",1,0),
        # Never measured: remain in source-photo opportunity denominator only.
        (13,101,"Northern species","o11","u11",35.,135.,"2026-01-05","white",1,0),
        (14,101,"Northern species","o12","u12",35.,135.,"2026-01-05","white",1,0),
    ]
    all_source=[]
    all_roi=[]
    for pid,tax,name,ob,owner,lat,lon,day,morph,white,yellow in entries:
        base={
            "inat_taxon_id":str(tax),"species":name,"photo_id":str(pid),
            "observation_id":ob,"observer_id":owner,
            "latitude":str(lat),"longitude":str(lon),"observed_on":day,
        }
        all_source.append({**base,"source_snapshot_date":"2026-09-27"})
        if pid<=12:
            passed=morph!="mixed_uncertain"
            all_roi.append({**base,
                "measurement_status":("classified_four_state_morph" if passed
                                      else "not_evaluable_roi_or_flip_gate"),
                "morph":morph,"flower_fraction_white":str(white),
                "flower_fraction_yellow":str(yellow),
                "global_classifiable":str(passed).lower()})
    put_csv(source,src_cols,all_source)
    put_csv(terminal,[*all_roi[0]],all_roi)
    source_receipt.write_text(json.dumps({
        "status":"COMPLETE_TARGETED_1000_MONTHLY_SOURCE_IDS_COUNT_MATCHED_NO_PHOTO_PIXELS",
        "source_snapshot_date":"2026-09-27",
        "selected_all_photo_IDs_csv_sha256":sha(source),
        "n_selected_original_photo_IDs":len(entries),
        "enrolled_outcome_blind_taxa":2,
        "historical_flower_colour_CSV_imported":False,
        "original_image_pixels_opened":False,
    }))
    roi_receipt.write_text(json.dumps({
        "status":"SYNTHETIC_ORIGINAL_GEOGRAPHIC_PAIR_TEST_NOT_BIOLOGY",
        "original_geographic_source_join_terminal_sha256":sha(terminal),
    }))
    return source,source_receipt,terminal,roi_receipt


def test_genuine_observed_day_is_required_and_hemispheres_use_local_seasons():
    assert str(source_day("2026-01-05"))=="2026-01-05"
    for raw in ("2026-01","2026-00-04","2026-02-30","2026-10-05",
                "2026/01/05","",None):
        assert source_day(raw) is None
    assert season(7,"north")=="local_summer"
    assert season(7,"south")=="local_winter"
    assert season(1,"south")=="local_summer"
    assert season(1,"north")=="local_winter"
    assert season(1,"tropical")=="tropical_month_only"


def test_postroi_monthly_white_denominator_is_independent_original_observations(tmp_path):
    args=fixtures(tmp_path)
    destination=tmp_path/"qc"
    receipt=analyze(*args,destination,_synthetic_only=True)
    assert receipt["status"]=="SYNTHETIC_TEST_NO_BIOLOGICAL_FINDING"
    assert receipt["source_taxa"]==2
    assert receipt["all_original_photos_in_supplied_source"]==14
    assert receipt["original_photos_with_valid_observation_day"]==12
    assert receipt["original_photos_with_invalid_or_unknown_observation_day"]==2
    assert receipt["new_roi_terminal_photos"]==12
    assert receipt["new_roi_photos_without_valid_observation_day"]==2
    assert receipt["new_roi_classifiable_flower_photos_with_valid_day"]==9
    assert receipt["same_roi_white_yellow_organ_mixed_excluded_from_primary"]==1
    assert receipt["all_93209351_photo_pixels_classified"] is False
    with (destination/"original_dated_month_site_white_qc.csv").open(newline="") as f:
        rows=list(csv.DictReader(f))
    assert len(rows)==3
    north=next(x for x in rows if x["inat_taxon_id"]=="101" and
               x["hemisphere"]=="north")
    assert north["year_month"]=="2026-01"
    assert north["season"]=="local_winter"
    assert north["dated_original_photo_opportunities"]=="9"
    assert north["new_ROI_photos_with_known_date"]=="7"
    assert north["classifiable_flower_ROIs"]=="6"
    assert north["same_ROI_white_yellow_organ_mixed"]=="1"
    assert north["unclassifiable_ROIs"]=="1"
    assert north["usable_flower_photo_labels"]=="5"
    assert north["independent_usable_observations"]=="2"
    assert north["white_observations"]=="1"
    assert north["nonwhite_observations"]=="1"
    assert north["contradictory_multiphoto_observations"]=="1"
    assert north["distinct_usable_observers"]=="2"
    assert float(north["white_rate_per_observation_descriptive"])==.5
    south=next(x for x in rows if x["inat_taxon_id"]=="104")
    assert south["season"]=="local_winter"
    assert south["calendar_month"]=="7"
    assert south["independent_usable_observations"]=="2"
    assert float(south["white_rate_per_observation_descriptive"])==.5
    trop=next(x for x in rows if x["hemisphere"]=="tropical")
    assert trop["season"]=="tropical_month_only"
    assert float(trop["white_rate_per_observation_descriptive"])==1.
    assert receipt["monthly_site_QC_csv_sha256"] == sha(
        destination/"original_dated_month_site_white_qc.csv")
    with pytest.raises(ValueError,match="overwrite"):
        analyze(*args,destination,_synthetic_only=True)


def test_modified_original_source_date_rejected_even_if_local_csv_has_same_photos(tmp_path):
    args=fixtures(tmp_path)
    src,receipt,measured,roi=args
    raw=src.read_text()
    src.write_text(raw.replace("2026-01-05","2026-01-06",1))
    with pytest.raises(ValueError,match="complete SHA-verified"):
        analyze(*args,tmp_path/"bad",_synthetic_only=True)


def test_modified_observer_is_blocked_after_roi_sha_reseal(tmp_path):
    args=fixtures(tmp_path)
    src,sreceipt,measured,rreceipt=args
    old=measured.read_text()
    measured.write_text(old.replace("u1","replaced_user",1))
    rec=json.loads(rreceipt.read_text())
    rec["original_geographic_source_join_terminal_sha256"]=sha(measured)
    rreceipt.write_text(json.dumps(rec))
    with pytest.raises(ValueError,match="source identity/date differs"):
        analyze(*args,tmp_path/"bad",_synthetic_only=True)


def test_orphan_roi_photo_is_rejected_even_with_new_roi_sha(tmp_path):
    args=fixtures(tmp_path)
    src,sreceipt,measured,rreceipt=args
    rows=list(csv.DictReader(measured.open(newline="")))
    rows[0]["photo_id"]="987654321"
    put_csv(measured,list(rows[0]),rows)
    rec=json.loads(rreceipt.read_text())
    rec["original_geographic_source_join_terminal_sha256"]=sha(measured)
    rreceipt.write_text(json.dumps(rec))
    with pytest.raises(ValueError,match="outside the exact original source"):
        analyze(*args,tmp_path/"bad",_synthetic_only=True)


def test_synthetic_receipts_cannot_be_used_by_production_path(tmp_path):
    args=fixtures(tmp_path)
    with pytest.raises(ValueError,match="genuinely remeasured"):
        analyze(*args,tmp_path/"bad",_synthetic_only=False)
