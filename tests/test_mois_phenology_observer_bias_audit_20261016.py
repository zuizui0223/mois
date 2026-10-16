"""Synthetic-only source-bias/phenology tests; no real annotations or pixels."""
import csv
import importlib.util
import json
from pathlib import Path

import pytest

ROOT=Path(__file__).resolve().parents[1]


def load(name,relative):
    spec=importlib.util.spec_from_file_location(name,ROOT/relative)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

collector_test=load("mois_synthetic_photo_collection","tests/test_mois_collect_sample_roi_white_20261016.py")
m=load("mois_obs_bias","scripts/analysis/mois_phenology_observer_bias_audit_20261016.py")


def synthetic(tmp_path):
    packets,terminal,_=collector_test.data(tmp_path)
    collected=tmp_path/"collected"
    collector_test.m.collect(packets,terminal,collected,synthetic=True)
    sample=collected/"private_joined_sample_flower_colours.csv"
    report=collected/"sample_collection_receipt.json"
    overlay=tmp_path/"annotation_overlay.csv"
    raw = [
      ["o1","12","13","LOOKED_UP","inat_api_v1","2026-10-16"],
      ["o2","12","21","LOOKED_UP","inat_api_v1","2026-10-16"],
      ["o9","12","","LOOKED_UP","inat_api_v1","2026-10-16"],
    ]
    with overlay.open("w",newline="",encoding="utf8") as fh:
        w=csv.writer(fh)
        w.writerow(m.OVERLAY_COLUMNS)
        w.writerows(raw)
    return sample,report,overlay,raw


def test_distinct_flowers_and_no_annotation_and_no_evidence(tmp_path):
    photo,receipt,overlay,rows=synthetic(tmp_path)
    private=tmp_path/"private.json"
    public=tmp_path/"aggregate.json"
    summary=m.audit(photo,receipt,overlay,m.digest(overlay),private,public,synthetic=True)
    assert summary["selected_observations"]==3
    assert summary["annotations_looked_up"]==3
    assert summary["annotations_NOT_verified"]==0
    assert summary["phenology_by_observation_status"]["FLOWERS"]["selected_observations"]==1
    assert summary["phenology_by_observation_status"]["NO_FLOWERS_FRUITS"]["selected_observations"]==1
    assert summary["phenology_by_observation_status"]["NO_ANNOTATION"]["selected_observations"]==1
    assert summary["phenology_by_observation_status"]["NO_ANNOTATION"]["design_weighted_white_given_ROI_classifiable"]==1.0
    assert summary["contradictions_explicit_no_flowers_but_ROI_classified"]==0
    assert summary["annotation_is_observation_not_photo_label"] is True
    assert summary["unsubmitted_plants_observation_bias_identified"] is False
    assert "o1" not in public.read_text() and "inat_taxon_id" not in public.read_text()
    assert "u1" not in public.read_text()
    assert private.exists()
    with pytest.raises(ValueError,match="overwrite"):
        m.audit(photo,receipt,overlay,m.digest(overlay),private,public,synthetic=True)


def test_not_looked_up_is_not_nonflower(tmp_path):
    photo,receipt,overlay,rows=synthetic(tmp_path)
    rows[2][3]="NOT_LOOKED_UP"
    with overlay.open("w",newline="") as f:
        w=csv.writer(f);w.writerow(m.OVERLAY_COLUMNS);w.writerows(rows)
    result=m.audit(photo,receipt,overlay,m.digest(overlay),
                   tmp_path/"private.json",tmp_path/"public.json",synthetic=True)
    assert result["annotations_NOT_verified"]==1
    assert result["phenology_by_observation_status"]["NOT_VERIFIED"]["selected_observations"]==1
    assert result["phenology_by_observation_status"]["NO_ANNOTATION"]["selected_observations"]==0


def test_negative_phenology_cannot_conflict_with_flowers(tmp_path):
    photo,receipt,overlay,rows=synthetic(tmp_path)
    rows[1][2]="13|21"
    with overlay.open("w",newline="") as f:
        w=csv.writer(f);w.writerow(m.OVERLAY_COLUMNS);w.writerows(rows)
    with pytest.raises(ValueError,match="Contradictory"):
        m.audit(photo,receipt,overlay,m.digest(overlay),
                tmp_path/"private.json",tmp_path/"public.json",synthetic=True)


def test_partial_or_forged_annotation_overlay_aborts(tmp_path):
    photo,receipt,overlay,rows=synthetic(tmp_path)
    with overlay.open("w",newline="") as f:
        w=csv.writer(f);w.writerow(m.OVERLAY_COLUMNS);w.writerows(rows[:-1])
    with pytest.raises(ValueError,match="annotation lookup row"):
        m.audit(photo,receipt,overlay,m.digest(overlay),
                tmp_path/"private.json",tmp_path/"public.json",synthetic=True)
    with pytest.raises(ValueError,match="SHA256"):
        m.audit(photo,receipt,overlay,"0"*64,
                tmp_path/"private.json",tmp_path/"public.json",synthetic=True)


def test_synthetic_roi_measurements_cannot_impersonate_real_source(tmp_path):
    photo,receipt,overlay,rows=synthetic(tmp_path)
    with pytest.raises(ValueError,match="sealed Mois"):
        m.audit(photo,receipt,overlay,m.digest(overlay),
                tmp_path/"private.json",tmp_path/"public.json",synthetic=False)


def test_missing_month_and_positive_roi_conflicting_annotation(tmp_path):
    photo,receipt,overlay,rows=synthetic(tmp_path)
    # Annotation "No Flowers" for the white-classified observation is a
    # disagreement that must be highlighted for manual review, not fixed away.
    rows[0][2]="21"
    with overlay.open("w",newline="") as f:
        w=csv.writer(f);w.writerow(m.OVERLAY_COLUMNS);w.writerows(rows)
    result=m.audit(photo,receipt,overlay,m.digest(overlay),
                   tmp_path/"private.json",tmp_path/"public.json",synthetic=True)
    assert result["contradictions_explicit_no_flowers_but_ROI_classified"]==1
    assert result["contradictions_require_blind_manual_review"] is True
