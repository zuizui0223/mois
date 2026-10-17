"""Frozen real 400-photo iNaturalist observation-only phenology source integrity.

Runs WITHOUT recontacting live iNaturalist or treating photos as botanical
flower specimens. Only the original source artifact from 2026-10-17 may
establish the frozen annotation sampling/denominator values.
"""
from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

DATA=Path("results/mois_400_live_observation_phenology_20261017")
COVERAGE=DATA/"source_siteyear_season_annotation_coverage.csv"
RECEIPT=DATA/"verified_source_400_annotation_coverage_receipt.json"
SOURCE_SHA="b961b50ee6efe9ec08a5575a8eb1e34b1edf7b7fcb426226fb11ee901cfd3329"
GROUP_SHA="34f6997b85d7d98c355991ea79e9d7efe7194322da53857062ed30e7fe3f52d0"
FIELDS=(
  "inat_taxon_id","species","source_site_lat0p25","source_site_lon0p25",
  "original_source_season_year","original_local_season",
  "original_photo_observation_opportunities","current_UUID_resolved",
  "current_photo_and_original_observed_on_both_match",
  "current_observation_flowering_NOT_photo",
  "current_observation_fruiting_NOT_photo",
  "current_observation_flower_budding_NOT_photo",
  "current_observation_no_evidence_of_flowering_NOT_photo",
  "current_observation_phenology_missing_or_unknown",
  "current_observation_other_or_conflicting_annotation",
  "original_picture_OPEN_petal_tepal_independently_reviewed",
)
COUNT=("current_observation_flowering_NOT_photo",
       "current_observation_fruiting_NOT_photo",
       "current_observation_flower_budding_NOT_photo",
       "current_observation_no_evidence_of_flowering_NOT_photo",
       "current_observation_phenology_missing_or_unknown",
       "current_observation_other_or_conflicting_annotation")


def source():
    assert hashlib.sha256(COVERAGE.read_bytes()).hexdigest()==GROUP_SHA
    d=json.loads(RECEIPT.read_text(encoding="utf-8"))
    with COVERAGE.open(newline="",encoding="utf-8") as f:
        q=csv.DictReader(f)
        assert tuple(q.fieldnames or ())==FIELDS
        rows=list(q)
    return d,rows


def test_live_annotation_qc_source_run_and_photo_ids_are_exact():
    r,rows=source()
    assert r["schema"]=="mois_400_original_photo_phenology_current_API_source_observation_quality_summary_v1"
    assert r["status"]=="VERIFIED_400_CURRENT_OBSERVATION_ANNOTATION_AVAILABILITY_NOT_FLOWER_PHOTO_STAGE"
    assert r["verified_github_run"]=="https://github.com/zuizui0223/mois/actions/runs/38989216871"
    assert r["original_photo_manifest_sha256"]==SOURCE_SHA
    assert r["group_csv_sha256"]==GROUP_SHA
    assert r["n_photo_opportunities"]==400
    assert r["n_observation_UUIDs_and_photo_ids_current_API_exact_match"]==400
    assert r["n_distinct_original_observers"]==400
    assert r["n_species"]==5
    assert r["n_taxon_siteyear_season_groups"]==10
    assert len(rows)==10


def test_species_season_photo_denominators_unchanged_after_live_lookup():
    r,rows=source()
    pairs=Counter(x["species"] for x in rows)
    assert all(n==2 for n in pairs.values()) and len(pairs)==5
    assert {int(x["original_photo_observation_opportunities"]) for x in rows}=={40}
    assert {int(x["current_UUID_resolved"]) for x in rows}=={40}
    assert {int(x["current_photo_and_original_observed_on_both_match"]) for x in rows}=={40}
    assert all(x["original_picture_OPEN_petal_tepal_independently_reviewed"]=="0" for x in rows)


def test_annotator_missingness_not_zero_flowering():
    r,rows=source()
    summed={k:sum(int(row[k]) for row in rows) for k in COUNT}
    assert summed=={
        "current_observation_flowering_NOT_photo":32,
        "current_observation_fruiting_NOT_photo":54,
        "current_observation_flower_budding_NOT_photo":6,
        "current_observation_no_evidence_of_flowering_NOT_photo":2,
        "current_observation_phenology_missing_or_unknown":280,
        "current_observation_other_or_conflicting_annotation":26,
    }
    assert sum(summed.values())==400
    assert r["n_source_observations_without_interpretable_plant_phenology_annotation"]==280
    assert r["n_source_observations_with_any_current_plant_phenology_annotation"]==120
    assert r["annotation_absence_does_not_imply_nonflowering"] is True
    assert r["several_phenology_values_on_one_observation_do_not_establish_photo_stage"] is True


def test_annotation_coverage_strongly_taxon_dependent_but_not_white_flower_rate():
    r,rows=source()
    by_tax=Counter()
    for x in rows:
        by_tax[x["species"]]+=(40-int(x["current_observation_phenology_missing_or_unknown"]))
    assert by_tax=={
        "Datura discolor":78,
        "Aralia californica":26,
        "Linaria vulgaris":9,
        "Asclepias syriaca":6,
        "Achillea millefolium":1,
    }
    assert r["n_botanically_reviewed_open_photo_petals"]==0
    assert r["original_source_photo_candidates_are_not_known_flowering_photos"] is True
    assert r["no_original_photo_colour_ROI_opened_from_new_400"] is True
    assert r["no_global_or_within_taxon_white_flowers_or_monthly_weather_effect_estimated"] is True


def test_current_api_annotation_date_not_equated_with_frozen_source_annotations():
    r,rows=source()
    assert r["recorded_at_utc"].startswith("2026-10-17T")
    assert r["public_current_observation_annotations_are_not_frozen_20260927"] is True
    assert r["several_phenology_values_on_one_observation_do_not_establish_photo_stage"] is True
    assert sorted(x["original_source_season_year"] for x in rows)==["2023"]*6+["2025"]*4
