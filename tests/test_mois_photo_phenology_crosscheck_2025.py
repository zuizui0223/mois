"""Fail-closed MOIS first64 photo-vs-observation phenology staging tests."""
from __future__ import annotations

import csv
import io
import zipfile
import json

import pytest

from scripts.analysis.mois_crosscheck_2025_photo_roi_current_phenology import (
    roi_label, annotation_label, cross, csv_read,
)


def photo(assessment="classified_four_state_morph",morph="white",mixed="False"):
    return {"measurement_status":assessment,"morph":morph,
            "same_roi_white_yellow_organ_mixed":mixed}


def live(status="current_observation_annotation_fruiting_NOT_photo",verified="False"):
    return {"plant_phenology_observation_annotation_status_NOT_PHOTO":status,
            "photo_open_flower_petal_tepal_independently_verified":verified}


def test_automatic_photo_roi_white_is_unreviewed_not_white_open_flower():
    assert roi_label(photo())=="white_automatic_ROI_UNREVIEWED"
    assert roi_label(photo(morph="yellow_orange"))=="chromatic_automatic_ROI_UNREVIEWED"


def test_white_yellow_same_roi_organs_cannot_be_called_white_morph():
    assert roi_label(photo(mixed="True"))=="same_ROI_white_yellow_floral_organs_mixed_not_morph"


def test_original_failed_roi_does_not_become_blooming_negative():
    assert roi_label(photo(assessment="not_evaluable_roi_or_flip_gate",
                           morph="mixed_uncertain"))=="ROI_failed_or_unclassifiable_UNKNOWN_flower_stage"
    with pytest.raises(ValueError,match="Failed original"):
        roi_label(photo(assessment="not_evaluable_roi_or_flip_gate",
                        morph="mixed_uncertain",mixed="True"))


def test_fruiting_observation_and_white_photo_roi_are_NOT_proof_of_plant_petal_colour():
    assert annotation_label(live())=="current_observation_annotation_fruiting_NOT_photo"
    assert roi_label(photo())=="white_automatic_ROI_UNREVIEWED"


def test_iNat_no_evidence_is_observation_not_photo_or_certified_nonflower():
    assert annotation_label(live("current_observation_annotation_no_evidence_of_flowering_NOT_photo")) == (
           "current_observation_annotation_no_evidence_of_flowering_NOT_photo")


def test_explicitly_reject_autobiographical_botanic_review_from_live_annotations():
    with pytest.raises(ValueError,match="falsely called botanical"):
        annotation_label(live("current_observation_annotation_flowering_NOT_photo",verified="True"))


def test_unrecognized_current_api_annotation_never_converts_to_known_true_colour():
    with pytest.raises(ValueError,match="unexpectedly recoded"):
        annotation_label(live("automatic_white_photo_blooms_inferred"))


def test_unsigned_original_or_live_source_never_produces_2025_colour_effect(tmp_path):
    p=tmp_path/"fcp.zip";l=tmp_path/"inat.zip";out=tmp_path/"output"
    with zipfile.ZipFile(p,"w") as z:
        z.writestr("postblind/first_original_source_matched_season_64_photo_roi_QC.csv",
                   "photo_id,measurement_status\n1,classified_four_state_morph\n")
        z.writestr("postblind/first_original_source_matched_season_64_roi_receipt.json",
                   json.dumps({"postblind_join_sha256":"notsha"}))
    with zipfile.ZipFile(l,"w") as z:
        z.writestr("mois_2025_live_plant_phenology_supplement/original_2025_64_live_iNaturalist_observation_plant_phenology_NOT_PHOTO.csv",
                   "photo_id,plant_phenology_observation_annotation_status_NOT_PHOTO\n1,fruiting\n")
        z.writestr("mois_2025_live_plant_phenology_supplement/2025_Asclepias_live_public_observation_phenology_receipt.json",
                   json.dumps({"status":"not verified"}))
    with pytest.raises(ValueError,match="SHA gate"):
        cross(p,l,out)


def test_csv_requires_readable_signed_schema():
    with pytest.raises(ValueError,match="no header"):
        csv_read(b"")
