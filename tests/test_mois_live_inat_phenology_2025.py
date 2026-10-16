"""Mois current iNaturalist plant-phenology codes are OBSERVATION-only.

Assertions deliberately do not elevate an observation annotation to a
visually confirmed petal/tepal phenotype of one original picture.
"""
from __future__ import annotations

import json
import zipfile

import pytest

from scripts.analysis.mois_audit_live_inat_observation_phenology_2025 import (
    TERM, VALUE, plant_status, load_source,
)


def a(term, val):
    return {"controlled_attribute_id":term,"controlled_value_id":val}


def test_official_2026_iNat_plant_annotation_codebook():
    assert TERM == 12
    assert VALUE == {
        13:"flowering",
        14:"fruiting",
        15:"flower_budding",
        21:"no_evidence_of_flowering",
    }


@pytest.mark.parametrize(("val","label"), [
    (13,"flowering"),
    (14,"fruiting"),
    (15,"flower_budding"),
    (21,"no_evidence_of_flowering"),
])
def test_live_annotation_labels_are_explicitly_observation_level(val,label):
    status, ids=plant_status([a(12,val)])
    assert status=="current_observation_annotation_"+label+"_NOT_photo"
    assert ids==[val]


def test_live_api_annotation_field_absence_is_not_nonflower():
    status,values=plant_status(None)
    assert status=="no_current_observation_phenology_field_unknown"
    assert values==[]
    status,values=plant_status([])
    assert status=="no_current_plant_phenology_annotation_unknown"
    assert values==[]


def test_irrelevant_terms_not_assumed_plant_phenology():
    status,values=plant_status([a(17,19),a(9,10)])
    assert status=="no_current_plant_phenology_annotation_unknown"
    assert values==[]


def test_conflict_and_duplicate_crowd_votes_do_not_imply_anthesis():
    status,values=plant_status([a(12,13),a(12,14)])
    assert status=="conflicting_current_observation_phenology_NOT_photo"
    assert values==[13,14]
    status,values=plant_status([a(12,21),a(12,21)])
    assert status=="current_observation_annotation_no_evidence_of_flowering_NOT_photo"
    assert values==[21]


def test_unknown_phenology_code_does_not_fall_into_negative_class():
    status,values=plant_status([a(12,999)])
    assert status=="unrecognized_current_observation_phenology_NOT_photo"
    assert values==[999]


def test_malformed_values_refuse_to_impute_nonflower():
    with pytest.raises(ValueError,match="phenology value"):
        plant_status([a(12,"not-a-code")])
    with pytest.raises(ValueError,match="annotation record"):
        plant_status(["not-a-dict"])


def test_unsigned_synthetic_photo_archive_never_admitted(tmp_path):
    bad=tmp_path/"fake.zip"
    with zipfile.ZipFile(bad,"w") as z:
        z.writestr("postblind/first_original_source_matched_season_64_photo_roi_QC.csv",
                   "photo_id,observation_id,observed_on\n1,fake,2025-07-01\n")
        z.writestr("postblind/first_original_source_matched_season_64_roi_receipt.json",
                   json.dumps({"postblind_join_sha256":"fake"}))
    with pytest.raises(ValueError,match="exact frozen FCP"):
        load_source(bad)
