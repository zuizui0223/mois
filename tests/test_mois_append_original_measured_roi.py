"""Protect source-only re-use of measured original FCP ROI photo labels."""
import pytest
from scripts.analysis.mois_append_fcp_six_real_roi_waves_20261017 import (
  normalize, SOURCES, HUES, PREVIOUS_PHOTO_COUNT, FINAL_PHOTO_COUNT,
  PREVIOUS_CLASSIFIED, FINAL_CLASSIFIED, EXPECTED_ADDED_BY_SOURCE, require,
)

def test_measured_original_photo_counts_make_valid_source_increment():
    assert len(SOURCES)==len(EXPECTED_ADDED_BY_SOURCE)==6
    assert PREVIOUS_PHOTO_COUNT==247032
    assert FINAL_PHOTO_COUNT==248312
    assert sum(x[0] for x in EXPECTED_ADDED_BY_SOURCE.values())==1280
    assert sum(x[1] for x in EXPECTED_ADDED_BY_SOURCE.values())==497
    assert PREVIOUS_CLASSIFIED==120310
    assert FINAL_CLASSIFIED==120807

def test_normalize_requires_provenance_and_never_infers_new_bloom():
    row={"photo_id":"123","inat_taxon_id":"44","morph":"mixed_uncertain",
         "species":"Test sp.","measurement_status":"not_evaluable_roi_or_flip_gate",
         "observed_on":"2025-06-17","photo_license":"cc-by"}
    x=normalize(row,"rgfca_monthly_first64")
    assert x["photo_id"]=="123"
    assert x["source_original_measured"]=="rgfca_monthly_first64"
    assert x["old_mois_or_new_fcp"]=="fcp_additional_ROI"
    assert x["observed_on"]=="2025-06-17"
    assert x["morph"] not in HUES
    assert x["time_observed_at"]=="" # an original date is not a local time

def test_abort_invalid_frozen_data_contract():
    with pytest.raises(ValueError):
        require(False,"must abort missing original source")
