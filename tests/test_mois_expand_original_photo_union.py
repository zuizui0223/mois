"""Check MOIS existing broad source plus FCP high-depth photo dedup conventions."""
import pytest
from scripts.analysis.mois_expand_original_100543_with_fcp_149900 import (
    normalized, COLORS, PHOTO_SHA, FCP_SHA, FCP_N,
)

def specimen(**overrides):
    row={"photo_id":"123","inat_taxon_id":"42","species":"Example plant",
         "morph":"white","global_classifiable":"True","latitude":"30",
         "longitude":"-120","observation_id":"o1","photo_license":"cc-by",
         "measurement_status":"classified_four_state_morph"}
    row.update(overrides)
    return row

def test_hue_sanitization_and_old_new_source():
    r=normalized(specimen(),"fcp_discovery","new_fcp")
    assert r["photo_id"]=="123" and r["morph"] in COLORS
    assert r["source_original_measured"]=="fcp_discovery"
    assert r["old_mois_or_new_fcp"]=="new_fcp"

def test_failed_colour_measurements_not_imputed():
    r=normalized(specimen(morph="mixed_uncertain",measurement_status="image_acquisition_failed"),
                 "fcp_prospective","new_fcp")
    assert r["morph"]=="mixed_uncertain"
    assert r["morph"] not in COLORS

def test_original_historical_pre_pixel_blank_is_preserved():
    r=normalized(specimen(morph="",measurement_status="transport_unresolved_pre_pixel"),
                 "historical_mois_20261014","historical")
    assert r["morph"]==""

@pytest.mark.parametrize("bad",["unknown", "light_pink", "background", "nonflower"])
def test_reject_nonfrozen_colour_states(bad):
    with pytest.raises(ValueError,match="classification"):
        normalized(specimen(morph=bad),"fcp_validation","new_fcp")

def test_immutability_3_original_fcp_cohorts():
    assert PHOTO_SHA == "2d6097ba1b0434632e13bbe79dd78180608a287ff5dabebc770b5c0a2cfb52a0"
    assert set(FCP_SHA)==set(FCP_N)=={"discovery","validation","prospective"}
    assert sum(FCP_N.values())==149900
