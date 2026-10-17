"""Synthetic guard: only original wave4 whole source photos, never FCP ecology."""
from scripts.analysis.mois_import_full_FCP_wave4_completed_original_photos import (
    EXPECTED_FCP_PHOTOS,EXPECTED_FCP_TAXA,HUES,SOURCE_SCHEMA,
    EXPECTED_SOURCE_SHA,BASE_COLS
)

def test_wave4_complete_roi_cannot_be_relabelled_as_source_only():
    assert EXPECTED_FCP_PHOTOS==49800
    assert EXPECTED_FCP_TAXA==498
    assert SOURCE_SCHEMA=="fcp_same_149900_method_wave4_real_photo_terminals_v1"
    assert len(EXPECTED_SOURCE_SHA)==64

def test_copy_photo_ids_not_FCP_polymorphism_measures():
    assert HUES==("white","yellow_orange","red_pink","blue_purple")
    assert {"photo_id","species","inat_taxon_id","morph","observation_id",
            "latitude","longitude","measurement_status"}.issubset(set(BASE_COLS))
    assert not {"D","W","Ppoly","Fst","genotype"}.intersection(BASE_COLS)
