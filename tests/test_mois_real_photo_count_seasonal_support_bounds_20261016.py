"""Source-bound real original FCP photo opportunities are NOT classifiable seasons."""
import importlib.util
import json
from pathlib import Path

import pytest

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location(
    "mois_real_source_season_ceiling",
    ROOT/"scripts/analysis/mois_real_photo_count_seasonal_support_bounds_20261016.py")
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
SOURCE=ROOT/"results/mois_real_photo_count_seasonal_support_bounds_20261016.json"


def fixture():
    return {
        "schema":"fcp_42111_source_photo_opportunity_scientific_receipt_v1",
        "status":"ORIGINAL_PROVIDER_42111_SOURCE_PHOTO_COUNTS_INDEPENDENTLY_VERIFIED_NO_COLOUR_MEASURED",
        "source_snapshot_date":"2026-09-27",
        "original_source_frame_sha256":m.PINNED_FRAME_SHA,
        "original_bundle_sha256":m.PINNED_ARCHIVE_SHA,
        "original_source_photo_count_ledger_sha256":m.PINNED_COUNT_SHA,
        "source_taxa":42111,"taxa_with_at_least_one_licensed_source_photo":42095,
        "taxa_without_licensed_source_photos":16,
        "source_photos_global_unique_licensed":93209351,
        "source_photo_count_run":38749681793,
        "independent_42111_identity_audit_run":38818718877,
        "data_scope_limitations":{
            "source_original_photo_ID_manifests_frozen":False,
            "original_image_pixels_fully_downloaded_and_measured":False
        },
        "species_by_minimum_source_photo_depth":{
            "1":42095,"10":39018,"20":36379,"40":32794,"100":26510
        }
    }


def test_replay_actual_original_fcp_photo_depth_and_upper_bound_restrictions():
    stored=json.loads(SOURCE.read_text())
    calculated=m.derive(fixture())
    assert stored==calculated
    assert stored["source_named_taxa"]==42111
    assert stored["two_seasons_ten_distinct_observations_each"][
        "at_least_species_ruled_out_by_insufficient_photo_counts"]==5732
    assert stored["four_seasons_ten_distinct_observations_each"][
        "at_least_species_ruled_out_by_insufficient_photo_counts"]==9317
    assert stored["two_seasons_ten_distinct_observations_each"][
        "upper_fraction_of_source_taxa"]==pytest.approx(.8639,abs=0.0001)
    assert stored["photographic_opportunities_are_not_unique_observations"] is True
    assert stored["same_species_seasons_actually_observed"] is False
    assert stored["ROI_four_colour_classified_photos_measured"]==0
    assert stored["photo_sampling_inclusion_probs_known_only_after_full_ID_freeze"] is True
    assert stored["not_evidence_of_seasonal_colour_or_weather_relation"] is True


def test_count_only_source_does_not_open_original_photo_id_gate():
    doc=fixture()
    doc["data_scope_limitations"]["source_original_photo_ID_manifests_frozen"]=True
    with pytest.raises(ValueError,match="real verified original"):
        m.derive(doc)
    doc=fixture()
    doc["species_by_minimum_source_photo_depth"]["20"]=36380
    with pytest.raises(ValueError,match="Real photo opportunity depth"):
        m.derive(doc)
