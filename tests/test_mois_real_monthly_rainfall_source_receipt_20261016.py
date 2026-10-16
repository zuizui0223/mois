"""Guard only the public real WorldClim source receipt; never fake biology.

The source ZIP digest and raster/manifest digest were observed in successful
original-provider-byte GitHub Actions #38834507200, artifact #12018706765.
This does not re-download the original bytes in ordinary Mois CI.
"""
from __future__ import annotations
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
PATH=ROOT/"results/mois_worldclim21_monthly_10m_source_receipt_20261016.json"

def test_actual_worldclim_source_is_real_and_no_biology_promoted():
    x=json.loads(PATH.read_text(encoding="utf-8"))
    assert x["schema"]=="mois_original_worldclim21_10m_monthly_rainfall_source_verification_v1"
    assert x["status"]=="ACTUAL_1970_2000_MONTHLY_WORLDCLIM_RASTERS_BYTE_CHECKED_NOT_PHOTOS"
    assert x["provider_archive_sha256"]=="1090f578f6402672a2b6438e252156f6cf867a25002a8ad6ffcbfa19e516ddc6"
    assert x["generated_monthly_manifest_sha256"]=="07f45ac2973589d8951ac3eb0fa1c81b700018c5ab7322e2346dcf7d6723ee23"
    assert x["verified_monthly_raster_count"]==12
    assert x["source_name"].startswith("WorldClim 2.1 historical 1970–2000")
    assert x["github_execution_url"]=="https://github.com/zuizui0223/mois/actions/runs/38834507200"
    assert x["github_artifact_id"]==12018706765
    assert x["provider_zip_contains_readme_txt"] is True
    assert x["original_FCP_42111_photo_ID_census_complete"] is False
    assert x["annual_BIO12_original_photo_site_reconciliation_complete"] is False
    assert x["new_classifiable_flower_photos"]==0
    assert x["season_by_flower_colour_new_coefficients_computed"] is False
    assert x["thirteen_variable_macro_models_rerun_on_new_source"] is False
    assert x["original_raster_bytes_deposited_in_GitHub"] is False
