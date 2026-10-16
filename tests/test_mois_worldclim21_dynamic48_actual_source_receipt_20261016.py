"""Actual 48 original WorldClim monthly T/S/P/vapour source receipt integrity."""
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
RECEIPT=ROOT/"results/mois_worldclim21_dynamic_48_real_source_receipt_20261016.json"

def test_actual_external_source_verified_four_variables_forty_eight_months():
    s=json.loads(RECEIPT.read_text(encoding="utf8"))
    assert s["schema"]=="mois_source_attested_worldclim21_original_dynamic_monthlies_v1"
    assert s["status"]=="VERIFIED_48_REAL_MONTHLY_GEOSPATIAL_CLIMATE_NORMALS_NO_PLANT_PHOTOS"
    assert s["original_execution_url"]=="https://github.com/zuizui0223/mois/actions/runs/38848205285"
    assert s["original_48_raster_manifest_sha256"]=="9d0cf37a64534542d7788debef28c13f106a15c3bf7bf91be0653f3ce106baae"
    assert s["n_source_variables"]==4
    assert s["n_original_monthly_rasters"]==48
    assert s["original_grid_n_columns"]==2160
    assert s["original_grid_n_rows"]==1080
    assert set(s["source_original_zip_sha256"])=={"prec","srad","tavg","vapr"}
    assert s["source_original_zip_sha256"]["prec"]=="1090f578f6402672a2b6438e252156f6cf867a25002a8ad6ffcbfa19e516ddc6"
    assert s["source_original_zip_sha256"]["srad"]=="c72ee7f4f9a0eb4b5f6cd7a003eddc05bda22a2e5666d968ed4e817fd36b9026"
    assert s["source_original_zip_sha256"]["tavg"]=="5e567dcfe94379b94229492849ce91078b1c6e5210aaf435fba449fae6b95405"
    assert s["source_original_zip_sha256"]["vapr"]=="0501b5b5cbbc6e2fc2a3801efa9a126a5f9669f0f9a00642d2519c2324ce726d"
    assert all(len(x)==64 for x in s["source_original_zip_sha256"].values())
    assert s["units"]=={"prec":"mm/month","srad":"kJ m-2 day-1","tavg":"degC","vapr":"kPa"}
    assert s["original_42111_taxon_photo_ids_complete"] is False
    assert s["photographed_flower_pixel_colours_joined"] is False
    assert s["complete_13_source_environment_reconciliation"] is False
    assert s["seasonal_white_climate_relation_estimated"] is False
    assert s["climate_is_1970_2000_normal_not_observed_photo_year_weather"] is True
