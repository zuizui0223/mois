"""Source-sealed synthetic month-by-month T/S/P/vapour and VPD photo join tests."""
from __future__ import annotations
import csv
import importlib.util
import json
import math
from pathlib import Path
import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[1]
def mod(name,path):
    spec=importlib.util.spec_from_file_location(name,ROOT/path)
    o=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(o)
    return o

fixture=mod("mois_48_fake_climate","tests/test_mois_pin_worldclim_dynamic_monthly_climate_20261016.py")
m=mod("mois_48_dynamic_join","scripts/analysis/mois_join_monthly_dynamic_climate_20261016.py")


def normal(tmp):
    orig,sha=fixture.fixture(tmp)
    root=tmp/"source"
    proof=fixture.m.validate_originals(orig,sha,root,synthetic=True)
    return m.DynamicNormals(root,proof["manifest_sha256"],synthetic=True)


def test_real_units_and_local_month_deviations_with_48_original_source_bands(tmp_path):
    a=normal(tmp_path)
    result=a.site_month("35.7","139.7","2021-08-12")
    assert result["source_climate_join_status"]=="VALID"
    assert result["original_month"]==8
    assert result["monthly_tavg_c"]==pytest.approx(18.6,abs=.01)
    assert result["monthly_prec_mm"]==33
    assert result["source_site_annual_prec_mm"]==pytest.approx(sum(25+x for x in range(1,13)))
    assert result["site_mean_prec_mm"]==pytest.approx(31.5)
    assert result["delta_prec_mm"]==pytest.approx(1.5)
    assert result["monthly_srad_kj_m2_day"]==12800
    assert result["monthly_vapr_kpa"]==pytest.approx(1.28,abs=.0001)
    expected=max(0,float(m.saturation_pressure_kpa(18.6))-1.28)
    assert result["monthly_vpd_kpa"]==pytest.approx(expected,abs=.0001)
    assert result["vpd_raw_negative_clipped_at_observation_month"] is False
    assert result["monthly_rainfall_phase"] in {"ASEASONAL","WET","DRY","TRANSITION"}
    assert len(a.month_sha)==48
    with pytest.raises(ValueError,match="Original measured WorldClim"):
        m.DynamicNormals(tmp_path/"source",a.manifest_sha)


def test_month_unknown_ocean_and_invalid_coordinates_do_not_guess_blooming(tmp_path):
    a=normal(tmp_path)
    missing=a.site_month("35.7","139.7","")
    assert missing["source_climate_join_status"]=="OBSERVATION_DATE_UNKNOWN"
    assert missing["monthly_tavg_c"]==""
    assert missing["site_mean_tavg_c"]!=""
    assert missing["monthly_rainfall_phase"]=="UNKNOWN_MONTH"
    ocean=a.site_month("-80","-175","2021-03-02")
    assert ocean["source_climate_join_status"]=="CLIMATE_UNAVAILABLE"
    assert ocean["monthly_vpd_kpa"]==""
    with pytest.raises(ValueError,match="Half-missing"):
        a.site_month("35.7","","2021-03-02")


def test_vpd_is_derived_with_correct_kilopascals_not_relative_humidity():
    assert m.saturation_pressure_kpa(20)==pytest.approx(2.3383,rel=.001)
    assert m.saturation_pressure_kpa(0)==pytest.approx(.6108)
    # Saturation-vapour difference can become negative due interpolated
    # monthly averages. Truncation happens only for VPD and is flagged.
    assert m.saturation_pressure_kpa(10)<2
    with pytest.raises(ValueError,match="Invalid"):
        m.saturation_pressure_kpa(-240)


def test_sealed_private_source_to_monthly_climate_and_no_biological_model(tmp_path):
    a=normal(tmp_path)
    source=tmp_path/"private_season_photo.csv"
    with source.open("w",newline="") as f:
        cols=["measurement_id","inat_taxon_id","genus","observation_id","observer_id",
              "latitude","longitude","observed_on","morph","measurement_status"]
        w=csv.DictWriter(f,fieldnames=cols);w.writeheader()
        w.writerows([
          {"measurement_id":"FCPR-A","inat_taxon_id":"1","genus":"Alpha",
           "observation_id":"obs1","observer_id":"observer-a","latitude":"35.7",
           "longitude":"139.7","observed_on":"2021-05-15","morph":"white",
           "measurement_status":"classified_four_state_morph"},
          {"measurement_id":"FCPR-B","inat_taxon_id":"1","genus":"Alpha",
           "observation_id":"obs2","observer_id":"observer-b","latitude":"35.7",
           "longitude":"139.7","observed_on":"2021-11-15",
           "morph":"mixed_uncertain","measurement_status":"not_evaluable_roi_or_flip_gate"},
        ])
    receipt=tmp_path/"private_support.json"
    receipt.write_text(json.dumps({
        "schema":"mois_private_photo_seasonal_trait_support_v1",
        "private_derived_photo_rows_sha256":m.sha256_file(source)
    }))
    result=m.join_original_sample(source,receipt,a,tmp_path/"private_join.csv",
                                  tmp_path/"public_receipt.json",synthetic=True)
    assert result["selected_original_photographs"]==2
    assert result["source_four_colour_classifiable_photos"]==1
    assert result["source_classifiable_white_photos"]==1
    assert result["n_successfully_joined_monthly_T_S_P_vapour"]==2
    assert result["seasonal_photo_white_environment_regression_fitted"] is False
    assert result["original_full13_photo_site_env_reconciled"] is False
    rows=list(csv.DictReader((tmp_path/"private_join.csv").open()))
    assert len(rows)==2
    assert float(rows[0]["monthly_prec_mm"])<float(rows[1]["monthly_prec_mm"])
    assert "photo_id" not in (tmp_path/"public_receipt.json").read_text()
    with pytest.raises(ValueError,match="never overwrite"):
        m.join_original_sample(source,receipt,a,tmp_path/"private_join.csv",
                               tmp_path/"public_receipt.json",synthetic=True)


def test_forged_source_identity_and_broken_photo_status_fail(tmp_path):
    a=normal(tmp_path)
    photos=tmp_path/"spoof.csv"
    photos.write_text("measurement_id,inat_taxon_id,genus,observation_id,latitude,longitude,observed_on,morph,measurement_status\nA,1,Alpha,o,35.7,139.7,2020-05-01,white,mixed_uncertain\n")
    receipt=tmp_path/"rec.json"
    receipt.write_text(json.dumps({
        "schema":"mois_private_photo_seasonal_trait_support_v1",
        "private_derived_photo_rows_sha256":m.sha256_file(photos)
    }))
    with pytest.raises(ValueError,match="classification mismatch"):
        m.join_original_sample(photos,receipt,a,tmp_path/"out.csv",
                               tmp_path/"report.json",synthetic=True)
