"""Synthetic full original photo/date-verified blooming to all 48 climate axes.

Proves original image-positive photograph identity is not dropped, multiple
images of ONE original iNaturalist observation give ONE independent observation
vote, and species-equal blooming composition remains the primary estimand.
This makes no claim that 93M source images have already been measured.
"""
from __future__ import annotations
import csv
import importlib.util
import json
from pathlib import Path

import pytest

ROOT=Path(__file__).resolve().parents[1]


def module(name,path):
    spec=importlib.util.spec_from_file_location(name,ROOT/path)
    x=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(x)
    return x


m=module("mois_full_observation_first_monthly_climate",
         "scripts/analysis/mois_full_bloom_monthly_TSPVPD_observation_species_FE_20261016.py")
source_test=module("mois_full_observed_flower_fixture",
              "tests/test_mois_original_full_fcp_dated_flower_ROI_20261016.py")
climate_test=module("mois_full_observed_climate_fixture",
               "tests/test_mois_join_monthly_dynamic_climate_20261016.py")


def prep(tmp):
    climate=climate_test.normal(tmp)
    photo,receipt,method,rows=source_test.fixture(tmp)
    positive=tmp/"full_original_classified_dated_blooms.ndjson"
    proof=tmp/"full_original_bloom_receipt.json"
    original=source_test.m.analyze(photo,receipt,method,positive,proof,synthetic=True)
    return climate,positive,proof,original


def test_end_to_end_real_FCP_format_positive_flower_to_observation_equal_site_month(tmp_path):
    climate,positive,proof,original=prep(tmp_path)
    private=tmp_path/"private_species_cell_month.csv"
    public=tmp_path/"public_monthly_white_TSPVPD.json"
    r=m.evaluate(positive,proof,climate,private,public,synthetic=True,min_repeat_cells=2)
    assert r["status"]==m.STATUS
    assert r["n_original_FCP_all_photo_opportunities_processed_with_ROI"]==40
    assert r["n_original_bloom_and_source_date_qualified_photo_IDs"]==original[
        "n_original_date_and_site_known_ROI_bloom_hues"]
    assert r["n_original_white_bloom_date_known_photos"]==original[
        "n_original_date_known_white_flower_photos"]
    assert r["n_distinct_original_blooming_photo_observation_UUIDs"]==r[
        "n_original_bloom_and_source_date_qualified_photo_IDs"]
    assert r["n_bloom_photos_monthly_worldclim_matched"]==r[
        "n_original_bloom_and_source_date_qualified_photo_IDs"]
    assert r["n_bloom_photos_no_monthly_worldclim"]==0
    assert r["n_original_named_taxa_with_one_or_more_source_date_known_bloom_photos"]==4
    assert r["n_original_genera_with_date_known_classifiable_bloom"]==4
    assert r["n_source_climate_cells_with_any_classified_original_bloom"]==2
    assert r["n_climate_cells_with_multiple_distinct_blooming_source_months"]==1
    assert len(r["three_geographic_month_effect_model_support_and_exploratory_estimates"])==3
    assert all(x["status"].startswith("NOT_ESTIMABLE") for x in r[
       "three_geographic_month_effect_model_support_and_exploratory_estimates"])
    assert r["within_observation_multiple_photos_dependent_and_aggregated_before_species"] is True
    assert r["source_photo_sampling_cap_per_species"] is None
    assert r["monthly_weather_is_1970_2000_climatological_not_photo_year"] is True
    summary=r["real_original_source_calendar_month_by_domain_photo_and_species_white"]
    # Source privacy suppression applies to the small 40-photo fake
    # cohort: every latitude-zone/month has fewer than ten flowers.
    # The total source-locked denominator remains complete privately.
    assert all(x["public_suppressed_below"]==10 for x in summary)
    assert r["n_original_bloom_and_source_date_qualified_photo_IDs"]==22
    rows=list(csv.DictReader(private.open()))
    assert rows
    assert all(int(z["n_original_independent_observations"])<=int(
        z["n_original_classifiable_photographs"]) for z in rows)
    assert "original_photo_id" not in public.read_text()
    assert "original_observation_uuid" not in public.read_text()
    assert "original_observer" not in public.read_text()
    with pytest.raises(ValueError,match="overwrite"):
        m.evaluate(positive,proof,climate,private,public,synthetic=True,min_repeat_cells=2)


def test_two_images_of_one_source_observation_only_count_one_independent_vote(tmp_path):
    climate,positive,receipt,_=prep(tmp_path)
    rows=[json.loads(x) for x in positive.read_text().splitlines()]
    # Change a fictitious second image's UUID, observer, taxon and date so
    # it is legitimately another image of the first source observation.
    r0=rows[0]
    r1=next(r for r in rows[1:] if (
        r["original_taxon_id"]==r0["original_taxon_id"]
        and r["original_observed_on"]==r0["original_observed_on"]
        and r["latitude"]==r0["latitude"]
        and r["longitude"]==r0["longitude"]))
    r1["original_observation_uuid"]=r0["original_observation_uuid"]
    r1["original_observer"]=r0["original_observer"]
    positive.write_text("\n".join(json.dumps(x) for x in rows)+"\n")
    d=json.loads(receipt.read_text())
    d["private_actual_bloom_photo_month_hue_manifest_sha256"]=m.sha(positive)
    receipt.write_text(json.dumps(d))
    result=m.evaluate(positive,receipt,climate,tmp_path/"private.csv",
                      tmp_path/"out.json",synthetic=True,min_repeat_cells=2)
    assert result["n_distinct_original_blooming_photo_observation_UUIDs"]==(
        len(rows)-1)
    assert result["n_original_bloom_and_source_date_qualified_photo_IDs"]==len(rows)
    private_rows=list(csv.DictReader((tmp_path/"private.csv").open()))
    assert sum(int(row["n_original_independent_observations"])
               for row in private_rows)==len(rows)-1
    assert sum(int(row["n_original_classifiable_photographs"])
               for row in private_rows)==len(rows)


def test_inconsistent_taxon_for_same_source_observation_uuid_refused(tmp_path):
    climate,positive,receipt,_=prep(tmp_path)
    rows=[json.loads(x) for x in positive.read_text().splitlines()]
    rows[1]["original_observation_uuid"]=rows[0]["original_observation_uuid"]
    # Original observation cannot have two unrelated source taxa / observers.
    positive.write_text("\n".join(json.dumps(x) for x in rows)+"\n")
    d=json.loads(receipt.read_text())
    d["private_actual_bloom_photo_month_hue_manifest_sha256"]=m.sha(positive)
    receipt.write_text(json.dumps(d))
    with pytest.raises(ValueError,match="changed source taxon/site/date"):
        m.evaluate(positive,receipt,climate,tmp_path/"private.csv",
                   tmp_path/"public.json",synthetic=True)
    assert not (tmp_path/"public.json").exists()


def test_strict_original_source_sha_prevents_renaming_nonflowers_as_white(tmp_path):
    climate,positive,receipt,_=prep(tmp_path)
    lines=positive.read_text().splitlines()
    first=json.loads(lines[0])
    first["original_ROI_v4_hue"]="white" if first["original_ROI_v4_hue"]!="white" else "blue_purple"
    lines[0]=json.dumps(first)
    positive.write_text("\n".join(lines)+"\n")
    with pytest.raises(ValueError,match="Unsealed all-photo"):
        m.evaluate(positive,receipt,climate,tmp_path/"priv.csv",
                   tmp_path/"pub.json",synthetic=True)
    assert not (tmp_path/"pub.json").exists()


def test_monthly_raster_LRU_does_not_expand_to_93million_photo_count(tmp_path):
    weather=climate_test.normal(tmp_path)
    weather.max_cache_cells=1
    a=weather.profile((1,1))
    assert a is not None and len(weather.cache)==1
    b=weather.profile((2,2))
    assert b is not None and len(weather.cache)==1
    assert (2,2) in weather.cache and (1,1) not in weather.cache
    assert weather.profile((1,1))["prec"].shape==(12,)
    assert len(weather.cache)==1
