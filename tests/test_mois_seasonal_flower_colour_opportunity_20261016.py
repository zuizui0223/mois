"""Synthetic-only tests for local-season floral hue and observation-ascertainment.

No original iNaturalist pixels, no independent biological result and no
13-variable climate model are run by this synthetic regression suite.
"""
from __future__ import annotations
import csv
import importlib.util
import json
from pathlib import Path
import pytest

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location(
    "mois_photo_seasons",ROOT/"scripts/analysis/mois_seasonal_flower_colour_opportunity_20261016.py")
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
COLS=[
    "inat_taxon_id","genus","measurement_id","observation_id","observer_id",
    "latitude","longitude","bio12_mm","region_stratum","precip_stratum",
    "observed_on","photo_design_weight","observation_inclusion_probability",
    "measurement_status","morph",
]

def fixture(tmp_path,items):
    joined=tmp_path/"private_joined_sample_flower_colours.csv"
    with joined.open("w",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=COLS)
        writer.writeheader()
        for j,(tax,genus,lat,when,morph) in enumerate(items):
            site={"latitude":str(lat),"longitude":"35","bio12_mm":"800"}
            group=m.classify_site(site)
            writer.writerow({
                "inat_taxon_id":str(tax),"genus":genus,
                "measurement_id":f"FCPR-TEST{j}","observation_id":f"UNIQUE_OBS_{j}",
                "observer_id":f"obs_{j%3}",
                **site,"region_stratum":group[0],"precip_stratum":group[1],
                "observed_on":when,"photo_design_weight":"2",
                "observation_inclusion_probability":"0.5",
                "measurement_status":"classified_four_state_morph" if morph else "image_acquisition_failed",
                "morph":morph or "mixed_uncertain",
            })
    receipt=tmp_path/"sample_collection_receipt.json"
    receipt.write_text(json.dumps({
        "schema":m.RECEIPT_SCHEMA,
        "status":m.RECEIPT_STATUS,
        "private_joined_photo_csv_sha256":m.sha256_file(joined),
        "original_photo_terminals_verified":len(items),
        "original_colour_four_state_classifiable":sum(bool(i[-1]) for i in items),
        "synthetic_only_test":True,
        "full_FCP_93m_photos_measured":False,
        "all13_photo_site_predictors_joined":False,
    }))
    return joined,receipt

def run(tmp,items,**kwargs):
    j,r=fixture(tmp,items)
    return m.audit(
        j,r,tmp/"private_long.csv",tmp/"private_support.json",
        tmp/"public_summary.json",
        min_classifiable_per_season=kwargs.get("min_classifiable",3),
        min_public_cell=kwargs.get("min_public_cell",1),synthetic=True
    )

def test_real_local_season_hemispheres_and_tropics():
    assert m.local_season("40","2020-03-15")[:2]==("NORTH_EXTRATROPICAL","SPRING")
    assert m.local_season("-40","2020-03-15")[:2]==("SOUTH_EXTRATROPICAL","AUTUMN")
    assert m.local_season("40","2020-06-15")[:2]==("NORTH_EXTRATROPICAL","SUMMER")
    assert m.local_season("-40","2020-06-15")[:2]==("SOUTH_EXTRATROPICAL","WINTER")
    assert m.local_season("-40","2020-12-16")[:2]==("SOUTH_EXTRATROPICAL","SUMMER")
    assert m.local_season("-40","2020-10-03")[:2]==("SOUTH_EXTRATROPICAL","SPRING")
    assert m.local_season("10","2020-06-15")[:2]==("TROPICAL","CALENDAR_Q2")
    assert m.local_season("10","2020-12-15")[:2]==("TROPICAL","CALENDAR_Q4")
    assert "WET" not in m.local_season("10","2020-06-15")[1]
    assert m.local_season("","2020-06-15")[1]=="UNKNOWN_SEASON"
    assert m.local_season("40","not-a-full-date")[1]=="UNKNOWN_SEASON"
    assert m.local_season("40","2020-08-15T11:00:00Z")[1]=="SUMMER"
    with pytest.raises(ValueError,match="exceeds immutable"):
        m.local_season("40","2026-10-01")

def test_same_species_seasonal_photo_white_with_fixed_region_and_bio12(tmp_path):
    items=[]
    # Northern species: white spring, yellow summer.
    items += [(1,"Alpha",40,"2020-04-12","white")]*3
    items += [(1,"Alpha",40,"2020-07-14","yellow_orange")]*3
    # Southern species: spring is October, summer is January.
    items += [(2,"Beta",-40,"2020-10-12","red_pink")]*3
    items += [(2,"Beta",-40,"2020-01-14","white")]*3
    report=run(tmp_path,items)
    assert report["selected_source_observations"]==12
    assert report["source_classifiable_photographs"]==12
    assert report["n_matched_species_season_pairs_before_public_suppression"]==2
    assert report["n_stricter_matched_species_region_precip_season_pairs_before_public_suppression"]==2
    pairs={(x["domain"],x["first_season"],x["second_season"]):x for x in report["same_species_season_matched"]}
    north=pairs[("NORTH_EXTRATROPICAL","SPRING","SUMMER")]
    south=pairs[("SOUTH_EXTRATROPICAL","SPRING","SUMMER")]
    assert north["equal_unit_mean_white_change_second_minus_first"]==-1
    assert south["equal_unit_mean_white_change_second_minus_first"]==1
    assert len(report["seasonal_observation_photo_white_conditioned"])==4
    assert report["all_13_photo_site_environment_predictors_available_here"] is False
    assert report["old_genus_13_environment_estimates_updated"] is False
    pub=(tmp_path/"public_summary.json").read_text()
    assert "UNIQUE_OBS_0" not in pub and "FCPR-TEST0" not in pub
    assert "obs_0" not in pub and '"inat_taxon_id"' not in pub
    prv=list(csv.DictReader((tmp_path/"private_long.csv").open()))
    assert len(prv)==12
    assert set(prv[0])==set(m.PRIVATE_FIELDS)
    with pytest.raises(ValueError,match="overwrite"):
        run(tmp_path,items)

def test_species_turnover_must_not_be_called_colour_change(tmp_path):
    items=(
        [(1,"Alpha",40,"2020-04-12","yellow_orange")]*3+
        [(2,"Beta",40,"2020-04-12","white")]*3+
        [(2,"Beta",40,"2020-07-12","white")]*3+
        [(3,"Gamma",40,"2020-07-12","white")]*3
    )
    report=run(tmp_path,items)
    groups={x["local_season"]:x for x in report["seasonal_observation_photo_white_conditioned"]}
    assert groups["SPRING"]["observation_weighted_white_given_classifiable"]==.5
    assert groups["SUMMER"]["observation_weighted_white_given_classifiable"]==1.0
    pairs=report["same_species_season_matched"]
    assert len(pairs)==1 and pairs[0]["n_matched_units"]==1
    assert pairs[0]["equal_unit_mean_white_change_second_minus_first"]==0.0
    assert report["n_stricter_matched_species_region_precip_season_pairs_before_public_suppression"]==1

def test_classifiability_unknown_dates_tropical_quarter_and_suppression(tmp_path):
    items=[
        (1,"Alpha",10,"2020-06-01","white"),
        (1,"Alpha",10,"2020-06-10",None),
        (1,"Alpha",10,"",None),
        (2,"Beta",40,"2020-05-17","blue_purple"),
    ]
    report=run(tmp_path,items,min_public_cell=5)
    assert report["source_classifiable_photographs"]==2
    assert report["n_selected_with_unknown_observation_date"]==1
    assert report["tropical_season_method"].startswith("calendar_month")
    assert len(report["seasonal_observation_photo_white_conditioned"])==3
    assert all("privacy_suppressed_count_below" in x for x in report["seasonal_observation_photo_white_conditioned"])
    assert report["n_matched_species_season_pairs_before_public_suppression"]==0

def test_modified_source_receipt_and_wrong_sample_count_abort(tmp_path):
    items=[(1,"Alpha",40,"2020-04-15","white")]
    joined,receipt=fixture(tmp_path,items)
    bad=json.loads(receipt.read_text())
    bad["private_joined_photo_csv_sha256"]="0"*64
    receipt.write_text(json.dumps(bad))
    with pytest.raises(ValueError,match="complete hash-bound"):
        m.audit(joined,receipt,tmp_path/"a.csv",tmp_path/"b.json",
                tmp_path/"c.json",synthetic=True)
    bad["private_joined_photo_csv_sha256"]=m.sha256_file(joined)
    bad["original_photo_terminals_verified"]=2
    receipt.write_text(json.dumps(bad))
    with pytest.raises(ValueError,match="differ"):
        m.audit(joined,receipt,tmp_path/"a.csv",tmp_path/"b.json",
                tmp_path/"c.json",synthetic=True)
    assert not (tmp_path/"a.csv").exists()

def test_missing_or_forged_annotation_does_not_create_flowering_claim(tmp_path):
    items=[(1,"Alpha",40,"2020-04-15","white")]
    joined,receipt=fixture(tmp_path,items)
    overlay=tmp_path/"annotation.csv"
    with overlay.open("w",newline="") as f:
        writer=csv.writer(f)
        writer.writerow(["observation_id","annotation_term_id","annotation_value_ids",
                          "lookup_status","annotation_source","lookup_snapshot_date"])
        writer.writerow(["UNIQUE_OBS_0","12","13","LOOKED_UP","inat_export","2026-10-16"])
    data=m.audit(joined,receipt,tmp_path/"a.csv",tmp_path/"b.json",
        tmp_path/"c.json",annotation_overlay=overlay,
        annotation_sha256=m.sha256_file(overlay),synthetic=True,min_public_cell=1)
    assert data["annotation_status_selected_counts"]=={"FLOWERS":1}
    assert data["synthetic_test_only"] is True
    assert data["causal_seasonal_phenotype_or_wild_flora_incidence_identified"] is False
