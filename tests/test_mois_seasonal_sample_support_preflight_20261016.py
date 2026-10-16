"""Synthetic-only pre-ROI support bounds, not a source flower-colour inference."""
import csv
import importlib.util
import json
from pathlib import Path

import pytest

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location(
    "mois_pre_roi_season_support",
    ROOT/"scripts/analysis/mois_seasonal_sample_support_preflight_20261016.py")
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def case(tmp):
    selected=tmp/"selected.csv"
    frames=[]
    for tid,genus,latitude,months,n in [
        ("1","Alpha","40",[(4,"SPRING"),(7,"SUMMER")],3),
        ("2","Beta","-40",[(10,"SPRING"),(1,"SUMMER")],3),
        ("3","Gamma","10",[(6,"CALENDAR_Q2")],1)
    ]:
        for month,_ in months:
            for j in range(n):
                row=dict.fromkeys(m.PRIVATE_FIELDS,"")
                i=len(frames)
                row.update({
                    "inat_taxon_id":tid,"genus":genus,
                    "photo_id":str(1000+i),"observation_id":f"uuid-{i}",
                    "observer_id":f"observer_{j}",
                    "observed_on":f"2020-{month:02d}-12",
                    "latitude":latitude,"longitude":"35","bio12_mm":"800",
                    "photo_license":"CC0","technical_pilot_flag":"0",
                    "observation_inclusion_probability":"1",
                    "photo_design_weight":"1",
                })
                region,precip=m.classify_site(row)
                row["region_stratum"]=region
                row["precip_stratum"]=precip
                frames.append(row)
    with selected.open("w",newline="",encoding="utf8") as f:
        w=csv.DictWriter(f,fieldnames=m.PRIVATE_FIELDS)
        w.writeheader()
        w.writerows(frames)
    receipt=tmp/"receipt.json"
    receipt.write_text(json.dumps({
        "schema":"mois_direct_fcp_source_stratified_photo_selection_v1",
        "status":"SOURCE_ADMITTED_SELECTED_PRIVATE_IDS_ONLY_NO_PHOTO_PIXELS_MEASURED",
        "selected_private_photo_csv_sha256":m.sha256_file(selected),
        "selection_photo_cap_per_taxon":40,
        "historical_colour_outcomes_seen":False,
        "pixels_opened_or_colours_measured":False,
        "original_13var_environment_join_complete":False,
        "source_species":4,
        "source_species_with_photos":3,
        "source_photo_ids_authenticated":19,
        "selected_original_photo_ids":len(frames),
    }))
    return selected,receipt,frames


def test_three_obs_per_local_season_are_only_preclassification_upper_bounds(tmp_path):
    sample,receipt,_=case(tmp_path)
    results=m.audit(sample,receipt,tmp_path/"private.json",tmp_path/"public.json",
                    synthetic=True)
    assert results["preclassification_selected_unique_observations"]==13
    assert results["selected_distinct_taxa"]==3
    assert results["n_source_taxon_with_zero_photos"]==1
    assert results["real_classifiable_photos_measured"]==0
    assert results["matched_counts_are_only_upper_bounds_before_flower_roi_classification"] is True
    pairs=results["matched_taxon_season_support_upper_bounds"]["3"]
    assert len(pairs)==2
    assert {item["domain"]:item["n_nominal_comparisons"] for item in pairs}=={
       "NORTH_EXTRATROPICAL":1,"SOUTH_EXTRATROPICAL":1
    }
    strict=results["matched_taxon_same_geographic_precip_stratum_support_upper_bounds"]["3"]
    assert sum(x["n_nominal_comparisons"] for x in strict)==2
    assert results["matched_taxon_season_support_upper_bounds"]["5"]==[]
    assert results["matched_taxon_season_support_upper_bounds"]["7"]==[]
    assert results["synthetic_test_only"] is True
    assert "uuid-" not in (tmp_path/"public.json").read_text()
    assert '"photo_id"' not in (tmp_path/"public.json").read_text()
    with pytest.raises(ValueError,match="overwrite"):
        m.audit(sample,receipt,tmp_path/"private.json",tmp_path/"public.json",synthetic=True)


def test_missing_month_does_not_count_as_eligible_season(tmp_path):
    selected,receipt,records=case(tmp_path)
    for row in records:
        if row["inat_taxon_id"]=="1" and row["observed_on"].startswith("2020-07"):
            row["observed_on"]=""
    with selected.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=m.PRIVATE_FIELDS);w.writeheader();w.writerows(records)
    doc=json.loads(receipt.read_text())
    doc["selected_private_photo_csv_sha256"]=m.sha256_file(selected)
    receipt.write_text(json.dumps(doc))
    x=m.audit(selected,receipt,tmp_path/"private.json",tmp_path/"public.json",
              synthetic=True)
    assert x["selected_missing_observation_month"]==3
    assert len(x["matched_taxon_season_support_upper_bounds"]["3"])==1
    assert x["matched_taxon_season_support_upper_bounds"]["3"][0]["domain"]=="SOUTH_EXTRATROPICAL"


def test_invalid_photo_identity_and_source_receipt_fail_closed(tmp_path):
    selected,receipt,records=case(tmp_path)
    records[-1]["photo_id"]=records[0]["photo_id"]
    with selected.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=m.PRIVATE_FIELDS);w.writeheader();w.writerows(records)
    doc=json.loads(receipt.read_text())
    doc["selected_private_photo_csv_sha256"]=m.sha256_file(selected)
    receipt.write_text(json.dumps(doc))
    with pytest.raises(ValueError,match="Duplicate/missing"):
        m.audit(selected,receipt,tmp_path/"a.json",tmp_path/"b.json",synthetic=True)
    assert not (tmp_path/"b.json").exists()


def test_synthetic_small_source_cannot_be_misrepresented_as_real(tmp_path):
    selected,receipt,_=case(tmp_path)
    with pytest.raises(ValueError,match="original real 42,111"):
        m.audit(selected,receipt,tmp_path/"a.json",tmp_path/"b.json",synthetic=False)
