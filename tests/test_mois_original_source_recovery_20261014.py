"""Verify source-recovered original ZIP audit without ever shipping raw photo-site data.

The aggregate JSON was computed directly from both exact original FCP Actions ZIPs
in an independent environment; these CI tests check its internal consistency and
agreement with the frozen 1,280 source coefficient ledger, not a new ZIP rerun.
"""
from __future__ import annotations
import csv
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
RECEIPT=ROOT/"results/mois_original_source_independent_recovery_audit_20261014.json"
COEFF=ROOT/"results/fcp_macro_direction_20261013/direction_ledger.csv"
SCRIPT=ROOT/"scripts/analysis/verify_mois_original_source_two_axis_wls_20261014.py"


def receipt():
    return json.loads(RECEIPT.read_text(encoding="utf8"))


def test_original_zip_recovery_and_sample_denominators():
    r=receipt()
    assert r["status"]=="PASS_INDEPENDENT_ORIGINAL_ZIP_JOIN_AND_EIGHT_WLS_COEFFICIENTS"
    assert r["photo_source_sha256"]=="2d6097ba1b0434632e13bbe79dd78180608a287ff5dabebc770b5c0a2cfb52a0"
    assert r["environment_source_sha256"]=="6e287294b3a05eaa383be3de644551825915bfecca9c4bde6950b4890f5e16ca"
    assert (r["source_photo_opportunities"],r["named_taxa"],r["classifiable_photo_hues"])==(100543,42111,45953)
    assert (r["climate26_classifiable_photo_rows"],r["climate26_source_species"])==(45836,25028)
    assert (r["soil36_classifiable_photo_rows"],r["soil36_source_species"])==(31925,20786)
    assert (r["genus_eligible_species"],r["genus_eligible_genera"])==(20694,1934)
    assert r["point_estimates_only"] and not r["bootstrap_CIs_recomputed"]


def test_eight_independent_wls_slopes_match_original_source_ledger():
    r=receipt()
    with COEFF.open(newline="",encoding="utf8") as f:
        ledger=list(csv.DictReader(f))
    assert len(ledger)==1280
    seen=0
    for axis,vals in r["eight_point_slopes_pp_per_source_sd"].items():
        assert axis in ("wc_bio12","wc_vapr_annual_mean")
        assert set(vals)=={"white","yellow_orange","red_pink","blue_purple"}
        assert abs(sum(vals.values()))<1e-10
        for hue,point in vals.items():
            q=[x for x in ledger if all((
                x["mode"]=="climate",x["quality"]=="all",x["region"]=="world",
                x["adjustment"]=="genus_geo",x["weight"]=="genus_equal",
                x["anchor"]==axis,x["hue"]==hue))]
            assert len(q)==1
            assert abs(float(q[0]["direction_pp_per_env_global_sd"])-point)<1e-9
            assert int(q[0]["n_species"])==20694
            seen+=1
    assert seen==8
    assert r["max_abs_difference_from_original_1280_ledger_pp"]<1e-9


def test_source_permissions_and_geography_are_aggregates_not_individual_records():
    r=receipt()
    licences=r["licence_counts"];site=r["source_geo_status_counts"]
    assert sum(licences.values())==100543
    assert sum(site.values())==100543
    assert licences["cc-by-nc"]==78367
    assert licences["UNSPECIFIED"]==2
    assert licences["cc-by-nc-nd"]==2
    assert site["VALID_PUBLIC_ORIGINAL_PHOTO_POINT"]==100348
    assert site["NO_UNOBSCURED_PUBLIC_COORDINATES"]==31
    assert not r["per_photo_ids_coordinates_and_attribution_released"]
    assert not r["all_original_photo_pixels_present"]
    assert r["source_rights_and_geoprivacy_human_review_required"] is True
    assert "raw_inputs/" in (ROOT/".gitignore").read_text(encoding="utf8")
    code=SCRIPT.read_text(encoding="utf8")
    assert "sha256" in code.lower() and "genus" in code.lower()
    assert "photo_path" in code and "environment_path" in code
