"""Source-aggregate safeguards for new photographed colour latitude/longitude/elevation analysis.

No original licensed photos, site coordinates, observation IDs or taxon names
are committed. Full original source replay is an owner-only separate check.
"""
from __future__ import annotations
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
REPORT=ROOT/"results/mois_photo_colour_geography_synthesis_20261015.json"
LAT=ROOT/"results/source_receipts/photo_white_latitudinal_null.json"
MACRO=ROOT/"results/fcp_macro_direction_20261013/decision_receipt.json"
SI=ROOT/"docs/MOIS_GEB_SUPPORTING_INFORMATION_20261014.md"


def data(p):return json.loads(p.read_text(encoding="utf8"))


def test_source_photo_scope_and_three_geography_bands():
    q=data(REPORT); l=data(LAT)
    assert q["schema"]=="mois_flower_white_original_geography_synthesis_v1"
    assert q["source_original_photo_sha256"]==l["source_sha256"]
    assert q["source_geotagged_all"]==100348
    assert q["source_geotagged_classified"]==l["n_original_classified_geotagged_photos"]==45879
    assert q["primary_adjusted_species"]==20694
    assert q["primary_adjusted_genera"]==1934
    assert q["primary_adjusted_classified_photos"]==38498
    raws=q["raw_geographic_classified_photo_white"]
    for category in ("absolute_latitude","longitude_90deg_sector"):
        rows=raws[category]
        assert sum(x["all_photos"] for x in rows)==100348
        assert sum(x["classifiable"] for x in rows)==45879
        assert sum(x["white"] for x in rows)==21550
        for r in rows:
            assert math.isclose(r["pct_white_of_classifiable"],100*r["white"]/r["classifiable"],abs_tol=.0001)
    assert [r["classifiable"] for r in raws["absolute_latitude"]]==l["original_latitude_band_n_photos"]
    assert raws["longitude_90deg_sector"][3]["pct_white_of_classifiable"]>49
    quartiles=raws["elevation_photo_opportunity_quartile"]
    assert sum(x["all_photos"] for x in quartiles)==100255
    assert quartiles[0]["pct_white_of_classifiable"]>quartiles[-1]["pct_white_of_classifiable"]
    assert [r["range_elevation_m"][0] for r in quartiles]==[-330,94,286,929]


def test_congeneric_geography_adjustment_reproduces_original_moisture_cohort():
    q=data(REPORT);macro=data(MACRO)
    assert abs(q["original_primary_white_BIO12_adjusted_pp_per_source_sd_reproduced"]-2.254545235123587)<1e-9
    assert macro["original_photo_opportunities"]==q["source_geotagged_all"]+195
    ms=q["genus_equal_20694_species_geographically_adjusted"]
    assert [x["model"] for x in ms]==[
        "geography_only","geography_plus_rain","geography_plus_vapour","geography_both_moisture"]
    assert q["source_bootstrap"]["reps"]==1999
    assert not q["source_bootstrap"]["across_models_or_research_program_adjustment"]
    assert ms[0]["white_vs_1000m_higher_elevation_pp"]<0
    assert ms[0]["within_model_geography4_maxT95_elevation_pp"][1]<0
    assert ms[1]["within_model_geography4_maxT95_elevation_pp"][0]<0<ms[1]["within_model_geography4_maxT95_elevation_pp"][1]
    assert ms[2]["within_model_geography4_maxT95_elevation_pp"][0]<0<ms[2]["within_model_geography4_maxT95_elevation_pp"][1]
    assert ms[2]["within_model_geography4_maxT95_lat_pp"][0]>0
    assert ms[3]["within_model_geography4_maxT95_lat_pp"][0]>0
    for m in ms:
        assert len(m["within_model_geography4_maxT95_lat_pp"])==2
        assert len(m["within_model_geography4_maxT95_elevation_pp"])==2
        assert 0<m["longitude_sin_cos_joint_approx_chisq2_p"]<1


def test_si_distinguishes_unadjusted_opportunity_and_climate_cohort():
    s=SI.read_text(encoding="utf8")
    for phrase in ("### S2.9 |","100,348","45,879","20,694","1,934",
                   "49.49%","42.96%","−1.903","−3.493","3.083","+0.600"):
        assert phrase in s,phrase
    assert "not causally linked" in s
    assert "not a monotone eastward ecological trend" in s
    assert "photographic classifiability" in s
    assert "not spatial-autocorrelation-corrected" in s
    report=REPORT.read_text(encoding="utf8")
    for sensitive in ('"photo_id"','"inat_taxon_id"','"latitude"','"longitude"'):
        assert sensitive not in report
