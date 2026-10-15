"""Source-locked ecological claim guards for soil-alone vs water/climate blocks.

The original photo/site ZIPs are private; public CI validates aggregate-only
source reproduction and synthetic equality/contrast semantics.
"""
from __future__ import annotations
import importlib.util,json,math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
R=ROOT/"results/mois_soil_crossdomain_block_audit_20261015.json"
F=ROOT/"results/mois_full26_climate_soil10_transfer_20261015.json"
E=ROOT/"results/mois_biogeography_climate_soil_elevation_full_audit_20261015.json"
SRC=ROOT/"scripts/analysis/mois_soil_vs_climate_blocks_20261015.py"
SI=ROOT/"docs/MOIS_GEB_SUPPORTING_INFORMATION_20261014.md"

def result(p):return json.loads(p.read_text(encoding="utf8"))

def test_original_source_identity_and_fixed_moisture_comparator():
    x=result(R);original=result(F);b=result(E)
    assert x["schema"]=="mois_soil_crossdomain_block_audit_v1"
    assert x["source_sha256"]["photo"]==original["photo_sha256"]
    assert x["source_sha256"]["environment"]==original["environment_sha256"]
    assert x["population"]=={"n_species":16801,"n_genera":1665,"n_original_classifiable_photos":26163}
    assert x["protocol"]["fold_species_counts"]==[3212,3505,3803,3306,2975]
    assert x["protocol"]["genus_holdout"]
    assert x["protocol"]["bootstrap_draws"]==1999
    a=x["by_alpha"]["0.01"]["models"]
    assert math.isclose(a["geography_water"]["gain_vs_geography"],
         b["heldout_genus_cv"]["models"]["geography_moisture"]["gain_vs_geography"],abs_tol=1e-9)
    for source_key,new_key in (
        ("geo_plus_all26_climate","geography_full26climate"),
        ("geo_plus_all26_climate_plus_soil5","geography_full26climate_soil5"),
        ("geo_plus_all26_climate_plus_soil10","geography_full26climate_soil10")):
        assert math.isclose(a[new_key]["loss"],
            original["model_results"][source_key]["heldout_species_photo_white_squared_error"],abs_tol=1e-9)

def test_key_predictions_comparable_on_identical_taxa_and_95_uncertainty():
    x=result(R)
    for alpha in ("0.001","0.01","0.1"):
        s=x["by_alpha"][alpha];m=s["models"];d=s["contrasts"]
        assert all(m[name]["gain_vs_geography"]>0 for name in
            ("geography_soil5","geography_soil10","geography_water",
             "geography_other24climate","geography_full26climate"))
        assert m["geography_water"]["loss"]<m["geography_soil5"]["loss"]
        assert m["geography_water"]["loss"]<m["geography_soil10"]["loss"]
        assert d["water_vs_soil5"]["gain"]>0
        assert d["water_vs_soil5"]["genus_bootstrap95"][0]>0
        assert d["water_vs_soil10"]["genus_bootstrap95"][0]>0
        other=d["increment_water_given_other24climate"]
        assert 0<other["gain"]<0.00005
        assert other["genus_bootstrap95"][0]<0<other["genus_bootstrap95"][1]
        for soil in ("increment_soil5_given_full26","increment_soil10_given_full26"):
            assert d[soil]["gain"]<0
            assert d[soil]["genus_bootstrap95"][1]<0
        assert math.isclose(d["water_vs_soil5"]["gain"],
            m["geography_soil5"]["loss"]-m["geography_water"]["loss"],abs_tol=1e-8)
        assert math.isclose(d["increment_water_given_other24climate"]["gain"],
            m["geography_other24climate"]["loss"]-m["geography_full26climate"]["loss"],abs_tol=1e-8)
    primary=x["by_alpha"]["0.01"]
    assert primary["models"]["geography_soil5"]["gain_ci95"][0]<0
    assert primary["models"]["geography_water"]["gain_ci95"][0]>0
    assert not x["no_raw_rows_exported"] is False
    assert "No original DEM neighbourhoods" in x["terrain"]

def test_source_replayer_synthetic_comparison_checker():
    spec=importlib.util.spec_from_file_location("soil_crossdomain_replayer",SRC)
    mod=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    mod.check_numeric({"a":1.0,"b":[1.0000000001,2]},{"a":1,"b":[1,2]})
    try:mod.check_numeric({"gain":0.1},{"gain":0.2})
    except AssertionError:pass
    else:raise AssertionError("Changed source-model numerical result silently accepted")
    assert tuple(mod.ALPHAS)==(.001,.01,.1)
    assert mod.DRAW_COUNT==1999

def test_supporting_information_preserves_noncausality_and_scale():
    s=SI.read_text(encoding="utf8")
    assert "### S2.14 |" in s
    for phrase in ("16,801","1,665","+0.001032","+0.002355",
        "+0.001323","+0.000040","−0.000404","water2","other 24"):
        assert phrase in s,phrase
    assert "do not evaluate site-specific soil moisture" in s.lower()
    assert "does not mean soil has no predictive association" in s.lower()
    assert '"inat_taxon_id"' not in R.read_text(encoding="utf8")
    assert '"photo_id"' not in R.read_text(encoding="utf8")
