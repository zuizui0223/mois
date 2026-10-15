"""Original-source ecosystem context and eight-source-matched comparisons.

The licensed original source photos and coordinates are NOT in GitHub;
the independent full ZIP replay is separate from this aggregate-only CI.
"""
from __future__ import annotations
import importlib.util, json, math
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
RECEIPT=ROOT/"results/mois_source_soil_moisture_context_eight_contrast_20261015.json"
FULL26=ROOT/"results/mois_full26_climate_soil10_transfer_20261015.json"
BLOCKS=ROOT/"results/mois_soil_crossdomain_block_audit_20261015.json"
SCRIPT=ROOT/"scripts/analysis/mois_soil_moisture_context_20261015.py"
SI=ROOT/"docs/MOIS_GEB_SUPPORTING_INFORMATION_20261014.md"
MEMO=ROOT/"docs/MOIS_SOIL_MOISTURE_CONTEXT_HETEROGENEITY_20261015.md"

def source():
    spec=importlib.util.spec_from_file_location("soil_moisture_context_src",SCRIPT)
    result=importlib.util.module_from_spec(spec);spec.loader.exec_module(result)
    return result

def data(path):return json.loads(path.read_text(encoding="utf8"))

def test_original_source_same_16801_genus_holdout_and_numeric_reconciliation():
    x=data(RECEIPT);m=source();older=data(FULL26)
    assert x["schema"]=="mois_source_soil_moisture_context_eight_contrast_v1"
    assert x["source_sha256"]=={"photo":older["photo_sha256"],"environment":older["environment_sha256"]}
    assert x["original_support"]=={"species":16801,"genera":1665,"photos_in_eligible_taxa":26163}
    assert x["models"]==list(m.MODELS)
    a=x["model_heldout_equal_genus_mse"]
    assert abs(a["climate26"]-older["model_results"]["geo_plus_all26_climate"]["heldout_species_photo_white_squared_error"])<1e-9
    assert abs(a["climate26_soil5"]-older["model_results"]["geo_plus_all26_climate_plus_soil5"]["heldout_species_photo_white_squared_error"])<1e-9
    assert abs(a["climate26_soil10"]-older["model_results"]["geo_plus_all26_climate_plus_soil10"]["heldout_species_photo_white_squared_error"])<1e-9
    assert abs(a["water2"]-0.22745535180555312)<1e-9
    assert abs(a["geography"]-0.22981006389209194)<1e-9
    assert abs(a["soil5"]-0.2287784849)<1e-8
    assert math.isclose(x["soil5_addition_to_full26"]["point_gain"],a["climate26"]-a["climate26_soil5"],abs_tol=1e-10)
    assert x["soil5_addition_to_full26"]["ordinary95"][1]<0

def test_context_exploratory_subsets_and_eightway_adjustment_no_effect_modification():
    o=data(RECEIPT);keys=list(o["context"])
    assert keys==["annual_rainfall","peak_temperature","elevation","latitude_abs"]
    u=o["uncertainty"]
    assert u["bootstrap_draws"]==1999
    assert u["half_seed_base"]==202610151755
    assert u["paired_seed"]==202610151801
    assert 2.7<u["eight_selected_contrast_maxT_critical"]<2.8
    assert u["maxT_only_eight_posthoc_contrasts"]
    assert u["no_study_wide_multiple_testing_or_nested_model_refitting"]
    assert math.isclose(o["context"]["annual_rainfall"]["source_median"],834.)
    assert math.isclose(o["context"]["elevation"]["source_median"],547.5)
    assert math.isclose(o["context"]["latitude_abs"]["source_median"],32.2769737052)
    for ctx,d in o["context"].items():
        g1=d["halves"]["lower"];g2=d["halves"]["higher"]
        assert g1["n_species"]+g2["n_species"]==16801
        assert 0<g1["n_genera_present"]<=1665 and 0<g2["n_genera_present"]<=1665
        assert g1["n_genera_present"]+g2["n_genera_present"]>1665  # same genus may have both climate classes
        for which in ("lower","higher"):
            a=d["halves"][which]["metrics"]
            assert a["water_gain_vs_geography"]["gain"]>0
            assert a["soil5_increment_vs_climate26"]["gain"]<0
            for estimate in a.values():
                lo,hi=estimate["ordinary95"];assert lo<=estimate["gain"]<=hi
        for target,m in d["high_minus_low"].items():
            assert target in ("water_gain_vs_geography","soil5_increment_vs_climate26")
            v=d["halves"]["higher"]["metrics"][target]["gain"]-d["halves"]["lower"]["metrics"][target]["gain"]
            assert math.isclose(v,m["higher_minus_lower"],abs_tol=1e-11)
            ci=m["eight_selected_contrast_maxT95"]
            assert ci[0]<0<ci[1]
            assert m["paired_genus_ordinary95"][0]<0<m["paired_genus_ordinary95"][1]

def test_synthetic_nested_aggregate_tolerance_and_estimator_source():
    m=source()
    assert len(m.MODELS["climate26"])==29
    assert len(m.MODELS["climate26_soil5"])==34
    assert len(m.MODELS["climate26_soil10"])==39
    assert len(m.MODELS["other24"])==27
    m.check_numbers({"point":[1.2,2.3],"text":"source","count":16801},
                    {"point":[1.20000001,2.30000001],"text":"source","count":16801})
    try:m.check_numbers({"num":1.0},{"num":2.0})
    except AssertionError:pass
    else:raise AssertionError("A 1-unit source numeric drift must fail")
    assert "mois_full_biogeography_climate_soil_audit_20261015.py" in SCRIPT.read_text(encoding="utf8")
    assert "original_loader" in SCRIPT.read_text(encoding="utf8")
    assert "training" in SCRIPT.read_text(encoding="utf8").lower()
    assert '"inat_taxon_id"' not in RECEIPT.read_text(encoding="utf8")
    assert '"photo_id"' not in RECEIPT.read_text(encoding="utf8")

def test_si_interpretation_and_real_terrain_still_not_claimed():
    si=SI.read_text(encoding="utf8");memo=MEMO.read_text(encoding="utf8")
    assert "### S2.15 |" in si
    for value in ("16,801","1,665","−0.000404461","−0.000509","−0.000480","2.7490","834mm","547.5m"):
        assert value in si or value in memo
    assert "All eight high-minus-low simultaneous 95% intervals include zero" in si
    assert "No original source photo-site terrain values have yet been extracted" in si
    assert "10.1594/PANGAEA.867115" in si
    assert "NO" not in memo or "NOT" in memo
