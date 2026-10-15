"""Original-source full climate26 + soil5/soil10 × ridge × genus-fold checks.

CI sees original safe numerical receipts and synthetic arrays, NOT licensed raw
photographs, named taxa, precise locations or individual source genus rows.
"""
from __future__ import annotations
import importlib.util, json, math
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
RESULT=ROOT/"results/mois_soil_full26_alpha_fold_robustness_20261015.json"
BASE=ROOT/"results/mois_full26_climate_soil10_transfer_20261015.json"
CODE=ROOT/"scripts/analysis/mois_soil_full26_alpha_fold_robustness_20261015.py"
SI=ROOT/"docs/MOIS_GEB_SUPPORTING_INFORMATION_20261014.md"
MEMO=ROOT/"docs/MOIS_BIOGEOGRAPHY_CLIMATE_SOIL_ELEVATION_TERRAIN_AUDIT_20261015.md"

def module():
    spec=importlib.util.spec_from_file_location("mois_source_soil_penalty_robustness",CODE)
    result=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result

def test_frozen_original_source_cohort_and_original_001_fits_match():
    data=json.loads(RESULT.read_text(encoding="utf8"))
    base=json.loads(BASE.read_text(encoding="utf8"))
    assert data["schema"]=="mois_soil_full26_alpha_fold_predictive_robustness_v1"
    assert (data["population"]["qualified_species"],
            data["population"]["qualified_genera"],
            data["population"]["qualified_source_photos"])==(16801,1665,26163)
    assert data["source_photo_sha256"]==base["photo_sha256"]
    assert data["source_environment_sha256"]==base["environment_sha256"]
    assert data["same_species_across_all_models"]
    assert data["ridge_alpha_grid"]==[.0003,.001,.003,.01,.03,.1,.3,1.]
    assert len(data["genus_split_salts"])==4
    orig=data["results"]["original"]["at_original_alpha_001"]
    source=base["model_results"]
    assert math.isclose(orig["climate26_loss"],source["geo_plus_all26_climate"]["heldout_species_photo_white_squared_error"],abs_tol=1e-9)
    for key,model in (
        ("soil5_increment","geo_plus_all26_climate_plus_soil5"),
        ("soil10_increment","geo_plus_all26_climate_plus_soil10"),
    ):
        gold=source[model]["incremental_soil_prediction_gain_over_geo_climate26"]
        assert math.isclose(orig[key],gold,abs_tol=1e-9)

def test_same_cohort_32_nonindependent_sensitivities_and_exception():
    d=json.loads(RESULT.read_text(encoding="utf8"))
    m=module()
    assert d["n_settings"]==32
    summary=d["summary"]
    assert (summary["soil5_negative_gain_settings"],summary["soil10_negative_gain_settings"])==(31,32)
    assert (summary["soil5_positive_gain_settings"],summary["soil10_positive_gain_settings"])==(1,0)
    assert set(d["results"])==set(d["genus_split_salts"])
    assert sum(v["soil5_negative"] for v in summary["by_alpha_signs"].values())==31
    assert sum(v["soil10_negative"] for v in summary["by_alpha_signs"].values())==32
    for k,val in d["results"].items():
        assert sum(val["source_species_per_heldout_fold"])==16801
        for prefix in ("soil5","soil10"):
            assert len(val["soil_block_increment_by_alpha"][prefix])==8
            assert math.isclose(
                val["soil_block_increment_by_alpha"][prefix][3],
                val["at_original_alpha_001"][prefix+"_increment"],abs_tol=1e-12)
        assert val["at_original_alpha_001"]["soil5_bootstrap95"][0]<0
        assert val["at_original_alpha_001"]["soil10_bootstrap95"][0]<0
    assert math.isclose(d["results"]["salt202610153"]["soil_block_increment_by_alpha"]["soil5"][-1],.000011947,abs_tol=1e-12)
    assert d["results"]["salt202610151"]["at_original_alpha_001"]["soil5_bootstrap95"][1]>0
    assert d["results"]["salt202610152"]["at_original_alpha_001"]["soil10_bootstrap95"][1]>0
    assert "not 32 independent" not in " ".join(d["limitations"]).lower() or "nonindependent" in " ".join(d["limitations"]).lower()
    for x in ('"photo_id"','"inat_taxon_id"','"genus"','"latitude"','"longitude"'):
        assert x not in RESULT.read_text(encoding="utf8")

def test_synthetic_fold_exclusion_and_prediction_loss_formula():
    m=module()
    x=np.arange(75,dtype=float).reshape(25,3)/25
    g=np.repeat(np.arange(5),5)
    size=np.bincount(g)
    w=1/size[g]
    y=np.full(25,.4)
    fold=g.astype(np.int8)
    means,genus=m.model_losses([x,x[:,[0,2]],x[:,[1]]],y,g,size,w,fold,.01)
    assert genus.shape==(5,3)
    np.testing.assert_allclose(means,0,atol=1e-23)
    np.testing.assert_allclose(genus,0,atol=1e-23)
    one=np.array(["A","B","C","D","E"])
    assert np.array_equal(m.split(one,lambda k:{"A":0,"B":1,"C":2,"D":3,"E":4}[k],None),np.arange(5))
    salted=m.split(one,None,202610151)
    assert len(salted)==5 and all(int(z) in range(5) for z in salted)
    assert np.array_equal(salted,m.split(one,None,202610151))

def test_source_and_si_cautions_not_claimed_selection():
    data=json.loads(RESULT.read_text(encoding="utf8"))
    caution=" ".join(data["limitations"]).lower()
    assert "nonindependent" in caution
    assert "not proof" in caution
    assert "soil36" in caution
    assert "fitness" in caution
    assert "terrain" in caution
    supp=SI.read_text(encoding="utf8")
    memo=MEMO.read_text(encoding="utf8")
    assert "### S2.13 |" in supp
    for z in ("31/32","32/32","16,801","1,665","0.000012"):
        assert z in supp and z in memo
