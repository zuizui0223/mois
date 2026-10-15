"""Joint rainfall/vapour source-aggregate and synthetic regression verification.

The original 34/39 MB licensed photo/environment ZIPs are intentionally absent
from public CI. Only full raw source owner-side replay checks these bytes.
"""
from __future__ import annotations
import importlib.util
import json
import math
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT/"scripts/analysis/mois_joint_rain_vapour_geography_20261015.py"
RESULT=ROOT/"results/mois_original_joint_rain_vapour_geography_20261015.json"
ORIGINAL=ROOT/"results/fcp_macro_direction_20261013/decision_receipt.json"
GEO=ROOT/"results/mois_geography_vapour_ovb_20261015.json"
SI=ROOT/"docs/MOIS_GEB_SUPPORTING_INFORMATION_20261014.md"


def module():
    spec=importlib.util.spec_from_file_location("joint_floral_moisture",SRC)
    loaded=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(loaded)
    return loaded


def test_synthetic_two_axis_rotation_and_wls_model_contrasts():
    mod=module()
    rng=np.random.default_rng(20261015)
    n=20
    ma=np.empty((n,6,6)); y=np.empty((n,6))
    beta=np.array([.3,-.25,.15,.22,.8,1.4])
    for i in range(n):
        z=rng.normal(size=(9,6))
        ma[i]=z.T@z+np.eye(6)
        y[i]=ma[i]@beta
    A,points,boots=mod.point_and_boot(ma,y,79,321)
    for name,idx in mod.MODELS.items():
        assert points[name].shape==(len(idx),)
        assert boots[name].shape==(79,len(idx))
        np.testing.assert_allclose(boots[name],
            np.tile(points[name],(79,1)),rtol=0,atol=1e-13)
    original=points["geo4_joint"][-2:]
    r=np.array([[1,1],[1,-1]])/np.sqrt(2)
    np.testing.assert_allclose(np.linalg.norm(r@original),np.linalg.norm(original),atol=1e-13)
    assert math.isclose((r@original)[0],sum(original)/np.sqrt(2),abs_tol=1e-13)
    assert -1<mod.partial_correlation_after_geography(A,[0,1,2,3])<1


def test_same_species_source_receipt_and_original_single_moisture_coefficients():
    q=json.loads(RESULT.read_text(encoding="utf8"))
    old=json.loads(ORIGINAL.read_text(encoding="utf8"))
    prior=json.loads(GEO.read_text(encoding="utf8"))
    assert q["schema"]=="mois_original_20694_joint_moisture_and_geography_v1"
    assert q["original_source_sha256"]["photo"]==old["source_source_SHA256"]["original_photo"]
    assert q["original_source_sha256"]["environment"]==old["source_source_SHA256"]["all36_abiotic"]
    assert q["source_support"]["main_eligible_species"]==20694
    assert q["source_support"]["main_eligible_genera"]==1934
    assert q["source_support"]["main_eligible_photos"]==38498
    assert q["bootstrap"]["draws"]==1999 and q["bootstrap"]["seed"]==202610151330
    assert not q["bootstrap"]["programme_wide_multiplicity_correction"]
    for key in ("geo3","geo4"):
        d=q["models"][key]
        assert (d["n_species"],d["n_genera"])==(20694,1934)
        assert d["water_axis_partial_correlation_given_geography"]>0
        coeff=d["both_axes_each_adjusted_for_other"]
        rain=coeff["annual_precipitation"];vap=coeff["atmospheric_vapour"]
        assert rain["point_pp"]>0 and vap["point_pp"]>0
        assert rain["two_contrast_simultaneous95_pp"][0]<0<rain["two_contrast_simultaneous95_pp"][1]
        assert vap["two_contrast_simultaneous95_pp"][1]>0
        rotated=d["rotated_equal_standardized_moisture_axes"]
        common=rotated["common_equal_increase"]
        diff=rotated["rain_minus_vapour_contrast"]
        assert common["two_contrast_simultaneous95_pp"][0]>0
        assert diff["two_contrast_simultaneous95_pp"][0]<0<diff["two_contrast_simultaneous95_pp"][1]
        assert math.isclose(common["point_pp"],(rain["point_pp"]+vap["point_pp"])/math.sqrt(2),abs_tol=1e-12)
        assert math.isclose(diff["point_pp"],(rain["point_pp"]-vap["point_pp"])/math.sqrt(2),abs_tol=1e-12)
        assert rotated["formula"].startswith("common=(")
    a=q["models"]["geo3"]["single_axis_pp_per_original_photo_site_sd"]
    assert math.isclose(a["annual_precipitation"],2.254545235123587,abs_tol=1e-10)
    assert math.isclose(a["atmospheric_vapour"],3.951850325152891,abs_tol=1e-10)
    b=q["models"]["geo4"]
    assert math.isclose(b["single_axis_pp_per_original_photo_site_sd"]["atmospheric_vapour"],
                        prior["identity_comparisons"]["add_vapour_to_geography"]["vapour_model_coefficient_pp_per_source_sd"],abs_tol=1e-9)
    assert math.isclose(b["both_axes_each_adjusted_for_other"]["atmospheric_vapour"]["point_pp"],
                        prior["identity_comparisons"]["add_vapour_to_geography_and_rain"]["vapour_model_coefficient_pp_per_source_sd"],abs_tol=1e-9)
    assert math.isclose(b["rotated_equal_standardized_moisture_axes"]["common_equal_increase"]["point_pp"],3.644471843122978,abs_tol=1e-12)
    assert b["both_axes_each_adjusted_for_other"]["atmospheric_vapour"]["two_contrast_simultaneous95_pp"][0]<0
    assert b["rotated_equal_standardized_moisture_axes"]["common_equal_increase"]["two_contrast_simultaneous95_pp"][0]>1.13
    raw=RESULT.read_text(encoding="utf8")
    for leak in ('"inat_taxon_id"','"photo_id"','"latitude"','"longitude"'):
        assert leak not in raw


def test_si_identifies_common_axis_not_individual_causal_fitness():
    si=SI.read_text(encoding="utf8")
    assert "### S2.11 |" in si
    for part in ("20,694","1,934","+1.383","−0.167","+3.771","−0.163","+3.644",
                 "+1.137","+6.152","1,999","maxT","rotated"):
        assert part in si,part
    assert "not a measured soil-water supply" in si
    assert "not an equivalence test" in si
    assert "not identify" in si.lower()
    assert "genetically" in si.lower() or "genetic" in si.lower()
