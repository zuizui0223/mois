"""Synthetic FWL two-way geographic/month controls and photographed-species means."""
from __future__ import annotations
import csv
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location(
    "mois_species_equal_month_FE",
    ROOT/"scripts/analysis/mois_monthly_white_species_equal_cell_FE_20261016.py")
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def clean_cell_month(rows):
    assert len(rows)>100
    assert all(0<=r["white_species_equal"]<=1 for r in rows)
    assert len(set(r["cell"] for r in rows))>=12
    return rows


def test_same_place_climate_vs_white_fixed_effects_recover_known_synthetic_slopes():
    result=[]
    for c in range(24):
        for month in range(1,13):
            t=2*math.pi*month/12
            x=[
                (1+c%7)*math.sin(t),
                (1+c%5)*math.cos(t)*1000,
                (1+c%11)*math.sin(2*t)*10,
                (1+c%3)*math.cos(3*t)*.10
            ]
            y=(.34+.01*(c-12)+.03*math.cos(t)+
               .020*x[0]+.000015*x[1]+.0005*x[2]-.10*x[3])
            result.append({"cell":f"C{c}","month":month,
                "domain":"NORTH_EXTRATROPICAL",
                "n_species":8,"white_species_equal":y,
                **dict(zip(m.DYNAMIC,x))})
    result=clean_cell_month(result)
    for month_effect in (False,True):
        fit=m.fit_fe(result,include_month=month_effect,min_cells=10)
        assert fit["status"]=="ESTIMATED_EXPLORATORY_POINT_WITH_SEPARATE_CELL_CLUSTER_GATE"
        assert fit["n_qualified_geographic_cells"]==24
        assert fit["joint_dynamic_4axis_effective_rank"]==4
        assert fit["stable_geographic_soil_elevation_terrain_cell_means_absorbed"] is True
        if month_effect:
            # Month fixed effects absorb the deliberately injected common
            # seasonal signal and recover all FOUR original climate axes.
            assert fit["within_fixed_effect_adjusted_variation_R2"]==pytest.approx(1,abs=1e-10)
            truth=(.020,.000015,.0005,-.10)
            for name,beta in zip(m.DYNAMIC,truth):
                out=fit["four_dynamic_partial_effect_coefficients"][name]
                assert out["slope_pp_per_original_monthly_climate_unit"]==pytest.approx(100*beta,abs=1e-8)
        else:
            # With only geographic-cell fixed effects, the unmodelled
            # calendar-month term remains, confounding seasonal solar T/P.
            assert .9<fit["within_fixed_effect_adjusted_variation_R2"]<1


def test_sparse_sites_refuse_joint_temperature_solar_water_slopes():
    items=[
      {"cell":"C1","month":1,"n_species":2,"white_species_equal":.2,
       **dict(zip(m.DYNAMIC,[1.,1000.,50.,.25]))},
      {"cell":"C1","month":2,"n_species":3,"white_species_equal":.7,
       **dict(zip(m.DYNAMIC,[3.,2000.,25.,.3]))},
      {"cell":"C2","month":1,"n_species":4,"white_species_equal":.6,
       **dict(zip(m.DYNAMIC,[0.,1100.,15.,.1]))},
      {"cell":"C2","month":2,"n_species":2,"white_species_equal":.3,
       **dict(zip(m.DYNAMIC,[3.,2500.,55.,.3]))}
    ]
    fit=m.fit_fe(items,include_month=True,min_cells=10)
    assert fit["status"]=="NOT_ESTIMABLE_TOO_FEW_REPEAT_MONTH_GEOGRAPHIC_CELLS"
    assert "four_dynamic_partial_effect_coefficients" not in fit


def test_source_bound_photo_classifiability_and_equal_species_cellmonth(tmp_path):
    path=tmp_path/"private_site_month.csv"
    rows=[]
    for month in (3,7):
        for tax,photo,morph in (
            (1,1,"white"),
            (1,2,"white"),
            (1,3,"white"),
            (2,4,"blue_purple"),
            (3,5,"mixed_uncertain"),
        ):
            row=dict.fromkeys(m.REQUIRED,"")
            row.update({
                "measurement_id":f"FCPR-{month}-{photo}",
                "observation_id":f"obs-{month}-{photo}",
                "inat_taxon_id":str(tax),"genus":"A" if tax==1 else "B",
                "original_month":str(month),
                "domain":"NORTH_EXTRATROPICAL",
                "climate_grid_row":"300","climate_grid_col":"1200",
                "source_climate_join_status":"VALID",
                "sample_observation_weight":"1",
                "morph":morph,
                "measurement_status":(
                    "classified_four_state_morph" if morph!="mixed_uncertain"
                    else "not_evaluable_roi_or_flip_gate"),
                **dict(zip(m.DYNAMIC,[month*1.,month*200.,month*10.,month*.05]))
            })
            rows.append(row)
    with path.open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=sorted(m.REQUIRED))
        writer.writeheader();writer.writerows(rows)
    receipt=tmp_path/"original_source.json"
    receipt.write_text(json.dumps({
        "schema":"mois_original_photo_all_monthly_dynamic_climate_site_join_v1",
        "status":"COMPLETE_EXACT_ORIGINAL_SOURCE_PHOTO_MONTHLY_T_S_P_VAPR_JOIN_NOT_13ENV_REGRESSION",
        "new_private_photo_monthly_climate_csv_sha256":m.sha256_file(path),
        "seasonal_photo_white_environment_regression_fitted":False,
        "selected_original_photographs":10,
        "source_four_colour_classifiable_photos":8,
        "source_classifiable_white_photos":6,
        "synthetic_test_only":True,
        "original_48_raster_manifest_sha256":"FAKE48"
    }))
    report=m.analyze(path,receipt,tmp_path/"public.json",synthetic=True)
    assert report["n_original_photos_in_full_source"]==10
    assert report["n_classifiable_four_hue"]==8
    assert report["n_classifiable_white"]==6
    assert report["n_month_source_climate_cells"]==2
    assert report["n_species_cell_month_rows"]==4
    assert report["n_climate_cells_with_at_least_two_distinct_classified_flower_months"]==1
    assert len(report["three_exploratory_fixed_effect_support_and_models"])==3
    assert all(x["status"].startswith("NOT_ESTIMABLE") for x in report[
        "three_exploratory_fixed_effect_support_and_models"])
    assert "measurement_id" not in (tmp_path/"public.json").read_text()
    assert report["target_is_photographed_blooming_species_turnover_not_individual_colour_plasticity"] is True
    with pytest.raises(ValueError,match="Refuse"):
        m.analyze(path,receipt,tmp_path/"public.json",synthetic=True)


def test_climate_field_variation_same_geographic_cell_month_rejected(tmp_path):
    items=[
       {"cell":"C1","month":3,"n_species":5,"white_species_equal":.5,
        **dict(zip(m.DYNAMIC,[0,0,0,0]))},
       {"cell":"C1","month":7,"n_species":5,"white_species_equal":.5,
        **dict(zip(m.DYNAMIC,[0,0,0,0]))},
    ]
    fitted=m.fit_fe(items,include_month=True,min_cells=2)
    assert fitted["status"]=="NOT_ESTIMABLE_TOO_FEW_REPEAT_MONTH_GEOGRAPHIC_CELLS"
    with pytest.raises(ValueError):
        m.fit_fe(items,include_month=True,min_cells=1)


def test_latitude_zone_calendar_month_fixed_effects_recover_known_dynamic_axes():
    points=[]
    for c in range(18):
        for month in range(1,13):
            t=month*2*math.pi/12
            x=[(1+c%6)*math.sin(t),
               (1+c%4)*math.cos(t)*900,
               (1+c%7)*math.sin(2*t)*16,
               (1+c%3)*math.cos(3*t)*.08]
            y=.35+.006*(c-9)+.02*math.cos(t)+.01*x[0]+.000012*x[1]               +.0003*x[2]-.06*x[3]
            points.append({
               "cell":f"source-cell-{c}","month":month,
               "domain":"NORTH_EXTRATROPICAL",
               "n_species":5,"white_species_equal":y,
               **dict(zip(m.DYNAMIC,x))
            })
    fit=m.fit_fe(points,include_month=True,zone_month=True,min_cells=10)
    assert fit["status"]=="ESTIMATED_EXPLORATORY_POINT_WITH_SEPARATE_CELL_CLUSTER_GATE"
    assert fit["fixed_effect_model"]=="cell_and_latitude_zone_month"
    for name,coef in zip(m.DYNAMIC,(.01,.000012,.0003,-.06)):
        assert fit["four_dynamic_partial_effect_coefficients"][name][
            "slope_pp_per_original_monthly_climate_unit"]==pytest.approx(
                coef*100,abs=1e-8)
    with pytest.raises(ValueError,match="requires month FE"):
        m.fit_fe(points,include_month=False,zone_month=True,min_cells=10)



def test_at_least_30_independent_geographic_cells_obtain_only_exploratory_cr1_cluster_intervals():
    rng=np.random.default_rng(20261016)
    points=[]
    for c in range(80):
        for month in range(1,13):
            t=2*math.pi*month/12
            x=[(1+c%4)*math.sin(t),
               (1+c%5)*math.cos(t)*800,
               (1+c%7)*math.sin(2*t)*8,
               (1+c%3)*math.cos(3*t)*.08]
            white=(.28+.0015*c+.015*math.cos(t)+
                   .004*x[0]+.000007*x[1]+.0002*x[2]-.03*x[3]+
                   rng.normal(0,.008))
            points.append({
                "cell":f"source_cell_{c:03}","domain":"NORTH_EXTRATROPICAL",
                "month":month,"white_species_equal":white,"n_species":4,
                **dict(zip(m.DYNAMIC,x))
            })
    for include_month,zone_month in ((False,False),(True,False),(True,True)):
        fit=m.fit_fe(points,include_month=include_month,
                     zone_month=zone_month,min_cells=30)
        assert fit["status"]=="ESTIMATED_EXPLORATORY_POINT_WITH_SEPARATE_CELL_CLUSTER_GATE"
        cl=fit["geographic_cell_cluster_robust_uncertainty"]
        assert cl["status"]=="EXPLORATORY_GEO_CELL_CR1_ASYMPTOTIC_AND_4AXIS_BONFERRONI"
        assert cl["n_independent_geographic_cell_clusters"]==80
        assert cl["n_climate_cell_clusters_with_nonzero_scores"]==80
        assert cl["conservative_residual_degrees_of_freedom"]>0
        assert cl["taxon_and_observer_dependence_not_resolved"] is True
        assert cl["photographic_missing_not_at_random_not_resolved"] is True
        assert cl["familywise_four_climate_axes_z_multiplier"]>cl["nominal_z_multiplier"]
        assert set(cl["intervals"])==set(m.DYNAMIC)
        for name in m.DYNAMIC:
            axis=cl["intervals"][name]
            lo,hi=axis["single_axis_asymptotic_95pct_interval_pp_per_original_unit"]
            flo,fhi=axis["four_axis_bonferroni_familywise_95pct_interval_pp_per_original_unit"]
            estimate=fit["four_dynamic_partial_effect_coefficients"][name][
                "slope_pp_per_original_monthly_climate_unit"]
            assert flo<=lo<estimate<hi<=fhi
            assert axis["cluster_se_pp_per_original_monthly_unit"]>0


def test_exactly_29_repeat_month_cells_may_have_a_point_but_must_not_get_cluster_intervals():
    items=[]
    for c in range(29):
        for month in range(1,13):
            t=2*math.pi*month/12
            x=[(c%5+1)*math.sin(t),
               (c%7+1)*math.cos(t)*1000,
               (c%3+1)*math.sin(2*t)*10,
               (c%4+1)*math.cos(3*t)*.1]
            items.append({
                "cell":f"c{c}","month":month,
                "domain":"NORTH_EXTRATROPICAL",
                "n_species":4,
                "white_species_equal":.3+.004*x[0]+.000006*x[1],
                **dict(zip(m.DYNAMIC,x))
            })
    fit=m.fit_fe(items,include_month=True,zone_month=True,min_cells=10)
    assert fit["status"]=="ESTIMATED_EXPLORATORY_POINT_WITH_SEPARATE_CELL_CLUSTER_GATE"
    assert fit["geographic_cell_cluster_robust_uncertainty"]["status"]==(
        "NOT_ESTIMABLE_FEWER_THAN_30_REPEAT_CLIMATE_CELLS")
    assert "intervals" not in fit["geographic_cell_cluster_robust_uncertainty"]
