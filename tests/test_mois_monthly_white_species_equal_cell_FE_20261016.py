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
                "n_species":8,"white_species_equal":y,
                **dict(zip(m.DYNAMIC,x))})
    result=clean_cell_month(result)
    for month_effect in (False,True):
        fit=m.fit_fe(result,include_month=month_effect,min_cells=10)
        assert fit["status"]=="ESTIMATED_EXPLORATORY_POINT_ONLY_NO_CLUSTER_UNCERTAINTY"
        assert fit["n_qualified_geographic_cells"]==24
        assert fit["joint_dynamic_4axis_effective_rank"]==4
        assert fit["within_fixed_effect_adjusted_variation_R2"]==pytest.approx(1,abs=1e-10)
        truth=(.020,.000015,.0005,-.10)
        for name,beta in zip(m.DYNAMIC,truth):
            out=fit["four_dynamic_partial_effect_coefficients"][name]
            assert out["slope_pp_per_original_monthly_climate_unit"]==pytest.approx(100*beta,abs=1e-8)
        assert fit["static_soil_elevation_latitude_terrain_cell_means_absorbed"] is True if "static_soil_elevation_latitude_terrain_cell_means_absorbed" in fit else fit["stable_geographic_soil_elevation_terrain_cell_means_absorbed"] is True


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
    assert all(x["status"].startswith("NOT_ESTIMABLE") for x in report[
        "two_exploratory_fixed_effect_support_and_models"])
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
