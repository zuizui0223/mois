"""Original-source biogeography/climate/soil/elevation and actual-DEM boundary.

Only safe aggregate receipts are used in public CI. The two original FCP ZIPs
and any private site coordinates/photography must stay outside public GitHub.
"""
from __future__ import annotations

import importlib.util
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
SAFE=ROOT/"results/mois_biogeography_climate_soil_elevation_full_audit_20261015.json"
CODE=ROOT/"scripts/analysis/mois_full_biogeography_climate_soil_audit_20261015.py"
DEM=ROOT/"scripts/analysis/mois_dem_terrain_metrics_gate_20261015.py"
OLD=ROOT/"results/fcp_macro_direction_20261013/decision_receipt.json"
SI=ROOT/"docs/MOIS_GEB_SUPPORTING_INFORMATION_20261014.md"
MEMO=ROOT/"docs/MOIS_BIOGEOGRAPHY_CLIMATE_SOIL_ELEVATION_TERRAIN_AUDIT_20261015.md"

def module(path,name):
 spec=importlib.util.spec_from_file_location(name,path)
 result=importlib.util.module_from_spec(spec)
 spec.loader.exec_module(result)
 return result

def data(p):return json.loads(p.read_text(encoding="utf8"))

def test_reproduce_original_soil_complete_population_and_primary_slope():
 a=data(SAFE);d=data(OLD)
 assert a["schema"]=="mois_biogeography_climate_soil_elevation_full_audit_safe_summary_v1"
 cohort=a["source_scope"]
 assert (cohort["photo_opportunities"],cohort["classifiable_photos"])==(100543,45953)
 assert (cohort["climate26_photos"],cohort["soil36_photos"])==(45836,31925)
 assert (cohort["soil36_species"],cohort["soil36_genera"],cohort["soil36_sampled_photos"])==(16801,1665,26163)
 assert (cohort["climate26_species"],cohort["climate26_genera"])==(20694,1934)
 assert cohort["source_sha256"]["photo"]==d["source_source_SHA256"]["original_photo"]
 assert cohort["source_sha256"]["environment"]==d["source_source_SHA256"]["all36_abiotic"]
 source=[q for q in d["key_white_directional_observations"] if
         q["mode"]=="soil" and q["quality"]=="all" and q["region"]=="world" and
         q["adjustment"]=="genus_geo" and q["axis"]=="wc_bio12" and q["hue"]=="white"]
 assert len(source)==1
 assert math.isclose(cohort["original_soil36_BIO12_genus_equal_coef_pp_per_photo_SD"],source[0]["point_pp"],abs_tol=1e-8)

def test_moisture_heldout_benefit_and_soil_increment_negative():
 q=data(SAFE);cv=q["heldout_genus_cv"]
 assert cv["folds"]==5 and cv["genus_holdout"]
 assert cv["boot_draws"]==1999 and cv["fold_species_counts"]==[3212,3505,3803,3306,2975]
 m=cv["models"];water=m["geography_moisture"]
 assert water["gain_vs_geography"]>0.0023
 assert water["gain_ci"][0]>0
 soil5=m["geography_climate_elevation_soil5"]
 soil10=m["geography_climate_elevation_soil10"]
 assert soil5["increment_vs_climate_elevation"]<0
 assert soil5["increment_ci"][1]<0
 assert soil10["increment_vs_climate_elevation"]<soil5["increment_vs_climate_elevation"]
 assert soil10["increment_ci"][1]<0
 assert m["geography_elevation"]["gain_ci"][0]<0<m["geography_elevation"]["gain_ci"][1]

def test_in_sample_soil_proxy_not_mistaken_for_transfer_or_measured_capacity():
 q=data(SAFE)
 z=q["global_same_soil36_taxa_OLS_congeneric"]["adding_soil5"]
 assert z["soil_available_water_proxy"]["point_pp_per_source_soil36_photo_site_SD"]>1.9
 assert z["soil_available_water_proxy"]["selected_within_model_maxT95_pp"][0]>0
 for name in ("soil_pH","soil_SOC","soil_N","soil_clay"):
  c=z[name]["selected_within_model_maxT95_pp"]
  assert c[0]<0<c[1]
 diag=q["global_same_soil36_taxa_OLS_congeneric"]["soil10_unpenalized_collinearity"]
 assert abs(diag["soil_texture_mean"]-1000)<.1
 assert diag["soil_texture_sd"]<.5
 assert diag["unpenalized_condition_number_about"]>1e6
 assert "proxy" in " ".join(q["scientific_limits"]).lower()

def test_approximate_continent_coverage_and_sector_heldout_nonuniversality():
 a=data(SAFE)
 b=a["biogeography"]
 assert b["climate26_species_polygon_matched"]+b["climate26_species_polygon_unresolved"]==20694
 assert b["climate26_species_polygon_unresolved"]==3158
 assert b["continent_white_species_equal_from_classified_photos"][0]["continent"]=="Africa"
 assert "not ecoregion" in b["continent_definition"]
 sec=b["longitude_sector_leaveout_exclude_target_genera"]
 assert len(sec)==4 and sum(x["test_taxa"] for x in sec)==16801
 assert sec[0]["moisture_gain_vs_geography"]>0
 assert sec[2]["moisture_gain_vs_geography"]<0
 assert sum(x["moisture_gain_vs_geography"]>0 for x in sec)==3

def test_soil_missingness_by_elevation_and_terrain_not_measured():
 q=data(SAFE)
 z=q["soil_missingness"]["altitude_photo_quartiles"]
 assert len(z)==4 and sum(x["photo_count"] for x in z)==45836
 assert sum(x["soil36_count"] for x in z)==31925
 assert z[0]["completeness"]<.45 and z[-1]["completeness"]>.87
 t=q["terrain"]
 assert not t["independent_dem_supplied"]
 for x in ("actual_slope","aspect","topographic_position_index","terrain_ruggedness_index"):
  assert t[x] is None
 assert "NO original terrain-geometry" in t["interpretation"]

def test_synthetic_genus_statistics_and_terrain_math():
 m=module(CODE,"mois_integrated")
 z=module(DEM,"mois_terrain")
 z.self_test()
 rng=np.random.default_rng(22)
 d=pd.DataFrame({"genus":np.repeat(["A","B","C","D"],5),
                 "white":np.tile(np.linspace(0.1,.9,5),4),
                 "abs_lat":rng.normal(size=20),
                 "wc_bio12":rng.normal(size=20)})
 vals,norm={},{"abs_lat":1.,"wc_bio12":1.}
 A,B,T=m.moments(d,norm,["abs_lat","wc_bio12"])
 assert A.shape==(4,2,2) and B.shape==(4,2)
 assert T>0 and np.linalg.cond(A.sum(axis=0))<100
 assert all(m.fold(g) in range(5) for g in ["A","B","C"])
 # Exactly constant photographic white within every genus -> zero numerator.
 d["white"]=np.repeat([0.,.25,.75,1.],5)
 _,rhs,_=m.moments(d,norm,["abs_lat","wc_bio12"])
 np.testing.assert_allclose(rhs,0,atol=1e-12)

def test_manuscript_si_marks_full_scope_and_real_dem_blocker():
 text=SI.read_text(encoding="utf8")
 memo=MEMO.read_text(encoding="utf8")
 assert "### S2.12 |" in text
 for phrase in ("16,801","1,665","0.002355","0.000435","44.6%","87.7%"):
  assert phrase in text and memo,phrase
 assert "terrain" in memo.lower() and "DEM" in memo
 assert "not" in memo.lower()
 assert "actual" in memo.lower()
 assert '"inat_taxon_id"' not in SAFE.read_text(encoding="utf8")
