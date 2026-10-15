"""Source-audited floral white climate-soil variance partition and transfer target.

Public CI does not carry original licensed photo or site ZIPs. Test synthetic
five-block Shapley math and source-derived anonymous aggregate against original
source-supported soil36 cohort and previously verified heldout predictions.
"""
from __future__ import annotations
import importlib.util,json,math
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
CODE=ROOT/"scripts/analysis/mois_environment_block_shapley_residual_cv_20261015.py"
RESULT=ROOT/"results/mois_env_shapley_and_genus_residual_transfer_20261015.json"
FULL=ROOT/"results/mois_biogeography_climate_soil_elevation_full_audit_20261015.json"
DOC=ROOT/"docs/MOIS_ENVIRONMENT_BLOCK_PARTITION_GENUS_RESIDUAL_TRANSFER_20261015.md"
SI=ROOT/"docs/MOIS_GEB_SUPPORTING_INFORMATION_20261014.md"

def load():
 spec=importlib.util.spec_from_file_location("mois_original_environmental_blocks",CODE)
 obj=importlib.util.module_from_spec(spec)
 spec.loader.exec_module(obj)
 return obj

def test_synthetic_block_shapley_monotone_and_sums_to_total():
 m=load()
 rng=np.random.default_rng(20261015)
 n_genus,n_per=15,10
 g=np.repeat(np.arange(n_genus),n_per)
 xi=rng.normal(size=(len(g),13))
 y=.2*xi[:,3]+.15*xi[:,6]+.3*xi[:,11]+rng.normal(scale=.5,size=len(g))
 X=xi.copy()
 for gen in range(n_genus):
  ix=g==gen
  X[ix]-=X[ix].mean(axis=0)
  y[ix]-=y[ix].mean()
 w=np.repeat(1/n_per,len(g))
 A,b,SST=m.moments(g,w,X,y,n_genus)
 assert A.shape==(15,13,13)
 keys=list(m.BLOCKS)
 names=[x for block in m.BLOCKS.values() for x in block]
 idx={k:np.array([names.index(t) for t in m.BLOCKS[k]]) for k in keys}
 subsets={mask:np.concatenate([idx[keys[j]] for j in range(5) if mask&(1<<j)])
          for mask in range(1,32)}
 F=m.fitted_r2_by_subset(A.sum(axis=0),b.sum(axis=0),SST.sum(),subsets)
 parts=m.shapley(F)
 assert 0<F[31]<1
 assert np.all(parts>=-1e-14)
 assert math.isclose(parts.sum(),F[31],abs_tol=1e-12)
 assert all(F[31]>=F[31^(1<<i)]-1e-12 for i in range(5))

def test_original_source_shapley_point_and_tiny_total_support():
 q=json.loads(RESULT.read_text(encoding="utf8"))
 old=json.loads(FULL.read_text(encoding="utf8"))
 assert q["schema"]=="mois_original_16801_environmental_block_partition_plus_residual_transfer_v1"
 assert (q["eligible_species"],q["eligible_genera"],q["source_photo_hue_observations_of_eligible_species"])==(16801,1665,26163)
 assert q["same_support_soil36_photo_opportunities"]==31925
 assert q["original_photo_sha256"]==old["source_scope"]["source_sha256"]["photo"]
 assert q["original_environment_sha256"]==old["source_scope"]["source_sha256"]["environment"]
 r=q["shapley_explanatory_in_sample"]
 assert r["shapley_bootstrap"]["draws"]==999 and r["shapley_bootstrap"]["seed"]==202610151610
 assert math.isclose(r["in_sample_total_R2"],.006202366064066511,abs_tol=1e-12)
 sh=r["block_attribution"]
 assert set(sh)==set(load().BLOCKS)
 assert math.isclose(sum(v["shapley_fraction_of_total_within_genus_white_variance"] for v in sh.values()),r["in_sample_total_R2"],abs_tol=1e-12)
 assert sh["Soil5"]["shapley_fraction_of_total_within_genus_white_variance"]>sh["Moisture"]["shapley_fraction_of_total_within_genus_white_variance"]
 assert sh["Geography"]["shapley_fraction_of_total_within_genus_white_variance"]>sh["Soil5"]["shapley_fraction_of_total_within_genus_white_variance"]
 for d in sh.values():
  assert d["partial_R2_given_all_other_four_blocks"]>=0
  assert d["shapley_genus_bootstrap95"][0]>=0
 assert "nonnegative by construction" in r["shapley_bootstrap"]["method"]

def test_matched_residual_cv_is_not_absolute_outofgenus_prediction():
 q=json.loads(RESULT.read_text(encoding="utf8"))
 old=json.loads(FULL.read_text(encoding="utf8"))
 cv=q["within_genus_residual_out_of_genus"]
 assert cv["fold_species_counts"]==old["heldout_genus_cv"]["fold_species_counts"]==[3212,3505,3803,3306,2975]
 assert cv["bootstrap"]["seed"]==202610151710
 assert cv["bootstrap"]["draws"]==1999
 assert "test-genus mean" in cv["target"]
 assert "NOT prospective" in cv["target"]
 m=cv["model_results"]
 original_cv=old["heldout_genus_cv"]["models"]
 # Same species and fold seed, DIFFERENT response grain/target: no equality.
 assert math.isclose(m["geography_moisture"]["gain_over_geography"],.00023194248948352203,abs_tol=1e-12)
 a=m["geography_moisture"]["gain_over_geography_genus_bootstrap95"]
 assert a[0]<0<a[1]
 soil=m["geography_moisture_temperature_elevation_soil5"]
 assert math.isclose(soil["increment_vs_geo_climate_elevation"],-.0000822771208470463,abs_tol=1e-12)
 assert soil["increment_vs_geo_climate_elevation_genus_bootstrap95"][0]<0<soil["increment_vs_geo_climate_elevation_genus_bootstrap95"][1]
 assert original_cv["geography_moisture"]["gain_vs_geography"]>.0023
 assert original_cv["geography_climate_elevation_soil5"]["increment_vs_climate_elevation"]<-.0004
 assert m["geography_moisture"]["heldout_residual_squared_error"]<.13
 assert m["geography_moisture"]["heldout_residual_squared_error"]>.12

def test_manuscript_scope_and_no_coordinates_leaked():
 txt=RESULT.read_text(encoding="utf8")
 for sensitive in ('"inat_taxon_id"','"photo_id"','"latitude"','"longitude"'):
  assert sensitive not in txt
 memo=DOC.read_text(encoding="utf8")
 si=SI.read_text(encoding="utf8")
 for v in ("16,801","1,665","0.62%","0.000232","−0.000082"):
  assert v in memo and v in si,v
 assert "### S2.16 |" in si
 assert "not" in memo.lower()
 assert "fitness" in memo.lower()
 assert "not" in q_limit()
def q_limit():
 return " ".join(json.loads(RESULT.read_text(encoding="utf8"))["limitations"]).lower()
