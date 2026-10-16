#!/usr/bin/env python3
"""Species-equal monthly photographed-white vs dynamic climate, controlling cell FE.

PRIMARY: for each original WorldClim 10m climate cell × source observation
calendar month, average the within-species observation-inclusion-weighted
classifiable photographic white fractions, one species vote each. This permits
seasonal TURNOVER of photographed flowering species, as user intended.
Cells absorb latitude/longitude, soil, elevation and terrain at the cell level.
Monthly climate T, solar, rainfall and derived VPD are source climatology,
NOT weather at historical observation year.

Two model diagnostics: cell fixed effects only; and cell PLUS calendar-month
fixed effects (different climate season profiles among the same geographic
cells). Weighted Frisch-Waugh-Lovell partialling handles large cell counts
without dense site dummy matrices. When rank/support is insufficient, refuse
to print four climate partial-effect coefficients (do NOT fill missing with 0).
CI, causal effects, genetic colour or true blooming population rates are NOT
established by this engineering module.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np

try:
    from scripts.acquisition.mois_stratified_photo_selection_20261016 import sha256_file
    from scripts.analysis.mois_join_monthly_dynamic_climate_20261016 import COLOUR,REAL_MANIFEST_SHA
except ModuleNotFoundError:
    import sys
    sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
    from scripts.acquisition.mois_stratified_photo_selection_20261016 import sha256_file
    from scripts.analysis.mois_join_monthly_dynamic_climate_20261016 import COLOUR,REAL_MANIFEST_SHA

SCHEMA="mois_species_equal_photo_white_monthly_same_cell_dynamic_FE_v1"
STATUS="SUPPORT_VERIFIED_MONTHLY_SOURCE_CELL_FE_NO_SOURCE_POPULATION_OR_CAUSAL_INFERENCE"
DYNAMIC=[
    "delta_tavg_c","delta_srad_kj_m2_day","delta_prec_mm","delta_vpd_kpa"
]
REQUIRED={
    "measurement_id","observation_id","inat_taxon_id","genus","domain","original_month",
    "climate_grid_row","climate_grid_col","source_climate_join_status",
    "sample_observation_weight","measurement_status","morph",*DYNAMIC
}


def partial_out(y_x:np.ndarray,w:np.ndarray,labels:list[np.ndarray],
                *,maximum_rounds=1500,tolerance=1e-10)->tuple[np.ndarray,int]:
    """Alternating exact group-weighted projections, Frisch-Waugh-Lovell."""
    if y_x.ndim!=2 or len(y_x)!=len(w) or np.any(w<=0):
        raise ValueError("Weighted design matrix geometry invalid")
    a=np.asarray(y_x,dtype=np.float64).copy()
    denoms=[np.bincount(ids,weights=w) for ids in labels]
    for rounds in range(1,maximum_rounds+1):
        for group,denom in zip(labels,denoms):
            for j in range(a.shape[1]):
                mean=np.bincount(group,weights=w*a[:,j],
                                minlength=len(denom))/np.maximum(denom,1e-15)
                a[:,j]-=mean[group]
        max_gradient=0.
        for group,denom in zip(labels,denoms):
            for j in range(a.shape[1]):
                means=np.bincount(group,weights=w*a[:,j],
                                 minlength=len(denom))/np.maximum(denom,1e-15)
                max_gradient=max(max_gradient,float(np.max(np.abs(means))))
        if max_gradient<tolerance:return a,rounds
    raise ValueError("Sparse cell×month fixed effect graph did not converge")


def fit_fe(cell_month:list[dict],*,include_month:bool,min_cells:int=10,
           zone_month:bool=False)->dict:
    if min_cells<2:raise ValueError("Require >=2 independent climate cells for model diagnostics")
    if zone_month and not include_month:
        raise ValueError("Latitude-zone × calendar-month FE requires month FE")
    months_by_cell=defaultdict(set)
    for row in cell_month:
        months_by_cell[row["cell"]].add(row["month"])
    repeated={cell for cell,months in months_by_cell.items() if len(months)>=2}
    a=[row for row in cell_month if row["cell"] in repeated]
    name=("cell_and_latitude_zone_month" if zone_month else
          "cell_and_month" if include_month else "cell_only")
    basic={
        "fixed_effect_model":name,
        "n_qualified_geographic_cells":len(repeated),
        "n_observed_cell_months_in_repeated_cells":len(a),
        "n_species_cell_month_votes":sum(row["n_species"] for row in a),
        "stable_geographic_soil_elevation_terrain_cell_means_absorbed":True,
        "separate_static_predictor_coefficients_not_identified":True,
        "model_not_a_selection_or_fitness_test":True,
    }
    if len(repeated)<min_cells or len(a)<min_cells+5:
        return {**basic,"status":"NOT_ESTIMABLE_TOO_FEW_REPEAT_MONTH_GEOGRAPHIC_CELLS"}
    cat_cell=np.unique([r["cell"] for r in a],return_inverse=True)[1]
    labels=[cat_cell]
    if include_month:
        if zone_month:
            if any(r.get("domain") not in (
                "NORTH_EXTRATROPICAL","SOUTH_EXTRATROPICAL","TROPICAL") for r in a):
                raise ValueError("Original observation hemisphere/tropics source domain missing")
            calendar=[str(r["domain"])+":"+str(r["month"]) for r in a]
        else:
            calendar=[r["month"] for r in a]
        labels.append(np.unique(calendar,return_inverse=True)[1])
    w=np.array([r["n_species"] for r in a],dtype=np.float64)
    y=np.array([r["white_species_equal"] for r in a],dtype=np.float64)
    x=np.array([[float(r[v]) for v in DYNAMIC] for r in a],dtype=np.float64)
    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)) \
       or np.any(y<0) or np.any(y>1):
        raise ValueError("Original photographed-white species fractions or climate incomplete")
    centered,iterations=partial_out(np.column_stack((y,x)),w,labels)
    yc,xc=centered[:,0],centered[:,1:]
    scale=np.sqrt(np.average(xc*xc,axis=0,weights=w))
    one={}
    for col,sd,var in zip(range(len(DYNAMIC)),scale,DYNAMIC):
        denom=float(np.dot(w*xc[:,col],xc[:,col]))
        one[var]={
            "status":"DESCRIPTIVE_UNIVARIATE_ONLY" if denom>1e-16 else "NO_INDEPENDENT_WITHIN_CELL_VARIATION",
            "slope_pp_per_original_monthly_unit":(
                float(100*np.dot(w*xc[:,col],yc)/denom) if denom>1e-16 else None
            ),
            "source_weighted_partialled_climate_SD":float(sd),
        }
    base={**basic,
          "partialling_converged_after_rounds":iterations,
          "fitted_white_response":"equal named species means among four-hue-classifiable original photos, then n-species weighted across climate cells",
          "monthly_X_are_site_month_minus_all_twelve_month_original_climatological_means":True,
          "univariate_axis_descriptive":one,
          "cluster_bootstrap_or_external_prediction_not_done":True}
    if any(v<=1e-12 for v in scale):
        return {**base,"status":"NOT_ESTIMABLE_SOME_DYNAMIC_CLIMATE_AXES_NO_VARIATION",
                "four_dynamic_partial_effect_coefficients":None}
    norm=xc/scale
    gram=(norm.T*w)@norm
    condition=float(np.linalg.cond(gram))
    rank=int(np.linalg.matrix_rank(norm*np.sqrt(w[:,None]),tol=1e-9))
    base["joint_dynamic_4axis_weighted_condition_number"]=condition
    base["joint_dynamic_4axis_effective_rank"]=rank
    if rank<4 or not math.isfinite(condition) or condition>1e7:
        return {**base,"status":"NOT_ESTIMABLE_DYNAMIC_AXES_COLLINEAR",
                "four_dynamic_partial_effect_coefficients":None}
    beta=np.linalg.solve(gram,norm.T@(w*yc))
    residual=yc-norm@beta
    null=float(np.dot(w*yc,yc))
    base["within_fixed_effect_adjusted_variation_R2"]=(
        float(1-np.dot(w*residual,residual)/null) if null>0 else None)
    base["four_dynamic_partial_effect_coefficients"]={
        name:{
            "slope_pp_per_original_monthly_climate_unit":float(100*b/sd),
            "slope_pp_per_weighted_within_cell_month_climate_SD":float(100*b),
            "weighted_partialled_climate_SD_original_units":float(sd)
        }
        for name,b,sd in zip(DYNAMIC,beta,scale)
    }
    return {**base,"status":"ESTIMATED_EXPLORATORY_POINT_ONLY_NO_CLUSTER_UNCERTAINTY"}


def analyze(climate_csv:Path,climate_receipt:Path,out:Path,*,
            synthetic=False,min_cells=10)->dict:
    if out.exists():raise ValueError("Refuse original source model receipt overwrite")
    source=json.loads(climate_receipt.read_text(encoding="utf8"))
    if (source.get("schema")!="mois_original_photo_all_monthly_dynamic_climate_site_join_v1"
        or source.get("status")!="COMPLETE_EXACT_ORIGINAL_SOURCE_PHOTO_MONTHLY_T_S_P_VAPR_JOIN_NOT_13ENV_REGRESSION"
        or source.get("new_private_photo_monthly_climate_csv_sha256")!=sha256_file(climate_csv)
        or source.get("synthetic_test_only") is not synthetic
        or source.get("seasonal_photo_white_environment_regression_fitted") is not False
        or not synthetic and source.get("original_48_raster_manifest_sha256")!=REAL_MANIFEST_SHA):
        raise ValueError("Not a full SHA-bound selected-original-photo 48-month climate measurement")
    original=Counter()
    taxa=defaultdict(lambda:[0.,0.,0])
    monthly_x={}
    monthly_domain={}
    seen=set();observations=set()
    with climate_csv.open(newline="",encoding="utf8") as f:
        reader=csv.DictReader(f)
        if not REQUIRED.issubset(reader.fieldnames or []):
            raise ValueError("Missing exact original source month/photo climate or probability")
        for row in reader:
            pid=row["measurement_id"];oid=row["observation_id"]
            if not pid or pid in seen or not oid or oid in observations:
                raise ValueError("Duplicated or unverified original independent-photo observation")
            seen.add(pid);observations.add(oid)
            original["all"]+=1
            four=row["measurement_status"]=="classified_four_state_morph"
            if four!=(row["morph"] in COLOUR):
                raise ValueError("Unclassifiable photograph cannot count as coloured or white")
            original["class"]+=int(four)
            if not four:continue
            original["white"]+=int(row["morph"]=="white")
            if row["source_climate_join_status"]!="VALID":
                original["class_photo_climate_missing"]+=1
                continue
            month=int(row["original_month"])
            if not 1<=month<=12:
                raise ValueError("Original observation month not one of twelve climate observations")
            cell=(row["climate_grid_row"],row["climate_grid_col"])
            if not all(x.isdecimal() for x in cell):
                raise ValueError("No source climate-cell identity")
            weight=float(row["sample_observation_weight"])
            if not math.isfinite(weight) or weight<1:
                raise ValueError("Invalid original independent observation inclusion probability")
            q=tuple(float(row[k]) for k in DYNAMIC)
            if not all(math.isfinite(v) for v in q):
                raise ValueError("Missing monthly source climate for classified photo")
            pair=(cell,month)
            if pair in monthly_x and any(abs(a-b)>1e-4 for a,b in zip(monthly_x[pair],q)):
                raise ValueError("Different climate values for same month and original WorldClim cell")
            monthly_x[pair]=q
            dom=row["domain"]
            if dom not in ("NORTH_EXTRATROPICAL","SOUTH_EXTRATROPICAL","TROPICAL"):
                raise ValueError("Original month-grid latitude domain missing")
            if pair in monthly_domain and monthly_domain[pair]!=dom:
                raise ValueError("One climate grid/month spans incompatible source latitude zones")
            monthly_domain[pair]=dom
            key=(cell,month,row["inat_taxon_id"])
            x=taxa[key]
            x[0]+=weight
            x[1]+=weight*int(row["morph"]=="white")
            x[2]+=1
    if original["all"]!=source["selected_original_photographs"] \
       or original["class"]!=source["source_four_colour_classifiable_photos"] \
       or original["white"]!=source["source_classifiable_white_photos"]:
        raise ValueError("Original photo source classifier denominator drift")
    units=defaultdict(list)
    for (cell,month,_taxon),(w,white,obs_n) in taxa.items():
        units[(cell,month)].append(white/w)
    aggregates=[]
    for (cell,month),fractions in sorted(units.items()):
        vals=monthly_x[(cell,month)]
        aggregates.append({
            "cell":f"{cell[0]}_{cell[1]}","month":month,
            "domain":monthly_domain[(cell,month)],
            "white_species_equal":sum(fractions)/len(fractions),
            "n_species":len(fractions),
            **dict(zip(DYNAMIC,vals))
        })
    months_by_cell=defaultdict(set)
    for row in aggregates:
        months_by_cell[row["cell"]].add(row["month"])
    enough=sum(len(v)>=2 for v in months_by_cell.values())
    model_a=fit_fe(aggregates,include_month=False,min_cells=min_cells)
    model_b=fit_fe(aggregates,include_month=True,min_cells=min_cells)
    model_c=fit_fe(aggregates,include_month=True,min_cells=min_cells,zone_month=True)
    public={
        "schema":SCHEMA,"status":STATUS,
        "verified_private_climate_photo_source_receipt_sha256":sha256_file(climate_receipt),
        "n_original_photos_in_full_source":original["all"],
        "n_classifiable_four_hue":original["class"],
        "n_classifiable_white":original["white"],
        "n_classifiable_with_missing_monthly_climate":original["class_photo_climate_missing"],
        "n_species_cell_month_rows":len(taxa),
        "n_month_source_climate_cells":len(aggregates),
        "n_distinct_climate_cells_with_classified_flowers":len(months_by_cell),
        "n_climate_cells_with_at_least_two_distinct_classified_flower_months":enough,
        "three_exploratory_fixed_effect_support_and_models":[model_a,model_b,model_c],
        "target_is_photographed_blooming_species_turnover_not_individual_colour_plasticity":True,
        "static_soil_elevation_latitude_terrain_coefficients_absorbed_by_exact_cell_FE":True,
        "unobserved_original_wild_flowering_population_estimated":False,
        "selection_bias_corrected_beyond_source_probability":False,
        "no_cluster_or_genus_uncertainty_in_this_point_estimator":True,
        "inference_needs_independent_geographic_and_genus_cluster_robustness":True,
        "synthetic_test_only":bool(synthetic)
    }
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(public,indent=2,sort_keys=True)+"\n",encoding="utf8")
    return public


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--private-monthly-original-photo-climate",required=True,type=Path)
    p.add_argument("--original-monthly-source-receipt",required=True,type=Path)
    p.add_argument("--out-public-supported-model",required=True,type=Path)
    p.add_argument("--min-qualified-repeat-climate-cells",type=int,default=10)
    a=p.parse_args()
    r=analyze(a.private_monthly_original_photo_climate,a.original_monthly_source_receipt,
         a.out_public_supported_model,min_cells=a.min_qualified_repeat_climate_cells)
    print(json.dumps({k:r[k] for k in (
        "status","n_species_cell_month_rows",
        "n_climate_cells_with_at_least_two_distinct_classified_flower_months"
    )},sort_keys=True))


if __name__=="__main__":main()
