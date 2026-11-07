#!/usr/bin/env python3
"""Observed same accepted taxon across mainland/island: sparse taxon fixed effects.

Analyzes ONLY original source-verified WCVP accepted taxon x photo-cell x year
strict9 white ROI measurements, and same-year WorldPop R2024B UC + VIIRS.
No source species is imputed at missing geographic sites. Prior photo-archive
taxon mean reference is NOT recycled: each observed source taxon's original
ROI9 value is the response, taxon FE identify its empirical variation.

Two cohorts: frozen strict 7,982 cell-years and measured two-accepted-taxa
20,543 cell-years. Each observed cell has total weight 1, divided equally
among accepted crossland taxa photographed there. Taxa must have genuine
original photographed ROI measurements on BOTH land classes; no inference of
true biological presence/absence from photographic nondetection.

Fit weighted OLS with species/year/region 5-degree fixed effects using
sparse weighted Frisch-Waugh-Lovell projection, geographic and accepted-taxon
cluster-robust sandwich SE, plus two-way cluster as qualified diagnostic.
Intercept is absorbed by fixed effects, not separately estimated.

This is not an IV, causal ALAN estimate, genetic species colour evolution,
MNAR missingness correction, or external independently verified floral trait.
"""
from __future__ import annotations
import argparse, json, math
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse.linalg import lsmr
from scipy.stats import t as student_t
from scripts.analysis import mois_island_ALAN_measured_photo_opportunity_tiers_20261107 as tier
from scripts.analysis import mois_island_ALAN_multiyear_sameyear_source_panel_20261107 as panel
from scripts.analysis import mois_island_ALAN_taxon_turnover_measured_photo_20261107 as taxa

TERM=("log_light","is_island","log_light:is_island","log_pop",
      "bio5","log_rain","elevation","log_attempts")
SPECIES=taxa.SPECIES_COL

def true_crossland_accepted_taxa(species, photo_cells):
    """No unobserved colour labels, no one-sided species allowed."""
    subset=species.loc[species.target_id.isin(set(photo_cells.target_id))].copy()
    group=subset.groupby(SPECIES,observed=True).is_island.agg(["min","max","size"])
    shared=set(group.index[(group["min"]==0)&(group["max"]==1)])
    both=subset.loc[subset[SPECIES].isin(shared)].copy()
    panel.require(both.target_id.isin(photo_cells.target_id).all(),
                  "Shared accepted taxa were not genuinely photographed")
    panel.require(both[SPECIES].nunique()==len(shared),
                  "No actual photographed species support for model")
    panel.require(not both.duplicated(["target_id",SPECIES]).any(),
                  "Taxon x original photographic cell-year not unique")
    panel.require(both.photo_white_fraction.between(0,1).all(),
                  "No source-verified strict9 photo colour for taxon")
    if len(both):
        n_in_cell=both.groupby("target_id",observed=True)[SPECIES].transform("size")
        both["photo_cell_equal_species_weight"]=1/n_in_cell
        panel.require(np.allclose(
            both.groupby("target_id").photo_cell_equal_species_weight.sum().values,1),
            "Original photograph cells not each equal weighted")
    return both,{
       "source_accepted_taxa_both_land_classes":len(shared),
       "actual_taxon_cell_year_measurements":len(both),
       "photo_cell_years_with_shared_accepted_taxa":int(both.target_id.nunique()),
       "source_island_taxon_cell_years":int(both.is_island.eq(1).sum()),
       "source_island_photo_cells":int(
           both.loc[both.is_island.eq(1)].target_id.nunique()),
       "source_geographic_5deg_blocks":int(both.region5.nunique()),
       "source_photo_years":sorted(set(int(x) for x in both.year)),
       "true_original_photo_WHITE_ROI_only":True,
       "not_detected_elsewhere_is_not_true_absence":True,
    }

def sparse_fixed_effect_design(frame,cols=("species","region5","year")):
    """Full one-hot FE with no intercept; minimum-norm sparse solver."""
    n=len(frame)
    rowid=np.arange(n,dtype=np.int32)
    arrays=[];levels={}
    for col in cols:
        cats,names=pd.factorize(frame[col],sort=True)
        panel.require(np.min(cats)>=0,"Source fixed effect contains NaN")
        ncol=len(names)
        z=sparse.coo_matrix((np.ones(n),(rowid,cats)),
                            shape=(n,ncol)).tocsr()
        arrays.append(z)
        levels[col]=ncol
    return sparse.hstack(arrays,format="csr"),levels

def residualize_against_FEs(z,weights,variables):
    """Return weighted outcomes and predictors after WLS-FWL residualization."""
    panel.require(np.isfinite(weights).all() and (weights>0).all(),
                  "Invalid original cell-equal weights")
    sw=np.sqrt(weights)
    weighted_Z=z.multiply(sw[:,None]).tocsr()
    v=np.asarray(variables,float)*sw[:,None]
    r=np.empty_like(v)
    convergence=[]
    for i in range(v.shape[1]):
        coeff,stop,iters,normr,*rest=lsmr(
            weighted_Z,v[:,i],atol=1.e-10,btol=1.e-10,maxiter=1800)
        panel.require(stop in (0,1,2,4,5) and np.isfinite(normr),
                      f"Sparse fixed-effect projection did not converge (variable {i}, stop {stop})")
        r[:,i]=v[:,i]-weighted_Z@coeff
        convergence.append({"stop_code":int(stop),"iterations":int(iters),
                            "weighted_projection_residual":float(normr)})
    return r,convergence

def cluster_meat(scores,categories):
    """Sandwich meat from original measured clusters, no fake independent photos."""
    code,uniques=pd.factorize(categories,sort=True)
    panel.require((code>=0).all(),"Cluster without actual source identity")
    pieces=np.column_stack([
        np.bincount(code,weights=scores[:,i],minlength=len(uniques))
        for i in range(scores.shape[1])])
    return pieces.T@pieces,len(uniques)

def cluster_uncertainty(beta,bread,scores,groups,n,p):
    meat,g=cluster_meat(scores,groups)
    panel.require(g>5,"Too few independently supported cluster units")
    correction=(g/(g-1))*((n-1)/(n-p))
    var=bread@meat@bread*correction
    var=(var+var.T)/2
    return var,g

def interaction_stats(beta,variance,clusters,n):
    """Nominal 95% t interval based on cluster df; exploratory p-value."""
    variance=float(variance)
    if not math.isfinite(variance) or variance<=0:
        return {"estimable":False,
                "reason":"nonpositive_or_nonfinite_cluster_sandwich_variance",
                "variance":variance,"clusters":int(clusters)}
    se=math.sqrt(variance)
    tstat=beta/se
    df=int(clusters)-1
    radius=float(student_t.ppf(.975,df)*se)
    return {"estimable":True,"beta":float(beta),"SE":se,
            "ci_low":float(beta-radius),"ci_high":float(beta+radius),
            "p":float(2*student_t.sf(abs(tstat),df)),
            "df":df,"clusters":int(clusters)}

def fit_sparse_taxon_FE(frame,effect_cols=("species","region5","year")):
    """Same genuine source rows; optionally add accepted-taxon FE to geography/year."""
    panel.require(set(effect_cols) in ({"species","region5","year"},{"region5","year"}),
                  "Only source-based nested taxon vs non-taxon fixed-effect models")
    panel.require(len(frame)>=300 and frame.region5.nunique()>=10 and
                  frame[SPECIES].nunique()>=20 and frame.is_island.sum()>=40,
                  "Too little actual same-taxonomy island/mainland support")
    d=frame.copy()
    d["species"]=d[SPECIES].astype(str)
    d["year"]=d.year.astype(str)
    d["interaction"]=d.log_light*d.is_island
    variables=np.column_stack([d.log_light,d.is_island,d.interaction,
                 d.log_pop,d.bio5,d.log_rain,d.elevation,d.log_attempts,
                 100*d.photo_white_fraction]).astype(float)
    panel.require(np.isfinite(variables).all(),
                  "True original WCVP photos have nonfinite predictors")
    weights=d.photo_cell_equal_species_weight.to_numpy(float)
    fe,levels=sparse_fixed_effect_design(d,cols=effect_cols)
    partial,solver_diagnostics=residualize_against_FEs(fe,weights,variables)
    x=partial[:,:-1]
    y=partial[:,-1]
    gram=x.T@x
    rank=int(np.linalg.matrix_rank(gram))
    panel.require(rank==len(TERM),"Taxon x land nightlight model not identified")
    bread=np.linalg.inv(gram)
    beta=bread@(x.T@y)
    e=y-x@beta
    scores=x*e[:,None]
    v_geo,n_geo=cluster_uncertainty(beta,bread,scores,d.region5,len(d),x.shape[1])
    v_taxon,n_taxa=cluster_uncertainty(beta,bread,scores,d.species,len(d),x.shape[1])
    v_intersection,n_pairs=cluster_uncertainty(
        beta,bread,scores,d.region5.astype(str)+"_"+d.species,len(d),x.shape[1])
    v_two=(v_geo+v_taxon-v_intersection)
    j=TERM.index("log_light:is_island")
    support={
       "observed_accepted_taxa":int(d.species.nunique()),
       "original_taxon_photo_cell_years":len(d),
       "original_photo_cell_years":int(d.target_id.nunique()),
       "island_taxon_photo_cell_years":int(d.is_island.sum()),
       "source_photo_years":sorted(set(int(a) for a in d.year)),
       "geographic_5degree_clusters":n_geo,
       "source_accepted_taxon_clusters":n_taxa,
       "observed_geographic_taxon_intersection_clusters":n_pairs,
       "fixed_effect_cardinalities":levels,
       "exact_fixed_effect_definition":list(effect_cols),
       "weighted_cell_sum_per_cell_is_one":True,
       "genuine_source_photo_roi_not_missingness_label":True,
       "no_external_genetic_or_true_native_species_phenotype_inference":True,
       "fwls_solver_iteration_max":max(x["iterations"] for x in solver_diagnostics),
       "fwls_max_projection_residual":max(x["weighted_projection_residual"] for x in solver_diagnostics),
       "nominal_after_exploratory_models":True,
    }
    return {
      "supported":True,
      "source":support,
      "interaction_geographic_spatial_cluster":interaction_stats(beta[j],v_geo[j,j],n_geo,len(d)),
      "interaction_accepted_taxon_cluster":interaction_stats(beta[j],v_taxon[j,j],n_taxa,len(d)),
      "interaction_two_way_spatial_plus_taxon_cluster":interaction_stats(
         beta[j],v_two[j,j],min(n_geo,n_taxa),len(d)),
      "fixed_effects_are_actual_accepted_taxa_region5_and_observation_year":("species" in effect_cols),
      "cell_weight_is_one_across_all_shared_taxa":True,
      "year_geography_may_still_confound_unmatched_distant_source_areas":True,
    }

def run(paths,out):
    out=Path(out);panel.require(not out.exists(),"Do not overwrite original source results")
    frames=tier.build_frozen_five_years(paths)
    cohorts=tier.by_year_tier(frames)
    strict=cohorts["strict_5p3r3s"]
    broader=cohorts["measured_3p2r2s"]
    panel.require(len(strict)==7982 and len(broader)==20543 and
                  set(strict.target_id).issubset(set(broader.target_id)),
                  "Original real measured photo source drift")
    canonical,wcvp=panel.load_frozen_photos(paths)
    species=taxa.real_species_cell_years(paths,canonical,wcvp,broader)
    panel.require(len(species)==66175,
                  "Previously SHA-verified original 66,175 source accepted-taxon photo-cell measurements changed")
    result={}
    for label,cohort in [("strict_original",strict),("measured_2species",broader)]:
        direct,original_support=true_crossland_accepted_taxa(species,cohort)
        fit=fit_sparse_taxon_FE(direct)
        without_taxon=fit_sparse_taxon_FE(direct,effect_cols=("region5","year"))
        original_beta=fit["interaction_geographic_spatial_cluster"]["beta"]
        no_taxon_beta=without_taxon["interaction_geographic_spatial_cluster"]["beta"]
        result[label]={
            "source_support":original_support,
            "source_same_taxon_FE_fit":fit,
            "identical_source_rows_without_taxon_FE":without_taxon,
            "taxon_FE_addition_point_estimate_change":float(original_beta-no_taxon_beta),
            "same_source_rows_weights_photographic_outcomes_all_covariates_year_and_region":True,
            "coefficient_difference_has_no_independent_standard_error_or_causal_meaning":True,
        }
    out.mkdir(parents=True)
    receipt={
      "schema":"MOIS_same_accepted_taxon_both_land_classes_sparse_FE_20261107_v1",
      "observed_years":[2019,2021,2023,2024,2025],
      "original_accepted_taxon_photo_cell_years":len(species),
      "original_strict_cell_years":len(strict),
      "original_measured_two_accepted_species_cell_years":len(broader),
      "all_original_cohorts_reused_without_photo_reinference":True,
      "real_strict9_roi_only":True,
      "observational_island_ALAN_not_causal":True,
      "source_accepted_taxon_is_not_verified_same_genotype":True,
      "missing_species_photo_has_not_been_imputed_or_called_absence":True,
      "no_MNAR_tipping_or_other_withdrawn_sensitivity":True,
      "source_input_SHA256":{k:panel.sha(v) for k,v in paths.items()},
      "results":result,
    }
    (out/"same_taxon_island_ALAN_FE_source_receipt.json").write_text(
         json.dumps(receipt,sort_keys=True,indent=2)+"\n")
    print("MOIS_TRUE_SAME_ACCEPTED_TAXON_SPARSE_FE",
          json.dumps(result,sort_keys=True),flush=True)
    return receipt

if __name__=="__main__":
    p=argparse.ArgumentParser()
    keys=("canonical","wcvp",
      "geo2019","viirs2019","pop2019",
      "geo2021","viirs2021","pop2021",
      "geo2023","viirs2023","pop2023",
      "geo2024","cn2024","uc2024","uc2024receipt",
      "geo2025","geo2025receipt","viirs2025")
    for key in keys:p.add_argument("--"+key.replace("_","-"),dest=key,type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    run({key:getattr(a,key) for key in keys},a.out)
