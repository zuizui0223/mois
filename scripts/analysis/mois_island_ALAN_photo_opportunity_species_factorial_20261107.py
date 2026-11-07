#!/usr/bin/env python3
"""MOIS original measured white ROI × ALAN: photographic opportunity factorial.

This is *not* a missingness correction or MNAR tipping analysis. All original
strict9-successful ROI measurements remain untouched. Compares an explicitly
fixed 2×2 observed eligibility design:
  strict: 5+ photos, 3+ measured accepted species (with >=3 strict9 successes)
  relaxed photos only: 3+ photos, 3+ measured species
  relaxed accepted-species count only: 5+ photos, 2+ measured species
  both relaxed: 3+ photos, 2+ measured species.

This separates photo-effort and observable accepted-species coverage from the
previously observed strict-to-moderate effect attenuation. All cohorts share
the same real 2019/2021/2023/2024/2025 VIIRS and WorldPop R2024B UC years,
source IDs, GSHHG. No 2022 VIIRS October substitution, 2020 Global1 mixing,
image downloads, inferred white species fractions, or technical unknown fills.

Careful: accepted_species_n <= strict9_success by construction. Hence the
3+ species / 3+ strict9 and 2+ species / 2+ strict9 constraints have redundant
ROI counts. Report the observed identity instead of falsely attributing
attenuation to an independently variable strict9 success gate.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scripts.analysis import mois_island_ALAN_measured_photo_opportunity_tiers_20261107 as tier
from scripts.analysis import mois_island_ALAN_multiyear_sameyear_source_panel_20261107 as panel
from scripts.analysis.mois_2024_sameyear_global_island_photo_white_20261107 import (
    cluster_fit,COVARS)

DESIGN={
  "strict_5photos_3species":(5,3,3),
  "photos_only_3photos_3species":(3,3,3),
  "species_only_5photos_2species":(5,2,2),
  "both_3photos_2species":(3,2,2)
}
STRICT_N=7982
BOTH_N=20543
BASE="white_pp ~ log_light * is_island + log_pop + "+COVARS

def cohort(frames,requirements):
    all_years=[]
    for year, (full, _) in frames.items():
        d=tier.verified_colour_frame(year,full,requirements)
        panel.require(d.target_id.is_unique and d.white_pp.between(0,100).all(),
                      "Original photographic ROI/taxon photo-cell source identity drift")
        all_years.append(d)
    x=pd.concat(all_years,ignore_index=True)
    panel.require(x.target_id.is_unique and x.year.nunique()==5,
                  "Year/source identities duplicated when photo gates relaxed")
    return x

def matched_year_regions(frame):
    shared,groups=tier.exact_year_block_support(frame)
    return shared,groups

def landmass_counts(frame):
    islands=frame.loc[frame.is_island.eq(1)]
    return int(islands.land_id.nunique())

def diagnostic_counts(frame):
    # Species require a successful genuine photographed ROI, so accepted
    # species ≤ successful unique photos even without gate.
    invalid=frame.loc[frame.accepted_species_n>frame.strict9_success]
    panel.require(len(invalid)==0,"Accepted WCVP species outnumber strict9 ROI successes")
    shared,groups=matched_year_regions(frame)
    return {
      "photo_cell_years":int(len(frame)),
      "island_cell_years":int(frame.is_island.sum()),
      "mainland_cell_years":int((1-frame.is_island).sum()),
      "independent_island_landmasses_with_valid_photos":landmass_counts(frame),
      "spatial_5degree_blocks":int(frame.region5.nunique()),
      "paired_sameyear_region5_groups":int(len(groups)),
      "paired_distinct_region5_blocks":int(shared.region5.nunique()),
      "paired_original_photo_cell_years":int(len(shared)),
      "paired_island_cell_years":int(shared.is_island.sum()),
      "all_accepted_species_count_le_original_strict9_success":True,
    }

def genuine_light_slope(frame):
    return cluster_fit(frame,BASE+" + C(region5) + C(year)",
                       "region5")["log_light:is_island"]

def direct_sameyear_slope(frame):
    shared, groups=matched_year_regions(frame)
    geo_n=int(shared.region5.nunique())
    isl_n=int(shared.is_island.sum())
    if geo_n<10 or isl_n<40 or len(shared)<100:
        return {"estimable":False, "cause":"inadequate original matched sameyear geography",
                "geographical_clusters":geo_n,"island_cells":isl_n,
                "original_photo_cells":int(len(shared))}
    frame=shared.copy()
    frame["year_region5"]=frame.year.astype(str)+"__"+frame.region5
    return {"estimable":True, "geographical_clusters":geo_n,
            "sameyear_region5_groups":len(groups),
            "island_cells":isl_n,"original_photo_cells":int(len(shared)),
            "interaction":cluster_fit(frame,BASE+" + C(year_region5)",
                                      "region5")["log_light:is_island"]}

def attempt_species_bins(frame):
    """Disjoint source-opportunity bins; no pooled pseudo-photos."""
    x=frame.copy()
    x["attempt_bin"]=pd.cut(x.photo_attempts.astype(int),
        bins=[2,4,9,np.inf],labels=["3-4","5-9","10+"],include_lowest=False)
    x["species_bin"]=pd.cut(x.accepted_species_n.astype(int),
        bins=[1,2,4,np.inf],labels=["2","3-4","5+"],include_lowest=False)
    panel.require(x.attempt_bin.notna().all() and x.species_bin.notna().all(),
                  "Only real measured ≥3-photo/≥2-species source eligible")
    x["match_stratum"]=(x.year.astype(str)+"|"+x.region5+"|"+
                        x.attempt_bin.astype(str)+"|"+x.species_bin.astype(str))
    return x

def exact_observed_opportunity_overlap(frame):
    """Year×georegion×photo opportunity×species coverage paired on real data."""
    x=attempt_species_bins(frame)
    grouped=x.groupby("match_stratum",observed=True).is_island.agg(["min","max","count"])
    support=set(grouped.index[(grouped["min"]==0)&(grouped["max"]==1)])
    shared=x.loc[x.match_stratum.isin(support)].copy()
    island_n=int(shared.is_island.sum())
    geo_n=int(shared.region5.nunique())
    by=shared.groupby(["match_stratum","is_island"],observed=True).size().unstack(fill_value=0)
    pairs=int(np.minimum(by[0],by[1]).sum()) if len(by) else 0
    return shared,{
       "sameyear_region_photo_species_paired_strata":len(support),
       "distinct_geographical_clusters":geo_n,
       "observed_photo_cell_years":len(shared),
       "observed_island_cell_years":island_n,
       "balanced_pairs_possible_from_real_strata":pairs,
       "geographically_paired_original_measurement_not_synthetic":True,
    }

def matched_opportunity_estimate(frame):
    matched,report=exact_observed_opportunity_overlap(frame)
    if (report["distinct_geographical_clusters"]<10 or
        report["observed_island_cell_years"]<35 or
        report["balanced_pairs_possible_from_real_strata"]<30):
        return {**report, "estimable":False,
                "reason":"original island/mainland year×5degree×effort×species strata have insufficient overlap"}
    # The stratum FE is based ONLY on source photo opportunity and year/region;
    # colour and VIIRS values never determine matching membership.
    f=BASE+" + C(match_stratum)"
    groups=matched.groupby(["match_stratum","is_island"],observed=True).size()
    # Within every supported stratum, both island and mainland class totals
    # get equal aggregate weight; these are descriptive standardization
    # weights, not true individual propensity scores or an MNAR model.
    matched["n_class"]=pd.MultiIndex.from_frame(
        matched[["match_stratum","is_island"]]).map(groups)
    min_counts=groups.unstack().min(axis=1)
    matched["n_pair"]=matched.match_stratum.map(min_counts)
    matched["balanced_weight"]=matched.n_pair / matched.n_class
    panel.require(matched.balanced_weight.gt(0).all() and
                  matched.balanced_weight.le(1).all(),
                  "Source-observed strata weights must be nonnegative and ≤1")
    fit=smf.wls(f,data=matched,weights=matched.balanced_weight).fit(
        cov_type="cluster",cov_kwds={"groups":matched.region5},use_t=True)
    term="log_light:is_island"
    panel.require(term in fit.params and np.isfinite(float(fit.params[term])),
                  "Matched original geographic photo opportunity model not identified")
    low,high=fit.conf_int().loc[term]
    return {**report,"estimable":True,
        "coefficient":{"beta":float(fit.params[term]),"SE":float(fit.bse[term]),
                       "ci_low":float(low),"ci_high":float(high),
                       "p":float(fit.pvalues[term])},
        "comparison":"source-only year×5degree×3bin photo attempts×3bin accepted species",
        "stratum_equal_total_island_mainland_weight":True,
        "photo_attempts_and_species_count_not_direct_nuisance_selection_model":True,
        "not_eligible_for_population_causal_inference":True}

def disjoint_original_measured_source_classes(cohorts):
    """Observed four-cell partition, not an invented missingness category.

    Each cell-year is in exactly one class:
    - strict (5+ photo; 3+ accepted measured species);
    - new solely by permitting 3-4 photos while maintaining ≥3 species;
    - new solely by permitting 2 species while maintaining ≥5 photos;
    - newly allowed only by simultaneous 3-4 photos AND exactly 2 species.

    This makes the observed nonadditive 2×2 cohort enlargement inspectable
    without modifying original ROI or inventing missing species colours.
    """
    name={"strict":"strict_5photos_3species",
          "photo":"photos_only_3photos_3species",
          "species":"species_only_5photos_2species",
          "both":"both_3photos_2species"}
    full=cohorts[name["both"]]
    a=set(cohorts[name["strict"]].target_id)
    b=set(cohorts[name["photo"]].target_id)
    c=set(cohorts[name["species"]].target_id)
    d=set(full.target_id)
    panel.require(a==b.intersection(c),"The observed 2×2 strict cohort is not photo∩species")
    partition={
       "strict_original":a,
       "photo_requirement_only_relaxed":b-a,
       "species_requirement_only_relaxed":c-a,
       "joint_photo_and_species_relaxation_only":d-(b|c),
    }
    all_ids=set().union(*partition.values())
    panel.require(all_ids==d and
                  sum(map(len,partition.values()))==len(d),
                  "Disjoint observational source cohort classes lost/duplicated photos")
    out={}
    for label,ids in partition.items():
        rows=full.loc[full.target_id.isin(ids)]
        is_island=rows.is_island.eq(1)
        def summarize_land(x):
            return {
              "cell_years":int(len(x)),
              "white_ROI_mean_pp":float(x.white_pp.mean()) if len(x) else None,
              "median_source_VIIRS":float(x.viirs_median.median()) if len(x) else None,
              "median_source_original_photo_attempts":float(x.photo_attempts.median()) if len(x) else None,
              "median_observed_accepted_species":float(x.accepted_species_n.median()) if len(x) else None,
            }
        out[label]={
          "cell_years":len(rows),
          "island_cell_years":int(is_island.sum()),
          "mainland_cell_years":int((~is_island).sum()),
          "source_years":{str(int(y)):int(n) for y,n in rows.groupby("year").size().items()},
          "island":summarize_land(rows.loc[is_island]),
          "mainland":summarize_land(rows.loc[~is_island]),
          "source_measured_only_true_original_strict9_colour":True,
        }
    return out

def run(paths,out):
    out=Path(out)
    panel.require(not out.exists(),"Cannot overwrite original measured source cohorts")
    frames=tier.build_frozen_five_years(paths)
    cohorts={key:cohort(frames,gate) for key,gate in DESIGN.items()}
    strict=cohorts["strict_5photos_3species"]
    full=cohorts["both_3photos_2species"]
    panel.require(len(strict)==STRICT_N and len(full)==BOTH_N,
                  "Frozen 5-year strict or broad real photo population drift")
    first=set(strict.target_id); last=set(full.target_id)
    for key,frame in cohorts.items():
        panel.require(first.issubset(set(frame.target_id)) and
                      set(frame.target_id).issubset(last),
                      "Factorial nested source cohort was not actually observed: "+key)
    counts={key:diagnostic_counts(frame) for key,frame in cohorts.items()}
    disjoint_classes=disjoint_original_measured_source_classes(cohorts)
    models={key:{"global_interaction":genuine_light_slope(frame),
                 "sameyear5degree_comparison":direct_sameyear_slope(frame)}
            for key,frame in cohorts.items()}
    baseline=models["strict_5photos_3species"]["global_interaction"]["beta"]
    panel.require(abs(baseline+2.4970445996320447)<1e-5,
                  "Strict primary signed original coefficient must remain unchanged")
    panel.require(counts["strict_5photos_3species"]["paired_island_cell_years"]==119 and
                  counts["both_3photos_2species"]["paired_island_cell_years"]==408,
                  "Previously validated geographic-year support was changed")
    effects={key:m["global_interaction"]["beta"] for key,m in models.items()}
    changes={
        "photo_opportunity_gate_relaxed_at_3species":effects["photos_only_3photos_3species"]-baseline,
        "accepted_species_count_gate_relaxed_at_5photos":effects["species_only_5photos_2species"]-baseline,
        "both_relaxed_vs_strict":effects["both_3photos_2species"]-baseline,
        "nonadditive_difference_of_observed_cohort_estimates":
            effects["both_3photos_2species"]-effects["photos_only_3photos_3species"]-
            effects["species_only_5photos_2species"]+baseline,
        "not_separate_causal_effects_nor_independent_model_tests":True,
    }
    exact={
        "strict":matched_opportunity_estimate(strict),
        "both_relaxed":matched_opportunity_estimate(full)
    }
    out.mkdir(parents=True)
    receipt={
       "schema":"MOIS_20261107_true_source_photographic_opportunity_species_factorial_v1",
       "source_original_observation_years":[2019,2021,2023,2024,2025],
       "design":{k:{"photo_attempts_min":v[0],
                    "true_strict9_success_min":v[1],
                    "independently_WCVP_accepted_species_min":v[2]}
                 for k,v in DESIGN.items()},
       "observed_count_only_no_estimated_missing_colour":True,
       "strict_primary_reference_unchanged":True,
       "source_quality_tier_counts":counts,
       "disjoint_source_photo_species_2x2_classes":disjoint_classes,
       "source_quality_tier_models":models,
       "observed_gate_relaxation_coefficient_differences":changes,
       "same_year_geographic_photo_species_strata_common_support":exact,
       "no_unknown_to_white_or_nonwhite_imputation":True,
       "no_MNAR_tipping_or_IPW_scenarios":True,
       "no_original_photo_RGB_downloads_or_reinference":True,
       "photographic_ROI_whiteness_not_species_white_frequency":True,
       "causal_ALAN_and_pollinator_mechanisms_unmeasured":True,
       "source_input_SHA256":{k:panel.sha(p) for k,p in paths.items()},
    }
    path=out/"factorial_photo_opportunity_species_source_support.csv"
    rows=[{"gate":gate,**values} for gate,values in counts.items()]
    pd.DataFrame(rows).to_csv(path,index=False)
    receipt["output_support_SHA256"]=panel.sha(path)
    (out/"factorial_source_receipt.json").write_text(
        json.dumps(receipt,sort_keys=True,indent=2)+"\n")
    print("MOIS_OBSERVED_PHOTO_SPECIES_OPPORTUNITY_FACTORIAL",
          json.dumps({
             "counts":counts,"models":models,"disjoint_classes":disjoint_classes,
             "gate_coefficient_differences":changes,
             "exact_photo_species_common_support":exact},sort_keys=True),flush=True)
    return receipt

if __name__=="__main__":
    p=argparse.ArgumentParser()
    keys=("canonical","wcvp",
      "geo2019","viirs2019","pop2019",
      "geo2021","viirs2021","pop2021",
      "geo2023","viirs2023","pop2023",
      "geo2024","cn2024","uc2024","uc2024receipt",
      "geo2025","geo2025receipt","viirs2025")
    for key in keys:
        p.add_argument("--"+key.replace("_","-"),dest=key,type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    run({key:getattr(a,key) for key in keys},a.out)
