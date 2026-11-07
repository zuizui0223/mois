#!/usr/bin/env python3
"""MOIS actual repeat-photo site-year island × ALAN temporal association.

One original true 0.05-degree cell + original GSHHG landmass = one site.
A photo-cell's presence in several observation years is NOT independent site
replication; site fixed effects remove all time-constant geographic attributes
and the observed photographic plant assemblage may still change over time.

Use ONLY the preexisting original photographic strict9/WCVP accepted-species
white ROI outcomes, actual 2019/2021/2023/2024/2025 same-year VIIRS and
WorldPop R2024B UC population. No new photographs, ROI inference, imputed
white or missing species occurrence, substituted 2022 October/provider values,
2020 Global1 population source or rejected MNAR/tipping calculations.

Reported coefficients are descriptive temporal associations of source
photographic colour, not causal ALAN, genetic species shifts, or true field
population trajectories. Keep distinct strict and measured 2-species arms.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np,pandas as pd
import statsmodels.api as sm
from scripts.analysis import mois_island_ALAN_measured_photo_opportunity_tiers_20261107 as tier
from scripts.analysis import mois_island_ALAN_multiyear_sameyear_source_panel_20261107 as panel

YEARS=(2019,2021,2023,2024,2025)
SOURCE_N={"strict_5p3r3s":7982,"measured_3p2r2s":20543}
MIN_ISLAND_SITES=25
MIN_ISLAND_LANDMASSES=8
MIN_GEO_BLOCKS=10

def authentic_site_identity(d):
    """A year-only repeated ID cannot make a new location."""
    p=d.copy()
    parts=p.target_id.astype(str).str.split("|",regex=False,expand=True)
    panel.require(parts.shape[1]==3 and parts[0].eq(p.year.astype(str)).all(),
                  "Unverified photo cell source year/target ID")
    panel.require(parts[1].str.contains("_",regex=False).all(),
                  "Missing original 0.05-degree photo cell ID")
    panel.require(parts[2].eq(p.land_id.astype(str)).all(),
                  "Original GSHHG landmass ID does not match target source")
    p["original_cell_id"]=parts[1].to_numpy()
    p["repeat_site_id"]=parts[1].astype(str)+"|"+parts[2].astype(str)
    panel.require(not p.duplicated(["repeat_site_id","year"]).any(),
                  "Original same site duplicated in one photo year")
    invariant=p.groupby("repeat_site_id",observed=True).agg(
        land_classes=("land_class","nunique"),
        landmasses=("land_id","nunique"),
        regions=("region5","nunique"),
        photo_cells=("original_cell_id","nunique"))
    panel.require((invariant==1).all().all(),
                  "Geographic/landmass membership varies for the same original site")
    return p

def repeat_site_frame(d):
    s=authentic_site_identity(d)
    grouped=s.groupby("repeat_site_id",observed=True)
    site_years=grouped.year.nunique()
    repeated=s.loc[s.repeat_site_id.isin(set(site_years.index[site_years>=2]))].copy()
    panel.require(not repeated.duplicated(["repeat_site_id","year"]).any(),
                  "Repeated original site photo year duplicated")
    panel.require(repeated.groupby("repeat_site_id").year.nunique().ge(2).all(),
                  "Single-year photo site accidentally included in within-site estimator")
    return s,repeated

def audit_support(allrows,repeat):
    def tally(q):
        island=q.loc[q.is_island.eq(1)]
        return {
           "photo_cell_years":int(len(q)),
           "independent_original_sites":int(q.repeat_site_id.nunique()),
           "geographic_5degree_blocks":int(q.region5.nunique()),
           "island_cell_years":len(island),
           "independent_original_island_sites":int(island.repeat_site_id.nunique()),
           "independent_island_landmasses":int(island.land_id.nunique()),
           "island_geographic_blocks":int(island.region5.nunique()),
           "mainland_cell_years":int((1-q.is_island).sum()),
        }
    support={"all_original_year_censused_photo_cells":tally(allrows),
             "same_geographic_cell_and_GSHHG_landmass_repeated_2plus_years":tally(repeat)}
    by_site=repeat.groupby("repeat_site_id",observed=True).agg(
        land_class=("land_class","first"),
        source_land_id=("land_id","first"),
        original_cell=("original_cell_id","first"),
        first_year=("year","min"),
        last_year=("year","max"),
        observed_years=("year","nunique"),
        year_light_span=("log_light",lambda v:float(v.max()-v.min())),
        year_white_span_pp=("white_pp",lambda v:float(v.max()-v.min())))
    out={}
    for cls in ("island","mainland"):
        s=by_site.loc[by_site.land_class.eq(cls)]
        spans=s.year_light_span
        out[cls]={
          "repeated_original_sites":len(s),
          "photo_year_count":{"2":int(s.observed_years.eq(2).sum()),
                              "3":int(s.observed_years.eq(3).sum()),
                              "4":int(s.observed_years.eq(4).sum()),
                              "5":int(s.observed_years.eq(5).sum())},
          "median_light_log_range":float(spans.median()) if len(s) else None,
          "sites_with_log_light_range_ge_0p05":int(spans.ge(.05).sum()),
          "sites_with_log_light_range_ge_0p20":int(spans.ge(.20).sum()),
          "median_observed_white_roi_span_pp":float(s.year_white_span_pp.median()) if len(s) else None,
        }
    support["measured_repeated_site_exposure_support"]=out
    support["no_assumption_that_same_site_photos_contain_same_plants"]=True
    support["not_same_individuals_or_verified_same_species_within_site"]=True
    return support

def site_demean_with_year_FE(frame):
    """Residualize true site intercept without dense thousands-of-site dummy matrix.

    Year indicators are also demeaned by site, so OLS is algebraically
    equivalent to site FE + calendar-year FE regression for the regressors.
    Constant BIO5, rainfall and elevation are absorbed by site effects and
    cannot be separately identified in a repeated-site fixed-effect design.
    """
    d=frame.copy().sort_values(["repeat_site_id","year"]).reset_index(drop=True)
    panel.require(d.repeat_site_id.nunique()>=2 and
                  d.groupby("repeat_site_id").year.nunique().ge(2).all(),
                  "Site fixed effects require real repeated observation years")
    d["island_loglight"]=d.log_light*d.is_island
    covars=["log_light","island_loglight","log_pop","log_attempts"]
    # Exact original years, not an interpolated annual series.
    indicators=[f"year_{year}" for year in YEARS[1:]]
    for year,col in zip(YEARS[1:],indicators):
        d[col]=d.year.eq(year).astype(float)
    cols=["white_pp",*covars,*indicators]
    panel.require(np.isfinite(d[cols].to_numpy(float)).all(),
                  "Original photo time-series contains nonfinite data")
    group=d.groupby("repeat_site_id",observed=True)
    adjusted=d[cols]-group[cols].transform("mean")
    x=adjusted[covars+indicators].to_numpy(float)
    y=adjusted.white_pp.to_numpy(float)
    return d,x,y,cols

def estimate_within_site(d):
    """Estimate genuine changes of light and photographic whiteness at same site.

    Inference clustered on the same 5-degree geographic block across years.
    Unbalanced repeated-site panel: no unobserved site-years are imputed.
    """
    _,r=repeat_site_frame(d)
    support=audit_support(authentic_site_identity(d),r)
    origin=support["same_geographic_cell_and_GSHHG_landmass_repeated_2plus_years"]
    requirements={
       "island_sites_ge_25":origin["independent_original_island_sites"]>=MIN_ISLAND_SITES,
       "island_landmasses_ge_8":origin["independent_island_landmasses"]>=MIN_ISLAND_LANDMASSES,
       "geographic_clusters_ge_10":origin["geographic_5degree_blocks"]>=MIN_GEO_BLOCKS,
       "both_land_classes_observed":0<origin["island_cell_years"]<origin["photo_cell_years"],
       "both_land_classes_have_exposure_change":all(
           support["measured_repeated_site_exposure_support"][t]["sites_with_log_light_range_ge_0p05"]>=5
           for t in ("island","mainland")),
    }
    if not all(requirements.values()):
        return {
          "estimable":False,
          "reason":"source_repeated_site_or_temporal_exposure_support_gate",
          "requirements":requirements,
          "source_support":support,
        }
    panel.require(r.repeat_site_id.nunique()>=50 and r.region5.nunique()>=10,
                  "Repeated source geometry below independent support gate")
    ordered,x,y,columns=site_demean_with_year_FE(r)
    names=["log_light","log_light:is_island","log_pop","log_attempts",
           *[f"year_{a}" for a in YEARS[1:]]]
    scale=np.linalg.norm(x,axis=0)
    # Drop only covariates that lack any site-within variation. Never create
    # a pseudo-estimable light slope by fitting a rank-deficient design.
    kept=[i for i,v in enumerate(scale) if v>1e-8]
    cols=[names[i] for i in kept]
    panel.require("log_light" in cols and "log_light:is_island" in cols,
                  "Not enough within-site original light changes to test island contrast")
    panel.require("log_pop" in cols or "log_attempts" in cols,
                  "No measured within-site variation in any nuisance covariate")
    xx=x[:,kept]
    rank=int(np.linalg.matrix_rank(xx))
    if rank<len(kept):
        return {
          "estimable":False,"reason":"within_site_site_and_year_FE_design_rank_deficient",
          "rank":rank,"predictor_columns":cols,
          "source_support":support
        }
    fit=sm.OLS(y,xx).fit(cov_type="cluster",
        cov_kwds={"groups":ordered.region5},
        use_t=True)
    terms={}
    for term in ("log_light","log_light:is_island"):
        i=cols.index(term)
        lo,hi=fit.conf_int()[i]
        terms[term]={"beta":float(fit.params[i]),"SE":float(fit.bse[i]),
                     "ci_low":float(lo),"ci_high":float(hi),"p":float(fit.pvalues[i])}
    return {
      "estimable":True,"source_support":support,
      "source_model":"white_pp ~ log_light + log_light:is_island + log_pop + log_attempts + site FE + observed-year FE; site-invariant BIO5, rainfall, elevation absorbed",
      "geographic_5degree_cluster_count":int(ordered.region5.nunique()),
      "repeated_source_original_photo_sites":int(ordered.repeat_site_id.nunique()),
      "model_rank":rank,"model_columns":cols,
      "terms":terms,
      "observational_temporal_not_causal_light_experiment":True,
      "photographed_source_taxa_may_change_even_when_site_does_not":True,
      "no_interpolation_or_composition_imputation_between_photo_years":True,
    }

def run(paths,out):
    out=Path(out);panel.require(not out.exists(),
        "Refuse overwrite of original source evidence")
    frames=tier.build_frozen_five_years(paths)
    cohorts=tier.by_year_tier(frames)
    main={}
    for name,expected in SOURCE_N.items():
        d=cohorts[name]
        panel.require(len(d)==expected and set(d.year)==set(YEARS),
                      "Original frozen photo yearly cohort/source count changed")
        main[name]=estimate_within_site(d)
    report={
      "schema":"MOIS_original_photo_same_site_repeatedyear_ALAN_white_v1",
      "true_observed_photo_years":list(YEARS),
      "same_year_R2024B_UC_and_original_monthly_VIIRS":True,
      "original_strict_source_photo_cell_years":7982,
      "original_measured_two_species_source_photo_cell_years":20543,
      "only_true_repeated_GSHHG_landmass_and_0p05degree_cells":True,
      "source_photo_ROI_white_not_repeated_individuals_nor_botanical_species_abundance":True,
      "no_missing_site_year_source_value_imputed":True,
      "no_new_original_photographs_or_ROI_inference":True,
      "no_MNAR_or_old_tipping_sensitivity":True,
      "results":main,
      "source_SHA256":{k:panel.sha(v) for k,v in paths.items()},
    }
    out.mkdir(parents=True)
    (out/"same_site_multiyear_photo_white_ALAN_receipt.json").write_text(
        json.dumps(report,indent=2,sort_keys=True)+"\n")
    print("MOIS_REAL_REPEATED_SAME_SITE_YEAR_NIGHTLIGHT_ISLAND",
          json.dumps(main,sort_keys=True),flush=True)
    return report

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
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
