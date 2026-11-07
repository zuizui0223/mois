#!/usr/bin/env python3
"""Original repeated MOIS photo sites: genuinely photographed accepted-taxon overlap.

Quantify what fraction of genuine WCVP accepted taxa photographed at the SAME
0.05-degree cell+GSHHG landmass overlap between its different source years.
Then measure photographic strict9 white ROI responses for the SAME accepted
taxon in the SAME site across ≥2 actually observed photo years, with
site×taxon and year fixed effects, and true-year VIIRS/WorldPop covariates.

This is an ecological observation-opportunity analysis of photographed taxa,
NOT true flora turnover: no photographed record != biological absence.
Likewise, individual plants/genotypes/flower anthesis and nocturnal service
are not tracked across years. No new original images or missing colours.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np,pandas as pd
import statsmodels.api as sm
from scripts.analysis import mois_island_ALAN_measured_photo_opportunity_tiers_20261107 as tier
from scripts.analysis import mois_island_ALAN_multiyear_sameyear_source_panel_20261107 as panel
from scripts.analysis import mois_island_ALAN_taxon_turnover_measured_photo_20261107 as tax
from scripts.analysis import mois_island_ALAN_repeated_sites_photo_year_20261108 as same_site

EXPECTED={"strict_5p3r3s":7982,"measured_3p2r2s":20543}
MIN_ISLAND_REPEAT_TAXON_SITES=30
MIN_ISLAND_LANDMASSES=8
MIN_SPATIAL_BLOCKS=10

def true_site_species_years(species,cohort):
    """Reconstruct original GSHHG+grid and WCVP taxa without faking absence."""
    source=cohort[["target_id","year","land_id","land_class","region5",
                   "white_pp","is_island","log_light","log_pop",
                   "log_attempts"]].copy()
    original=same_site.authentic_site_identity(source)
    d=species.loc[species.target_id.isin(set(cohort.target_id))].copy()
    panel.require(d.target_id.nunique()==len(cohort),
                  "Genuine taxon ROI records do not span all original measured photo cells")
    d=d.merge(original[["target_id","original_cell_id","repeat_site_id"]],
              on="target_id",how="inner",validate="many_to_one")
    panel.require(not d.duplicated(["target_id",tax.SPECIES_COL]).any(),
                  "One original accepted species repeated in same photographed cell")
    panel.require(d.photo_white_fraction.between(0,1).all(),
                  "Unobserved original strict9 white ROI entered species-site chronology")
    d["species_repeat_site_id"]=d.repeat_site_id+"|"+d[tax.SPECIES_COL].astype(str)
    panel.require(not d.duplicated(["species_repeat_site_id","year"]).any(),
                  "Same accepted taxon/photograph site has duplicate calendar year")
    return original,d

def observed_taxon_overlap(repeated_cells,species):
    """Between consecutive ACTUALLY PHOTOGRAPHED years in same real site."""
    d=species.loc[species.target_id.isin(set(repeated_cells.target_id))].copy()
    observed=d.groupby(["repeat_site_id","year"],observed=True)[tax.SPECIES_COL].agg(
        lambda a:frozenset(str(v) for v in a)).reset_index(name="taxa")
    site_years=repeated_cells.groupby("repeat_site_id").year.nunique()
    panel.require(len(observed)==len(repeated_cells),
                  "Species photos cannot reconstruct repeated site-year source")
    bysite={site:g.sort_values("year") for site,g in observed.groupby("repeat_site_id",sort=False)}
    rows=[]
    mapping=repeated_cells.drop_duplicates("repeat_site_id").set_index("repeat_site_id")
    for site,group in bysite.items():
        panel.require(len(group)>=2 and len(group)==site_years.loc[site],
                      "Observed species photo years do not match original repeated site")
        vals=list(group[["year","taxa"]].itertuples(index=False,name=None))
        for (y0,t0),(y1,t1) in zip(vals[:-1],vals[1:]):
            panel.require(y1>y0,"Original photographed observation years not ordered")
            shared=len(t0&t1)
            union=len(t0|t1)
            panel.require(union>0,"Source accepted photos have empty species lists")
            rows.append({
               "land_class":mapping.loc[site,"land_class"],
               "year_gap":int(y1-y0),
               "observed_taxa_first":len(t0),
               "observed_taxa_last":len(t1),
               "observed_taxa_shared":shared,
               "observed_taxa_union":union,
               "jaccard":shared/union,
               "no_same_taxon_both_photographic_years":shared==0,
            })
    return pd.DataFrame(rows)

def species_site_repeats(d):
    counts=d.groupby("species_repeat_site_id",observed=True).year.nunique()
    good=set(counts.index[counts>=2])
    rep=d.loc[d.species_repeat_site_id.isin(good)].copy()
    panel.require(not rep.duplicated(["species_repeat_site_id","year"]).any(),
                  "Source species-site-year never allowed to duplicate")
    if len(rep):
        panel.require(rep.groupby("species_repeat_site_id").year.nunique().ge(2).all(),
                      "Nonrepeat source accepted taxon entered within-taxon estimand")
    rep["white_pp"]=100*rep.photo_white_fraction
    return rep

def support_for_species_site_repeats(all_species,repeat_cells,rep,pairs):
    island=rep.loc[rep.is_island.eq(1)]
    all_island=repeat_cells.loc[repeat_cells.is_island.eq(1)]
    out={
       "genuine_repeated_photo_cell_years":int(len(repeat_cells)),
       "genuine_repeated_photo_island_cell_years":int(len(all_island)),
       "observed_taxon_year_pair_comparisons":int(len(pairs)),
       "actual_accepted_taxon_site_years_at_repeated_photographed_sites":int(len(all_species.loc[
           all_species.target_id.isin(set(repeat_cells.target_id))])),
       "genuine_repeated_accepted_taxon_site_years":len(rep),
       "genuine_repeated_accepted_taxon_sites":int(rep.species_repeat_site_id.nunique()),
       "island_repeated_accepted_taxon_site_years":len(island),
       "island_repeated_accepted_taxon_sites":int(island.species_repeat_site_id.nunique()),
       "distinct_island_landmasses_with_repeat_taxa":int(island.land_id.nunique()),
       "source_5degree_spatial_blocks":int(rep.region5.nunique()),
       "source_accepted_taxa_repeated_in_identical_locations":int(rep[tax.SPECIES_COL].nunique()),
       "not_detected_next_year_means_photo_unknown_not_biological_absence":True,
    }
    for label in ("island","mainland"):
        x=pairs.loc[pairs.land_class.eq(label)]
        out[label]={
          "consecutive_actual_photo_year_pairs":len(x),
          "pairs_with_no_same_accepted_taxon_recorded_both_years":int(
              x.no_same_taxon_both_photographic_years.sum()),
          "fraction_zero_observed_taxon_overlap":float(
              x.no_same_taxon_both_photographic_years.mean()) if len(x) else None,
          "median_photographed_taxon_jaccard":float(x.jaccard.median()) if len(x) else None,
          "mean_photographed_taxon_jaccard":float(x.jaccard.mean()) if len(x) else None,
          "median_calendar_year_gap_between_photographs":float(x.year_gap.median()) if len(x) else None,
        }
    return out

def fit_within_real_species_site(d):
    """Actual source white pixel ROI within the SAME accepted species at SAME site."""
    frame=d.copy()
    site=frame.groupby("species_repeat_site_id",observed=True)
    y_counts=site.year.nunique()
    panel.require(y_counts.ge(2).all(),
                  "No single-year pseudo species-site source units")
    if len(frame):
        lightspans=site.log_light.agg(lambda x:float(x.max()-x.min()))
    else:
        lightspans=pd.Series(dtype=float)
    islands=frame.loc[frame.is_island.eq(1)]
    requirements={
       "island_species_sites_ge30":islands.species_repeat_site_id.nunique()>=MIN_ISLAND_REPEAT_TAXON_SITES,
       "island_landmasses_ge8":islands.land_id.nunique()>=MIN_ISLAND_LANDMASSES,
       "spatial_blocks_ge10":frame.region5.nunique()>=MIN_SPATIAL_BLOCKS,
       "both_land_categories":int(frame.is_island.sum())>0 and
           int(frame.is_island.sum())<len(frame),
       "at_least_10_taxon_sites_with_real_light_range_ge0p05":int(lightspans.ge(.05).sum())>=10,
    }
    if not all(requirements.values()):
        return {"estimable":False,"requirements":requirements,
                "reason":"true_same_site_same_taxon_multi_year_support_below_identification_gate"}
    frame["repeat_site_id"]=frame.species_repeat_site_id
    # The accepted taxon and original site are both held constant here.
    ordered,x,y,columns=same_site.site_demean_with_year_FE(frame)
    names=["log_light","log_light:is_island","log_pop","log_attempts",
           *[f"year_{a}" for a in same_site.YEARS[1:]]]
    scales=np.linalg.norm(x,axis=0)
    kept=[i for i,v in enumerate(scales) if v>1e-8]
    predictors=[names[i] for i in kept]
    if not all(n in predictors for n in ("log_light","log_light:is_island")):
        return {"estimable":False,
                "reason":"source_same_taxon_same_site_light_exposure_without_identifiable_variance"}
    xx=x[:,kept]
    rank=int(np.linalg.matrix_rank(xx))
    if rank!=len(kept):
        return {"estimable":False,
                "reason":"source_same_taxon_year_site_FE_collinear",
                "rank":rank,"cols":predictors}
    fit=sm.OLS(y,xx).fit(cov_type="cluster",
        cov_kwds={"groups":ordered.region5},use_t=True)
    i=predictors.index("log_light:is_island")
    lo,hi=fit.conf_int()[i]
    return {
      "estimable":True,
      "interaction":{
         "beta":float(fit.params[i]),"SE":float(fit.bse[i]),
         "ci_low":float(lo),"ci_high":float(hi),"p":float(fit.pvalues[i])},
      "source": {
         "geographic_5degree_cluster_count":int(ordered.region5.nunique()),
         "genuine_accepted_taxon_same_site_repeated_years":int(len(ordered)),
         "genuine_accepted_taxon_same_sites":int(ordered.repeat_site_id.nunique()),
         "model_rank":rank,
         "model_columns":predictors,
         "within_species_site_FE_removes_constant_taxon_identity_and_geography":True,
         "genuine_photo_year_fixed_effects":True,
         "only_true_source_photo_ROI_and_sameyear_VIIRS_population":True,
      },
      "same_site_same_taxon_not_same_photographed_individual":True,
      "possible_photo_acquisition_and_phenology_changes_not_removed":True,
      "not_causal_ALAN_or_stable_species_colour_evolution":True,
    }

def run(paths,out):
    out=Path(out);panel.require(not out.exists(),
                                "Do not overwrite source photographic ecology results")
    frames=tier.build_frozen_five_years(paths)
    cohorts=tier.by_year_tier(frames)
    original,wcvp=panel.load_frozen_photos(paths)
    broad=cohorts["measured_3p2r2s"]
    panel.require(len(broad)==20543,"Original measured two accepted species photo source changed")
    actual=tax.real_species_cell_years(paths,original,wcvp,broad)
    panel.require(len(actual)==66175,"Original photographed taxon×site×year results drift")
    res={}
    for key,n in EXPECTED.items():
        frame=cohorts[key]
        panel.require(len(frame)==n,"Original strict/expanded photo source changed")
        source,species=true_site_species_years(actual,frame)
        all_source,repeat=same_site.repeat_site_frame(source)
        pairs=observed_taxon_overlap(repeat,species)
        repeated_taxon=species_site_repeats(species.loc[
            species.target_id.isin(set(repeat.target_id))])
        info=support_for_species_site_repeats(species,repeat,repeated_taxon,pairs)
        fit=fit_within_real_species_site(repeated_taxon)
        res[key]={"source_support":info,"same_accepted_taxon_same_site_year_FE":fit}
    receipt={
        "schema":"MOIS_original_repeat_photo_site_and_accepted_taxon_year_ALAN_white_v1",
        "actual_source_years":list(same_site.YEARS),
        "strict_original_photo_cell_years":7982,
        "expanded_original_photo_cell_years":20543,
        "original_WCVP_taxon_photo_cell_years":66175,
        "observation_taxon_Jaccard_is_not_field_plant_turnover":True,
        "same_source_photographed_taxon_in_same_site_not_same_individual":True,
        "no_2020_or_2022_substituted_population_or_radiance":True,
        "no_new_photographs_or_ROI_inference":True,
        "no_technical_unknown_colour_label_imputation":True,
        "no_MNAR_tipping_or_IPW_scenarios":True,
        "true_original_input_SHA256":{k:panel.sha(v) for k,v in paths.items()},
        "results":res
    }
    out.mkdir(parents=True)
    (out/"repeated_photo_taxon_site_year_light_receipt.json").write_text(
        json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print("MOIS_REAL_REPEATED_SAME_TAXON_SAME_SITE_YEAR",
          json.dumps(res,sort_keys=True),flush=True)
    return receipt

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
