#!/usr/bin/env python3
"""MOIS original photographic flower-white ROI: taxon turnover versus taxon-relative colour.

Genuine original strict9 ROI and WCVP *accepted species IDs* are joined to the
previously SHA-locked 2019, 2021, 2023, 2024, 2025 same-year VIIRS/WorldPop
R2024B UC records and original GSHHG island/mainland geometry.

The community-level photographed white ROI component can be decomposed exactly
into (a) mean photographically observed source taxon white ROI among the taxa
actually found at a cell and (b) cell-specific deviations from those source
taxon means. This is descriptive not genetic/species functional inference:
source taxon means are estimated from the SAME photo archive, singletons have
zero residuals, and species replacement and photo detection remain confounded.

The original strict 5photos/3accepted spp regression and expanded 3photos/
2accepted spp regression remain unchanged. No missing-colour imputation, image
downloads, night moth causal claim or rejected MNAR tipping scenario.
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd
from scripts.analysis import mois_island_ALAN_measured_photo_opportunity_tiers_20261107 as tier
from scripts.analysis import mois_island_ALAN_multiyear_sameyear_source_panel_20261107 as panel
from scripts.analysis.mois_2024_sameyear_global_island_photo_white_20261107 import (
    cell_id,cluster_fit,COVARS)

YEARS=(2019,2021,2023,2024,2025)
SPECIES_COL="wcvp_accepted_plant_name_id"
STRICT_N=7982
BROAD_N=20543

def real_species_cell_years(paths,canonical,wcvp,broad):
    """Return original (year, 0.05° photographed cell, real WCVP accepted taxon).

    Independent photo->cell reconstruction with strict single-GSHHG-land
    evidence, and verify reconstruction exactly matches historical
    pre-aggregated species counts and ROI values at all broad cells.
    """
    p=canonical.loc[
        canonical.source_subarm.eq("historical_main") &
        canonical.strict9 &
        canonical.strict9_white.notna() &
        canonical.latitude.between(-90,90) &
        canonical.longitude.between(-180,180)].copy()
    p["photo_year"]=pd.to_numeric(p.observed_on.astype(str).str[:4],errors="coerce")
    p=p.loc[p.photo_year.isin(YEARS)].copy()
    lookup=wcvp.loc[
        wcvp.wcvp_match_state.eq("unique_exact_name_to_accepted_id"),
        ["inat_taxon_id",SPECIES_COL]].copy()
    panel.require(lookup.inat_taxon_id.is_unique,"WCVP source taxa duplicate")
    p=p.merge(lookup,on="inat_taxon_id",validate="many_to_one",how="inner")
    p=p.loc[p[SPECIES_COL].notna()].copy()
    panel.require(p.photo_id.is_unique,"An original photograph was counted twice")
    p["cell_id"]=cell_id(p.longitude,p.latitude)
    species_parts=[]
    for year in YEARS:
        g=pd.read_csv(paths[f"geo{year}"],
                      usecols=["cell_id","target_id","land_class"],
                      dtype={"cell_id":str,"target_id":str,"land_class":str})
        # Same geographic ambiguity gate as the old cell-outcome estimation.
        ambiguity=g.cell_id.duplicated(keep=False)
        pure=g.loc[(~ambiguity)&g.land_class.isin(("island","mainland")),
                   ["cell_id","target_id"]].copy()
        panel.require(pure.cell_id.is_unique,"Unresolved original GSHHG cell ambiguity")
        valid=broad.loc[broad.year.eq(year),["target_id"]].copy()
        pure=pure.merge(valid,on="target_id",how="inner",validate="one_to_one")
        py=p.loc[p.photo_year.eq(year),
                 ["photo_id","cell_id",SPECIES_COL,"strict9_white"]].copy()
        py=py.merge(pure,on="cell_id",how="inner",validate="many_to_one")
        panel.require(py.photo_id.is_unique,"Photo duplicated by source land assignment")
        spp=py.groupby(["target_id",SPECIES_COL],observed=True).agg(
            photo_white_fraction=("strict9_white","mean"),
            actual_ROI_photos=("photo_id","nunique")).reset_index()
        spp["year"]=year
        species_parts.append(spp)
    all_species=pd.concat(species_parts,ignore_index=True)
    source_broad=broad[["target_id","year","region5","land_id","is_island",
        "white_pp","log_light","log_pop","bio5","log_rain","elevation",
        "log_attempts"]].copy()
    species=all_species.merge(source_broad,on=["target_id","year"],
                             validate="many_to_one",how="inner")
    panel.require(len(species)==len(all_species),"Original accepted species-cell rows lost")
    panel.require(species.groupby("target_id")[SPECIES_COL].nunique().eq(
        broad.set_index("target_id").accepted_species_n.astype(int)
    ).all(),"Independent source WCVP accepted-species counts do not match")
    mean=species.groupby("target_id").photo_white_fraction.mean().mul(100)
    pd.testing.assert_index_equal(mean.index.sort_values(),
        broad.set_index("target_id").white_pp.index.sort_values(),
        check_names=False)
    aligned=broad.set_index("target_id").loc[mean.index,"white_pp"]
    panel.require(np.allclose(mean.to_numpy(),aligned.to_numpy(),rtol=0,atol=1e-8),
                  "Independent original photo white ROI values do not match cell colours")
    return species

def accepted_taxon_photo_support(species,eligible,shared):
    """Counts are source-observed opportunities; missing means unknown, not absent."""
    frame=species.loc[species.target_id.isin(set(eligible.target_id))]
    cls=frame.groupby(SPECIES_COL).is_island.agg(["min","max","count"])
    both=set(cls.index[(cls["min"]==0)&(cls["max"]==1)])
    all_ids=set(frame[SPECIES_COL])
    common=frame.loc[frame[SPECIES_COL].isin(both)]
    return {
        "actual_accepted_taxa_observed":len(all_ids),
        "accepted_taxa_observed_in_both_island_and_mainland":len(both),
        "observed_taxon_cell_years":len(frame),
        "taxon_cell_years_with_both_land_class_observations":len(common),
        "cells_with_any_taxon_observed_in_both_land_classes":int(
           frame.groupby("target_id")[SPECIES_COL].apply(
               lambda s:bool(set(s)&both)).sum()),
        "island_photo_cells_with_shared_taxon":int(
           common.loc[common.is_island.eq(1)].target_id.nunique()),
        "original_island_photo_cell_years":int(eligible.is_island.sum()),
        "shared_taxa_across_land_classes_are_not_verified_same_genotype":True,
    }

def source_taxon_reference_white_components(species):
    """Exact algebraic decomposition of observed taxon-equal cell white ROI.

    A species' reference is its unweighted mean across eligible photo
    cell-years (NOT a species' true genetically-fixed colour).
    """
    taxon=species.groupby(SPECIES_COL,observed=True).agg(
        source_mean_photo_white_fraction=("photo_white_fraction","mean"),
        measured_cell_years=("target_id","nunique"),
        in_island=("is_island","max")).reset_index()
    d=species.merge(taxon[[SPECIES_COL,"source_mean_photo_white_fraction",
                          "measured_cell_years"]],
                    on=SPECIES_COL,validate="many_to_one")
    ref=d.groupby("target_id",observed=True).agg(
         taxon_average_white_fraction=("source_mean_photo_white_fraction","mean"),
         source_singleton_taxa=("measured_cell_years",lambda s:int((s==1).sum())),
         taxon_species_count=(SPECIES_COL,"nunique")).reset_index()
    return taxon,ref

def source_composition_decomposition_fit(frame,taxonref):
    m=frame.merge(taxonref,on="target_id",how="left",validate="one_to_one")
    panel.require(m.taxon_average_white_fraction.notna().all() and
                  m.taxon_species_count.eq(m.accepted_species_n).all(),
                  "Original photo-source taxon composition reference incomplete")
    m["taxon_source_mean_white_pp"]=100*m.taxon_average_white_fraction
    m["photo_deviation_from_taxon_source_mean_pp"]=(
        m.white_pp-m.taxon_source_mean_white_pp)
    panel.require(np.allclose(
        m.white_pp,
        m.taxon_source_mean_white_pp+m.photo_deviation_from_taxon_source_mean_pp,
        rtol=0,atol=1e-9),"Original photographed white cannot be decomposed exactly")
    formula=" ~ log_light * is_island + log_pop + "+COVARS+" + C(region5) + C(year)"
    observed=cluster_fit(m,"white_pp"+formula,"region5")["log_light:is_island"]
    community=cluster_fit(m,"taxon_source_mean_white_pp"+formula,
                          "region5")["log_light:is_island"]
    photo_deviation=cluster_fit(m,"photo_deviation_from_taxon_source_mean_pp"+formula,
                                "region5")["log_light:is_island"]
    panel.require(abs(observed["beta"]-community["beta"]-
                      photo_deviation["beta"])<0.00001,
                  "Original photo colour interaction decomposition failed")
    taxon_n=m.taxon_species_count.sum()
    singleton_n=m.source_singleton_taxa.sum()
    return {
        "observed_white_ROI":observed,
        "taxon_reference_component":community,
        "taxon_relative_photographic_component":photo_deviation,
        "cell_years":len(m),
        "observed_accepted_taxon_cell_years":int(taxon_n),
        "observed_source_singleton_taxon_cell_years":int(singleton_n),
        "taxon_mean_plus_deviation_beta_equals_original_exactly":True,
        "taxon_reference_is_from_same_photo_archive_not_external_trait":True,
        "taxon_relative_deviation_for_singleton_is_zero_by_construction":True,
        "this_is_descriptive_association_decomposition_not_causal_mediation":True
    }

def genuinely_crossland_taxon_5deg_year_support(species,eligible):
    """Can one same accepted species be directly seen on BOTH land classes
    inside the same observation year AND the same 5-degree region?"""
    x=species.loc[species.target_id.isin(set(eligible.target_id))].copy()
    keys=["year","region5",SPECIES_COL]
    g=x.groupby(keys,observed=True).is_island.agg(["min","max","size"])
    paired=g.loc[(g["min"]==0)&(g["max"]==1)]
    # Pair original observed species-cell-years only, no taxon imputation.
    tuples=set(paired.index)
    idx=pd.MultiIndex.from_frame(x[keys])
    match=x.loc[idx.isin(tuples)].copy()
    n_islands=int(match.loc[match.is_island.eq(1)].target_id.nunique())
    n_regions=int(match.region5.nunique())
    out={
       "paired_same_taxon_year_region_groups":len(paired),
       "independent_5degree_blocks":n_regions,
       "taxa_directly_seen_on_both_island_mainland_within_year_region":int(
           match[SPECIES_COL].nunique()),
       "actual_photo_species_cell_years":len(match),
       "actual_island_photo_species_cell_years":int(match.is_island.sum()),
       "island_photo_cells_with_matched_same_taxon":n_islands,
       "not_recorded_elsewhere_does_not_mean_species_absent":True,
    }
    # This is an opportunistic feasibility gate: not a global species-FE causal
    # model. With too few regions the coefficient would extrapolate wildly.
    if n_regions<10 or len(paired)<20 or match.is_island.sum()<35 or len(match)<100:
        return {**out,"estimable":False,
                "reason":"insufficient independent source-accepted species within same-year geographic island/mainland pairs"}
    import statsmodels.formula.api as smf
    match["real_pair"]=match.year.astype(str)+"_"+match.region5+"_"+match[SPECIES_COL].astype(str)
    match["taxon_photo_white_pp"]=100*match.photo_white_fraction
    formula="taxon_photo_white_pp ~ log_light * is_island + log_pop + "+COVARS+" + C(real_pair)"
    fit=smf.ols(formula,data=match).fit(
        cov_type="cluster",cov_kwds={"groups":match.region5},use_t=True)
    t="log_light:is_island"
    panel.require(t in fit.params and np.isfinite(fit.params[t]),
                  "No actual same-species island×light interaction identifiable")
    low,high=fit.conf_int().loc[t]
    return {**out,"estimable":True,
      "interaction":{"beta":float(fit.params[t]),"SE":float(fit.bse[t]),
                     "ci_low":float(low),"ci_high":float(high),
                     "p":float(fit.pvalues[t])},
      "year_region_accepted_species_fixed_effects":True,
      "geographic_5degree_clustered_errors":True,
      "not_observationally_independent_same_species_genotypes":True}

def run(paths,out):
    out=Path(out)
    panel.require(not out.exists(),"Do not overwrite published source receipts")
    frames=tier.build_frozen_five_years(paths)
    strict=tier.by_year_tier(frames)["strict_5p3r3s"]
    broad=tier.by_year_tier(frames)["measured_3p2r2s"]
    panel.require(len(strict)==STRICT_N and len(broad)==BROAD_N and
                  set(strict.target_id).issubset(set(broad.target_id)),
                  "Original 5-year strict/broader source cohorts changed")
    canonical,wcvp=panel.load_frozen_photos(paths)
    species=real_species_cell_years(paths,canonical,wcvp,broad)
    taxon,taxonref=source_taxon_reference_white_components(species)
    main={}
    for label,frame in (("strict_original",strict),("measured_two_species",broad)):
        main[label]={
          "taxonomic_photo_support":accepted_taxon_photo_support(
              species,frame,True),
          "global_cell_white_ROI_decomposition":
              source_composition_decomposition_fit(frame,taxonref),
          "same_accepted_species_sameyear_5degree_island_mainland":
              genuinely_crossland_taxon_5deg_year_support(species,frame)
        }
    panel.require(abs(
        main["strict_original"]["global_cell_white_ROI_decomposition"]["observed_white_ROI"]["beta"]
        - (-2.4970445996320447))<1e-5,"Strict original island×light slope not unchanged")
    panel.require(abs(
        main["measured_two_species"]["global_cell_white_ROI_decomposition"]["observed_white_ROI"]["beta"]
        - (-1.204806413294285))<1e-5,"Broader genuine colour cohort changed")
    source_taxa_shared=species.groupby(SPECIES_COL).is_island.nunique()
    reported={
        "schema":"MOIS_original_photo_taxon_turnover_ALAN_island_sameyear_20261107_v1",
        "source_photo_years":list(YEARS),
        "strict_source_cell_years":STRICT_N,
        "broader_2species_source_cell_years":BROAD_N,
        "actual_original_source_accepted_taxa":int(len(taxon)),
        "original_source_taxa_seen_both_island_and_mainland":int(source_taxa_shared.ge(2).sum()),
        "actual_species_photo_cell_years":len(species),
        "results":main,
        "exact_original_ROI_and_WCVP_white_values_reconstructed":True,
        "no_new_RGB_download_or_ROI_reinference":True,
        "no_photo_unknown_imputation_or_MNAR_tipping":True,
        "no_species_loss_nocturnal_pollination_or_ALAN_causality_inference":True,
        "species_reference_trait_is_observed_photographic_average_only":True,
        "source_input_SHA256":{k:panel.sha(v) for k,v in paths.items()}
    }
    out.mkdir(parents=True)
    (out/"taxon_turnover_and_within_taxon_photo_association_receipt.json").write_text(
        json.dumps(reported,sort_keys=True,indent=2)+"\n")
    print("MOIS_REAL_TAXON_TURNOVER_ISLAND_ALAN_WHITE",json.dumps({
        k:reported[k] for k in ("actual_original_source_accepted_taxa",
        "original_source_taxa_seen_both_island_and_mainland",
        "actual_species_photo_cell_years","results")},sort_keys=True),flush=True)
    return reported
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
