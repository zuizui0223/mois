#!/usr/bin/env python3
"""Increase genuine same-year island/mainland photographed colour common support.

Separate *source-measured* cohorts; no missing hue/ROI status is imputed.
  strict_core: original >=5 opportunities, >=3 actual strict9 successes,
               >=3 uniquely accepted species (frozen 7,982 cell-years).
  measured_2species: >=3 opportunities, >=2 strict9, >=2 accepted species.
  measured_1species: >=1 opportunity, >=1 strict9, >=1 accepted species;
                     SUPPORT COUNTS ONLY, not a confirmatory refit.

All same-year original 2019/21/23/24/25 landmasses, VIIRS and WP UC are
rechecked by the pre-existing source SHA-pinned full-panel loaders.
No 2022 month substitution, no 2020 Global1 mixing, no image downloads.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import pandas as pd
import numpy as np
from scripts.analysis import mois_island_ALAN_multiyear_sameyear_source_panel_20261107 as panel
from scripts.analysis.mois_2024_sameyear_global_island_photo_white_20261107 import (
    COVARS, cluster_fit)

TIERS={
 "strict_5p3r3s":(5,3,3),
 "measured_3p2r2s":(3,2,2),
 "measured_1p1r1s_support_only":(1,1,1),
}

def exact_year_block_support(df,region="region5"):
    """Require actual island and mainland same year AND geographic block."""
    a=df.groupby(["year",region],observed=True).is_island.agg(["min","max"])
    ids=set(a.index[(a["min"]==0)&(a["max"]==1)])
    idx=pd.MultiIndex.from_frame(df[["year",region]])
    return df.loc[idx.isin(ids)].copy(),ids

def verified_colour_frame(year, full, tier):
    """One real photo ROI per accepted species; never treat unknown as white/nonwhite."""
    attempts, successes, species=tier
    p=full["population_status"] if "population_status" in full else full["population_source_status"]
    density=(full["population_density_per_km2"] if "population_density_per_km2" in full
             else full["population_density"])
    eligible=(
        full.photo_attempts.ge(attempts)&
        full.strict9_success.ge(successes)&
        full.accepted_species_n.ge(species)&
        full.viirs_valid_months.ge(3)&full.viirs_median.notna()&
        p.eq("OK")&density.notna()&
        full.land_class.isin(("island","mainland"))&
        full.land_id.notna())
    x=full.loc[eligible].copy()
    x["year"]=year
    x["is_island"]=x.land_class.eq("island").astype(int)
    x["log_light"]=np.log1p(x.viirs_median.clip(lower=0))
    x["log_pop"]=np.log1p(density.loc[x.index].clip(lower=0))
    x["log_rain"]=np.log1p(x.bio12.clip(lower=0))
    x["log_attempts"]=np.log1p(x.photo_attempts.clip(lower=0))
    x["white_pp"]=100*x.white_equal
    x["region5"]=(np.floor((x.query_longitude+180)/5).astype(int).astype(str)+
                  "_"+np.floor((x.query_latitude+90)/5).astype(int).astype(str))
    x=x.dropna(subset=["white_pp","log_light","log_pop","log_rain","bio5",
                       "elevation","log_attempts","region5"])
    panel.require(x.target_id.is_unique and x.white_pp.between(0,100).all(),
                  "New apparent photo coverage contains duplicate/imputed source colour")
    panel.require(x.target_id.str.startswith(str(year)+"|").all(),
                  "Photo-year exposure switched to another year")
    return x

def build_frozen_five_years(paths):
    photo,wcvp=panel.load_frozen_photos(paths)
    frames={}
    for year in (2019,2021,2023):
        full,strict=panel.year_historic(year,paths,photo,wcvp,return_all=True)
        frames[year]=(full,strict)
    for year,fn in ((2024,panel.year2024),(2025,panel.year2025)):
        full,strict=fn(paths,photo,wcvp,return_all=True)
        frames[year]=(full,strict)
    return frames

def by_year_tier(frames):
    cohorts={}
    original_strict_ids={}
    for year,(full,strict) in frames.items():
        panel.require(len(strict)==panel.YEAR_COUNTS[year],
                      f"Original frozen strict {year} cell count was altered")
        original_strict_ids[year]=set(strict.target_id)
        for key,tier in TIERS.items():
            d=verified_colour_frame(year,full,tier)
            if key=="strict_5p3r3s":
                panel.require(set(d.target_id)==original_strict_ids[year],
                    f"New tier changes existing strict/source eligibility in {year}")
            cohorts.setdefault(key,[]).append(d)
    return {key:pd.concat(dfs,ignore_index=True) for key,dfs in cohorts.items()}

def support_audit(tier,df):
    sub, groups=exact_year_block_support(df)
    row={
        "tier":tier,
        "source_photo_cells_years":len(df),
        "island_cell_years":int(df.is_island.sum()),
        "mainland_cell_years":int((1-df.is_island).sum()),
        "source_island_landmasses":int(df.loc[df.is_island.eq(1),"land_id"].nunique()),
        "source_geoblocks":int(df.region5.nunique()),
        "both_land_classes_sameyear_block_groups":len(groups),
        "sameyear_distinct_supported_5deg_blocks":int(sub.region5.nunique()),
        "sameyear_shared_photo_cell_years":len(sub),
        "sameyear_shared_island_cell_years":int(sub.is_island.sum()),
        "true_photo_strict9_and_wcvp_only":True,
    }
    return row,sub

def fit_tier(df,shared):
    """Prespecify fixed model forms, no p-hacking over photo quality tier."""
    f="white_pp ~ log_light * is_island + log_pop + "+COVARS
    # Full original 5-year "different light slope on islands" model.
    full=cluster_fit(df,f+" + C(region5) + C(year)","region5")
    # Matched same YEAR+region block: enforce year:block fixed effects,
    # not merely separate year FE and region FE. Cluster across all years
    # by same region5, to prevent pseudo-replicated geographic groups.
    shared=shared.copy()
    shared["year_region5"]=(shared.year.astype(str)+"__"+shared.region5)
    groups=int(shared.year_region5.nunique())
    if shared.region5.nunique()<10 or groups<10 or shared.is_island.sum()<40:
        matched={"estimable":False,
                 "reason":"fewer than 10 matched year-geographical units or 40 island cell-years",
                 "region5_blocks":int(shared.region5.nunique()),
                 "year_region5_groups":groups,
                 "island_cell_years":int(shared.is_island.sum())}
    else:
        result=cluster_fit(shared,f+" + C(year_region5)","region5")
        matched={"estimable":True,"region5_blocks":int(shared.region5.nunique()),
                 "year_region5_groups":groups,
                 "island_cell_years":int(shared.is_island.sum()),
                 "interaction":result["log_light:is_island"],
                 "geographic_cluster_repeated_across_years":True,
                 "year_by_region_fixed_effects":True}
    return {
        "global_region_year_adjusted_interaction":full["log_light:is_island"],
        "within_sameyear_5deg_region_matched_interaction":matched,
        "tier_is_exploratory_quality_expansion_not_primary_confirmatory":True,
    }

def provenance_of_expanded_co_support(strict,moderate):
    """Where do additional true-ROI pairs arise? No observation-level export.

    Keep old strict paired year×region units fixed, then ask whether adding
    measured but lower-opportunity photos changes the interaction WITHIN
    the same paired units. Separately report the newly paired year×regions.
    The two are exploratory support diagnostics, not a multiple-testing
    search for a favourable significance threshold.
    """
    old_shared,oldgroups=exact_year_block_support(strict)
    new_shared,newgroups=exact_year_block_support(moderate)
    panel.require(oldgroups.issubset(newgroups),"New measured data lost strict original paired groups")
    coords=pd.MultiIndex.from_frame(moderate[["year","region5"]])
    old_area=moderate.loc[coords.isin(oldgroups)].copy()
    new_area=moderate.loc[coords.isin(newgroups-oldgroups)].copy()
    ids=set(strict.target_id)
    newly=moderate.loc[~moderate.target_id.isin(ids)]
    report={
      "old_strict_paired_year_region_groups":len(oldgroups),
      "new_measured_paired_year_region_groups":len(newgroups),
      "newly_paired_year_region_groups":len(newgroups-oldgroups),
      "new_original_strict9_measured_added_cell_years":int(len(newly)),
      "new_original_strict9_measured_added_island_cells":int(newly.is_island.sum()),
      "expanded_measured_cells_in_previously_paired_groups":int(len(old_area)),
      "newly_added_cells_in_previously_paired_groups":int(
          (~old_area.target_id.isin(ids)).sum()),
      "island_cells_in_previously_paired_groups":int(old_area.is_island.sum()),
      "expanded_measured_cells_in_newly_paired_groups":int(len(new_area)),
      "island_cells_in_newly_paired_groups":int(new_area.is_island.sum()),
    }
    f="white_pp ~ log_light * is_island + log_pop + "+COVARS+" + C(year_region5)"
    def restricted(frame, name):
        regions=int(frame.region5.nunique())
        n_island=int(frame.is_island.sum())
        if regions<10 or n_island<40 or len(frame)<100:
            return {"estimable":False,"reason":"insufficient geographic/year direct-pair support",
                    "clusters":regions,"island_cells":n_island,"cells":len(frame)}
        frame=frame.copy()
        frame["year_region5"]=frame.year.astype(str)+"__"+frame.region5
        return {"estimable":True,"clusters":regions,"island_cells":n_island,
                "cells":len(frame),
                "interaction":cluster_fit(frame,f,"region5")["log_light:is_island"]}
    report["measured_2species_within_strict_existing_paired_blocks"]=restricted(old_area,"existing")
    report["measured_2species_only_new_paired_blocks"]=restricted(new_area,"new")
    return report

def run(paths,out):
    out=Path(out)
    panel.require(not out.exists(),"Refuse to overwrite source-verified photo cohort")
    frames=build_frozen_five_years(paths)
    tiers=by_year_tier(frames)
    keys=list(TIERS)
    for left,right in zip(keys,keys[1:]):
        panel.require(set(tiers[left].target_id).issubset(set(tiers[right].target_id)),
                      f"Measured tier is not superset of strict original: {left}")
    summary=[]
    fits={}
    for name,df in tiers.items():
        row,matched=support_audit(name,df)
        summary.append(row)
        if name in ("strict_5p3r3s","measured_3p2r2s"):
            fits[name]=fit_tier(df,matched)
    provenance=provenance_of_expanded_co_support(
        tiers["strict_5p3r3s"],tiers["measured_3p2r2s"])
    baseline=fits["strict_5p3r3s"]["global_region_year_adjusted_interaction"]
    panel.require(abs(baseline["beta"]-(-2.4970445996320447))<0.00001 and
                  len(tiers["strict_5p3r3s"])==7982,
                  "Previously signed strict-cohort estimate must remain unchanged")
    frame=pd.DataFrame(summary)
    panel.require(set(frame.tier)==set(TIERS),"Missing quality tier")
    out.mkdir(parents=True)
    csv_path=out/"source_measured_flower_photo_quality_support_counts.csv"
    frame.to_csv(csv_path,index=False)
    receipt={
       "schema":"MOIS_island_mainland_actual_measured_photo_opportunity_tier_v1",
       "year_specific_original_photo_source":[2019,2021,2023,2024,2025],
       "no_2022_provider_400_month_replacement":True,
       "no_2020_global1_population_splicing":True,
       "white_outcome":"100 * species-equal original strict9 ROI white component",
       "source_population":"R2024B unconstrained sameyear 1km pixel density",
       "source_VIIRS":"same photo year at least three authentic monthly values",
       "island_area_cutoff":None,
       "unknown_colour_imputed_as_white_or_nonwhite":0,
       "new_original_image_downloads":0,
       "new_original_ROI_predictions":0,
       "tier_gates":{"photo_attempts_strict9_accepted_species":{k:list(v) for k,v in TIERS.items()}},
       "tier_support":summary,
       "tier_models":fits,
       "source_measured_co_support_provenance":provenance,
       "strict_primary_preserved":True,
       "measured_1_species_is_support_counts_only":True,
       "observational_not_causal_nor_biologically_white_species_fraction":True,
       "5degree_co_support_requires_simultaneous_island_mainland_same_year":True,
       "source_input_SHA256":{k:panel.sha(v) for k,v in paths.items()},
       "output_support_counts_SHA256":panel.sha(csv_path),
    }
    (out/"source_receipt.json").write_text(json.dumps(receipt,sort_keys=True,indent=2)+"\n")
    print("MOIS_REAL_MEASURED_OPPORTUNITY_TIERS_ISLAND_ALAN",
        json.dumps({"counts":summary,"models":fits,"support_provenance":provenance},
                   sort_keys=True),flush=True)
    return receipt

if __name__=="__main__":
    p=argparse.ArgumentParser()
    source_keys=("canonical","wcvp",
       "geo2019","viirs2019","pop2019",
       "geo2021","viirs2021","pop2021",
       "geo2023","viirs2023","pop2023",
       "geo2024","cn2024","uc2024","uc2024receipt",
       "geo2025","geo2025receipt","viirs2025")
    for key in source_keys:
        p.add_argument("--"+key.replace("_","-"),dest=key,type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    run({k:getattr(a,k) for k in source_keys},a.out)
