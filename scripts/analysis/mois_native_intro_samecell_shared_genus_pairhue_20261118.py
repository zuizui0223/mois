#!/usr/bin/env python3
"""MOIS genuine ecological test: within the SAME GBIF 100km cell and SAME
genus, compare taxonomic white/9-hue assortativity in regional-native vs
regional-introduced species recorded by accepted WCVP TDWG3.

Ecological predeclared estimand: for genus g and species pool status s with
n>=2 distinct accepted species in a grid, unbiased probability that TWO
DISTINCT sampled species have the same photo9 hue label:
  P_same(g,s) = sum_k n_k*(n_k-1) / [n*(n-1)].
It avoids mechanical small-sample downward bias of plug-in Shannon.
Within each cell+genus, compare introduced minus native P_same, without
ever changing geographic cell or genus (and without calling same taxa
independent across cells).

Core stratum: WCVP STRICT assigned grid with >=30 native AND >=30 introduced
species in total; >=3 shared genera with >=2 species in each status.
Sensitivity: >=3/5 species per shared genus, >=3/5 shared genera, moderate
geographic status, minimum 100 native+introduced total per same grid,
min(species) weighting. 20k 500km block bootstrap; within-block highland
minus lowland 500km geographical comparisons when 12+ shared blocks.
Photo hues are species-level modal iNat ROI, not independently verified petal
colour; WCVP status is regional TDWG3, not demonstrated local wild occurrence.
"""
from __future__ import annotations
import argparse,hashlib,io,json
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pandas as pd

from scripts.analysis.mois_global_hue9_community_assembly_genus_null_20261118 import (
    signed_original, GBIF_SHA,ENV_SHA,check,sha
)
from scripts.analysis.mois_wcvp_native_hue9_geographical_null_20261118 import (
    get_wcvp_source, match_gbif_to_wcvp, load_palaeoclimate_grid,
    assign_cells_broad_region,regional_species_cell_rows,load_geog,
    LGM_SHA,TDWG_GEO_SHA,GEO_SHA
)

HUES=("white","yellow","orange","red","pink","magenta","purple","blue","bronze")
MIN_NATIVE_INTRO_TOTAL=30
N_BOOT=20000
SEED=20261118

def signed_wcvp_source(gbif,env,lgm,wcvp_dir,tdwg,geo_zip):
    pairs,spp=signed_original(gbif,env)
    geo=load_geog(geo_zip)
    check(len(geo)==18577 and
          not geo[["grid_x","grid_y"]].duplicated().any(),
          "Original 100km geographical reference changed")
    lm=load_palaeoclimate_grid(lgm)
    regions=assign_cells_broad_region(lm,tdwg)
    check(len(regions)==18577 and int(regions.region_strict.sum())==9288
          and int(regions.region_moderate.sum())==10391,
          "Original strict/moderate WCVP geographical support changed")
    nm,di,receipt=get_wcvp_source(wcvp_dir)
    resolved=match_gbif_to_wcvp(spp,nm)
    check(int(resolved.wcvp_match_exact_binomial.sum())==36762,
          "Original exact WCVP accepted binomial source coverage changed")
    ids=regions.loc[regions.region_moderate,["grid_x","grid_y"]]
    source=regional_species_cell_rows(pairs,spp,nm,di,regions,ids)
    check(len(source)==3127248 and
          not source[["specieskey","grid_x","grid_y"]].duplicated().any(),
          "Accepted species WCVP same-GBIF-cell source count changed")
    src=source.merge(spp[["key","source_hue9"]],left_on="specieskey",
                     right_on="key",how="left",validate="many_to_one")
    check(src.source_hue9.isin(HUES).all() and
          src.white.eq(src.source_hue9.eq("white")).all(),
          "Authentic nine-hue species source mismatch at WCVP regional level")
    regional_source_counts={}
    for mode in ("strict","moderate"):
        subset=src.loc[src["region_"+mode]]
        for status in ("native","introduced"):
            n=int(subset.status.eq(status+"_at_TDWG3").sum())
            regional_source_counts[mode+"_"+status]=n
    check(regional_source_counts=={
        "strict_native":2271852,"strict_introduced":436422,
        "moderate_native":2604870,"moderate_introduced":522378
    },"Source WCVP native and introduced class counts changed")
    return src,geo,regional_source_counts

def genus_pairs(src,mode):
    selected=src.loc[src["region_"+mode]].copy()
    # Original Gini-Simpson across *accepted species*, never weighted by GBIF
    # number of photograph/occurrence records or repeated species per cell.
    color=selected.groupby(["grid_x","grid_y","genus","status",
                            "source_hue9"],sort=False,observed=True).size()
    col=color.rename("n_hue").reset_index()
    col["ordered_pairs_same_hue"]=col.n_hue*(col.n_hue-1)
    col["n_white"]=np.where(col.source_hue9.eq("white"),col.n_hue,0)
    totals=col.groupby(["grid_x","grid_y","genus","status"],sort=False,observed=True).agg(
        n_species=("n_hue","sum"),
        ordered_pairs_same_hue=("ordered_pairs_same_hue","sum"),
        white_species=("n_white","sum")
    ).reset_index()
    a=totals.loc[totals.status.eq("native_at_TDWG3")].drop(columns="status")
    b=totals.loc[totals.status.eq("introduced_at_TDWG3")].drop(columns="status")
    paired=a.merge(b,on=["grid_x","grid_y","genus"],how="inner",
                   validate="one_to_one",suffixes=("_native","_introduced"))
    check(len(paired)>10000,"Too few genuine shared genera at original same cells")
    paired["n_min"]=paired[["n_species_native",
                           "n_species_introduced"]].min(axis=1)
    paired["pair_prob_same_native"]=np.divide(
        paired.ordered_pairs_same_hue_native,
        paired.n_species_native*(paired.n_species_native-1),
        out=np.full(len(paired),np.nan),
        where=paired.n_species_native.ge(2))
    paired["pair_prob_same_introduced"]=np.divide(
        paired.ordered_pairs_same_hue_introduced,
        paired.n_species_introduced*(paired.n_species_introduced-1),
        out=np.full(len(paired),np.nan),
        where=paired.n_species_introduced.ge(2))
    paired["introduced_minus_native_same_photo_hue_prob"]=(
        paired.pair_prob_same_introduced-paired.pair_prob_same_native)
    paired["intro_minus_native_white_genus_fraction"]=(
        paired.white_species_introduced/paired.n_species_introduced-
        paired.white_species_native/paired.n_species_native)
    check(paired.loc[paired.n_min.ge(2),
          "introduced_minus_native_same_photo_hue_prob"].between(-1,1).all(),
          "Unbiased pairwise hue same-class probability invalid")
    check(paired.n_min.ge(1).all(),"Genus occurrence supports absent")
    return paired,totals

def group_cell_support(totals,geo,mode):
    whole=totals.groupby(["grid_x","grid_y","status"],sort=False,observed=True).n_species.sum().unstack("status",fill_value=0)
    check("native_at_TDWG3" in whole and "introduced_at_TDWG3" in whole,
          "WCVP native introduced total species source unavailable")
    whole=whole.rename(columns={
        "native_at_TDWG3":"total_native",
        "introduced_at_TDWG3":"total_introduced"
    }).reset_index()
    d=geo.merge(whole,on=["grid_x","grid_y"],how="left",validate="one_to_one")
    d[["total_native","total_introduced"]]=d[
        ["total_native","total_introduced"]].fillna(0).astype(int)
    n=(d.total_native.ge(30)&d.total_introduced.ge(30)).sum()
    check(int(n)=={"strict":3106,"moderate":3602}[mode],
          "Original exact WCVP matched 100km grid count changed")
    d["original_cohort_both_status30"]=(d.total_native.ge(30)&
                                        d.total_introduced.ge(30))
    d["original_cohort_both_status100"]=(d.total_native.ge(100)&
                                         d.total_introduced.ge(100))
    d["geographic_500km_block"]=(d.grid_x.floordiv(5).astype(str)+"|"+
                                   d.grid_y.floordiv(5).astype(str))
    return d

def build_cell_genus_comparison(paired,geo,mode,per_genus_min=2,
                                min_shared_genera=3,status_min_species=30):
    check(per_genus_min in (2,3,5) and
          min_shared_genera in (3,5) and
          status_min_species in (30,100),"Non-prespecified genus support threshold")
    eligible=paired.loc[paired.n_min.ge(per_genus_min)].copy()
    eligible=eligible.dropna(subset=["introduced_minus_native_same_photo_hue_prob"])
    eligible["delta_equal"]=eligible.introduced_minus_native_same_photo_hue_prob
    eligible["delta_support_num"]=(eligible.delta_equal*eligible.n_min)
    eligible["delta_white"]=eligible.intro_minus_native_white_genus_fraction
    group=eligible.groupby(["grid_x","grid_y"],sort=False,observed=True).agg(
       n_shared_genera=("genus","size"),
       native_species_in_shared_genera=("n_species_native","sum"),
       intro_species_in_shared_genera=("n_species_introduced","sum"),
       total_min_genus_n=("n_min","sum"),
       avg_genus_probability=("delta_equal","mean"),
       support_weighted_num=("delta_support_num","sum"),
       avg_genus_white_diff=("delta_white","mean")
    ).reset_index()
    group["weighted_genus_probability"]=(
        group.support_weighted_num/group.total_min_genus_n)
    joined=geo.merge(group,on=["grid_x","grid_y"],how="inner",
                     validate="one_to_one")
    keep=joined.n_shared_genera.ge(min_shared_genera)
    keep &= (joined.total_native.ge(status_min_species)&
             joined.total_introduced.ge(status_min_species))
    d=joined.loc[keep].copy().reset_index(drop=True)
    d["min_whole_status_species"]=d[["total_native","total_introduced"]].min(axis=1)
    check(not d[["grid_x","grid_y"]].duplicated().any(),
          "Paired genus hue geometry is not one row per 100km GBIF cell")
    if per_genus_min==2 and min_shared_genera==3 and status_min_species==30:
        check(len(d)>=500 and d.geographic_500km_block.nunique()>100,
              "Insufficient source within-genus native introduced paired floras")
    d["highland"]=d.sampled_land_median_elevation_m.ge(1500)
    d["lowland"]=(d.sampled_land_median_elevation_m.ge(-430)&
                   d.sampled_land_median_elevation_m.lt(500))
    d["north"]=d.center_latitude.ge(0)
    d["south"]=d.center_latitude.lt(0)
    d["island"]=d.gshhg_land_class.eq("island")
    return d

def block_bootstrap(d,var="avg_genus_probability",weighting="equal",
                    reps=N_BOOT,seed=SEED):
    if len(d)<25:
        return {"status":"few_cells","n_cells":int(len(d))}
    key="geographic_500km_block"
    nblock=d[key].nunique()
    if nblock<8:return {"status":"few_blocks","n_cells":len(d),"n_blocks":nblock}
    w=(np.ones(len(d)) if weighting=="equal"
       else d.min_whole_status_species.to_numpy(float))
    temp=pd.DataFrame({"b":d[key].to_numpy(str),"w":w,
                       "wy":w*d[var].to_numpy(float)})
    g=temp.groupby("b",sort=False)[["w","wy"]].sum().to_numpy(float)
    rng=np.random.default_rng(seed+len(d))
    out=np.empty(reps)
    for start in range(0,reps,250):
        n=min(250,reps-start)
        ids=rng.integers(0,len(g),size=(n,len(g)))
        ag=g[ids].sum(axis=1)
        out[start:start+n]=ag[:,1]/ag[:,0]
    val=np.average(d[var].to_numpy(float),weights=w)
    return {
        "status":"ok","n_cells":len(d),"n_500km_blocks":nblock,
        "estimate":float(val),"CI95":[float(x) for x in np.percentile(out,[2.5,97.5])],
        "n_boot":reps
    }

def within_common_block_hilo(d,variable,weighting,seed):
    key="geographic_500km_block"
    a=d.loc[d.highland];b=d.loc[d.lowland]
    blocks=pd.Index(a[key].unique()).intersection(b[key].unique())
    if len(blocks)<12:
        return {"status":"few_shared_high_low_500km_blocks",
                "n_common_blocks":len(blocks)}
    def mean_by_block(df):
        df=df[df[key].isin(blocks)]
        w=(np.ones(len(df)) if weighting=="equal"
           else df.min_whole_status_species.to_numpy(float))
        t=pd.DataFrame({"block":df[key],"w":w,
                        "wy":w*df[variable].to_numpy(float)})
        sums=t.groupby("block")[["w","wy"]].sum().reindex(blocks)
        return (sums.wy/sums.w).to_numpy(float)
    delta=mean_by_block(a)-mean_by_block(b)
    rng=np.random.default_rng(seed+len(blocks)+4000)
    samples=np.empty(N_BOOT)
    for i in range(0,N_BOOT,250):
        take=min(250,N_BOOT-i)
        ids=rng.integers(0,len(blocks),size=(take,len(blocks)))
        samples[i:i+take]=delta[ids].mean(axis=1)
    return {
        "status":"ok","n_common_500km_blocks":len(blocks),
        "high_minus_low_estimate":float(np.mean(delta)),
        "CI95":[float(x) for x in np.percentile(samples,[2.5,97.5])],
        "n_boot":N_BOOT
    }

def regions(d):
    lat=d.center_latitude
    return {
      "global":np.ones(len(d),dtype=bool),
      "highland_1500mplus":d.highland.to_numpy(bool),
      "lowland_under500m":d.lowland.to_numpy(bool),
      "island_GSHHG_majority7":d.island.to_numpy(bool),
      "mainland_GSHHG_majority7":d.gshhg_land_class.eq("mainland").to_numpy(bool),
      "arctic_N66p562plus":lat.ge(66.562).to_numpy(bool),
      "north_hemisphere":d.north.to_numpy(bool),
      "south_hemisphere":d.south.to_numpy(bool),
      "north_midlat_30to55":(lat.ge(30)&lat.lt(55)).to_numpy(bool),
      "south_midlat_minus55to_minus30":(lat.le(-30)&lat.gt(-55)).to_numpy(bool)
    }

def model(mode,paired,totals,geo):
    geometry=group_cell_support(totals,geo,mode)
    candidate={}
    flat=[]
    for pergenus,mingenera,minstatus in [
        (2,3,30),(3,3,30),(5,3,30),(2,5,30),(2,3,100)
    ]:
        key=f"min{pergenus}species_each_same_genus_min{mingenera}genera_wholestatus{minstatus}"
        d=build_cell_genus_comparison(paired,geometry,mode,pergenus,
                                       mingenera,minstatus)
        reg=regions(d)
        estimates={}
        for weight in ("equal","min_whole_status_species"):
            selected={}
            for region,mask in reg.items():
                sub=d.loc[mask]
                for var in ("avg_genus_probability","weighted_genus_probability"):
                    res=block_bootstrap(sub,var,weight,
                         seed=SEED+(0 if weight=="equal" else 10000)+
                              (0 if var=="avg_genus_probability" else 20000))
                    selected[f"{region}__{var}"]=res
                    if res["status"]=="ok":
                        flat.append({
                            "WCVP_mode":mode,
                            "genus_filter":key,
                            "weighting":weight,"geo":region,
                            "measure":var,**res})
            hilow={}
            for variable in ("avg_genus_probability",
                             "weighted_genus_probability"):
                hilow[variable]=within_common_block_hilo(
                    d,variable,weight,SEED+
                    (1 if weight=="equal" else 9000))
            estimates[weight]={"subregions":selected,
                               "nearby500km_high_minus_low":hilow}
        candidate[key]={
            "n_cells":len(d),
            "n_actual_common_genera_per_cell_median":float(
                d.n_shared_genera.median()) if len(d)>0 else None,
            "n_shared_500km_blocks":int(
                d.geographic_500km_block.nunique()),
            "estimates":estimates
        }
        if key=="min2species_each_same_genus_min3genera_wholestatus30":
            print("MOIS_TRUE_COMMON_GENUS_SAMECELL_HUE9_PRIMARY",json.dumps({
               "mode":mode,"n_cells":len(d),
               "median_common_genera":candidate[key][
                    "n_actual_common_genera_per_cell_median"],
               "global_equal":estimates["equal"]["subregions"].get(
                    "global__avg_genus_probability"),
               "nearby_high_minus_low":estimates["equal"][
                    "nearby500km_high_minus_low"]
            },allow_nan=False),flush=True)
    return candidate,flat

def run(args):
    root=Path(args.out)
    check(not root.exists(),"Never replace frozen real ecological results")
    src,geo,original=signed_wcvp_source(
        args.GBIF,args.ENV,args.LGM,args.WCVP,args.TDWG,args.GEO)
    results={}
    rows=[]
    for mode in ("strict","moderate"):
        pairs,totals=genus_pairs(src,mode)
        result,table=model(mode,pairs,totals,geo)
        results[mode]={
            "n_genus_cell_status_pairs":len(pairs),
            "n_unique_same_cell_shared_genera":int(pairs[
                ["grid_x","grid_y","genus"]].drop_duplicates().shape[0]),
            "statistics":result
        }
        rows.extend(table)
    root.mkdir(parents=True)
    outcsv=root/"true_shared_genus_hue9_same_species_pair_probability_500km.csv"
    pd.DataFrame(rows).to_csv(outcsv,index=False)
    output={
       "study":"9class same-photo-colour probability of 2 distinct accepted plant species within each shared genus and exact 100km grid, intro-minus-native",
       "source_accepted_species":37104,
       "source_GBIF_species_cell_pairs":5493931,
       "source_WCVP_resolved_native_introduced_GBIF_pairs":3127248,
       "source_status_gbif_counts":original,
       "source_sha256":{
         "GBIF":GBIF_SHA,"ENV":ENV_SHA,"LGM":LGM_SHA,
         "GEO":GEO_SHA,"TDWG3":TDWG_GEO_SHA
       },
       "unbiased_same_photo9_pair_estimator":True,
       "species_photo_colour_not_vetted_true_petal":True,
       "WCVP_status_TDWG3_not_individual_native":True,
       "within_cell_within_genus_restricts_taxonomy_but_not_species_causality":True,
       "spatial_500km_bootstrap_keeps_cell_statuses_paired":True,
       "n_500km_bootstraps":N_BOOT,
       "strata_not_independent_populations":True,
       "results":results
    }
    report=root/"true_common_genus_hue9_ecological_assortment.json"
    report.write_text(json.dumps(output,indent=2,allow_nan=False)+"\n")
    output["numeric_output_sha256"]={p.name:sha(p) for p in (outcsv,report)}
    (root/"receipt.json").write_text(json.dumps(output,indent=2,allow_nan=False)+"\n")
    print("MOIS_TRUE_COMMON_GENUS_9HUE_SAMECELL_REAL_SUCCESS",json.dumps({
         "n_original_source_pairs":len(src),
         "strict_primary":results["strict"]["statistics"][
              "min2species_each_same_genus_min3genera_wholestatus30"],
         "moderate_primary":results["moderate"]["statistics"][
              "min2species_each_same_genus_min3genera_wholestatus30"]
    },allow_nan=False)[:16000],flush=True)

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    for flag in ("GBIF","ENV","LGM","WCVP","TDWG","GEO","out"):
        parser.add_argument("--"+flag,required=True)
    run(parser.parse_args())
