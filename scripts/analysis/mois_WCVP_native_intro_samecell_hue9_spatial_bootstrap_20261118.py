#!/usr/bin/env python3
"""MOIS same-cell WCVP-native vs introduced photo9 assembly anomaly.

The original global WCVP spatial source was already authenticated against
37,104 accepted species, 5,493,931 species×GBIF100km cells, the original
Kew TDWG3 native/introduced source, GSHHG 2.3.7 coastline/WorldClim elevation
and the *same 199 genus-by-occupancy-breadth globally shuffled colour labels*.

Here we only perform a *spatial cluster bootstrap* of signed 100km cell
contrasts between regional native and introduced species assemblages, keeping
both groups in the EXACT same cell (when each has >=30 recorded species).
Each group's Shannon anomaly is observed minus genus×range-width expected.
WCVP "introduced" is regional source status, not proof of recent wild
invasion, natural selection on flower petals, or post-introduction colour
change. Latitude/altitude/island contrasts compare distinct taxa and localities
and are not causal tests.
"""
from __future__ import annotations
import argparse,hashlib,io,json
from pathlib import Path
from zipfile import ZipFile
import numpy as np
import pandas as pd

ORIGINAL_SHA="cbe1eef381e3e1a0062dbf5134b047f947338b7c499d2ecc0cd22ee8855909fd"
CELL_FILE="global_WCVP_TDWG3_native_only_photo9_island_alpine_polar_south_cells.csv.gz"
N_BOOT=20000
SEED=20261118
ARCTIC=66.562

def check(ok,msg):
    if not bool(ok):raise ValueError("MOIS_SAMECELL_NATIVE_INTRO_HUE9: "+msg)

def sha(v):
    if isinstance(v,(Path,str)):v=Path(v).read_bytes()
    return hashlib.sha256(v).hexdigest()

def load(path):
    check(sha(path)==ORIGINAL_SHA,"Original native status nine-hue ZIP signature changed")
    with ZipFile(path) as z:
        receipt=json.loads(z.read("receipt.json"))
        for name,digest in receipt["output_sha256"].items():
            check(sha(z.read(name))==digest,
                  "Original WCVP nine hue geographic file internal SHA changed")
        d=pd.read_csv(io.BytesIO(z.read(CELL_FILE)),compression="gzip")
    check(receipt["source_full_GBIF_species_cell_pairs"]==5493931 and
          receipt["source_37104_accepted_photo_hue_taxa"]==37104 and
          receipt["source_null_199_same_global_genus_by_width_colour_label_permutations"],
          "Original verified accepted colour species scope or permutation null changed")
    check(len(d)==18577 and not d[["grid_x","grid_y"]].duplicated().any(),
          "Original source 18577 GBIF cells changed")
    return d

def paired(d,mode):
    s=d.copy()
    n=f"{mode}_native_"
    i=f"{mode}_introduced_"
    require_cols=[n+"n_photo9_species",i+"n_photo9_species",
                  n+"white_species_fraction",i+"white_species_fraction",
                  n+"Shannon_observed",i+"Shannon_observed",
                  n+"Shannon_null",i+"Shannon_null",
                  n+"Shannon_anomaly",i+"Shannon_anomaly"]
    check(set(require_cols).issubset(s),"WCVP native/intro required columns absent")
    keep=s[n+"n_photo9_species"].ge(30)&s[i+"n_photo9_species"].ge(30)
    sub=s.loc[keep].copy().reset_index(drop=True)
    check(len(sub)==({"strict":3106,"moderate":3602}[mode]),
          "Original WCVP matched exact same-cell source cohort drifted")
    sub["status_support_min"]=sub[[n+"n_photo9_species",
                                   i+"n_photo9_species"]].min(axis=1)
    sub["native_minus_introduced_anomaly"]=(
        sub[n+"Shannon_anomaly"]-sub[i+"Shannon_anomaly"])
    sub["native_minus_introduced_9hue_observed"]=(
        sub[n+"Shannon_observed"]-sub[i+"Shannon_observed"])
    sub["native_minus_introduced_genus_null"]=(
        sub[n+"Shannon_null"]-sub[i+"Shannon_null"])
    sub["native_minus_introduced_white_fraction"]=(
        sub[n+"white_species_fraction"]-sub[i+"white_species_fraction"])
    check(np.allclose(sub.native_minus_introduced_anomaly,
          sub.native_minus_introduced_9hue_observed-
          sub.native_minus_introduced_genus_null,atol=1e-10,rtol=0),
          "WCVP same-cell Shannon native-intro anomaly identity broken")
    sub["block_500km"]=(sub.grid_x.floordiv(5).astype(str)+"|"+
                         sub.grid_y.floordiv(5).astype(str))
    return sub

def masks(d):
    elev=d.sampled_land_median_elevation_m
    lat=d.center_latitude
    return {
      "all_same_100km_cells":np.ones(len(d),dtype=bool),
      "GSHHG_islands_majority7":d.gshhg_land_class.eq("island").to_numpy(),
      "GSHHG_mainlands_majority7":d.gshhg_land_class.eq("mainland").to_numpy(),
      "high_elevation1500plus":elev.ge(1500).to_numpy(),
      "low_elevation_below500":(elev.ge(-430)&elev.lt(500)).to_numpy(),
      "Arctic_circle_N66p562plus":lat.ge(ARCTIC).to_numpy(),
      "north_boreal_N50to66p562":(lat.ge(50)&lat.lt(ARCTIC)).to_numpy(),
      "southern_hemisphere":lat.lt(0).to_numpy(),
      "northern_hemisphere":lat.ge(0).to_numpy(),
      "southern_midlatitude_S30to55":(lat.le(-30)&lat.gt(-55)).to_numpy(),
      "northern_midlatitude_N30to55":(lat.ge(30)&lat.lt(55)).to_numpy(),
    }

def all_stat_500km_bootstrap(d,weight):
    grouped_masks=masks(d)
    k=[name for name,arr in grouped_masks.items() if arr.sum()>=30]
    check("all_same_100km_cells" in k,"Original same-cell source absent")
    # Each group/metric has signed sum and positive denominator. Every
    # bootstrap draw resamples the SAME 500km block IDs across all groups,
    # so group contrasts preserve geographic clustering and covariance.
    variables=(
      "native_minus_introduced_anomaly",
      "native_minus_introduced_white_fraction",
      "native_minus_introduced_9hue_observed",
      "native_minus_introduced_genus_null",
    )
    blocks,bidx=np.unique(d.block_500km.to_numpy(str),return_inverse=True)
    nb=len(blocks)
    check(nb>=100,"Too few spatially distinct WCVP same cells for block bootstrap")
    w=(np.ones(len(d)) if weight=="equal_100km_cells"
       else d.status_support_min.to_numpy(float))
    cells=[]
    for gname in k:
        selector=np.asarray(grouped_masks[gname],bool)
        support=np.bincount(bidx,weights=w*selector,minlength=nb)
        block_vals=[support]
        for v in variables:
            block_vals.append(np.bincount(
                bidx,weights=d[v].to_numpy(float)*w*selector,minlength=nb))
        cells.append(np.column_stack(block_vals))
    stack=np.concatenate(cells,axis=1)
    nvariables=len(variables)
    # Shape: boot chunks × groups × (weight, weighted responses)
    rng=np.random.default_rng(SEED+(0 if weight=="equal_100km_cells" else 9000)
                              +len(d))
    dist={name:{v:np.empty(N_BOOT) for v in variables} for name in k}
    width=nvariables+1
    for pos in range(0,N_BOOT,200):
        take=min(200,N_BOOT-pos)
        draw=rng.integers(0,nb,size=(take,nb))
        sums=stack[draw].sum(axis=1)
        for gi,name in enumerate(k):
            kslice=sums[:,gi*width:(gi+1)*width]
            den=kslice[:,0]
            check((den>0).all(),"Bootstrap occasionally dropped entire geo stratum: "+name)
            for j,v in enumerate(variables):
                dist[name][v][pos:pos+take]=kslice[:,j+1]/den
    results={}
    for name in k:
        m=np.asarray(grouped_masks[name],bool)
        a=d.loc[m]
        ww=w[m]
        result={
            "n_samecell_100km":int(m.sum()),
            "n_500km_blocks":int(a.block_500km.nunique()),
            "n_bootstrap":N_BOOT,
        }
        for v in variables:
            result[v]={
                "mean":float(np.average(a[v],weights=ww)),
                "CI95":[float(z) for z in np.percentile(dist[name][v],[2.5,97.5])]
            }
        results[name]=result
    contrasts={}
    pairs={
        "island_minus_mainland":("GSHHG_islands_majority7","GSHHG_mainlands_majority7"),
        "alpine_minus_lowland":("high_elevation1500plus","low_elevation_below500"),
        "Arctic_minus_boreal":("Arctic_circle_N66p562plus","north_boreal_N50to66p562"),
        "southern_minus_northern_hemisphere":("southern_hemisphere",
                                               "northern_hemisphere"),
        "southern_minus_northern_midlatitudes":("southern_midlatitude_S30to55",
                                                 "northern_midlatitude_N30to55"),
    }
    for name,(A,B) in pairs.items():
        if A not in k or B not in k:continue
        contrasts[name]={"n_A":results[A]["n_samecell_100km"],
                         "n_B":results[B]["n_samecell_100km"]}
        for v in variables:
            boot=dist[A][v]-dist[B][v]
            contrasts[name][v]={
                "difference":float(results[A][v]["mean"]-results[B][v]["mean"]),
                "CI95":[float(z) for z in np.percentile(boot,[2.5,97.5])]
            }
    return results,contrasts

def run(path,out):
    root=Path(out)
    check(not root.exists(),"Refuse overwriting same-cell spatial bootstrap result")
    src=load(path)
    output={}
    combined=[]
    for mode in ("strict","moderate"):
        same=paired(src,mode)
        cols=["grid_x","grid_y","center_latitude","center_longitude",
              "gshhg_land_class","sampled_land_median_elevation_m",
              "block_500km","status_support_min",
              "native_minus_introduced_anomaly",
              "native_minus_introduced_white_fraction",
              "native_minus_introduced_9hue_observed",
              "native_minus_introduced_genus_null"]
        snap=same[cols].copy()
        snap["region_assignment_mode"]=mode
        combined.append(snap)
        output[mode]={}
        for weighting in ("equal_100km_cells","min_of_two_status_species"):
            res,diff=all_stat_500km_bootstrap(same,weighting)
            output[mode][weighting]={
                "geographic_groups":res,
                "group_differences":diff
            }
            print("MOIS_SAMECELL_WCVP_NATIVE_INTRO_HUE9_BOOT_SUCCESS",json.dumps({
                "mode":mode,"weighting":weighting,
                "n_cells":len(same),
                "mean_anomaly_diff":
                   res["all_same_100km_cells"][
                       "native_minus_introduced_anomaly"],
                "geo_differences":{
                  k:v["native_minus_introduced_anomaly"]
                  for k,v in diff.items()}
            },allow_nan=False),flush=True)
    root.mkdir(parents=True)
    f=root/"WCVP_native_intro_same_100km_cell_Shannon_anomaly_geography.csv.gz"
    pd.concat(combined,ignore_index=True).to_csv(f,index=False,compression="gzip")
    summary=root/"WCVP_native_intro_paired_500km_cluster_bootstrap.json"
    report={
      "signed_original_WCVP_source_SHA256":ORIGINAL_SHA,
      "n_accepted_species_universe":37104,
      "n_true_source_species_grid_pairs":5493931,
      "WCVP_native_introduced_only_broad_TDWG3_range_status":True,
      "same_100km_grid_comparison_with_both_statuses_min30_species":True,
      "n_500km_cluster_bootstrap":N_BOOT,
      "strict_moderate_same_underlying_species_not_independent":True,
      "does_not_identify_introduction_fitness_selection_or_flowers":True,
      "contrast_is_permuted_photo9_entropy_anomaly_not_true_flower_colour":True,
      "results":output,
    }
    summary.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    report["files_sha256"]={p.name:sha(p) for p in (f,summary)}
    (root/"receipt.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--source",required=True)
    p.add_argument("--out",required=True)
    a=p.parse_args()
    run(a.source,a.out)
