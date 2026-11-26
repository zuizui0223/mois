#!/usr/bin/env python3
"""MOIS ecological stress-test: did shared-genus photo9 pair assortativity
at high altitude or in the Southern Hemisphere simply arise from a few
repeated genera and/or a few geographic 500km blocks?

Read the source-signed 2026-11-18 primary WCVP within-100km-cell,
within-genus (2+ native and 2+ introduced species, 3+ genera) result.
Perform reproducible crossed independent Exponential(1) multiplier weights
on 500km blocks AND accepted genera, applied coherently across all cells
where that genus occurs. Estimate equal-cell mean probability of two
different species matching photo9 hues, intro minus native, with its exact
white-white and matching nonwhite components.

For comparison report "space only" and "genus only" multiplier intervals.
These are descriptive *multiplier sensitivity intervals*, NOT established
nominal frequentist confidence intervals under all geographic/phylogenetic
dependence. Leave one and top N geographically frequent genera out as
non-independent, deliberately harsh influence diagnostics, always report
lost cell support and recompute original estimate on the retained support.

Never present an unstable single hue/intro relation as evidence of floral
selection or photographed true local petal identity.
"""
from __future__ import annotations
import argparse,hashlib,io,json
from pathlib import Path
from zipfile import ZipFile
import numpy as np
import pandas as pd

PARENT_SHA="f352fd89958c89a0ccb45136b555b97a13ad6fab2374496ae11b1000f83896c0"
CELL="same_grid_shared_genus_white_chromatic_pair_components.csv.gz"
GENUS="strict_same_grid_per_genus_colour_pair_evidence.csv.gz"
OUTCOMES=("introduced_minus_native_same_photo_hue_prob",
          "delta_white_white","delta_same_nonwhite")
SEED=20261118
N_MULT=4000

def ensure(x,why):
    if not bool(x):raise ValueError("MOIS_GENUS_SPATIAL_INFLUENCE: "+str(why))
def sha(v):
    if isinstance(v,(Path,str)):v=Path(v).read_bytes()
    return hashlib.sha256(v).hexdigest()

def read_source(path):
    ensure(sha(path)==PARENT_SHA,"Signed two-dimensional original ecology source ZIP drifted")
    with ZipFile(path) as z:
        receipt=json.loads(z.read("receipt.json"))
        for file,digest in receipt["output_sha256"].items():
            ensure(sha(z.read(file))==digest,
                   "Original numeric per-genus/per-grid source SHA mismatch "+file)
        grid=pd.read_csv(io.BytesIO(z.read(CELL)),compression="gzip")
        genus=pd.read_csv(io.BytesIO(z.read(GENUS)),compression="gzip")
    strict=grid.loc[grid.WCVP_mode.eq("strict")].copy()
    ensure(len(strict)==1363 and
           int(grid.WCVP_mode.eq("moderate").sum())==1640,
           "Exact primary supported grid source denominator changed")
    ensure(len(genus)==14639 and genus.genus.nunique()==236,
           "Signed original 236-genera matched species pair rows changed")
    ensure(not genus[["grid_x","grid_y","genus"]].duplicated().any(),
           "Duplicated observed genus per original cell")
    ensure(np.allclose(genus[OUTCOMES[0]],
           genus[OUTCOMES[1]]+genus[OUTCOMES[2]],atol=1e-10,rtol=0),
           "Exact white/nonwhite pair count identity broken")
    cols=["grid_x","grid_y","center_latitude",
          "sampled_land_median_elevation_m","geographic_500km_block"]
    joined=genus.merge(strict[cols],on=["grid_x","grid_y"],
                       how="inner",validate="many_to_one")
    ensure(len(joined)==len(genus),"Signed strict cell/genus keys altered")
    return joined,strict

def regions(t):
    lat=t.center_latitude.to_numpy(float)
    alt=t.sampled_land_median_elevation_m.to_numpy(float)
    land=t.geographic_500km_block
    return {
        "all":np.ones(len(t),bool),
        "highland":alt>=1500,
        "lowland":(alt>=-430)&(alt<500),
        "southern_hemisphere":lat<0,
        "northern_hemisphere":lat>=0,
        "highland_north":(alt>=1500)&(lat>=0),
        "highland_south":(alt>=1500)&(lat<0),
        "lowland_north":(alt>=-430)&(alt<500)&(lat>=0),
        "lowland_south":(alt>=-430)&(alt<500)&(lat<0),
        "southern_midlatitude":(lat<=-30)&(lat>-55),
        "northern_midlatitude":(lat>=30)&(lat<55)
    }

def multiplier(g,spatial,seed=SEED,nmult=N_MULT):
    if len(spatial)<25 or len(g)<50:
        return {"status":"insufficient_observed_source_support",
                "n_cells":len(spatial)}
    cid=pd.MultiIndex.from_frame(g[["grid_x","grid_y"]])
    orig_idx=pd.MultiIndex.from_frame(spatial[["grid_x","grid_y"]])
    cell_idx=orig_idx.get_indexer(cid)
    ensure((cell_idx>=0).all(),"Genus-cell row outside frozen spatial sample")
    genus_index,genera=pd.factorize(g.genus,sort=True)
    block_index,blocks=pd.factorize(spatial.geographic_500km_block,sort=True)
    ncell=len(spatial)
    ngen=len(genera);nblock=len(blocks)
    if min(nblock,ngen)<8:
        return {"status":"insufficient_genus_or_spatial_clusters",
                "n_cells":ncell,"n_genera":ngen,"n_blocks":nblock}
    y=g[list(OUTCOMES)].to_numpy(float)
    cc=np.bincount(cell_idx,minlength=ncell)
    ensure((cc>=3).all(),"Frozen >=3 shared genus per cell support changed")
    baseline=np.column_stack([
        np.bincount(cell_idx,weights=y[:,j],minlength=ncell)/cc
        for j in range(len(OUTCOMES))
    ]).mean(axis=0)
    ensure(abs(baseline[0]-baseline[1]-baseline[2])<1e-10,
           "Same-colour baseline does not equal white and nonwhite")
    rng=np.random.default_rng(seed+ncell)
    records={}
    for flavor in ("spatial_only","genus_only","genus_x_spatial"):
        arr=np.empty((nmult,len(OUTCOMES)))
        for i in range(nmult):
            wg=(np.ones(ngen,float) if flavor=="spatial_only"
                 else rng.exponential(1,size=ngen))
            wb=(np.ones(nblock,float) if flavor=="genus_only"
                 else rng.exponential(1,size=nblock))
            row_weight=wg[genus_index]
            den=np.bincount(cell_idx,weights=row_weight,minlength=ncell)
            ensure(den.min()>0,
                   "Zero species-photo hue support under strictly positive weights")
            block_weight=wb[block_index]
            for j in range(len(OUTCOMES)):
                weighted_sum=np.bincount(cell_idx,weights=row_weight*y[:,j],
                                         minlength=ncell)
                cell_value=weighted_sum/den
                arr[i,j]=float(np.average(cell_value,weights=block_weight))
        ensure(np.max(abs(arr[:,0]-arr[:,1]-arr[:,2]))<1e-10,
               "Two-way genus × spatial photo nine-hue arithmetic split failed")
        q=np.percentile(arr,[2.5,97.5],axis=0)
        records[flavor]={
           "N_multiplier_replicates":nmult,
           "variation":"independent positive Exponential1 weights on "+
                     ("blocks" if flavor=="spatial_only" else
                      "genera" if flavor=="genus_only" else "genera and blocks"),
           "effects":{name:{
                "original_equal_cell_estimate":float(baseline[j]),
                "sensitivity_95pct_interval":[float(q[0,j]),float(q[1,j])],
                "multiplier_replicate_mean":float(arr[:,j].mean())}
                for j,name in enumerate(OUTCOMES)}
        }
    return {
        "status":"ok","n_cells":ncell,"n_genera":ngen,
        "n_distinct_500km_blocks":nblock,
        "original_equal_cell_estimates":{
           name:float(baseline[j]) for j,name in enumerate(OUTCOMES)},
        "multiplier_methods":records
    }

def leave_top_n(g,grid):
    counts=g.groupby("genus").size().sort_values(ascending=False)
    group_names=counts.index.to_list()
    out=[]
    for k in (0,1,3,5,10,15,20):
        removed=group_names[:k]
        remain=g.loc[~g.genus.isin(removed)]
        cells=remain.groupby(["grid_x","grid_y"],sort=False).agg(
            n_after=("genus","size"),
            after=(OUTCOMES[0],"mean"),
            after_white=(OUTCOMES[1],"mean"),
        ).reset_index()
        cells=cells.loc[cells.n_after.ge(3)]
        joined=cells.merge(grid[["grid_x","grid_y",
                    "center_latitude","sampled_land_median_elevation_m",
                    "avg_genus_probability"]],on=["grid_x","grid_y"],
                    validate="one_to_one")
        for name,mask in (
            ("all",np.ones(len(joined),bool)),
            ("highland",joined.sampled_land_median_elevation_m.ge(1500)),
            ("south",joined.center_latitude.lt(0)),
            ("north",joined.center_latitude.ge(0)),
        ):
            d=joined.loc[mask]
            out.append({
                "top_genuscell_frequency_removed":k,
                "analysis_region":name,
                "n_eligible_cells_after_deleting_genera":int(len(d)),
                "original_estimate_on_same_retained_cells":float(
                    d.avg_genus_probability.mean()) if len(d) else None,
                "recomputed_estimate_without_removed_genera":float(
                    d.after.mean()) if len(d) else None,
                "recomputed_white_white_component":float(
                    d.after_white.mean()) if len(d) else None,
                "removed_genera_by_global_cell_genus_frequency":removed
            })
    return out,group_names[:20]

def run(path,out):
    root=Path(out)
    ensure(not root.exists(),"Do not overwrite existing ecological evidence")
    all_rows,grid=read_source(path)
    masks=regions(grid)
    result={}
    for i,(name,sel) in enumerate(masks.items()):
        subgrid=grid.loc[sel]
        subgen=all_rows.merge(subgrid[["grid_x","grid_y"]],
                 on=["grid_x","grid_y"],how="inner",validate="many_to_one")
        result[name]=multiplier(subgen,subgrid,
                                  seed=SEED+i*1000,nmult=N_MULT)
        print("MOIS_GENUS_SPATIAL_2WAY_MULTIPLIER",name,json.dumps({
            "n":result[name].get("n_cells"),
            "n_gen":result[name].get("n_genera"),
            "n_block":result[name].get("n_distinct_500km_blocks"),
            "space_only":result[name].get("multiplier_methods",{}).get("spatial_only",{}).get("effects",{}).get(OUTCOMES[0]),
            "both":result[name].get("multiplier_methods",{}).get("genus_x_spatial",{}).get("effects",{}).get(OUTCOMES[0])},allow_nan=False),flush=True)
    influ,top20=leave_top_n(all_rows,grid)
    root.mkdir(parents=True)
    f=root/"samecell_shared_genus_top20_influence.csv"
    pd.DataFrame([{k:v for k,v in row.items()
                   if k!="removed_genera_by_global_cell_genus_frequency"}
                  for row in influ]).to_csv(f,index=False)
    report={
       "analysis":"within-genus species photo9 same hue probability, WCVP intro−native, 500km geographic×genus two-way multiplier weights",
       "source_parent_sha256":PARENT_SHA,
       "original_37104_accepted_species":True,
       "n_original_source_species_by_100km_grid":5493931,
       "n_original_WCVP_matched_cells":1363,
       "n_original_shared_genus_cell_rows":14639,
       "n_original_unique_genera":236,
       "n_replicates_per_multiplier_method":N_MULT,
       "spatial_only_versus_crossed_are_sensitivity_intervals_not_causal_effects":True,
       "photo9_modal_species_color_not_true_petal_validation":True,
       "WCVP_TDWG3_status_not_direct_wild_introduction_success":True,
       "genus_multiplier_samples_are_not_an_explicit_species_phylogeny":True,
       "posthoc_influence_diagnostics_not_pre_registered":True,
       "top20_geographically_most_repeated_genera":top20,
       "geographical_multiplier_results":result,
       "leave_topN_genera_influence":influ,
    }
    j=root/"genus_and_spatial_independent_multiplier_stress.json"
    j.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    report["output_sha256"]={p.name:sha(p) for p in (f,j)}
    (root/"receipt.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print("MOIS_SHARED_GENUS_SPATIAL_CROSSED_STRESS_REAL_SUCCESS",
         json.dumps({
            "highland":result["highland"],
            "south":result["southern_hemisphere"],
            "topgenus_influence_high_south":[r for r in influ
                   if r["analysis_region"] in ("highland","south") and
                      r["top_genuscell_frequency_removed"] in (0,10,20)]
         },allow_nan=False)[:14000],flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--source",required=True)
    p.add_argument("--out",required=True)
    a=p.parse_args()
    run(a.source,a.out)
