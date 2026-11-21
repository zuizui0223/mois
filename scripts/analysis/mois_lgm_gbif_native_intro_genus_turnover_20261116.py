#!/usr/bin/env python3
"""GBIF ONLY taxonomic outcome: overlap of native/introduced genera per 100km cell.

Source species presence: 2026-11-01 GBIF snapshot, 37104 accepted taxa.
Broad-region introduced/native status: Kew WCVP TDWG3.
Environmental exposure: same-period WorldClim14 paired LGM/current + AGB.
No image or species floral colour used to CALCULATE genus overlap/Jaccard.
Photo-hue is carried only in previously signed comparative reference models.
Neither genus overlap nor LGM climate change reconstructs past plant ranges.
"""
from __future__ import annotations

import argparse
import io
import json
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pandas as pd
import statsmodels.api as sm

from scripts.analysis.mois_lgm_climate_history_species_20261116 import (
    read_signed_sources,sha,require
)
from scripts.analysis.mois_lgm_wcvp_native_regional_block_20261116 import (
    assign_cells_broad_region, get_wcvp_source
)
from scripts.analysis.mois_lgm_native_intro_common_genus_20261116 import (
    regional_species_cell_rows
)
from scripts.analysis.mois_lgm_native_introduced_matched_cells_20261116 import (
    signed_tables,same_cell_cohort,matrix
)

LGM_SHA="db3b5b4fa2133194963be9233411bb6bb50bc86393d1dc50d996cc4cc78b25ab"
LGM_CELL="WorldClim1p4_LGM_threeGCM_18577_cells.parquet"
COMMON_SHA="ba725090afbd732c9dacc8e098161aad2bf88042994f5c683db465519969356a"
COMMON_CELL="samecell_native_intro_shared_genera_100km_summary.parquet"
TARGET={"strict":3106,"moderate":3602}

def genus_overlap(source):
    """Counts from GBIF botanical-country species pairs only, no hues."""
    require(not source[["specieskey","grid_x","grid_y"]].duplicated().any(),
            "GBIF species presence duplicated within a cell")
    g=source[["grid_x","grid_y","genus","status"]].drop_duplicates()
    require(g.status.isin(["native_at_TDWG3","introduced_at_TDWG3"]).all(),
            "only WCVP accepted regional native/introduced statuses")
    g["native"]=g.status.eq("native_at_TDWG3")
    g["introduced"]=g.status.eq("introduced_at_TDWG3")
    bygenus=g.groupby(["grid_x","grid_y","genus"],sort=False).agg(
        native=("native","any"),introduced=("introduced","any")).reset_index()
    bygenus["both"]=bygenus.native & bygenus.introduced
    bygenus["native_exclusive"]=bygenus.native & ~bygenus.introduced
    bygenus["introduced_exclusive"]=bygenus.introduced & ~bygenus.native
    z=bygenus.groupby(["grid_x","grid_y"],sort=False).agg(
        native_genus_n=("native","sum"),
        introduced_genus_n=("introduced","sum"),
        shared_genus_n=("both","sum"),
        native_exclusive_genus_n=("native_exclusive","sum"),
        introduced_exclusive_genus_n=("introduced_exclusive","sum")).reset_index()
    union=(z.native_genus_n+z.introduced_genus_n-z.shared_genus_n)
    require((union>0).all(),"unexpected empty genus union")
    require((z.native_exclusive_genus_n==z.native_genus_n-z.shared_genus_n).all()
            and (z.introduced_exclusive_genus_n==z.introduced_genus_n-z.shared_genus_n).all(),
            "GBIF species pool arithmetic not consistent")
    z["jaccard_native_introduced_genus"]=z.shared_genus_n/union
    z["genus_turnover_one_minus_jaccard"]=1-z.jaccard_native_introduced_genus
    z["native_exclusive_genus_fraction"]=(
        z.native_exclusive_genus_n/z.native_genus_n.replace(0,np.nan))
    z["introduced_exclusive_genus_fraction"]=(
        z.introduced_exclusive_genus_n/z.introduced_genus_n.replace(0,np.nan))
    require(z.jaccard_native_introduced_genus.between(0,1).all()
            and z.genus_turnover_one_minus_jaccard.between(0,1).all(),
            "genus turnover outside 0-1")
    return z

def model(d,response,block="500",add_genus_overlap=False):
    x,sd=matrix(d)
    if add_genus_overlap:
        term=d.genus_turnover_one_minus_jaccard.astype(float)
        term_sd=term.std()
        require(np.isfinite(term_sd) and term_sd>1e-12,
                "genus overlap term does not vary")
        x["genus_turnover"]=(term-term.mean())/term_sd
    y=d[response].astype(float)
    b=d.block_500km.astype(str) if block=="500" else d.block_1000km.astype(str)
    xx=x-x.groupby(b).transform("mean")
    yy=y-y.groupby(b).transform("mean")
    useful=xx.abs().sum(axis=1)>0
    xx=xx.loc[useful]; yy=yy.loc[useful]
    cl=d.loc[useful,"block_500km"].astype(str)
    require(len(yy)>200 and cl.nunique()>25,"too little spatial replication")
    m=sm.OLS(yy,xx).fit(cov_type="cluster",
                        cov_kwds={"groups":cl},use_t=False)
    cols=["LGM_warming"]+(["genus_turnover"] if add_genus_overlap else [])
    effects={}
    for col in cols:
        lo,hi=m.conf_int().loc[col].tolist()
        effects[col]={"beta_per_SD":float(m.params[col]),
                      "95CI":[float(lo),float(hi)],
                      "p":float(m.pvalues[col])}
    return {"n_cells_model":len(yy),"n_clusters_500km":int(cl.nunique()),
            "outcome":response,"block_km":int(block),
            "adjusted_for_genus_turnover":bool(add_genus_overlap),
            "effects":effects,"r_squared":float(m.rsquared),
            "global_SD_LGM_warming_C":sd}

def signed_context(gbif,env,lgm,native,spatial,common):
    import pyarrow.parquet as pq
    pairs,species=read_signed_sources(gbif,env)
    require(sha(lgm)==LGM_SHA,"LGM 3 GCM ZIP source signature mismatch")
    with ZipFile(lgm) as z:
        v=json.loads(z.read("receipt.json"))
        b=z.read(LGM_CELL)
    require(sha(b)==v["output_SHA256"][LGM_CELL],"LGM cell Parquet changed")
    climate=pq.read_table(io.BytesIO(b)).to_pandas()
    signed_native,signed_spatial=signed_tables(native,spatial)
    require(sha(common)==COMMON_SHA,"common genus artifact source signature changed")
    with ZipFile(common) as z:
        v=json.loads(z.read("receipt.json"))
        raw=z.read(COMMON_CELL)
    require(sha(raw)==v["output_sha256"][COMMON_CELL],
            "original genus shared count table signed SHA mismatch")
    shared=pq.read_table(io.BytesIO(raw)).to_pandas()[[
         "grid_x","grid_y","n_shared_genera"]].rename(
         columns={"n_shared_genera":"previous_signed_shared_genera"})
    return pairs,species,climate,signed_native,signed_spatial,shared

def run(gbif,env,lgm,native,spatial,common,wcvp,tdwg,out):
    out=Path(out)
    require(not out.exists(),"refuse overwrite")
    pairs,species,climate,native_s,spatial_s,shared=signed_context(
        gbif,env,lgm,native,spatial,common)
    names,distribution,receipt=get_wcvp_source(wcvp)
    regions=assign_cells_broad_region(climate,tdwg)
    matched={}
    ids=[]
    for mode in TARGET:
        d=same_cell_cohort(native_s,spatial_s,mode,30)
        require(len(d)==TARGET[mode],"WCVP GBIF matched 100km cohort drift")
        matched[mode]=d
        ids.append(d[["grid_x","grid_y"]])
    ids=pd.concat(ids,ignore_index=True).drop_duplicates()
    # Source uses only GBIF species presence and WCVP region statuses.
    source=regional_species_cell_rows(
        pairs,species,names,distribution,regions,ids)
    overlap=genus_overlap(source)
    z=overlap.merge(shared,on=["grid_x","grid_y"],how="left",
                    validate="one_to_one")
    both=z.previous_signed_shared_genera.notna()
    require(np.array_equal(
        z.loc[both,"shared_genus_n"].to_numpy(int),
        z.loc[both,"previous_signed_shared_genera"].to_numpy(int)),
        "independent native vs introduced genus sets do not reproduce signed common genera")
    result={"GBIF_accepted_species":37104,
            "GBIF_current_species_x_100km_cell_pairs":5493931,
            "native_introduced_status_broad_WCVP_TDWG3_only":True,
            "genus_turnover_computed_from_GBIF_WCVP_no_flower_colour":True,
            "LGM_source_measures_historical_climate_at_current_occupied_cells":True,
            "LGM_refugia_inferred":False,"regional_models":{}}
    tables=[]
    for mode,d in matched.items():
        u=d.merge(z.drop(columns=["previous_signed_shared_genera"]),
                  on=["grid_x","grid_y"],how="left",validate="one_to_one")
        require(len(u)==TARGET[mode] and u.genus_turnover_one_minus_jaccard.notna().all(),
                "missing true GBIF taxonomic community contrast")
        result["regional_models"][mode]={
            "n_matched_cells":len(u),
            "mean_jaccard":float(u.jaccard_native_introduced_genus.mean()),
            "median_jaccard":float(u.jaccard_native_introduced_genus.median()),
            "mean_genus_turnover":float(u.genus_turnover_one_minus_jaccard.mean()),
            "mean_native_genus_count":float(u.native_genus_n.mean()),
            "mean_introduced_genus_count":float(u.introduced_genus_n.mean()),
            "mean_shared_genus_count":float(u.shared_genus_n.mean()),
            "models":{}}
        for outcome in ("genus_turnover_one_minus_jaccard",
                        "native_exclusive_genus_fraction",
                        "introduced_exclusive_genus_fraction",
                        "native_minus_intro_white"):
            result["regional_models"][mode]["models"][outcome]={
                "500km":model(u,outcome,"500"),
                "1000km":model(u,outcome,"1000")}
        result["regional_models"][mode]["models"]["photo_white_difference_with_turnover_covariate"]={
            "500km":model(u,"native_minus_intro_white","500",True),
            "1000km":model(u,"native_minus_intro_white","1000",True)}
        u["tdwg_assignment_mode"]=mode
        columns=["grid_x","grid_y","tdwg_assignment_mode",
                 "native_genus_n","introduced_genus_n",
                 "native_exclusive_genus_n","introduced_exclusive_genus_n",
                 "shared_genus_n","jaccard_native_introduced_genus",
                 "genus_turnover_one_minus_jaccard",
                 "native_exclusive_genus_fraction",
                 "introduced_exclusive_genus_fraction",
                 "threeGCM_delta_temperature_C_median","lat",
                 "native_minus_intro_white","native_species_n","introduced_species_n"]
        tables.append(u[columns].copy())
        m=result["regional_models"][mode]["models"]
        print("MOIS_GBIF_ONLY_GENUS_TURNOVER",json.dumps({
            "mode":mode,"n":len(u),
            "genus_turnover_500km_beta":m["genus_turnover_one_minus_jaccard"]["500km"]["effects"]["LGM_warming"]["beta_per_SD"],
            "genus_turnover_500km_p":m["genus_turnover_one_minus_jaccard"]["500km"]["effects"]["LGM_warming"]["p"],
            "hue_difference_500km_beta":m["native_minus_intro_white"]["500km"]["effects"]["LGM_warming"]["beta_per_SD"],
            "hue_adjusted_for_genus_turnover_500km_beta":m["photo_white_difference_with_turnover_covariate"]["500km"]["effects"]["LGM_warming"]["beta_per_SD"]
        },sort_keys=True),flush=True)
    out.mkdir(parents=True)
    f=out/"GBIF100km_WCVP_native_introduced_genus_overlap_no_photo_response.parquet"
    pd.concat(tables,ignore_index=True).to_parquet(f,index=False)
    j=out/"models.json"
    j.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    audit={"schema":"mois_GBIF_WCVP_real_current_genus_turnover_vs_LGM_history",
           "source_sha256":{"GBIF":sha(gbif),"ENV":sha(env),
                "LGM":sha(lgm),"WCVP_native":sha(native),
                "spatial":sha(spatial),"common":sha(common),
                "WCVP_names":sha(Path(wcvp)/"mois350399_wcvp_source_taxon_exact_name_resolution.csv.gz"),
                "TDWG3":sha(tdwg)},
           "trait_free_taxonomic_response":True,
           "cannot_infer_genus_change_since_LGM_from_current_data":True,
           "result":result,"output_sha256":{f.name:sha(f),j.name:sha(j)}}
    (out/"receipt.json").write_text(json.dumps(audit,indent=2,sort_keys=True)+"\n")
    return audit

if __name__=="__main__":
    p=argparse.ArgumentParser()
    for name in ("gbif","env","lgm","native","spatial","common","wcvp","tdwg","output"):
        p.add_argument("--"+name,required=True)
    a=p.parse_args()
    run(a.gbif,a.env,a.lgm,a.native,a.spatial,a.common,
        a.wcvp,a.tdwg,a.output)
