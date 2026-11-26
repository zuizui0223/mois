#!/usr/bin/env python3
"""MOIS: decompose 9-hue same-genus native/introduced taxonomic assortativity
into WHITE-WHITE vs same NONWHITE photo hue pairs, and inspect geographical
confounding of high elevation and Southern Hemisphere contrasts.

Use same original signed accepted species GBIF/WCVP photo9 inputs and
the identical previously validated primary within-100km-cell+genus
n_native>=2, n_introduced>=2, >=3 eligible shared genera and >=30
source species of each status per cell. Compare 500km geographical
block bootstrap of 3 additive outcomes; 500km FE and 500km-clustered
environment adjustments using signed 2026-11-18 climate/vegetation grids.

Our outcome is P(two different accepted species share photo9 hue).
This is a compositional *species assemblage* proxy, neither a directly
photographed individual within-site hue nor an experimental mechanism of
introgression, flower/pollinator selection or introduction success.
"""
from __future__ import annotations
import argparse, hashlib, io, json
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scripts.analysis.mois_global_hue9_community_assembly_genus_null_20261118 import (
    GBIF_SHA,ENV_SHA,check,sha
)
from scripts.analysis.mois_wcvp_native_hue9_geographical_null_20261118 import (
    LGM_SHA,GEO_SHA,TDWG_GEO_SHA
)
from scripts.analysis.mois_native_intro_samecell_shared_genus_pairhue_20261118 import (
    signed_wcvp_source,genus_pairs,group_cell_support,build_cell_genus_comparison
)

SIGNED_GRID_ECOLOGY_SHA="4f5d8fc3a308482fc3645d796a779d1a4725c51ed2afbf6a06e61fd6213f3ee6"
SIGNED_GRID_ECOLOGY_MEMBER="all_current_GBIF_100km_cells_genus_breadth_entropy_ecology.csv.gz"
NPP="NPP_2001_2023_mean_kgC_m2_yr"
N_BOOT=20000
SEED=20261118
COMP=("all_same_hue","white_white","same_nonwhite_hue")
REGIONS=("all","highland","lowland","north","south","highland_north",
         "highland_south","lowland_north","lowland_south",
         "north_midlat","south_midlat")
MODELS={
 "geo_only":("highland","southern","abs_latitude",
             "log_native_species","log_intro_species","log_common_genera","island"),
 "modern_environment":("highland","southern","abs_latitude",
            "log_native_species","log_intro_species","log_common_genera","island",
            "temperature","log_precip","NPP","NDVI","log_AGB"),
 "modern_plus_LGM":("highland","southern","abs_latitude",
            "log_native_species","log_intro_species","log_common_genera","island",
            "temperature","log_precip","NPP","NDVI","log_AGB",
            "LGM_temperature_change","LGM_precip_change"),
 "modern_plus_LGM_highland_south_interaction":(
            "highland","southern","highland_x_southern","abs_latitude",
            "log_native_species","log_intro_species","log_common_genera","island",
            "temperature","log_precip","NPP","NDVI","log_AGB",
            "LGM_temperature_change","LGM_precip_change"),
}

def ecology_source(path):
    check(sha(path)==SIGNED_GRID_ECOLOGY_SHA,
          "Source-signed regional ecological grid ZIP changed")
    with ZipFile(path) as z:
        receipt=json.loads(z.read("receipt.json"))
        for filename,checksum in receipt["numeric_outputs_sha256"].items():
            check(sha(z.read(filename))==checksum,
                  "Source ecology numerical CSV/JSON SHA changed")
        frame=pd.read_csv(io.BytesIO(z.read(SIGNED_GRID_ECOLOGY_MEMBER)),
                          compression="gzip")
    check(len(frame)==12295 and
          not frame[["grid_x","grid_y"]].duplicated().any(),
          "12,295 authentic complete case GBIF ecological grid changed")
    return frame

def enrich_pair_columns(paired):
    df=paired.loc[paired.n_min.ge(2)].copy()
    n=df.n_species_native.to_numpy(float)
    m=df.n_species_introduced.to_numpy(float)
    wn=df.white_species_native.to_numpy(float)
    wi=df.white_species_introduced.to_numpy(float)
    # Exactly partition the unbiased 9-hue pair statistic:
    # P(same hue) = P(white-white) + P(same non-white hue).
    pww_n=wn*(wn-1)/(n*(n-1))
    pww_i=wi*(wi-1)/(m*(m-1))
    ptotal_n=df.pair_prob_same_native.to_numpy(float)
    ptotal_i=df.pair_prob_same_introduced.to_numpy(float)
    df["delta_white_white"]=pww_i-pww_n
    df["delta_same_nonwhite"]=(ptotal_i-pww_i)-(ptotal_n-pww_n)
    delta=df.introduced_minus_native_same_photo_hue_prob.to_numpy(float)
    check(np.allclose(
        delta,df.delta_white_white+df.delta_same_nonwhite,atol=1e-12,rtol=0),
        "Original two-species same-hue probability not exactly partitioned")
    return df

def source_matched_primary(paired,totals,geo,mode):
    d=build_cell_genus_comparison(paired,
          group_cell_support(totals,geo,mode),mode,
          per_genus_min=2,min_shared_genera=3,status_min_species=30)
    expect={"strict":1363,"moderate":1640}[mode]
    check(len(d)==expect,"Original same-cell 9-hue common-genus primary denominator changed")
    sub=enrich_pair_columns(paired)
    sub=sub.merge(d[["grid_x","grid_y"]],
                  on=["grid_x","grid_y"],how="inner",validate="many_to_one")
    sums=sub.groupby(["grid_x","grid_y"],sort=False).agg(
        mean_white_white=("delta_white_white","mean"),
        mean_same_nonwhite=("delta_same_nonwhite","mean"),
        n_shared_genera_replay=("genus","size")
    ).reset_index()
    d=d.merge(sums,on=["grid_x","grid_y"],validate="one_to_one")
    check((d.n_shared_genera==d.n_shared_genera_replay).all(),
          "Original per-cell matched genus count not reproduced")
    d["delta_all_same_hue"]=d.avg_genus_probability
    d["delta_white_white"]=d.mean_white_white
    d["delta_same_nonwhite_hue"]=d.mean_same_nonwhite
    check(np.allclose(d.delta_all_same_hue,
          d.delta_white_white+d.delta_same_nonwhite_hue,
          atol=1e-12,rtol=0),"Cell total probability must exactly decompose")
    return d,sub

def masks(d):
    e=d.sampled_land_median_elevation_m
    lat=d.center_latitude
    high=e.ge(1500)
    low=e.ge(-430)&e.lt(500)
    south=lat.lt(0)
    north=lat.ge(0)
    return {
      "all":np.ones(len(d),bool),
      "highland":high.to_numpy(bool),"lowland":low.to_numpy(bool),
      "north":north.to_numpy(bool),"south":south.to_numpy(bool),
      "highland_north":(high&north).to_numpy(bool),
      "highland_south":(high&south).to_numpy(bool),
      "lowland_north":(low&north).to_numpy(bool),
      "lowland_south":(low&south).to_numpy(bool),
      "north_midlat":(lat.ge(30)&lat.lt(55)).to_numpy(bool),
      "south_midlat":(lat.le(-30)&lat.gt(-55)).to_numpy(bool),
    }

def blockboot(d,subset,mode):
    df=d.loc[subset]
    if len(df)<25:
        return {"status":"insufficient_25_cells","n_cells":len(df)}
    ids=df.geographic_500km_block.to_numpy(str)
    cols=["delta_all_same_hue","delta_white_white","delta_same_nonwhite_hue"]
    wg=np.ones(len(df),float)
    t=pd.DataFrame({"b":ids,"w":wg})
    for c in cols:t[c]=wg*df[c].to_numpy(float)
    grouped=t.groupby("b")[["w",*cols]].sum().to_numpy(float)
    if len(grouped)<8:
        return {"status":"insufficient_8_spatial_blocks","n_cells":len(df),
                "n_blocks":len(grouped)}
    rng=np.random.default_rng(SEED+len(df)+(0 if mode=="strict" else 6000))
    draws=np.empty((N_BOOT,len(cols)))
    for k in range(0,N_BOOT,250):
        n=min(250,N_BOOT-k)
        iid=rng.integers(0,len(grouped),size=(n,len(grouped)))
        batch=grouped[iid].sum(axis=1)
        draws[k:k+n]=batch[:,1:]/batch[:,:1]
    res={
        "n_cells":len(df),"n_500km_blocks":len(grouped),
        "status":"ok","n_boot":N_BOOT,
        "effects":{c:{"mean":float(df[c].mean()),
                      "CI95":[float(z) for z in
                            np.percentile(draws[:,i],[2.5,97.5])]}
                   for i,c in enumerate(cols)}
    }
    check(abs(sum([res["effects"][c]["mean"] for c in cols[1:]])-
              res["effects"][cols[0]]["mean"])<1e-10,
          "Bootstrap population mean components don't add up")
    return res

def environmental_matrix(d):
    e=d.sampled_land_median_elevation_m.astype(float)
    lat=d.center_latitude.astype(float)
    raw=pd.DataFrame({
       "highland":e.ge(1500).astype(float),
       "southern":lat.lt(0).astype(float),
       "highland_x_southern":(e.ge(1500)&lat.lt(0)).astype(float),
       "abs_latitude":lat.abs(),
       "log_native_species":np.log1p(d.total_native.astype(float)),
       "log_intro_species":np.log1p(d.total_introduced.astype(float)),
       "log_common_genera":np.log1p(d.n_shared_genera.astype(float)),
       "island":d.gshhg_land_class.eq("island").astype(float),
       "temperature":d.present_temperature_C.astype(float),
       "log_precip":np.log1p(d.present_precip_mm.astype(float)),
       "NPP":d[NPP].astype(float),
       "NDVI":d.MODIS_NDVI_2001_2021_mean.astype(float),
       "log_AGB":np.log1p(d.ESA_AGB_woody_Mg_ha_2020.astype(float)),
       "LGM_temperature_change":d.threeGCM_delta_temperature_C_median.astype(float),
       "LGM_precip_change":d.threeGCM_delta_log_precip_median.astype(float),
    },index=d.index)
    check(np.isfinite(raw.to_numpy()).all(),"Nonfinite independent ecological exposure")
    scale={}
    for col in raw:
        if col in ("highland","southern","highland_x_southern","island"):
            continue
        mean=float(raw[col].mean());sd=float(raw[col].std(ddof=1))
        check(sd>1e-9,"No ecological variance: "+col)
        raw[col]=(raw[col]-mean)/sd
        scale[col]={"mean":mean,"SD":sd}
    return raw,scale

def regress(d,field,terms,spatial=False):
    X,scales=environmental_matrix(d)
    X=X[list(terms)].copy()
    y=d[field].astype(float)
    blocks=d.geographic_500km_block.astype(str)
    if spatial:
        X=X-X.groupby(blocks).transform("mean")
        y=y-y.groupby(blocks).transform("mean")
        keep=X.abs().sum(axis=1)>1e-10
        X=X.loc[keep]
        y=y.loc[keep]
        blocks=blocks.loc[keep]
        # Indicators that cannot vary in any 500km block cannot be identified.
        X=X.loc[:,X.std().gt(1e-9)]
        if "highland" not in X:
            return {"status":"highland_has_no_within_500km_block_support"}
    else:
        X=sm.add_constant(X,has_constant="add")
    rank=int(np.linalg.matrix_rank(X.to_numpy()))
    check(rank==len(X.columns),"Design has no full rank in same-cell geographic model")
    model=sm.OLS(y,X).fit(cov_type="cluster",
                          cov_kwds={"groups":blocks},use_t=False)
    return {
        "status":"ok","n_cells":len(y),
        "n_blocks":int(blocks.nunique()),
        "cluster":"500km",
        "within_500km_block_FE":spatial,
        "R2_in_sample":float(model.rsquared),
        "terms":{term:{
             "beta":float(model.params[term]),
             "CI95":[float(z) for z in model.conf_int().loc[term]],
             "p_cluster":float(model.pvalues[term])
         } for term in model.params.index}
    }

def analyze(src,geo,env,mode):
    genus,totals=genus_pairs(src,mode)
    d,paired=source_matched_primary(genus,totals,geo,mode)
    d=d.merge(env[["grid_x","grid_y","present_temperature_C",
        "present_precip_mm",NPP,"MODIS_NDVI_2001_2021_mean",
        "ESA_AGB_woody_Mg_ha_2020",
        "threeGCM_delta_temperature_C_median",
        "threeGCM_delta_log_precip_median"]],
        on=["grid_x","grid_y"],how="left",validate="one_to_one")
    usable=d.dropna(subset=["present_temperature_C","present_precip_mm",
      NPP,"MODIS_NDVI_2001_2021_mean","ESA_AGB_woody_Mg_ha_2020",
      "threeGCM_delta_temperature_C_median",
      "threeGCM_delta_log_precip_median"]).reset_index(drop=True)
    check(len(usable)>1200 and
          d.geographic_500km_block.nunique()>250,
          "Insufficient original WCVP samecell paired colour and environment")
    m=masks(d)
    b={k:blockboot(d,which,mode) for k,which in m.items()}
    # Global regression with 500km cluster errors; cannot identify the
    # 'southern' exposure under 500km FE because hemispheres do not share
    # 500km equal-area land blocks. This is stated, not concealed.
    fits={}
    for name,terms in MODELS.items():
        fits[name]=regress(usable,"delta_all_same_hue",terms)
    full=MODELS["modern_plus_LGM"]
    fits["500km_withinblock_modern_LGM"]=regress(
        usable,"delta_all_same_hue",full,spatial=True)
    for outcome in ("delta_white_white","delta_same_nonwhite_hue"):
        fits[outcome+"_full_climate"]=regress(
            usable,outcome,full)
    source_n=len(d)
    evidence={
      "n_eligible_samecell_common_genus_assemblages":source_n,
      "n_complete_case_climate_cells":len(usable),
      "n_500km_original_blocks":int(d.geographic_500km_block.nunique()),
      "n_highland_both_hemispheres":int((d.sampled_land_median_elevation_m.ge(1500)&
          d.center_latitude.ge(0)).sum()),
      "n_highland_southern_hemisphere":int((d.sampled_land_median_elevation_m.ge(1500)&
          d.center_latitude.lt(0)).sum()),
      "region_components":b,
      "environment_models":fits,
    }
    return d,paired,evidence

def run(args):
    root=Path(args.out)
    check(not root.exists(),"Never overwrite a previously verified biological result")
    src,geo,original=signed_wcvp_source(
        args.GBIF,args.ENV,args.LGM,args.WCVP,args.TDWG,args.GEO)
    eco=ecology_source(args.ECO)
    outcomes={}
    rows=[]
    genera=[]
    for mode in ("strict","moderate"):
        d,g,result=analyze(src,geo,eco,mode)
        print("MOIS_SHARED_GENUS_WHITEPAIR_CHROMATIC_CLIMATE_REAL_SUCCESS",
              json.dumps({
               "mode":mode,
               "n_cells":len(d),
               "n_high_north":result["n_highland_both_hemispheres"],
               "n_high_south":result["n_highland_southern_hemisphere"],
               "highland":result["region_components"]["highland"],
               "southern":result["region_components"]["south"],
               "full_climate":{
                 term:result["environment_models"]["modern_plus_LGM"]["terms"][term]
                 for term in ("highland","southern","temperature","log_precip","NPP")},
               "spatial_FE_highland":result["environment_models"][
                 "500km_withinblock_modern_LGM"]["terms"].get("highland")
              },allow_nan=False),flush=True)
        outcomes[mode]=result
        export=["grid_x","grid_y","center_latitude",
                "sampled_land_median_elevation_m","gshhg_land_class",
                "geographic_500km_block","n_shared_genera",
                "total_native","total_introduced",
                "avg_genus_probability",
                "delta_white_white","delta_same_nonwhite_hue",
                "avg_genus_white_diff","present_temperature_C","present_precip_mm",
                NPP]
        z=d[export].copy()
        z["WCVP_mode"]=mode
        rows.append(z)
        if mode=="strict":
            p=g[["grid_x","grid_y","genus","n_species_native",
                 "n_species_introduced","n_min",
                 "introduced_minus_native_same_photo_hue_prob",
                 "delta_white_white","delta_same_nonwhite"]].copy()
            genera.append(p)
    root.mkdir(parents=True)
    grid_file=root/"same_grid_shared_genus_white_chromatic_pair_components.csv.gz"
    pd.concat(rows,ignore_index=True).to_csv(grid_file,compression="gzip",index=False)
    genus_file=root/"strict_same_grid_per_genus_colour_pair_evidence.csv.gz"
    pd.concat(genera,ignore_index=True).to_csv(genus_file,compression="gzip",index=False)
    result={
      "source_photo9_accepted_species":37104,
      "source_actual_GBIF_species_grid_pairs":5493931,
      "original_WCVP_native_introduced_GBIF_pairs":3127248,
      "source_sha256":{"GBIF":GBIF_SHA,"ENV":ENV_SHA,"GEO":GEO_SHA,
                       "LGM":LGM_SHA,"TDWG3":TDWG_GEO_SHA,
                       "same_cell_ECO":SIGNED_GRID_ECOLOGY_SHA},
      "n_block_bootstraps":N_BOOT,
      "actual_national_regional_status_not_wild_introduced_individual":True,
      "photo9_hue_is_not_independently_vetted_flower_petal":True,
      "two_species_same_hue_9class_exact_decomposition":True,
      "geographical_climate_covariates_not_causal_elevation_selection":True,
      "global_south_north_cluster_errors_do_not_establish_realm_independence":True,
      "results":outcomes
    }
    f=root/"samecell_genus_photo9_whitepair_chromatic_climate.json"
    f.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    result["output_sha256"]={p.name:sha(p) for p in (grid_file,genus_file,f)}
    (root/"receipt.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")

if __name__=="__main__":
    p=argparse.ArgumentParser()
    for name in ("GBIF","ENV","LGM","GEO","WCVP","TDWG","ECO","out"):
        p.add_argument("--"+name,required=True)
    run(p.parse_args())
