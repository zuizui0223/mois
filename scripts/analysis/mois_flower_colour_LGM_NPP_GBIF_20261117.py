#!/usr/bin/env python3
"""MOIS primary FLOWER-COLOUR response: NPP / woody AGB / NDVI / LGM.

Primary response is species-level WHITE versus NONWHITE photographic modal
9-hue label, linked by accepted plant key to genuine GBIF occupied 100km
cells. Secondary response is the proportion of these photo-labelled species
in actual GBIF native/introduced region-wise 100km assemblages.

- LGM (~22ka) paired vs current WorldClim1.4 environmental exposure is
  measured at *current occupied cells*, not evidence of refugial presence.
- Native/introduced is WCVP broad TDWG3 status, not wild individuals/islands.
- Photo9 labels may describe nonpetal image content and are not proven
  intraspecific genetic fixed floral colours.
- MOD17 C061 NPP 2001-2023 is annual kg C/m2/year FLUX, not floral biomass.
- Species-level comparison: one row/species, genus FE and genus-cluster SE.
- Assemblage comparison: 100km GBIF cells, 500km/1000km spatial FE and SE.
Do not present as a causal effect or independent paleobotanical range evidence.
"""
from __future__ import annotations
import argparse,io,json,math
from pathlib import Path
from zipfile import ZipFile
import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from scripts.analysis.mois_lgm_climate_history_species_20261116 import (
    read_signed_sources,genus_model,require,sha,HUES
)
from scripts.analysis.mois_gbif_family_excess_biomass_ndvi_20261116 import (
    signed_parquet,
    NATIVE_SHA,SPATIAL_SHA,NDVI_SHA,
    NATIVE_MEMBER,SPATIAL_MEMBER,NDVI_MEMBER
)
from scripts.analysis.mois_gbif_family_genus_npp_AGB_NDVI_20261116 import (
    NPP_SHA,NPP_MEMBER,NPP_FIELD
)

LGM_SHA="db3b5b4fa2133194963be9233411bb6bb50bc86393d1dc50d996cc4cc78b25ab"
HISTORY_MEMBER="GBIF37104_photo9_LGM_climate_species.parquet"
NPP_NAME="NPP_mean_GBIF_2001_2023_kgC_m2_yr"
NPP_COVER="NPP_fraction_occupied_cells_18plus_years"
LGM_TEMP="threeGCM_delta_temperature_C_median"
LGM_RAIN="threeGCM_delta_log_precip_median"
CTRL=["BIO01","BIO12","MODIS_NDVI_2001_2021_mean",
      "ESA_AGB_woody_Mg_ha_2020","mean_abs_latitude_over_occupied_cells",
      "records_passed_spatial_QA"]

def load_species_history(lgm_archive):
    import pyarrow.parquet as pq
    require(sha(lgm_archive)==LGM_SHA,"signed original true LGM ZIP changed")
    with ZipFile(lgm_archive) as z:
        raw=z.read(HISTORY_MEMBER)
        rec=json.loads(z.read("receipt.json"))
    require(sha(raw)==rec["output_SHA256"][HISTORY_MEMBER]
            and rec["gbif_accepted_species"]==37104
            and rec["species_x_grid_pairs"]==5493931,
            "previously confirmed true GBIF photographic-colour universe changed")
    d=pq.read_table(io.BytesIO(raw)).to_pandas()
    d["key"]=d.key.astype(str)
    require(len(d)==37104 and d.key.is_unique
            and d.source_hue9.isin(HUES).all()
            and d.genus.notna().all(),"species photograph label/genus identity changed")
    return d

def make_species_NPP(pairs,npp_cells,historical):
    require(len(pairs)==5493931 and len(npp_cells)==18577,
            "true source GBIF species / NPP cell count drift")
    require(not npp_cells[["grid_x","grid_y"]].duplicated().any(),
            "NPP cell duplicated")
    cols=["grid_x","grid_y",NPP_FIELD,"NPP_valid_years_2001_2023",
          "NPP_2020_kgC_m2_yr"]
    p=pairs[["specieskey","grid_x","grid_y"]].merge(
        npp_cells[cols],on=["grid_x","grid_y"],how="left",validate="many_to_one")
    require(len(p)==5493931,"GBIF occurrence distribution not preserved")
    p["NPP_good"]=p["NPP_valid_years_2001_2023"].ge(18) & p[NPP_FIELD].notna()
    p["npp_usable"]=p[NPP_FIELD].where(p.NPP_good,np.nan)
    p["npp_2020_valid"]=p["NPP_2020_kgC_m2_yr"].where(
        p.NPP_good,np.nan)
    grouped=p.groupby("specieskey",sort=False).agg(
        occupied_cells=("grid_x","size"),
        NPP_occupied_cells_valid=("NPP_good","sum"),
        NPP_mean_GBIF_2001_2023_kgC_m2_yr=("npp_usable","mean"),
        NPP_2020_GBIF_kgC_m2_yr=("npp_2020_valid","mean"),
        NPP_sd_across_occupied_cells=("npp_usable","std"))
    grouped[NPP_COVER]=grouped.NPP_occupied_cells_valid/grouped.occupied_cells
    s=historical.merge(grouped,left_on="key",right_index=True,how="left",
                       validate="one_to_one")
    require(len(s)==37104 and s.key.is_unique,"joined species no longer unique")
    require((s.occupied_cells.fillna(0).astype(int)
             ==s.lgm_all_occupied_cells.fillna(0).astype(int)).all(),
            "species occupation denominator differs from signed LGM source")
    s["white_photo_label"]=s.source_hue9.eq("white").astype("int8")
    return s

def prepare_species(s,minimum_LGM=.75,minimum_NPP=.75,
                    high_photo_confidence=False,complete_23years=False):
    cols=[*CTRL,LGM_TEMP,LGM_RAIN,NPP_NAME,"lgm_climate_coverage",NPP_COVER]
    ok=s[cols].notna().all(axis=1)
    ok &= (s.lgm_climate_coverage>=minimum_LGM) & (s[NPP_COVER]>=minimum_NPP)
    ok &= (s.BIO12>=0)&(s.ESA_AGB_woody_Mg_ha_2020>=0)
    ok &= s.hue_status.notna()
    if high_photo_confidence:
        ok &= s.hue_status.eq("sampling_consensus_single")
    if complete_23years:
        # Full 23-years at *all* included modern occupied cells.
        # GBIF species-level NPP_n years are summarized by cell coverage.
        # Current data have NPP valid cell count, not separate count of 23-year
        # cells; avoid claiming that >=18y support implies complete 23y.
        raise ValueError("full 23-year species sensitivity unavailable without 23-year cell numerator")
    return s.loc[ok].copy()

def species_model_advanced(d,y,modes=("species_equal","genus_equal"),
                           extra=("NPP",),year2020=False):
    """Species within-genus hue comparison on exactly SAME complete rows."""
    response=d.source_hue9.eq(y).astype(float)
    grp=d.genus.astype(str)
    keep=response.groupby(grp).transform("nunique").eq(2)
    z=d.loc[keep].copy()
    yv=response.loc[keep].copy()
    if len(z)<150 or z.genus.nunique()<25:
        return {"status":"insufficient_within_genus_photographic_colour_contrast",
                "n_species":len(z),"n_genera":z.genus.nunique()}
    G=z.genus.astype(str)
    nppcol="NPP_2020_GBIF_kgC_m2_yr" if year2020 else NPP_NAME
    col_values={
        "LGM_warming":z[LGM_TEMP].astype(float),
        "LGM_rainfall_change":z[LGM_RAIN].astype(float),
        "NPP":z[nppcol].astype(float),
        "current_temperature":z.BIO01.astype(float),
        "log_current_precip":np.log1p(z.BIO12.astype(float)),
        "NDVI":z.MODIS_NDVI_2001_2021_mean.astype(float),
        "log_woody_AGB":np.log1p(z.ESA_AGB_woody_Mg_ha_2020.astype(float)),
        "absolute_lat":z.mean_abs_latitude_over_occupied_cells.astype(float),
        "log_GBIF_records":np.log1p(z.records_passed_spatial_QA.astype(float)),
    }
    X=pd.DataFrame(col_values,index=z.index)
    for status in sorted(z.hue_status.astype(str).unique())[1:]:
        X["status_"+status]=(z.hue_status.astype(str)==status).astype(float)
    require(np.isfinite(X.to_numpy()).all(),"nonfinite species model covariates")
    standard=["LGM_warming","LGM_rainfall_change","NPP",
              "current_temperature","log_current_precip","NDVI",
              "log_woody_AGB","absolute_lat","log_GBIF_records"]
    std=X[standard].std()
    require((std>1e-9).all(),"species predictor lacks variance")
    X[standard]=(X[standard]-X[standard].mean())/std
    XC=X-X.groupby(G).transform("mean")
    YC=yv-yv.groupby(G).transform("mean")
    XC=XC.loc[:,XC.std()>1e-9]
    results={"status":"fitted","n_species":len(z),"n_genera":G.nunique(),
             "fraction_white_on_full_sample":float(yv.mean()),
             "year2020_only":bool(year2020),"models":{}}
    for weighted in modes:
        weights=(1/G.map(G.value_counts()).astype(float)
                 if weighted=="genus_equal" else None)
        results["models"][weighted]={}
        for npp_included in (True,False):
            keep_cols=XC.columns if npp_included else XC.columns.drop("NPP")
            x=XC[keep_cols]
            if weights is None:
                model=sm.OLS(YC,x).fit(cov_type="cluster",
                                        cov_kwds={"groups":G},use_t=False)
            else:
                model=sm.WLS(YC,x,weights=weights).fit(
                    cov_type="cluster",cov_kwds={"groups":G},use_t=False)
            eff={}
            for name in ("LGM_warming","LGM_rainfall_change","NPP",
                         "NDVI","log_woody_AGB"):
                if name not in x:continue
                low,hi=model.conf_int().loc[name].tolist()
                eff[name]={"beta_probability_points_per_sd":float(model.params[name]),
                           "95CI":[float(low),float(hi)],
                           "p":float(model.pvalues[name])}
            results["models"][weighted]["with_NPP" if npp_included else "without_NPP"]={
                "effects":eff,"R2_in_sample":float(model.rsquared),
                "n_species":len(z),"n_genera":G.nunique(),
                "model":"within_genus_linear_probability_clustered_genus"}
        withN=results["models"][weighted]["with_NPP"]
        without=results["models"][weighted]["without_NPP"]
        results["models"][weighted]["delta_R2_add_NPP"]=(withN["R2_in_sample"]
                                                           -without["R2_in_sample"])
    return results

def run_species_models(species):
    d=prepare_species(species)
    require(len(d)>20000,"photographic species and NPP severely undercovered")
    out={"n_complete_species_before_mixed_genus_filter":len(d),
         "n_full_photo_species":len(species),
         "hue_counts_all":species.source_hue9.value_counts().to_dict(),
         "hue_counts_complete":d.source_hue9.value_counts().to_dict(),
         "photo_status_counts_complete":d.hue_status.value_counts().to_dict(),
         "n_median_NPP_kgC_m2yr":float(d[NPP_NAME].median()),
         "white":species_model_advanced(d,"white"),
         "white_photo_consensus_only":species_model_advanced(
             prepare_species(species,high_photo_confidence=True),"white"),
         "white_NPP_Coverage_100pct":species_model_advanced(
             prepare_species(species,minimum_NPP=1.0),"white"),
         "nine_hues_exploratory":{}}
    for hue in HUES:
        out["nine_hues_exploratory"][hue]=species_model_advanced(
            d,hue,modes=("genus_equal",))
    fitted=[k for k,v in out["nine_hues_exploratory"].items()
            if v["status"]=="fitted"]
    for field in ("NPP","LGM_warming"):
        p=[out["nine_hues_exploratory"][k]["models"]["genus_equal"][
               "with_NPP"]["effects"][field]["p"] for k in fitted]
        q=multipletests(p,method="fdr_bh")[1]
        for k,qq in zip(fitted,q):
            out["nine_hues_exploratory"][k].setdefault("q_BH_nine_hues",{})[
                field]=float(qq)
    require(out["white"]["status"]=="fitted","cannot test white photo colour")
    return out,d

def signed_cell_inputs(native,spatial,ndvi,npp):
    n=signed_parquet(native,NATIVE_SHA,NATIVE_MEMBER,"output_sha256")
    s=signed_parquet(spatial,SPATIAL_SHA,SPATIAL_MEMBER,"output_sha256")
    v=signed_parquet(ndvi,NDVI_SHA,NDVI_MEMBER,"output_SHA256")
    p=signed_parquet(npp,NPP_SHA,NPP_MEMBER,"output_SHA256")
    require(len(p)==len(v)==18577 and len(s)==12560,
            "true original 100km cell source population changed")
    cols=["grid_x","grid_y","block_500km","block_1000km",
          "lat","lon","present_temperature_C","present_precip_mm",
          "threeGCM_delta_temperature_C_median",
          "threeGCM_delta_log_precip_median","ESA_AGB_woody_Mg_ha_2020"]
    d=n.merge(s[cols],on=["grid_x","grid_y"],how="inner",validate="many_to_one")
    d=d.merge(v[["grid_x","grid_y","MODIS_NDVI_2001_2021_mean"]],
              on=["grid_x","grid_y"],validate="many_to_one")
    d=d.merge(p[["grid_x","grid_y",NPP_FIELD,"NPP_valid_years_2001_2023",
                 "NPP_2020_kgC_m2_yr"]],
              on=["grid_x","grid_y"],validate="many_to_one")
    require(not d[["tdwg3_assignment_mode","grid_x","grid_y"]].duplicated().any(),
            "region+cell no longer unique")
    return d

def cell_cohort(d,region_mode,response,min_support=30):
    g=d[d.tdwg3_assignment_mode.eq(region_mode)].copy()
    if response=="native":
        n,w=g.native_species_n,g.native_white_n
    elif response=="introduced":
        n,w=g.introduced_species_n,g.introduced_white_n
    elif response=="all":
        n,w=g.photo_label_species_n,g.white_photo_label_species_n
    else:raise ValueError("unknown photo-colour cohort")
    g["colour_species_n"]=n
    g["white_colour_species_n"]=w
    g["white_photo_species_fraction"]=w/n.replace(0,np.nan)
    needed=["white_photo_species_fraction",
            "threeGCM_delta_temperature_C_median",
            "threeGCM_delta_log_precip_median",
            "present_temperature_C","present_precip_mm",
            "MODIS_NDVI_2001_2021_mean",
            "ESA_AGB_woody_Mg_ha_2020",NPP_FIELD,
            "GBIF_occurrence_count","lat"]
    mask=(g.colour_species_n.ge(min_support)
          & g[needed].notna().all(axis=1)
          & g.NPP_valid_years_2001_2023.ge(18)
          & g.present_precip_mm.ge(0)
          & g.ESA_AGB_woody_Mg_ha_2020.ge(0))
    return g.loc[mask].copy()

def cell_white_model(d,block="500",weight=False,npp=True):
    require(len(d)>350 and d.block_500km.nunique()>35,
            "insufficient independent GBIF 100km flower-colour cells")
    X=pd.DataFrame({
       "NPP":d[NPP_FIELD].astype(float),
       "LGM_warming":d[LGM_TEMP].astype(float),
       "LGM_rain":d[LGM_RAIN].astype(float),
       "current_temp":d.present_temperature_C.astype(float),
       "log_current_precip":np.log1p(d.present_precip_mm.astype(float)),
       "log_woody_AGB":np.log1p(d.ESA_AGB_woody_Mg_ha_2020.astype(float)),
       "NDVI":d.MODIS_NDVI_2001_2021_mean.astype(float),
       "absolute_lat":d.lat.abs().astype(float),
       "log_GBIF_records":np.log1p(d.GBIF_occurrence_count.astype(float)),
       "log_species_count":np.log1p(d.colour_species_n.astype(float))
    },index=d.index)
    require(X.notna().all().all(),"missing in cell model features")
    standard=X.std()
    require((standard>1e-9).all(),"no variation in geographic exposure features")
    X=(X-X.mean())/standard
    Y=d.white_photo_species_fraction.astype(float)
    B=d.block_500km.astype(str) if block=="500" else d.block_1000km.astype(str)
    X=X-X.groupby(B).transform("mean")
    Y=Y-Y.groupby(B).transform("mean")
    good=X.abs().sum(axis=1)>0
    X=X.loc[good];Y=Y.loc[good]
    cluster=d.loc[good,"block_500km"].astype(str)
    if not npp:X=X.drop(columns=["NPP"])
    weights=d.loc[good,"colour_species_n"].astype(float)
    if weight:
        model=sm.WLS(Y,X,weights=weights/weights.median()).fit(
            cov_type="cluster",cov_kwds={"groups":cluster},use_t=False)
    else:
        model=sm.OLS(Y,X).fit(cov_type="cluster",
                             cov_kwds={"groups":cluster},use_t=False)
    eff={}
    for k in ("NPP","LGM_warming","LGM_rain","log_woody_AGB","NDVI"):
        if k not in X:continue
        lo,hi=model.conf_int().loc[k].tolist()
        eff[k]={"beta_fraction_per_SD":float(model.params[k]),
                "95CI":[float(lo),float(hi)],
                "p":float(model.pvalues[k])}
    return {"n_cells":len(Y),"n_spatial_blocks_500km":int(cluster.nunique()),
            "effects":eff,"R2_in_sample":float(model.rsquared),
            "outcome":"GBIF_species_share_with_photo9_white_proxy",
            "block_km":int(block),"weights":"n_photo_label_species" if weight else "equal_100km_cells",
            "no_local_flower_colour_measured":True}

def run_grid_models(d):
    out={}
    tables=[]
    for reg in ("strict","moderate"):
        out[reg]={}
        for label in ("native","introduced","all"):
            p=cell_cohort(d,reg,label)
            entry={"n_eligible_100km_cells":len(p),"models":{}}
            if len(p)>350 and p.block_500km.nunique()>35:
                entry["models"]["500km_full"]=cell_white_model(p)
                entry["models"]["500km_no_NPP"]=cell_white_model(p,npp=False)
                entry["models"]["1000km_full"]=cell_white_model(p,block="1000")
                entry["models"]["500km_species_support_weighted"]=cell_white_model(p,weight=True)
            else:entry["status"]="insufficient_region_cells"
            out[reg][label]=entry
            p["region_mode"]=reg
            p["source_group"]=label
            tables.append(p[["grid_x","grid_y","region_mode","source_group",
                "colour_species_n","white_photo_species_fraction",NPP_FIELD,
                LGM_TEMP,"block_500km"]].copy())
            print("MOIS_FLOWER_COLOUR_GBIF_NATIVE_NPP_LGM",json.dumps({
               "region_mode":reg,"source_group":label,"cells":len(p),
               "beta_NPP":entry["models"].get("500km_full",{}).get(
                   "effects",{}).get("NPP",{}).get("beta_fraction_per_SD"),
               "p_NPP":entry["models"].get("500km_full",{}).get(
                   "effects",{}).get("NPP",{}).get("p"),
               "beta_LGM":entry["models"].get("500km_full",{}).get(
                   "effects",{}).get("LGM_warming",{}).get("beta_fraction_per_SD"),
               "p_LGM":entry["models"].get("500km_full",{}).get(
                   "effects",{}).get("LGM_warming",{}).get("p")
            },sort_keys=True),flush=True)
    return out,pd.concat(tables,ignore_index=True)

def run(gbif,env,lgm,npp,native,spatial,ndvi,out):
    out=Path(out)
    require(not out.exists(),"never overwrite true-source floral-colour report")
    p,ignored=read_signed_sources(gbif,env)
    hist=load_species_history(lgm)
    v=signed_parquet(npp,NPP_SHA,NPP_MEMBER,"output_SHA256")
    full=make_species_NPP(p,v,hist)
    main,complete=run_species_models(full)
    cell=signed_cell_inputs(native,spatial,ndvi,npp)
    cellresult,celltab=run_grid_models(cell)
    out.mkdir(parents=True)
    f1=out/"GBIF37104_photo9_species_NPP23years_LGM_AGB_NDVI.parquet"
    f2=out/"GBIF100km_native_introduced_photo_white_NPP_LGM.parquet"
    f3=out/"flower_colour_NPP_LGM_models.json"
    full.to_parquet(f1,index=False)
    celltab.to_parquet(f2,index=False)
    result={
       "source_species_level":main,
       "source_current_native_vs_introduced_cells":cellresult,
       "n_accepted_species":37104,
       "n_source_GBIF_species_100km_pairs":5493931,
       "source_white_is_photo_species_proxy_not_local_petal_measurement":True,
       "genetic_fixed_flower_colour_established":False,
       "historical_LGM_refugia_established":False,
       "WCVP_native_means_only_broad_TDWG3":True,
       "NPP_kgC_m2_yr_not_flower_biomass":True,
       "model_500km_cells_are_not_independent_GBIF_species_samples":True,
       "sources_sha256":{"GBIF":sha(gbif),"ENV":sha(env),"LGM":sha(lgm),
          "NPP":sha(npp),"WCVP_native":sha(native),
          "spatial_history":sha(spatial),"NDVI":sha(ndvi)}
    }
    f3.write_text(json.dumps(result,indent=2,sort_keys=True,default=str)+"\n")
    receipt={
      "source_37k_species":37104,"source_5493931_GBIF_species_cells":5493931,
      "primary_goal_is_flower_colour":True,
      "source_HUE9_is_observational_photo_proxy":True,
      "species_level_one_row_per_species":True,
      "grid_level_one_outcome_per_GBIF_cell":True,
      "native_intro_is_broad_region_not_individual_native":True,
      "NPP_annual_flux_not_floral_biomass":True,
      "LGM_climate_history_not_refugia_reconstruction":True,
      "source_zip_SHA256":result["sources_sha256"],
      "n_complete_species":len(complete),
      "n_eligible_native_strict_cells":cellresult["strict"]["native"]["n_eligible_100km_cells"],
      "results":result,
      "output_SHA256":{i.name:sha(i) for i in (f1,f2,f3)}
    }
    (out/"receipt.json").write_text(json.dumps(receipt,indent=2,sort_keys=True,default=str)+"\n")
    w=main["white"]
    print("MOIS_PRIMARY_FLOWER_COLOUR_NPP_LGM_RESULTS",json.dumps({
        "n_complete_species":len(complete),
        "white_n_species":w.get("n_species"),
        "white_n_genera":w.get("n_genera"),
        "species_equal":w.get("models",{}).get("species_equal",{}).get("with_NPP",{}).get("effects"),
        "genus_equal":w.get("models",{}).get("genus_equal",{}).get("with_NPP",{}).get("effects"),
        "n_strict_native_cells":cellresult["strict"]["native"]["n_eligible_100km_cells"]},
        sort_keys=True),flush=True)
    return receipt

if __name__=="__main__":
    args=argparse.ArgumentParser()
    for k in ("gbif","env","lgm","npp","native","spatial","ndvi","out"):
        args.add_argument("--"+k,required=True)
    a=args.parse_args()
    run(a.gbif,a.env,a.lgm,a.npp,a.native,a.spatial,a.ndvi,a.out)
