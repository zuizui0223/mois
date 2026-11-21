#!/usr/bin/env python3
"""GBIF/WCVP nonrandom genus segregation versus current woody AGB and NDVI.

Response is observed minus EXACT family-conditioned genus Jaccard expectation.
NO photo hue is used in constructing or regressing this taxonomic response.
ESA CCI woody AGB (Mg/ha) is stock, and MODIS NDVI is greenness: neither
equals NPP, floral allocation or historical biomass. This is observational.
"""
from __future__ import annotations
import argparse,io,json
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scripts.analysis.mois_lgm_climate_history_species_20261116 import sha,require

FAMILY_SHA="0259071f76ffad271e5113b792c7c21f27c0c1d10ef99a215ff3a7779ce6a27b"
SPATIAL_SHA="e209257adbc185b4855f3b8d76aeaa50d2e2afe8a1feeb0f76013366e57c4d16"
NATIVE_SHA="b72b7358ed5382d625583ad61b90952e4e77377296fe90ae677bf9481aa80a7f"
NDVI_SHA="e08e81d1c2afd1a2ee5312f3ebc4aeb5ab7678988c2453a2a9c63d28c82eb17a"
FAMILY_MEMBER="GBIF100km_WCVP_family_conditioned_native_intro_genus_null.parquet"
SPATIAL_MEMBER="LGM_white_photo9_species_assemblages_100km_cells.parquet"
NATIVE_MEMBER="WCVP_TDWG3_native_introduced_photo9_LGM_100km_cell_summary.parquet"
NDVI_MEMBER="WorldClim19_weather_NDVI_18577_global_100km_cells.parquet"

def signed_parquet(archive,zip_sha,member,meta_sha_key):
    import pyarrow.parquet as pq
    require(sha(archive)==zip_sha, f"signed source ZIP changed: {member}")
    with ZipFile(archive) as z:
        b=z.read(member)
        receipt=json.loads(z.read("receipt.json"))
    require(sha(b)==receipt[meta_sha_key][member],
            "signed internal Parquet digest changed "+member)
    return pq.read_table(io.BytesIO(b)).to_pandas()

def build_comparison(family,spatial,native,ndvi):
    require(len(family)>6500 and len(spatial)==12560
            and len(ndvi)==18577, "source population sizes not aligned")
    require(not family[["grid_x","grid_y","region_mode"]].duplicated().any(),
            "family conditional double-counted source grid")
    keep=["grid_x","grid_y","block_500km","block_1000km","lon",
          "present_temperature_C","present_precip_mm",
          "threeGCM_delta_log_precip_median","ESA_AGB_woody_Mg_ha_2020"]
    d=family.merge(spatial[keep],on=["grid_x","grid_y"],validate="many_to_one")
    ni=ndvi[["grid_x","grid_y","MODIS_NDVI_2001_2021_mean",
              "MODIS_NDVI_2001_2021_sd"]].copy()
    require(not ni[["grid_x","grid_y"]].duplicated().any(),
            "NDVI global grid should have 18577 unique cells")
    d=d.merge(ni,on=["grid_x","grid_y"],validate="many_to_one")
    n=native[["grid_x","grid_y","tdwg3_assignment_mode",
              "GBIF_occurrence_count","native_species_n",
              "introduced_species_n"]].rename(
                  columns={"tdwg3_assignment_mode":"region_mode"}).copy()
    require(not n[["grid_x","grid_y","region_mode"]].duplicated().any(),
            "WCVP region source has duplicate cells")
    d=d.merge(n,on=["grid_x","grid_y","region_mode"],validate="one_to_one")
    require(len(d)==len(family),"lost family-null GBIF cells")
    require((d.native_species_n>=30).all() &
            (d.introduced_species_n>=30).all(),
            "original native intro matched sample changed")
    # Exclude missing satellite climate, not equal to zero woody biomass.
    required=["residual_after_family_conditioning","observed_jaccard",
              "family_conditioned_expected_jaccard",
              "threeGCM_delta_temperature_C_median",
              "threeGCM_delta_log_precip_median",
              "ESA_AGB_woody_Mg_ha_2020","MODIS_NDVI_2001_2021_mean",
              "present_temperature_C","present_precip_mm",
              "GBIF_occurrence_count","lat","lon"]
    good=d[required].notna().all(axis=1)
    good &= d.ESA_AGB_woody_Mg_ha_2020.ge(0) & d.present_precip_mm.ge(0)
    d["complete_woodyNDVI"]=good
    return d

def normalized_model_matrix(d):
    x=pd.DataFrame({
        "LGM_warming":d.threeGCM_delta_temperature_C_median.astype(float),
        "LGM_rainfall_change":d.threeGCM_delta_log_precip_median.astype(float),
        "current_temperature":d.present_temperature_C.astype(float),
        "current_log_precip":np.log1p(d.present_precip_mm.astype(float)),
        "log_woody_AGB":np.log1p(d.ESA_AGB_woody_Mg_ha_2020.astype(float)),
        "NDVI_mean":d.MODIS_NDVI_2001_2021_mean.astype(float),
        "mean_abs_lat":d.lat.abs().astype(float),
        "log_GBIF_records":np.log1p(d.GBIF_occurrence_count.astype(float)),
        "log_native_species_n":np.log1p(d.native_species_n.astype(float)),
        "log_intro_species_n":np.log1p(d.introduced_species_n.astype(float)),
        "log_family_n":np.log1p(d.family_union.astype(float)),
        "log_genus_n":np.log1p(d.genus_union.astype(float))
    },index=d.index)
    require(x.notna().all().all(),"model input contains missing")
    sd=x.std()
    require((sd>1e-10).all(),"nonvarying model predictor")
    return (x-x.mean())/sd,{k:float(v) for k,v in sd.items()}

def fit(d,response="residual_after_family_conditioning",block="500",
        include_agb=True,include_ndvi=True,weighted=False):
    require(len(d)>400 and d.block_500km.nunique()>30,
            "too few independent clusters for family-conditional test")
    x,std=normalized_model_matrix(d)
    if not include_agb:
        x=x.drop(columns=["log_woody_AGB"])
    if not include_ndvi:
        x=x.drop(columns=["NDVI_mean"])
    y=d[response].astype(float)
    b=d.block_500km.astype(str) if block=="500" else d.block_1000km.astype(str)
    x=x-x.groupby(b).transform("mean")
    y=y-y.groupby(b).transform("mean")
    useful=x.abs().sum(axis=1)>0
    x=x.loc[useful];y=y.loc[useful]
    groups=d.loc[useful,"block_500km"].astype(str)
    if weighted:
        support=d.loc[useful,["native_species_n",
                              "introduced_species_n"]].min(axis=1).astype(float)
        m=sm.WLS(y,x,weights=support/support.median()).fit(
            cov_type="cluster",cov_kwds={"groups":groups},use_t=False)
    else:
        m=sm.OLS(y,x).fit(cov_type="cluster",
            cov_kwds={"groups":groups},use_t=False)
    effect={}
    for name in ("LGM_warming","log_woody_AGB","NDVI_mean"):
        if name in x:
            lo,hi=m.conf_int().loc[name].tolist()
            effect[name]={"beta_per_global_SD":float(m.params[name]),
                          "95CI":[float(lo),float(hi)],
                          "p":float(m.pvalues[name])}
    return {
        "n_cells":len(y),"n_500km_clusters":int(groups.nunique()),
        "response":response,"block_km":int(block),
        "weights":"min_native_intro_species_count" if weighted else "equal_cells",
        "AGB_present":bool(include_agb),"NDVI_present":bool(include_ndvi),
        "R2_descriptive_not_holdout":float(m.rsquared),"effect":effect,
        "raw_global_SD":{"LGM_warming_C":std["LGM_warming"],
                         "log_woody_AGB":std["log_woody_AGB"],
                         "NDVI_mean":std["NDVI_mean"]},
        "woody_AGB_is_NPP_or_flower_biomass":False,
        "observational_500km_cluster_not_causal":True
    }

def run(familyzip,spatialzip,nativezip,ndvizip,out):
    out=Path(out)
    require(not out.exists(),"refuse overwrite previous source-verified results")
    f=signed_parquet(familyzip,FAMILY_SHA,FAMILY_MEMBER,"output_SHA256")
    s=signed_parquet(spatialzip,SPATIAL_SHA,SPATIAL_MEMBER,"output_sha256")
    n=signed_parquet(nativezip,NATIVE_SHA,NATIVE_MEMBER,"output_sha256")
    v=signed_parquet(ndvizip,NDVI_SHA,NDVI_MEMBER,"output_SHA256")
    allrows=build_comparison(f,s,n,v)
    result={"strict":{},"moderate":{}}
    snapshots=[]
    for mode in ("strict","moderate"):
        d=allrows[allrows.region_mode.eq(mode)].copy()
        complete=d[d.complete_woodyNDVI].copy()
        require(len(d)==({"strict":3106,"moderate":3602}[mode])
                and len(complete)>2500,"source cohort or NDVI coverage drift")
        models={}
        for response in ("residual_after_family_conditioning",
                         "observed_jaccard"):
            models[response]={
                "500km_full_AGB_and_NDVI":fit(complete,response,"500"),
                "500km_no_AGB":fit(complete,response,"500",False,True),
                "500km_no_NDVI":fit(complete,response,"500",True,False),
                "500km_both_excluded":fit(complete,response,"500",False,False),
                "1000km_full_AGB_and_NDVI":fit(complete,response,"1000"),
                "500km_support_weighted":fit(complete,response,"500",True,True,True)
            }
        result[mode]={"n_raw":len(d),"n_complete":len(complete),
            "mean_family_excess":float(
                complete.residual_after_family_conditioning.mean()),
            "AGB2020_Mg_ha_median":float(complete.ESA_AGB_woody_Mg_ha_2020.median()),
            "NDVI_mean_median":float(complete.MODIS_NDVI_2001_2021_mean.median()),
            "AGB_NDVI_Spearman":float(
                complete.ESA_AGB_woody_Mg_ha_2020.corr(
                    complete.MODIS_NDVI_2001_2021_mean,method="spearman")),
            "models":models}
        snapshots.append(complete[["grid_x","grid_y","region_mode",
            "residual_after_family_conditioning",
            "observed_jaccard","family_conditioned_expected_jaccard",
            "ESA_AGB_woody_Mg_ha_2020","MODIS_NDVI_2001_2021_mean",
            "threeGCM_delta_temperature_C_median",
            "block_500km"]].copy())
        v0=models["residual_after_family_conditioning"]["500km_full_AGB_and_NDVI"]
        print("MOIS_GBIF_FAMILY_NULL_BIOMASS_NDVI",json.dumps({
            "mode":mode,"n":len(complete),
            "AGB_beta":v0["effect"]["log_woody_AGB"]["beta_per_global_SD"],
            "AGB_p":v0["effect"]["log_woody_AGB"]["p"],
            "NDVI_beta":v0["effect"]["NDVI_mean"]["beta_per_global_SD"],
            "NDVI_p":v0["effect"]["NDVI_mean"]["p"],
            "LGM_beta":v0["effect"]["LGM_warming"]["beta_per_global_SD"],
            "LGM_p":v0["effect"]["LGM_warming"]["p"],
            "AGB_NDVI_rho":result[mode]["AGB_NDVI_Spearman"]},
            sort_keys=True),flush=True)
    out.mkdir(parents=True)
    fp=out/"GBIF100km_family_null_genus_segregation_AGB_NDVI_samecells.parquet"
    pd.concat(snapshots,ignore_index=True).to_parquet(fp,index=False)
    modelpath=out/"models.json"
    modelpath.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    receipt={"source_GBIF_accepted_species":37104,
      "original_GBIF_species_100km_pairs":5493931,
      "taxon_response_without_photo_hue":True,
      "woody_AGB_not_NPP_or_individual_flower_allocation":True,
      "NDVI_not_biomass":True,
      "native_status_broad_WCVP_TDWG3_only":True,
      "LGM_present_cell_climate_not_refugia":True,
      "sources_zip_sha256":{
         "family_genus_null":FAMILY_SHA,
         "spatial_LGM":SPATIAL_SHA,
         "native_WCVP":NATIVE_SHA,
         "NDVI":NDVI_SHA},
      "models":result,
      "output_sha256":{fp.name:sha(fp),modelpath.name:sha(modelpath)}}
    (out/"receipt.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    return receipt

if __name__=="__main__":
    p=argparse.ArgumentParser()
    for a in ("family","spatial","native","ndvi","output"):
        p.add_argument("--"+a,required=True)
    a=p.parse_args()
    run(a.family,a.spatial,a.native,a.ndvi,a.output)
