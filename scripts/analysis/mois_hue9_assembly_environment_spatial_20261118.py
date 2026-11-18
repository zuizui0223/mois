#!/usr/bin/env python3
"""MOIS ecological drivers of genus × occupancy-breadth residual photo-hue entropy.

One original accepted species × 100-km GBIF presence contributes one count
per grid; 199 global species-label permutations within each genus × range-size
bin define the independent null. We analyze observed minus permutation-
expected 9-class Shannon entropy at those SAME cells.

This study estimates descriptive environmental gradients, not natural
selection, flower pigmentation, actual petal colour or pollinator mechanism.

Real signed source GBIF x 9 hues: 37,104 accepted taxa, 5,493,931 species-
cell pairs, 18,577 cells (36,758 taxa with source grid records). Current
GBIF, WorldClim, MOD17 2001-2023, MODIS NDVI, ESA woody AGB and LGM grids
are not contemporaneous direct organism-level phenotypes.
"""
from __future__ import annotations
import argparse,hashlib,io,json
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge
from sklearn.metrics import r2_score,mean_squared_error

SOURCES={
 "null":{"sha":"596166fa2fe8bfb28653c325a3a8f1a1495c0792688a7fdb781014cfedaa88c7"},
 "spatial":{"sha":"e209257adbc185b4855f3b8d76aeaa50d2e2afe8a1feeb0f76013366e57c4d16",
            "member":"LGM_white_photo9_species_assemblages_100km_cells.parquet",
            "digest_key":"output_sha256","rows":12560},
 "NDVI":{"sha":"e08e81d1c2afd1a2ee5312f3ebc4aeb5ab7678988c2453a2a9c63d28c82eb17a",
         "member":"WorldClim19_weather_NDVI_18577_global_100km_cells.parquet",
         "digest_key":"output_SHA256","rows":18577},
 "NPP":{"sha":"68b2cfa1aeac9dc5d02cd310c9b2f45553a5be4db3319a1e40897612ad9b1bef",
        "member":"MOD17_C061_NPP_2001_2023_18577_true_GBIF100km_cells.parquet",
        "digest_key":"output_SHA256","rows":18577}
}
ROOT="mois_global_9hue_community_assembly_breadth_20261118"
GENUS="mois_global_9hue_community_assembly_20261118"
GRID="genuine_global_100km_photo9_diversity_vs_genus_null.csv.gz"
NPP="NPP_2001_2023_mean_kgC_m2_yr"
NPP_YEARS="NPP_valid_years_2001_2023"
LGM_TEMP="threeGCM_delta_temperature_C_median"
LGM_RAIN="threeGCM_delta_log_precip_median"
CLIMATE=("temperature","log_precip","LGM_temp_delta","LGM_precip_delta")
VEGETATION=("NPP","NDVI","log_woody_AGB")
ECOLOGY=CLIMATE+VEGETATION
BASE=("log_species_richness","expected_null_Shannon",
      "log_GBIF_records","abs_latitude","longitude_sine","longitude_cosine")
SEED=20261118

def require(cond,msg):
    if not bool(cond): raise ValueError("MOIS_HUE9_ECOLOGY_DRIVERS: "+str(msg))

def sha(x):
    if isinstance(x,(Path,str)):x=Path(x).read_bytes()
    return hashlib.sha256(x).hexdigest()

def table(zip_path,name):
    import pyarrow.parquet as pq
    spec=SOURCES[name]
    require(sha(zip_path)==spec["sha"],name+" outer original archive SHA changed")
    with ZipFile(zip_path) as z:
        rec=json.loads(z.read("receipt.json"))
        b=z.read(spec["member"])
    require(sha(b)==rec[spec["digest_key"]][spec["member"]],
            name+" original internal parquet SHA differs")
    d=pq.read_table(io.BytesIO(b)).to_pandas()
    require(len(d)==spec["rows"] and
            not d[["grid_x","grid_y"]].duplicated().any(),
            name+" grid keys or source denominator changed")
    return d

def null_tables(zip_path):
    require(sha(zip_path)==SOURCES["null"]["sha"],
            "Signed two-null global hue9 archive altered")
    out={}
    with ZipFile(zip_path) as z:
        for root,key in ((ROOT,"genus_breadth"),(GENUS,"genus")):
            receipt=json.loads(z.read(root+"/receipt.json"))
            require(receipt["null_model"]==key
                    and receipt["source_species"]==37104
                    and receipt["source_species_grid_pairs"]==5493931
                    and receipt["num_global_grid_cells"]==18577
                    and receipt["genus_preserving_global_species_hue_label_permutations"]==199,
                    "Signed source photo9 permuted-null model metadata drifted")
            for fn,digest in receipt["output_sha256"].items():
                require(sha(z.read(root+"/"+fn))==digest,
                        "Photographic hue9 null internal numeric artifact changed")
            d=pd.read_csv(io.BytesIO(z.read(root+"/"+GRID)),compression="gzip")
            require(len(d)==18577 and not d[["grid_x","grid_y"]].duplicated().any(),
                    "Source colour-grid uniqueness changed")
            out[key]=d
    a=out["genus_breadth"];b=out["genus"]
    key=["grid_x","grid_y"]
    require(a[key].equals(b[key]),"Strict and simple genus grid row order differs")
    require(np.allclose(a.observed_hue9_Shannon,b.observed_hue9_Shannon,atol=1e-12),
            "Actual observed floral photo hue Shannon changed between nulls")
    require((a.n_accepted_species==b.n_accepted_species).all(),
            "Actual richness changed between nulls")
    return a,b

def join_sources(null,genus,spatial,ndvi,npp):
    original_n=len(null)
    need_spatial=["grid_x","grid_y","block_500km","block_1000km",
                  "lat","lon","present_temperature_C","present_precip_mm",
                  "ESA_AGB_woody_Mg_ha_2020",
                  LGM_TEMP,LGM_RAIN,"GBIF_occurrence_count"]
    need_ndvi=["grid_x","grid_y","MODIS_NDVI_2001_2021_mean"]
    need_npp=["grid_x","grid_y",NPP,NPP_YEARS]
    for cols,df,name in ((need_spatial,spatial,"SPATIAL"),
                         (need_ndvi,ndvi,"NDVI"),(need_npp,npp,"NPP")):
        require(set(cols).issubset(df),name+" required columns unavailable")
    extra=genus[["grid_x","grid_y","observed_minus_genus_null_Shannon"]].rename(
        columns={"observed_minus_genus_null_Shannon":"genus_only_Shannon_anomaly"})
    d=null.merge(extra,on=["grid_x","grid_y"],validate="one_to_one")
    d=d.merge(spatial[need_spatial],on=["grid_x","grid_y"],how="left",
              validate="one_to_one")
    d=d.merge(ndvi[need_ndvi],on=["grid_x","grid_y"],how="left",
              validate="one_to_one")
    d=d.merge(npp[need_npp],on=["grid_x","grid_y"],how="left",
              validate="one_to_one")
    require(len(d)==original_n and not d[["grid_x","grid_y"]].duplicated().any(),
            "True source GBIF grid cells duplicated by environment join")
    fields=["observed_minus_genus_null_Shannon",
            "expected_hue9_Shannon_given_genus",
            "genus_only_Shannon_anomaly",
            "present_temperature_C","present_precip_mm",
            "MODIS_NDVI_2001_2021_mean",
            "ESA_AGB_woody_Mg_ha_2020",NPP,LGM_TEMP,LGM_RAIN,
            "lat","lon","GBIF_occurrence_count",
            "block_500km","block_1000km"]
    raw={
      "source_genuine_grid_cells":len(d),
      "richness_30plus_global":int(d.n_accepted_species.ge(30).sum()),
      "richness_100plus_global":int(d.n_accepted_species.ge(100).sum()),
      "climate_block_grid_rows":int(d.block_500km.notna().sum()),
      "npp_18plus_years_grid_rows":int(d[NPP_YEARS].ge(18).sum())
    }
    ok=d.n_accepted_species.ge(30)&d[fields].notna().all(axis=1)
    ok&=d[NPP_YEARS].ge(18)&d.present_precip_mm.ge(0)
    ok&=d.ESA_AGB_woody_Mg_ha_2020.ge(0)
    ok&=d.GBIF_occurrence_count.ge(0)
    result=d.loc[ok].copy().reset_index(drop=True)
    require(len(result)>=7000
            and result.block_500km.nunique()>500
            and result.block_1000km.nunique()>150,
            "Global shared NPP+AGB+strict null sample too small")
    result["white_species_fraction"]=(result.n_white_photo_species/
                                      result.n_accepted_species)
    result["entropy_anomaly"]=result.observed_minus_genus_null_Shannon
    result["genus_only_anomaly"]=result.genus_only_Shannon_anomaly
    raw.update({
      "strict_complete_30plus_cells":len(result),
      "strict_complete_100plus_cells":int(result.n_accepted_species.ge(100).sum()),
      "valid_500km_blocks":int(result.block_500km.nunique()),
      "valid_1000km_blocks":int(result.block_1000km.nunique()),
      "global_mean_entropy_anomaly_in_complete_cells":
          float(result.entropy_anomaly.mean()),
      "white_fraction_median":float(result.white_species_fraction.median()),
    })
    return result,raw

def matrix(d,include_white=False):
    lat=d.lat.astype(float);lon=d.lon.astype(float)
    cov=pd.DataFrame({
      "temperature":d.present_temperature_C.astype(float),
      "log_precip":np.log1p(d.present_precip_mm.astype(float)),
      "LGM_temp_delta":d[LGM_TEMP].astype(float),
      "LGM_precip_delta":d[LGM_RAIN].astype(float),
      "NPP":d[NPP].astype(float),
      "NDVI":d.MODIS_NDVI_2001_2021_mean.astype(float),
      "log_woody_AGB":np.log1p(d.ESA_AGB_woody_Mg_ha_2020.astype(float)),
      "log_species_richness":np.log1p(d.n_accepted_species.astype(float)),
      "expected_null_Shannon":d.expected_hue9_Shannon_given_genus.astype(float),
      "log_GBIF_records":np.log1p(d.GBIF_occurrence_count.astype(float)),
      "abs_latitude":lat.abs(),
      "longitude_sine":np.sin(np.deg2rad(lon)),
      "longitude_cosine":np.cos(np.deg2rad(lon))
    },index=d.index)
    if include_white:
        cov["white_species_fraction"]=d.white_species_fraction.astype(float)
    require(np.isfinite(cov.to_numpy()).all(),"Nonfinite ecological predictors")
    sd=cov.std(ddof=1)
    require((sd>1e-10).all(),"No predictor variation")
    return (cov-cov.mean())/sd,{k:{"mean":float(cov[k].mean()),
                                   "SD":float(sd[k])} for k in cov.columns}

def spatial_fit(d,response="entropy_anomaly",block_km=500,
                include_white=False):
    require(len(d)>2500,"Insufficient actual global climatic assemblages")
    X,scales=matrix(d,include_white)
    key="block_500km" if block_km==500 else "block_1000km"
    b=d[key].astype(str)
    xc=X-X.groupby(b).transform("mean")
    y=d[response].astype(float)
    yc=y-y.groupby(b).transform("mean")
    keep=(xc.abs().sum(axis=1)>1e-11)
    xc=xc.loc[keep]
    yc=yc.loc[keep]
    block=b.loc[keep]
    require(len(yc)>2000 and block.nunique()>100,
            "Too few independent spatial blocks")
    rank=int(np.linalg.matrix_rank(xc.to_numpy()))
    require(rank==len(xc.columns),"Severe collinearity / singular fixed effect matrix")
    fit=sm.OLS(yc,xc).fit(cov_type="cluster",cov_kwds={"groups":block},
                          use_t=False)
    out={}
    for k in X.columns:
        ci=fit.conf_int().loc[k]
        out[k]={
            "beta_Shannon_per_global_SD":float(fit.params[k]),
            "CI95":[float(ci.iloc[0]),float(ci.iloc[1])],
            "p_cluster":float(fit.pvalues[k])
        }
    q=multipletests([out[k]["p_cluster"] for k in ECOLOGY],
                    method="fdr_bh")[1]
    for k,v in zip(ECOLOGY,q):
        out[k]["q_BH_seven_ecological_exposures"]=float(v)
    # Joint group Wald tests with exact same rows and covariates.
    wald={}
    for label,grp in (("current_temp_precip",("temperature","log_precip")),
                      ("LGM_paleoclimate",("LGM_temp_delta","LGM_precip_delta")),
                      ("productivity_vegetation",VEGETATION)):
        R=np.zeros((len(grp),len(fit.params)))
        for i,k in enumerate(grp):R[i,fit.params.index.get_loc(k)]=1
        stat=fit.wald_test(R,scalar=True)
        wald[label]={"chi_square":float(stat.statistic),
                     "df":len(grp),"p_cluster":float(stat.pvalue)}
    result={
        "response":response,"spatial_FE_km":block_km,
        "n_source_cells":len(d),"n_fit_cells":len(yc),
        "n_spatial_clusters":int(block.nunique()),
        "n_covariates":len(xc.columns),
        "R2_within_block_in_sample":float(fit.rsquared),
        "scales":scales,"coefficients":out,"joint_Wald":wald,
        "not_causal_environmental_effect":True
    }
    return result

def spatial_oof(d):
    X,_=matrix(d)
    y=d.entropy_anomaly.to_numpy(float)
    blocks=d.block_1000km.astype(str).to_numpy()
    variants={
      "geographic_richness_null_baseline":list(BASE),
      "baseline_plus_current_and_LGM_climate":list(BASE+CLIMATE),
      "baseline_plus_vegetation_only":list(BASE+VEGETATION),
      "baseline_plus_all_environment":list(BASE+ECOLOGY)
    }
    allpred={name:np.full(len(d),np.nan) for name in variants}
    groups=GroupKFold(n_splits=5)
    foldrows=[]
    for i,(tr,te) in enumerate(groups.split(X,y,blocks)):
        require(set(blocks[tr]).isdisjoint(set(blocks[te])),
                "1000km spatial holdout leaked a block")
        foldrecord={"fold":i,"n_test_cells":len(te),
                    "n_holdout_1000km_blocks":len(set(blocks[te]))}
        for name,columns in variants.items():
            train=X.iloc[tr][columns].to_numpy(float)
            val=X.iloc[te][columns].to_numpy(float)
            # All centering, scaling, regularisation learned from training only.
            scaler=StandardScaler().fit(train)
            reg=Ridge(alpha=100)
            reg.fit(scaler.transform(train),y[tr])
            yp=reg.predict(scaler.transform(val))
            allpred[name][te]=yp
            foldrecord[name+"_mse"]=float(mean_squared_error(y[te],yp))
        foldrows.append(foldrecord)
    require(all(np.isfinite(pred).all() for pred in allpred.values()),
            "Spatial 1000km OOF predictions incomplete")
    model_metrics={
      name:{
        "R2_out_of_fold":float(r2_score(y,pred)),
        "MSE":float(mean_squared_error(y,pred)),
        "RMSE":float(np.sqrt(mean_squared_error(y,pred))),
        "MAE":float(np.mean(np.abs(y-pred))),
      } for name,pred in allpred.items()}
    base=allpred["geographic_richness_null_baseline"]
    ids,gidx=np.unique(blocks,return_inverse=True)
    sizes=np.bincount(gidx,minlength=len(ids))
    rng=np.random.default_rng(SEED+700)
    changes=[]
    for name,pred in allpred.items():
        if name=="geographic_richness_null_baseline":continue
        delta=(y-pred)**2-(y-base)**2
        sums=np.bincount(gidx,weights=delta,minlength=len(ids))
        boot=[]
        for _ in range(2000):
            counts=np.bincount(rng.integers(0,len(ids),size=len(ids)),
                               minlength=len(ids))
            boot.append(float(np.dot(counts,sums)/np.dot(counts,sizes)))
        changes.append({
          "model_vs_baseline":name,
          "MSE_increment_added_minus_baseline":float(delta.mean()),
          "CI95_1000km_block_bootstrap":
             [float(v) for v in np.percentile(boot,[2.5,97.5])],
          "n_1000km_blocks":int(len(ids)),
          "n_bootstrap":2000
        })
    return {"n_cells":len(y),"n_holdout_1000km_blocks":int(len(ids)),
            "n_folds":5,"models":model_metrics,
            "paired_delta_MSE":changes,
            "folds":foldrows,
            "leave_one_1000km_block_group_out_not_contiguous_realm_holdout":True}

def run(args):
    root=Path(args.out)
    require(not root.exists(),"Cannot overwrite prior source-backed ecological run")
    strict,genus=null_tables(args.null)
    spatial=table(args.spatial,"spatial")
    ndvi=table(args.NDVI,"NDVI")
    npp=table(args.NPP,"NPP")
    d,coverage=join_sources(strict,genus,spatial,ndvi,npp)
    primary=spatial_fit(d)
    richer=spatial_fit(d[d.n_accepted_species.ge(100)].reset_index(drop=True))
    broad=spatial_fit(d,block_km=1000)
    white_control=spatial_fit(d,include_white=True)
    genus_only=spatial_fit(d,response="genus_only_anomaly")
    cv=spatial_oof(d)
    root.mkdir(parents=True)
    filegrid=root/"all_current_GBIF_100km_cells_genus_breadth_entropy_ecology.csv.gz"
    fields=["grid_x","grid_y","n_accepted_species","lat","lon",
            "block_500km","block_1000km",
            "observed_hue9_Shannon","expected_hue9_Shannon_given_genus",
            "entropy_anomaly","genus_only_anomaly","SES_genus_conditioned",
            "white_species_fraction",NPP,NPP_YEARS,
            "MODIS_NDVI_2001_2021_mean","ESA_AGB_woody_Mg_ha_2020",
            "present_temperature_C","present_precip_mm",LGM_TEMP,LGM_RAIN]
    d[fields].to_csv(filegrid,index=False,compression="gzip")
    outputs={
        "source_cohort":coverage,
        "models":{
          "primary_500km_richness30plus":primary,
          "100plus_richness_500km":richer,
          "1000km_FE_richness30plus":broad,
          "white_fraction_adjusted_diagnostic":white_control,
          "simple_genus_only_null_diagnostic":genus_only
        },
        "1000km_block_heldout_predictions":cv,
        "source_zip_sha256":{key:val["sha"] for key,val in SOURCES.items()},
        "source_199_global_perm_null_not_genuine_petal_colour":True,
        "photo_hue_not_local_flower_petal_validation":True,
        "native_only_not_filtered":True,
        "sampling_density_not_actual_flower_abundance":True,
        "MOD17_NPP_annual_carbon_flux_not_floral_biomass":True,
        "ESA_AGB_woody_stock_not_herbaceous_floral_biomass":True,
        "LGM_current_site_paleoclimate_not_refugium":True,
        "spatial_FE_and_1000km_holdout_are_not_causal_identification":True
    }
    filejson=root/"ecological_environmental_effects.json"
    filejson.write_text(json.dumps(outputs,indent=2,allow_nan=False)+"\n")
    rows=[]
    for model,info in outputs["models"].items():
        for term,c in info["coefficients"].items():
            rows.append({"model":model,"n_cells":info["n_fit_cells"],
                         "spatial_FE_km":info["spatial_FE_km"],
                         "term":term,"beta_per_SD":c["beta_Shannon_per_global_SD"],
                         "ci_lo":c["CI95"][0],"ci_hi":c["CI95"][1],
                         "p":c["p_cluster"],
                         "q_ecology_BH":c.get("q_BH_seven_ecological_exposures")})
    filecsv=root/"environment_effects_5_sensitivities.csv"
    pd.DataFrame(rows).to_csv(filecsv,index=False)
    outputs["numeric_outputs_sha256"]={f.name:sha(f) for f in
                                      (filegrid,filejson,filecsv)}
    (root/"receipt.json").write_text(json.dumps(outputs,indent=2,allow_nan=False)+"\n")
    print("MOIS_HUE9_ECOLOGY_GLOBAL_ENVIRONMENT_REAL_SUCCESS",json.dumps({
       "coverage":coverage,
       "primary":{
          "n":primary["n_fit_cells"],
          "n_500km_blocks":primary["n_spatial_clusters"],
          "exposures":{k:primary["coefficients"][k] for k in ECOLOGY},
          "joint_tests":primary["joint_Wald"]
       },
       "sensitivity_100km_richness_ge100":{k:richer["coefficients"][k] for k in ECOLOGY},
       "sensitivity_1000km":{k:broad["coefficients"][k] for k in ECOLOGY},
       "white_control":{k:white_control["coefficients"][k] for k in ECOLOGY},
       "spatial_OOF":cv["models"],
       "spatial_OOF_vs_baseline":cv["paired_delta_MSE"]
       },allow_nan=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--null",required=True)
    p.add_argument("--spatial",required=True)
    p.add_argument("--NDVI",required=True)
    p.add_argument("--NPP",required=True)
    p.add_argument("--out",required=True)
    run(p.parse_args())
