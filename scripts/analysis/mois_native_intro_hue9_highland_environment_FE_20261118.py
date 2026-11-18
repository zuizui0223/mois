#!/usr/bin/env python3
"""MOIS highland native vs introduced hue9 community contrast: environmental confounding.

Original observational sources:
* WCVP TDWG3 strict/moderate photo9 colour Shannon at each 100km grid, 18,577.
* Global 100km current WorldClim, 2001-2023 NPP, NDVI, ESA woody AGB,
  3GCM LGM-vs-current climatic differences at 12,295 data-complete grids.
* Same 37,104 accepted iNaturalist photographic hue labels and signed
  genus × GBIF occupancy-breadth global 199-permutation expected Shannon.

Response = [H_native−E_native]-[H_introduced−E_introduced] at each
SAME grid with >=30 (or >=100) recorded species in BOTH WCVP statuses.
Elev positive means more introduced clustering relative to native, but not
an effect of elevation independent of climate or introduction history.

500km fixed blocks + 500km cluster robust SE; repeat at 1000km, equal-cell
and smaller-status-species support-weighted. Adjust observation effort
through species numbers, and preserve exact genotype and current ecology
source denominator. This model does not identify elevation/biotic causality.
"""
from __future__ import annotations
import argparse,hashlib,io,json
from pathlib import Path
from zipfile import ZipFile
import numpy as np
import pandas as pd
import statsmodels.api as sm

GEO_SHA="cbe1eef381e3e1a0062dbf5134b047f947338b7c499d2ecc0cd22ee8855909fd"
ENV_SHA="4f5d8fc3a308482fc3645d796a779d1a4725c51ed2afbf6a06e61fd6213f3ee6"
GEO_CSV="global_WCVP_TDWG3_native_only_photo9_island_alpine_polar_south_cells.csv.gz"
ENV_CSV="all_current_GBIF_100km_cells_genus_breadth_entropy_ecology.csv.gz"
NPP="NPP_2001_2023_mean_kgC_m2_yr"
SPEC_COLS={
 "baseline":["elev_km","abs_lat_z","log_native_species_z",
             "log_intro_species_z","island"],
 "modern_climate_vegetation":["elev_km","abs_lat_z",
           "log_native_species_z","log_intro_species_z","island",
           "temp_z","log_rain_z","NPP_z","NDVI_z","AGB_z"],
 "full_present_LGM":["elev_km","abs_lat_z",
           "log_native_species_z","log_intro_species_z","island",
           "temp_z","log_rain_z","NPP_z","NDVI_z","AGB_z",
           "LGM_temp_z","LGM_rain_z"],
 "full_high_altitude_indicator":["high_1500m","abs_lat_z",
           "log_native_species_z","log_intro_species_z","island",
           "temp_z","log_rain_z","NPP_z","NDVI_z","AGB_z",
           "LGM_temp_z","LGM_rain_z"],
 "full_plus_white_fraction_diagnostic":["elev_km","abs_lat_z",
           "log_native_species_z","log_intro_species_z","island",
           "temp_z","log_rain_z","NPP_z","NDVI_z","AGB_z",
           "LGM_temp_z","LGM_rain_z","white_fraction_difference_z"],
}

def check(x,msg):
    if not bool(x):raise ValueError("MOIS_HIGHLAND_HUE9_ECOLOGY: "+msg)

def sha(data):
    if isinstance(data,(str,Path)):data=Path(data).read_bytes()
    return hashlib.sha256(data).hexdigest()

def load_signed(zippath,expected,member,receipt_sha_key):
    check(sha(zippath)==expected,"Original 9hue geography/environment ZIP SHA drift")
    with ZipFile(zippath) as z:
        rec=json.loads(z.read("receipt.json"))
        for fname,digest in rec[receipt_sha_key].items():
            check(sha(z.read(fname))==digest,
                  "Previous real-data numerical CSV/JSON SHA changed "+fname)
        d=pd.read_csv(io.BytesIO(z.read(member)),compression="gzip")
    check(not d[["grid_x","grid_y"]].duplicated().any(),"Source grid duplicates")
    return d,rec

def assemble(geo,environment,mode,min_species):
    check(mode in ("strict","moderate") and min_species in (30,100),
          "Unspecified ecological sensitivity")
    n=f"{mode}_native";i=f"{mode}_introduced"
    columns=[
        "grid_x","grid_y","sampled_land_median_elevation_m",
        "center_latitude","gshhg_land_class",
        f"{n}_n_photo9_species",f"{i}_n_photo9_species",
        f"{n}_Shannon_anomaly",f"{i}_Shannon_anomaly",
        f"{n}_white_species_fraction",f"{i}_white_species_fraction"
    ]
    check(set(columns).issubset(geo),"Signed WCVP native/intro geographical fields missing")
    d=geo[columns].merge(environment[["grid_x","grid_y",
         "present_temperature_C","present_precip_mm",
         NPP,"MODIS_NDVI_2001_2021_mean",
         "ESA_AGB_woody_Mg_ha_2020",
         "threeGCM_delta_temperature_C_median",
         "threeGCM_delta_log_precip_median"]],
         on=["grid_x","grid_y"],how="inner",validate="one_to_one")
    d=d[d[f"{n}_n_photo9_species"].ge(min_species)&
        d[f"{i}_n_photo9_species"].ge(min_species)].reset_index(drop=True)
    if min_species==30:
        check(len(d)=={"strict":3106,"moderate":3602}[mode],
              "Original exact native/intro same 100km cell cohort drift")
    else:
        check(len(d)>=1000,"Few paired rich WCVP source grids")
    check(d.sampled_land_median_elevation_m.notna().all(),"Missing signed true land elevation")
    check(d.present_precip_mm.ge(0).all() and
          d.ESA_AGB_woody_Mg_ha_2020.ge(0).all(),"Unphysical source environment")
    check(d[NPP].notna().all(),"Signed NPP missing in selected samecell cohort")
    d["source_n_native"]=d[f"{n}_n_photo9_species"]
    d["source_n_intro"]=d[f"{i}_n_photo9_species"]
    d["min_source_species"]=d[["source_n_native","source_n_intro"]].min(axis=1)
    d["delta_Shannon_native_minus_intro"]=(
          d[f"{n}_Shannon_anomaly"]-d[f"{i}_Shannon_anomaly"])
    d["white_fraction_difference_pp"]=100*(
        d[f"{n}_white_species_fraction"]-d[f"{i}_white_species_fraction"])
    d["block_500km"]=(d.grid_x.floordiv(5).astype(str)+"|"+
                      d.grid_y.floordiv(5).astype(str))
    d["block_1000km"]=(d.grid_x.floordiv(10).astype(str)+"|"+
                       d.grid_y.floordiv(10).astype(str))
    check(np.isfinite(d[["delta_Shannon_native_minus_intro",
         "white_fraction_difference_pp"]].to_numpy()).all(),
          "Source native/intro delta invalid")
    return d

def standardized_covariates(d):
    raw=pd.DataFrame({
       "elev_km":d.sampled_land_median_elevation_m.astype(float)/1000.,
       "high_1500m":d.sampled_land_median_elevation_m.ge(1500).astype(float),
       "abs_lat_z":d.center_latitude.abs().astype(float),
       "log_native_species_z":np.log1p(d.source_n_native),
       "log_intro_species_z":np.log1p(d.source_n_intro),
       "island":d.gshhg_land_class.eq("island").astype(float),
       "temp_z":d.present_temperature_C.astype(float),
       "log_rain_z":np.log1p(d.present_precip_mm),
       "NPP_z":d[NPP].astype(float),
       "NDVI_z":d.MODIS_NDVI_2001_2021_mean.astype(float),
       "AGB_z":np.log1p(d.ESA_AGB_woody_Mg_ha_2020),
       "LGM_temp_z":d.threeGCM_delta_temperature_C_median.astype(float),
       "LGM_rain_z":d.threeGCM_delta_log_precip_median.astype(float),
       "white_fraction_difference_z":d.white_fraction_difference_pp.astype(float)
    },index=d.index)
    scale={}
    for col in raw.columns:
        if col in ("elev_km","high_1500m","island"):continue
        mean=float(raw[col].mean());sd=float(raw[col].std(ddof=1))
        check(sd>1e-9,"No variance in ecology covariate "+col)
        scale[col]={"mean":mean,"SD":sd}
        raw[col]=(raw[col]-mean)/sd
    check(np.isfinite(raw.to_numpy()).all(),"Nonfinite native-introduced predictors")
    return raw,scale

def fit_one(d,cols,block_km,weighting):
    check(block_km in (500,1000),"Unsupported spatial block")
    b=(d.block_500km if block_km==500 else d.block_1000km).astype(str)
    X,scales=standardized_covariates(d)
    x=X[cols].copy()
    w=(np.ones(len(d),float)
       if weighting=="equal_cells" else d.min_source_species.to_numpy(float))
    wseries=pd.Series(w,index=d.index)
    total_w=wseries.groupby(b).transform("sum")
    xavg=x.mul(w,axis=0).groupby(b).transform("sum").div(total_w,axis=0)
    xc=x-xavg
    y=d.delta_Shannon_native_minus_intro.astype(float)
    yavg=(y*w).groupby(b).transform("sum")/total_w
    yc=y-yavg
    valid=xc.abs().sum(axis=1)>1e-12
    xv=xc.loc[valid];yv=yc.loc[valid]
    groups=b.loc[valid]
    check(len(yv)>700 and groups.nunique()>100,
          "Not enough geographically informative WCVP cells")
    rank=int(np.linalg.matrix_rank(xv.to_numpy()))
    check(rank==len(xv.columns),"Severe source covariate collinearity")
    fit=sm.WLS(yv,xv,weights=w[valid.to_numpy()]).fit(
         cov_type="cluster",cov_kwds={"groups":groups},use_t=False)
    key="high_1500m" if "high_1500m" in cols else "elev_km"
    ci=fit.conf_int().loc[key]
    return {
       "exposure":key,
       "beta_per_1km_or_highland_flag":float(fit.params[key]),
       "CI95":[float(ci.iloc[0]),float(ci.iloc[1])],
       "p_cluster":float(fit.pvalues[key]),
       "n_fit_cells":int(len(yv)),
       "n_fit_spatial_blocks":int(groups.nunique()),
       "block_km":block_km,
       "cell_weights":weighting,
       "n_total_X":len(cols),
       "R2_in_sample":float(fit.rsquared),
       "environment_SD":scales,
       "regression_effects":{s:{
         "beta":float(fit.params[s]),
         "95CI":[float(x) for x in fit.conf_int().loc[s]],
         "p":float(fit.pvalues[s])}
         for s in fit.params.index},
       "not_causal_effect":True,
    }

def run(path_geography,path_environment,out):
    root=Path(out)
    check(not root.exists(),"Refuse overwriting authentic ecology estimates")
    geo,rg=load_signed(path_geography,GEO_SHA,GEO_CSV,"output_sha256")
    env,re=load_signed(path_environment,ENV_SHA,ENV_CSV,
                       "numeric_outputs_sha256")
    check(len(geo)==18577 and len(env)==12295,
          "Signed global cell denominator changed")
    results={}
    for mode in ("strict","moderate"):
        results[mode]={}
        for min_n in (30,100):
            d=assemble(geo,env,mode,min_n)
            label=f"min_both_status_species_{min_n}"
            model={}
            for weighting in ("equal_cells","min_native_intro_species"):
                for block_km in (500,1000):
                    for spec in SPEC_COLS:
                        if (spec.startswith("full_plus_white") and
                            (weighting!="equal_cells" or block_km!=500)):
                            continue
                        keys=SPEC_COLS[spec]
                        model[f"{spec}_{block_km}km_{weighting}"]=fit_one(
                            d,keys,block_km,weighting)
            results[mode][label]={
                "n_matched_same_GBIF_cells":len(d),
                "n_distinct_500km_blocks":int(d.block_500km.nunique()),
                "n_altitude1500plus_cells":int(
                    d.sampled_land_median_elevation_m.ge(1500).sum()),
                "n_lowland_below500_cells":int(
                    d.sampled_land_median_elevation_m.lt(500).sum()),
                "models":model,
            }
            main=model["full_present_LGM_500km_equal_cells"]
            print("MOIS_SAMECELL_ALTITUDE_CLIMATE_CONFOUNDING_REAL_RESULT",
                  json.dumps({
                     "mode":mode,"min_each_status":min_n,
                     "n_cells":len(d),
                     "baseline":model["baseline_500km_equal_cells"][
                         "regression_effects"]["elev_km"],
                     "modern":model["modern_climate_vegetation_500km_equal_cells"][
                         "regression_effects"]["elev_km"],
                     "full":main["regression_effects"]["elev_km"],
                     "full_weighted":model[
                         "full_present_LGM_500km_min_native_intro_species"][
                         "regression_effects"]["elev_km"],
                     "highflag":model["full_high_altitude_indicator_500km_equal_cells"][
                         "regression_effects"]["high_1500m"]
                  },allow_nan=False),flush=True)
    root.mkdir(parents=True)
    csv=root/"native_intro_photo9_highland_500km_environment_effects.csv"
    rows=[]
    for mode,cohorts in results.items():
        for label,z in cohorts.items():
            for name,model in z["models"].items():
                for term,e in model["regression_effects"].items():
                    rows.append({
                      "WCVP_geo_mode":mode,"cohort":label,"model":name,
                      "n_fit_cells":model["n_fit_cells"],
                      "n_fit_blocks":model["n_fit_spatial_blocks"],
                      "term":term,
                      "beta":e["beta"],"CI95_low":e["95CI"][0],
                      "CI95_high":e["95CI"][1],"cluster_p":e["p"]})
    pd.DataFrame(rows).to_csv(csv,index=False)
    report={
      "source_zip_SHA256":{"native_WCVP_geography":GEO_SHA,
                           "environment":ENV_SHA},
      "source_accepted_species":37104,
      "source_GBIF_species_x_100km_pairs":5493931,
      "n_original_199_global_genus_width_permutations":199,
      "same_actual_GBIF_cells_native_introduced":True,
      "spatial_500_and_1000km_block_fixed_effects":True,
      "elevation_raster_sampled_9points_at_100km_not_flower_site_height":True,
      "botanical_native_TDWG3_not_native_individual":True,
      "photo9_labels_not_independently_validated_flowers":True,
      "covariate_adjustment_does_not_identify_direct_elevation_causality":True,
      "results":results,
    }
    model=root/"full_highland_environment_models.json"
    model.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    report["output_sha256"]={p.name:sha(p) for p in (csv,model)}
    (root/"receipt.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--geography",required=True)
    p.add_argument("--environment",required=True)
    p.add_argument("--out",required=True)
    a=p.parse_args()
    run(a.geography,a.environment,a.out)
