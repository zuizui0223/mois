#!/usr/bin/env python3
"""MOIS source-authenticated continuous 9-colour macroecology on latest 1m photos.

No new image download or 4-class forced assignments. Photo-level palette
white fraction (0..1) and Shannon 9-palette entropy (0..1) are responses,
not biological pigments or calibrated CIELAB. Source strict ROI only; explicit
ambiguous palette is recorded but not silently admitted into the primary fit.
"""
from __future__ import annotations
import argparse,json,re
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from scipy.optimize import linear_sum_assignment
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scripts.analysis.mois_real1043930_continuous9_source_lake import sha, PALETTE
from scripts.analysis.mois_latest1043930_water_elevation import (
    PRED, SRC, aggregate, simple_slopes,
)

SCHEMA="mois_real1043930_original_roi_continuous9_macroecology_exploratory_v1"
def require(ok,msg):
    if not ok:raise ValueError("MOIS_REAL1043930_PALETTE_ECOLOGY: "+msg)

def frame(source,source_receipt,climate,geo):
    source,source_receipt,climate,geo=map(Path,
        (source,source_receipt,climate,geo))
    p=json.loads(source_receipt.read_text(encoding="utf-8"))
    require(p["schema"]=="mois_real1043930_original9_with_FCP447500_independent_source_reuse_v1"
        and p["latest_real_photo_ID_attempts"]==1043930
        and p["four_hue_original_labels_unchanged"]==459235
        and p["output_photo_lake_SHA256"]==sha(source),
        "Not the fully recovered and signed 1m original nine-colour ROI data")
    names=["photo_id","inat_taxon_id","species","arm","date","latitude",
           "longitude","continuous9_strict","continuous9_exploratory_ambiguous",
           "continuous9_white","continuous9_entropy","continuous9_bronze"]
    ph=pd.read_parquet(source,columns=names)
    env=pd.read_parquet(climate)
    island=pd.read_parquet(geo,columns=["photo_id","land_type","landmass","area_km2"])
    require(len(ph)==len(env)==len(island)==1043930
            and ph.photo_id.tolist()==env.photo_id.tolist()==island.photo_id.tolist(),
            "Annual environment and GSHHG rows are not aligned by measured photo ID")
    require(int(ph.continuous9_strict.sum())==390541
            and int(ph.continuous9_exploratory_ambiguous.sum())==36224,
            "Original source classes unexpectedly lost nine-colour coverage")
    a=ph[["photo_id","inat_taxon_id","species","arm","date","latitude",
          "longitude","continuous9_strict","continuous9_exploratory_ambiguous",
          "continuous9_white","continuous9_entropy","continuous9_bronze"]].copy()
    a["land_type"]=island.land_type.to_numpy()
    a["landmass"]=island.landmass.to_numpy()
    a["area_km2"]=pd.to_numeric(island.area_km2,errors="coerce").to_numpy(dtype=float)
    for key,col in SRC.items():
        a[key]=pd.to_numeric(env[col],errors="coerce").to_numpy(float)
    return a,p,sha(climate),sha(geo)

def get_analyzable(a,trait):
    d=a[a.continuous9_strict].copy()
    d=d[np.isfinite(d[trait]) & d.latitude.between(-90,90) &
        d.longitude.between(-180,180)].copy()
    for col in PRED:
        d=d[np.isfinite(d[col])].copy()
    d["inat_taxon_id"]=d.inat_taxon_id.fillna("").astype(str)
    d=d[d.inat_taxon_id.str.fullmatch("[0-9]+")]
    d["genus"]=d.species.fillna("").astype(str).str.split().str[0]
    d=d[d.genus.str.fullmatch("[A-Z][A-Za-z-]+",na=False)].copy()
    d["cohort"]=np.where(d.arm.str.startswith("FCP5_"),"FCP5",d.arm).astype(str)
    d["log1p_annual_rainfall"]=np.log1p(d.log1p_annual_rainfall.clip(lower=0))
    d["white_float"]=d[trait].astype(float)
    require(len(d)>50000,"Too few actual original continuous colour locations")
    return d

def within_genus_environment(a,traits):
    outputs=[]
    for trait in traits:
        d=get_analyzable(a,trait)
        for degrees,minimum in ((5,2),(10,3)):
            t=aggregate(d,degrees)
            slopes,_=simple_slopes(t,degrees,min_taxa=minimum,draws=249)
            for r in slopes:
                if r["predictor"] not in ("log1p_annual_rainfall",
                                          "mean_atmospheric_vapour",
                                          "elevation"):
                    continue
                outputs.append({
                    "outcome":trait,"grid_degrees":degrees,
                    "source_strict_photo_classified_env_eligible":len(d),
                    "source_taxon_geocell_comparison_rows":int(
                        r["original_taxon_geocells"]),
                    "source_genera_for_bootstrap":r["genera_for_bootstrap"],
                    "predictor":r["predictor"],
                    "continuous_fraction_pp_per_1full_SD":
                        r["white_pp_per_original_taxon_cell_SD"],
                    "nominal_genus_bootstrap95_pp":
                        r["nominal_95_genus_cluster_bootstrap_pp"],
                    "not_independent_causal_adaptation":True
                })
    return outputs

def island_environment(a,trait):
    d=get_analyzable(a,trait)
    d=d[d.land_type.isin(("island","mainland"))&d.landmass.notna()].copy()
    d["ilat"]=np.floor(d.latitude).astype(int)
    d["ilon"]=np.floor(d.longitude).astype(int)
    for dest,col in [("heat","warmest_month_temperature"),
                     ("rain","log1p_annual_rainfall"),
                     ("solar","annual_solar"),
                     ("elevation","elevation")]:
        d[dest]=d[col].astype(float)
    group=["landmass","land_type","ilat","ilon"]
    tax=d.groupby(group+["inat_taxon_id"],observed=True).agg(
        p_trait=(trait,"mean"),latitude=("latitude","mean"),
        longitude=("longitude","mean"),heat=("heat","mean"),
        rain=("rain","mean"),solar=("solar","mean"),
        elevation=("elevation","mean"),photo_id=("photo_id","size")).reset_index()
    cell=tax.groupby(group,observed=True).agg(
        p_trait=("p_trait","mean"),latitude=("latitude","mean"),
        longitude=("longitude","mean"),heat=("heat","mean"),
        rain=("rain","mean"),solar=("solar","mean"),
        elevation=("elevation","mean"),
        n_photos=("photo_id","sum"),n_taxa=("inat_taxon_id","size")).reset_index()
    cell=cell[cell.n_taxa.ge(5)&cell.n_photos.ge(10)].copy()
    cell["region"]=(np.floor((cell.latitude+90)/10).astype(int).astype(str)
                    +"_"+np.floor((cell.longitude+180)/30).astype(int).astype(str))
    cell["block"]=(np.floor((cell.latitude+90)/10).astype(int).astype(str)
                    +"_"+np.floor((cell.longitude+180)/10).astype(int).astype(str))
    cell["island"]=cell.land_type.eq("island").astype(int)
    both=cell.groupby("region").island.nunique()
    cell=cell[cell.region.isin(both[both.eq(2)].index)].copy()
    require(len(cell)>=100 and cell.island.sum()>=20,
            "Insufficient botanically classifiable original photo island-region cells")
    for name in ("heat","solar","elevation","rain"):
        z=cell[name]
        sd=z.std()
        require(sd>0 and np.isfinite(sd),"No island adjusted environment SD")
        cell["z_"+name]=(z-z.mean())/sd
    fit=smf.ols(
        "p_trait ~ island + z_heat + z_rain + z_solar + z_elevation + C(region)",
        cell).fit(cov_type="cluster",cov_kwds={"groups":cell.block},use_t=True)
    adjusted=fit.conf_int().loc["island"]
    pairs=[]
    predictors=["z_heat","z_rain","z_solar","z_elevation"]
    for region,part in cell.groupby("region"):
        island=part[part.island.eq(1)]
        main=part[part.island.eq(0)]
        if island.empty or main.empty:continue
        dis=cdist(island[predictors].to_numpy(float),main[predictors].to_numpy(float))
        # 1-unit multivariate Euclidean caliper, before any colour outcome.
        cost=np.concatenate((np.where(dis<=1,dis,10000.),
                             np.full((len(island),len(island)),100.)),axis=1)
        ii,jj=linear_sum_assignment(cost)
        for i,j in zip(ii,jj):
            if j<len(main) and dis[i,j]<=1:
                a1=island.iloc[i];b1=main.iloc[j]
                pairs.append(dict(region=str(region),island_landmass=str(a1.landmass),
                                  diff=float(a1.p_trait-b1.p_trait)))
    paired=pd.DataFrame(pairs)
    require(len(paired)>20 and paired.region.nunique()>5,
            "Colour-blind island matches unavailable for strict palette photos")
    model=sm.OLS(paired["diff"],np.ones((len(paired),1))).fit(
        cov_type="cluster",cov_kwds={"groups":paired.region},use_t=True)
    ci=model.conf_int().iloc[0]
    return {
        "outcome":trait,
        "strict_palette_photo_rows_island_mainland":len(d),
        "geographic_source_taxon_equal_cells":int(len(cell)),
        "island_cells":int(cell.island.sum()),
        "matched_environment_pairs":len(paired),
        "matched_regions":int(paired.region.nunique()),
        "unadjusted_strict_palette_photo_island_minus_mainland_pp":
            100*(float(d.loc[d.land_type.eq("island"),trait].mean())-
                 float(d.loc[d.land_type.eq("mainland"),trait].mean())),
        "adjusted_island_minus_mainland_pp":100*float(fit.params["island"]),
        "adjusted_nominal_spatial_cluster95_pp":[100*float(z) for z in adjusted],
        "matched_island_minus_mainland_pp":100*float(model.params.iloc[0]),
        "matched_nominal_region_cluster95_pp":[100*float(z) for z in ci],
        "not_a_same_species_or_native_island_difference":True,
    }

def seasonal(a,trait):
    d=a[a.continuous9_strict].copy()
    d=d[np.isfinite(d[trait])&d.latitude.between(-90,90)].copy()
    d["day"]=pd.to_datetime(d.date,errors="coerce")
    month=d.day.dt.month
    n=d.latitude.ge(23.5)
    s=d.latitude.le(-23.5)
    summer=((n&month.isin((6,7,8)))|(s&month.isin((12,1,2))))
    winter=((n&month.isin((12,1,2)))|(s&month.isin((6,7,8))))
    d=d[summer|winter].copy()
    summer=summer.loc[d.index]
    d["hemisphere"]=np.where(d.latitude.ge(0),"north","south")
    d["season"]=np.where(summer,"summer","winter")
    ans={}
    for hemi,group in d.groupby("hemisphere"):
        z={}
        for season in ("summer","winter"):
            x=group[group.season.eq(season)]
            z[season]={"original_strict_palette_photo_count":int(len(x)),
                       "mean_continuous_value":float(x[trait].mean()),
                       "source_nominal_taxa":int(x.inat_taxon_id.nunique())}
        z["summer_minus_winter_continuous_fraction_pp"] = (
            100*(z["summer"]["mean_continuous_value"]-
                 z["winter"]["mean_continuous_value"]))
        ans[hemi]=z
    return ans

def run(trait,proof,climate,geo,out):
    trait,proof,climate,geo,out=map(Path,(trait,proof,climate,geo,out))
    require(not out.exists(),"Original continuous source ecology cannot overwrite previous result")
    d,r,env_hash,geo_hash=frame(trait,proof,climate,geo)
    traits=("continuous9_white","continuous9_entropy")
    env=within_genus_environment(d,traits)
    island={target:island_environment(d,target) for target in traits}
    season={target:seasonal(d,target) for target in traits}
    report={
        "schema":SCHEMA,
        "source_real_original_photo_attempts":1043930,
        "original_automated_four_hue_photo_ROIs_unchanged":459235,
        "strict_9palette_photo_vectors":390541,
        "separate_ambiguous_9palette_photo_vectors_not_fitted":36224,
        "source_original_9palette_output_SHA256":sha(trait),
        "climate_photo_SHA256":env_hash,"landmass_photo_SHA256":geo_hash,
        "within_genus_geographic_source_cohort_slopes":env,
        "island_and_mainland_comparisons":island,
        "seasonal_photo_composition":season,
        "continuous9_white_is_visible_ROI_colour_fraction_not_anthocyanin":True,
        "continuous9_entropy_is_palette_diversity_not_species_diversity":True,
        "no_unknown_roi_photo_was_assigned_a_colour":True,
        "no_image_reinference_or_download":True,
        "technical_classifiability_selection_bias_remains":True,
        "no_causal_selection_or_botanical_flower_opening_verified":True,
    }
    out.mkdir(parents=True)
    (out/"mois1043930_continuous9_ecological_gradients.json").write_text(
        json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+"\n")
    pd.DataFrame(env).to_csv(out/"mois1043930_continuous9_genus_climate_slopes.csv",index=False)
    with (out/"mois1043930_continuous9_primary_results_receipt.json").open("w") as f:
        json.dump({"source":1043930,"strict":390541,
                   "ecology_source_SHA256":sha(out/"mois1043930_continuous9_ecological_gradients.json"),
                   "ecology_model_slopes_SHA256":sha(out/"mois1043930_continuous9_genus_climate_slopes.csv"),
                   "new_image_fetches":0},f,indent=2,sort_keys=True)
        f.write("\n")
    print("MOIS_1043930_REAL_CONTINUOUS9_MACROECOLOGY",json.dumps({
        "strict_palette_photos":390541,
        "vapour_white_climate_5deg":[r for r in env if r["grid_degrees"]==5
             and r["outcome"]=="continuous9_white"
             and r["predictor"] in ("mean_atmospheric_vapour","elevation","log1p_annual_rainfall")],
        "island_9white":island["continuous9_white"],
        "season_9white":season["continuous9_white"],
        "entropy_5deg":[r for r in env if r["grid_degrees"]==5 and r["outcome"]=="continuous9_entropy"]
    },sort_keys=True),flush=True)
    return report

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("trait","proof","climate","geo","output"):
        p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args()
    run(a.trait,a.proof,a.climate,a.geo,a.output)
if __name__=="__main__":main()
