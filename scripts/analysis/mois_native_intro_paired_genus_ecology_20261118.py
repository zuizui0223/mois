#!/usr/bin/env python3
"""Matched native vs introduced plant assemblages: environmental effects on
photo-white species share, decomposed into genus-pool and within-genus
species-composition components at identical genuine GBIF 100-km cells.

Primary ecological question: do contemporary rainfall/NPP gradients associate
with the *native-minus-introduced* white-photo-species-share contrast, and
is any association attributable to genus-pool composition or residual
within-genus species turnover? WCVP is only TDWG3 regional status. Colours
are MOIS photographic proxy, not independently annotated petals.

Data sources are immutable GitHub Actions archives verified twice:
outer ZIP sha256, inner Parquet sha256. No new GBIF/images collected.
"""
from __future__ import annotations
import argparse
import hashlib
import io
import json
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

SOURCE = {
    "decomp": {
        "sha":"e014f500c806242b06ae7251ac399f8849a006f8a5a2836000bdad59c83c691f",
        "member":"GBIF_native_introduced_genus_expected_photo_white_vs_NPP_LGM.parquet",
        "sha_key":"output_SHA256",
        "n":22036,
    },
    "spatial": {
        "sha":"e209257adbc185b4855f3b8d76aeaa50d2e2afe8a1feeb0f76013366e57c4d16",
        "member":"LGM_white_photo9_species_assemblages_100km_cells.parquet",
        "sha_key":"output_sha256",
        "n":12560,
    },
    "NDVI": {
        "sha":"e08e81d1c2afd1a2ee5312f3ebc4aeb5ab7678988c2453a2a9c63d28c82eb17a",
        "member":"WorldClim19_weather_NDVI_18577_global_100km_cells.parquet",
        "sha_key":"output_SHA256",
        "n":18577,
    },
}
NPP="NPP_2001_2023_mean_kgC_m2_yr"
LGM_T="threeGCM_delta_temperature_C_median"
LGM_P="threeGCM_delta_log_precip_median"
PAIR_FIELDS=("observed_white","expected_white_from_genus_composition",
             "within_genus_colour_residual","colour_species_n")
OUTCOMES={
    "observed_native_minus_introduced":"observed_white",
    "genus_pool_native_minus_introduced":"expected_white_from_genus_composition",
    "within_genus_native_minus_introduced":"within_genus_colour_residual",
}
PRIORITY=("log_current_precip","NPP")
SEED=20261118

def require(test,message):
    if not bool(test):raise ValueError("MOIS_NATIVE_PAIRED_ECOLOGY: "+str(message))

def sha(inp):
    if isinstance(inp,(str,Path)):inp=Path(inp).read_bytes()
    return hashlib.sha256(inp).hexdigest()

def read_source(zip_path,which):
    import pyarrow.parquet as pq
    spec=SOURCE[which]
    require(sha(zip_path)==spec["sha"],which+" signed ZIP SHA mismatch")
    with ZipFile(zip_path) as z:
        member=z.read(spec["member"])
        receipt=json.loads(z.read("receipt.json"))
    require(sha(member)==receipt[spec["sha_key"]][spec["member"]],
            which+" internal Parquet SHA changed")
    d=pq.read_table(io.BytesIO(member)).to_pandas()
    require(len(d)==spec["n"],which+" source denominator changed")
    return d

def prepare(decomp,spatial,ndvi):
    keys=["WCVP_geo_mode","WCVP_status","grid_x","grid_y"]
    require(not decomp.duplicated(keys).any(),"Native/genus input duplicate cells")
    require(set(decomp.WCVP_geo_mode)=={"strict","moderate"} and
            set(decomp.WCVP_status)=={"native","introduced"},
            "Original status or geography groups drift")
    require((decomp.colour_species_n>=30).all(),"Source originally required 30 accepted species")
    require(decomp["observed_white"].between(0,1).all() and
            decomp["expected_white_from_genus_composition"].between(0,1).all(),
            "Original source white fraction invalid")
    require(np.allclose(
        decomp.observed_white,
        decomp.expected_white_from_genus_composition+
            decomp.within_genus_colour_residual,atol=1e-12,rtol=0),
        "Original exact genus-pool decomposition drift")
    require(len(spatial)==12560 and not spatial[["grid_x","grid_y"]].duplicated().any(),
            "Original spatial source duplicate")
    require(len(ndvi)==18577 and not ndvi[["grid_x","grid_y"]].duplicated().any(),
            "Original NDVI source duplicate")
    group=["WCVP_geo_mode","grid_x","grid_y"]
    # Every environment value must agree within matched native/introduced cell.
    for v in (NPP,LGM_T):
        spread=decomp.groupby(group)[v].agg(lambda x:float(x.max()-x.min()))
        require((spread.fillna(0).abs()<1e-8).all(),
                "Conflicting environmental exposure between native/introduced")
    pivot=decomp.pivot(index=group,columns="WCVP_status",
                        values=list(PAIR_FIELDS))
    require(pivot.columns.nlevels==2,"Pivot failed")
    flat=pivot.copy()
    flat.columns=[f"{a}__{b}" for a,b in pivot.columns.to_flat_index()]
    flat=flat.reset_index()
    # This is a matched SAME-CELL comparison, not all native/intro singletons.
    for field in PAIR_FIELDS:
        require(f"{field}__native" in flat and
                f"{field}__introduced" in flat,"Missing paired status column")
    wanted=[f"{field}__{status}" for field in PAIR_FIELDS
            for status in ("native","introduced")]
    complete=flat[wanted].notna().all(axis=1)
    matched=flat.loc[complete].copy()
    # From decomposition model source, 2424/2839 were eligible in prior
    # native-minus-introduced matched analyses; do not hardcode those exact
    # counts because this newer complete-case definition can differ slightly.
    require(len(matched)>4500 and
            matched.groupby("WCVP_geo_mode").size().min()>2000,
            "Insufficient actual native-introduced matched cells")
    aux=decomp.drop_duplicates(group)[group+[NPP,LGM_T,LGM_P,
                                "block_500km"]].copy()
    require(not aux.duplicated(group).any(),"Duplicate environmental cell keys")
    matched=matched.merge(aux,on=group,how="left",validate="one_to_one")

    s_cols=["grid_x","grid_y","lat","lon","block_500km","block_1000km",
            "present_temperature_C","present_precip_mm",
            "ESA_AGB_woody_Mg_ha_2020"]
    require(set(s_cols).issubset(spatial.columns),
            "Original WorldClim1.4 present climate or AGB not accessible")
    matched=matched.merge(spatial[s_cols],on=["grid_x","grid_y"],how="left",
                          validate="many_to_one",suffixes=("_original","_spatial"))
    require((matched.block_500km_original.astype(str)==
             matched.block_500km_spatial.astype(str)).all(),
            "Spatial 500 km block source changed")
    matched=matched.rename(columns={"block_500km_spatial":"block_500km"})
    matched=matched.drop(columns=["block_500km_original"])
    n_cols=["grid_x","grid_y","MODIS_NDVI_2001_2021_mean"]
    matched=matched.merge(ndvi[n_cols],on=["grid_x","grid_y"],how="left",
                          validate="many_to_one")
    require(len(matched)==complete.sum(),"Environment join damaged same cells")
    numeric=[NPP,LGM_T,LGM_P,"lat","present_temperature_C",
             "present_precip_mm","ESA_AGB_woody_Mg_ha_2020",
             "MODIS_NDVI_2001_2021_mean"]
    eligible=matched[numeric].notna().all(axis=1)
    eligible &= matched.present_precip_mm.ge(0)
    eligible &= matched.ESA_AGB_woody_Mg_ha_2020.ge(0)
    # Original source has already complete NPP; enforce before fitting.
    require(eligible.all(),"Signed original matched cells lost crucial environment data")
    for name,col in OUTCOMES.items():
        matched[name]=(matched[f"{col}__native"]
                       -matched[f"{col}__introduced"]).astype(float)
    require(np.allclose(
        matched["observed_native_minus_introduced"],
        matched["genus_pool_native_minus_introduced"]+
        matched["within_genus_native_minus_introduced"],atol=1e-12,rtol=0),
        "Paired genus decomposition does not add up")
    require(matched.block_500km.notna().all() and
            matched.block_1000km.notna().all(),"Missing spatial clusters")
    matched["minimum_native_introduced_species_count"]=matched[[
        "colour_species_n__native","colour_species_n__introduced"]].min(axis=1)
    return matched

def design(d):
    a=pd.DataFrame({
        "NPP":d[NPP].astype(float),
        "log_current_precip":np.log1p(d.present_precip_mm.astype(float)),
        "current_temp":d.present_temperature_C.astype(float),
        "LGM_warming":d[LGM_T].astype(float),
        "LGM_rain":d[LGM_P].astype(float),
        "NDVI":d.MODIS_NDVI_2001_2021_mean.astype(float),
        "log_woody_AGB":np.log1p(d.ESA_AGB_woody_Mg_ha_2020.astype(float)),
        "abs_lat":d.lat.abs().astype(float),
        "log_native_species":np.log1p(d["colour_species_n__native"].astype(float)),
        "log_introduced_species":np.log1p(d["colour_species_n__introduced"].astype(float)),
    },index=d.index)
    require(np.isfinite(a.to_numpy()).all(),"Nonfinite model predictors")
    sd=a.std(ddof=1)
    require((sd>1e-10).all(),"No predictor variation")
    x=(a-a.mean())/sd
    return x,{col:{"mean":float(a[col].mean()),"SD":float(sd[col])} for col in a.columns}

def fit_cohort(d,block_km=500,weighted=False,drop_npp=False):
    require(len(d)>500 and d.block_500km.nunique()>30,
            "Insufficient matched native/intro spatial data")
    x,scales=design(d)
    if drop_npp:x=x.drop(columns="NPP")
    b=(d.block_500km if block_km==500 else d.block_1000km).astype(str)
    clusters=b
    weights=(d.minimum_native_introduced_species_count.astype(float).to_numpy()
             if weighted else np.ones(len(d),dtype=float))
    require(np.isfinite(weights).all() and weights.min()>0,"Invalid species support")
    # True weighted within-block transform (unlike unweighted demeaning then WLS).
    gr=pd.Series(weights,index=d.index).groupby(b,sort=False).transform("sum")
    weighted_x=x.mul(weights,axis=0)
    group_x=weighted_x.groupby(b,sort=False).transform("sum").div(gr,axis=0)
    xc=x-group_x
    good=xc.abs().sum(axis=1)>1e-12
    require(good.sum()>500 and b[good].nunique()>25,"No spatial within-block support")
    xc=xc.loc[good]
    bx=b.loc[good]
    wg=weights[good.to_numpy()]
    y={}
    for field in OUTCOMES:
        original=d[field].astype(float)
        ymean=(original*weights).groupby(b,sort=False).transform("sum")/gr
        y[field]=(original-ymean).loc[good]
    require(np.allclose(
        y["observed_native_minus_introduced"],
        y["genus_pool_native_minus_introduced"]+
        y["within_genus_native_minus_introduced"],atol=1e-11,rtol=0),
        "Within-block paired responses must add up")
    fitted={}
    for field,target in y.items():
        result=sm.WLS(target,xc,weights=wg).fit(
            cov_type="cluster",cov_kwds={"groups":bx},use_t=False)
        eff={}
        for name in xc.columns:
            ci=result.conf_int().loc[name]
            eff[name]={
                "coef_pp_per_global_cohort_SD":float(result.params[name]*100),
                "CI95_pp":[float(ci.iloc[0]*100),float(ci.iloc[1]*100)],
                "p_cluster":float(result.pvalues[name])
            }
        fitted[field]={
            "n_cells":int(len(target)),
            "n_spatial_clusters":int(bx.nunique()),
            "environmental_effects":eff,
            "R2_within_block":float(result.rsquared),
            "coefficient_vector":{key:float(v) for key,v in result.params.items()},
        }
    deviation={
        k:(fitted["observed_native_minus_introduced"]["coefficient_vector"][k]-
           fitted["genus_pool_native_minus_introduced"]["coefficient_vector"][k]-
           fitted["within_genus_native_minus_introduced"]["coefficient_vector"][k])
        for k in xc.columns
    }
    require(max(abs(v) for v in deviation.values())<1e-10,
            "Fitted same-cell environmental beta must decompose exactly")
    return {
        "n_source_matched_cells":int(len(d)),
        "n_fit_cells":int(good.sum()),
        "n_spatial_clusters":int(bx.nunique()),
        "block_km":int(block_km),
        "weighting":"min_both_status_species" if weighted else "equal_100km_cells",
        "NPP_in_model":not drop_npp,
        "scales":scales,
        "components":fitted,
        "decomposition_max_abs_parameter_error":float(max(abs(v) for v in deviation.values())),
    }

def run(source_paths,out):
    dst=Path(out)
    require(not dst.exists(),"Refuse prior output overwrite")
    inputs={key:read_source(source_paths[key],key) for key in SOURCE}
    d=prepare(inputs["decomp"],inputs["spatial"],inputs["NDVI"])
    models={}
    for mode in ("strict","moderate"):
        t=d[d.WCVP_geo_mode.eq(mode)].copy().reset_index(drop=True)
        require(len(t)>2000,"Insufficient WCVP mode cells")
        models[mode]={}
        for name,block,weighted,drop_npp in [
            ("primary_500km_equal_cell",500,False,False),
            ("500km_min_status_support_weighted",500,True,False),
            ("1000km_equal_cell",1000,False,False),
            ("500km_equal_cell_without_NPP",500,False,True),
        ]:
            models[mode][name]=fit_cohort(t,block,weighted,drop_npp)
            print("MOIS_NATIVE_INTRO_GENUINE_MATCHED",mode,name,
                  "n",len(t),flush=True)
        primary=models[mode]["primary_500km_equal_cell"]["components"]
        # Four predeclared tests per mode, correction within that mode:
        # precipitation/NPP on raw contrast and within-genus residual.
        tests=[(f,k) for f in (
            "observed_native_minus_introduced",
            "within_genus_native_minus_introduced"
        ) for k in PRIORITY]
        q=multipletests([primary[f]["environmental_effects"][k]["p_cluster"]
                         for f,k in tests],method="fdr_bh")[1]
        for (f,k),v in zip(tests,q):
            primary[f]["environmental_effects"][k]["q_BH_four_preregistered"]=float(v)
        models[mode]["cell_mean_native_photo_white_fraction"]=float(
            t["observed_white__native"].mean())
        models[mode]["cell_mean_introduced_photo_white_fraction"]=float(
            t["observed_white__introduced"].mean())
        models[mode]["cell_mean_paired_native_minus_introduced_white_pp"]=float(
            100*t["observed_native_minus_introduced"].mean())
        models[mode]["cell_mean_paired_genus_pool_component_pp"]=float(
            100*t["genus_pool_native_minus_introduced"].mean())
        models[mode]["cell_mean_paired_within_genus_component_pp"]=float(
            100*t["within_genus_native_minus_introduced"].mean())
    dst.mkdir(parents=True)
    outputs=[]
    for mode in models:
        for name,v in models[mode].items():
            if not isinstance(v,dict) or "components" not in v:continue
            for field,model in v["components"].items():
                for predictor,effect in model["environmental_effects"].items():
                    outputs.append({
                        "WCVP_mode":mode,"model":name,"response_component":field,
                        "predictor":predictor,"n_matched_cells":v["n_fit_cells"],
                        "n_spatial_clusters":v["n_spatial_clusters"],
                        "beta_pp_per_SD":effect["coef_pp_per_global_cohort_SD"],
                        "CI95_low_pp":effect["CI95_pp"][0],
                        "CI95_high_pp":effect["CI95_pp"][1],
                        "p_cluster":effect["p_cluster"],
                        "q_BH_predeclared":effect.get("q_BH_four_preregistered"),
                    })
    fcsv=dst/"native_minus_introduced_genus_pool_vs_within_genus_coefficients.csv"
    pd.DataFrame(outputs).to_csv(fcsv,index=False)
    fgrid=dst/"WCVP_native_introduced_matched_100km_white_genus_components.csv.gz"
    export=["WCVP_geo_mode","grid_x","grid_y","lat","lon","block_500km","block_1000km",
            "observed_white__native","observed_white__introduced",
            "colour_species_n__native","colour_species_n__introduced",
            "observed_native_minus_introduced","genus_pool_native_minus_introduced",
            "within_genus_native_minus_introduced",NPP,LGM_T,LGM_P,
            "present_temperature_C","present_precip_mm",
            "MODIS_NDVI_2001_2021_mean","ESA_AGB_woody_Mg_ha_2020"]
    d[export].to_csv(fgrid,index=False,compression="gzip")
    (dst/"full_models.json").write_text(json.dumps(models,indent=2,allow_nan=False)+"\n")
    summarize={
        "sources_sha256":{name:SOURCE[name]["sha"] for name in SOURCE},
        "source_units":"accepted photo-labelled species per real GBIF 100km cell",
        "n_true_accepted_species_from_source":37104,
        "n_original_GBIF_species_cell_pairs":5493931,
        "strict_matched_cells":int(d.WCVP_geo_mode.eq("strict").sum()),
        "moderate_matched_cells":int(d.WCVP_geo_mode.eq("moderate").sum()),
        "model":"identical-current-grid WCVP native minus introduced, genus pool plus within genus exact decomposition, environment-adjusted spatial FE",
        "primary_predictors":list(PRIORITY),
        "WCVP_native_is_TDWG3_status_not_wild_individual":True,
        "not_measured_local_flower_petals":True,
        "current_LGM_climate_not_historical_range_or_fossil_refugia":True,
        "NPP_is_annual_flux_not_flower_allocation":True,
        "genus_reference_from_same_photo_catalogue_not_external_crossfit":True,
        "findings":models
    }
    (dst/"summary.json").write_text(json.dumps(summarize,indent=2,allow_nan=False)+"\n")
    summarize["output_sha256"]={p.name:sha(p) for p in sorted(dst.iterdir())}
    (dst/"receipt.json").write_text(json.dumps(summarize,indent=2,allow_nan=False)+"\n")
    for mode in ("strict","moderate"):
        prim=models[mode]["primary_500km_equal_cell"]["components"]
        print("MOIS_NATIVE_PAIRED_ECOLOGY_SUCCESS",json.dumps({
           "mode":mode,"n":models[mode]["primary_500km_equal_cell"]["n_fit_cells"],
           "white_native_minus_intro_pp":
             models[mode]["cell_mean_paired_native_minus_introduced_white_pp"],
           "effects":{y:{k:prim[y]["environmental_effects"][k]
                            for k in ("NPP","log_current_precip","LGM_warming","NDVI")}
                   for y in OUTCOMES},
           "max_decomp_error":models[mode]["primary_500km_equal_cell"][
                   "decomposition_max_abs_parameter_error"],
           "weight_sensitivity":{y:{
                   k:models[mode]["500km_min_status_support_weighted"]["components"][y][
                      "environmental_effects"][k]["coef_pp_per_global_cohort_SD"]
                   for k in PRIORITY} for y in OUTCOMES}
        },allow_nan=False),flush=True)

def smoke():
    # Independent synthetic test, never present as ecology. Includes exact
    # additivity of each cell and each spatial-FE slope in WLS/OLS.
    rng=np.random.default_rng(SEED)
    n=3000
    d=pd.DataFrame({
       "block_500km":[f"b{i//12}" for i in range(n)],
       "block_1000km":[f"z{i//24}" for i in range(n)],
       NPP:rng.uniform(.1,3,n),LGM_T:rng.uniform(1,15,n),
       LGM_P:rng.normal(0,.5,n),"present_precip_mm":rng.uniform(80,4500,n),
       "present_temperature_C":rng.normal(10,8,n),
       "ESA_AGB_woody_Mg_ha_2020":rng.uniform(0,180,n),
       "MODIS_NDVI_2001_2021_mean":rng.uniform(.15,.8,n),
       "lat":rng.uniform(-65,65,n),
       "colour_species_n__native":rng.integers(30,120,n),
       "colour_species_n__introduced":rng.integers(30,90,n)
    })
    a=rng.uniform(-.2,.2,n);b=rng.uniform(-.1,.1,n)
    d["genus_pool_native_minus_introduced"]=a
    d["within_genus_native_minus_introduced"]=b
    d["observed_native_minus_introduced"]=a+b
    d["minimum_native_introduced_species_count"]=d[
        ["colour_species_n__native","colour_species_n__introduced"]].min(axis=1)
    for weighted in (False,True):
        m=fit_cohort(d,500,weighted)
        require(m["decomposition_max_abs_parameter_error"]<1e-9,
                "Additivity failed")
    print("MOIS_NATIVE_INTRO_ECOLOGICAL_MATCHED_SYNTHETIC_PASS",flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--decomp");p.add_argument("--spatial");p.add_argument("--NDVI")
    p.add_argument("--out");p.add_argument("--selftest",action="store_true")
    a=p.parse_args()
    if a.selftest:smoke()
    else:
        require(a.decomp and a.spatial and a.NDVI and a.out,
                "Three signed inputs and output required")
        run({"decomp":a.decomp,"spatial":a.spatial,"NDVI":a.NDVI},a.out)
