#!/usr/bin/env python3
"""MOIS 2026-11-18: test whether native-introduced photographic white-share
differences persist when comparing ONLY genera shared within the same real
GBIF 100-km cell, and whether current moisture or NPP predict the contrast.

Input 1: SHA-locked real species-level WCVP+GBIF shared-genus summary.
Input 2: SHA-locked real same-cell native/introduced colour+climate+NPP.
No new pictures, GBIF observations, or new native classifications.
The two inputs independently originate from the same 37,104 species universe.

The global 9-hue photo label is NOT validated botanically as petal colour.
Native is WCVP TDWG3 regional status, not locally documented wild occurrence.
"""
from __future__ import annotations
import argparse,hashlib,io,json
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

COMMON_SHA="ba725090afbd732c9dacc8e098161aad2bf88042994f5c683db465519969356a"
COMMON_MEMBER="samecell_native_intro_shared_genera_100km_summary.parquet"
PAIRED_SHA="3efaf6ee359aca4a49713d55f5bbadd0c86467ccd4d27c300f1e06936245759a"
PAIRED_MEMBER="WCVP_native_introduced_matched_100km_white_genus_components.csv.gz"
NPP="NPP_2001_2023_mean_kgC_m2_yr"
LGM="threeGCM_delta_temperature_C_median"
RAIN="threeGCM_delta_log_precip_median"
SEED=20261118
N_BOOT=20000

TARGETS={
    "all_genera_native_minus_intro":"observed_native_minus_introduced",
    "shared_genera_equal_genus_native_minus_intro":"genus_equal_native_minus_intro",
    "shared_genera_min_support_native_minus_intro":"matched_genera_support_weighted_difference",
    "shared_genera_pooled_species_native_minus_intro":"matched_genera_species_pooled_difference",
    "raw_minus_shared_equal_gap":"raw_minus_common_genus_equal_gap",
}

def require(cond,message):
    if not bool(cond):
        raise ValueError("MOIS_SAMECELL_SHARED_GENUS_ECOLOGY: "+message)

def sha(x):
    if isinstance(x,(Path,str)):x=Path(x).read_bytes()
    return hashlib.sha256(x).hexdigest()

def load_signed(common,paired):
    import pyarrow.parquet as pq
    require(sha(common)==COMMON_SHA,"Signed old real shared-genera source changed")
    with ZipFile(common) as z:
        receipt=json.loads(z.read("receipt.json"))
        for name,digest in receipt["output_sha256"].items():
            require(sha(z.read(name))==digest,"Old shared-genus inner SHA: "+name)
        old=pq.read_table(io.BytesIO(z.read(COMMON_MEMBER))).to_pandas()
        past=json.loads(z.read("within_genus_native_intro_500km_models.json"))
    require(past["native_introduced_groups"]["strict"][
        "all_matched_cells_before_genus_filter"]==3106,
        "Old strict cohort changed")
    require(past["native_introduced_groups"]["moderate"][
        "all_matched_cells_before_genus_filter"]==3602,
        "Old moderate cohort changed")
    require(not old[["grid_x","grid_y"]].duplicated().any(),
            "Old shared-genus table duplicated geographic cells")
    require(sha(paired)==PAIRED_SHA,"Signed 2026-11-18 paired grid ZIP changed")
    with ZipFile(paired) as z:
        rec=json.loads(z.read("receipt.json"))
        for name,digest in rec["output_sha256"].items():
            require(sha(z.read(name))==digest,"Recent native paired inner SHA: "+name)
        current=pd.read_csv(io.BytesIO(z.read(PAIRED_MEMBER)),
                            compression="gzip",low_memory=False)
    require(len(current)==6708,"Current matched native/intro cell denominator changed")
    require(current.WCVP_geo_mode.value_counts().to_dict()==
            {"strict":3106,"moderate":3602},
            "Current strict/moderate samples changed")
    require(not current[["WCVP_geo_mode","grid_x","grid_y"]].duplicated().any(),
            "Repeated geographical cells within same regional mode")
    return old,current,past

def assemble(old,current,past):
    # Old controls derived from the EXACT same current GBIF/WCVP source.
    keep=["grid_x","grid_y","n_shared_genera",
          "shared_native_species_n","shared_intro_species_n",
          "shared_native_white_n","shared_intro_white_n",
          "genus_equal_native_minus_intro",
          "matched_genera_support_weighted_difference",
          "matched_genera_species_pooled_difference"]
    require(set(keep).issubset(old.columns),"Missing real shared-genus columns")
    d=current.merge(old[keep],on=["grid_x","grid_y"],how="left",
                    validate="many_to_one")
    require(len(d)==len(current),"Merge duplicated actual same-cell source rows")
    d["n_shared_genera"]=d.n_shared_genera.fillna(0).astype(int)
    for mode in ("strict","moderate"):
        t=d[d.WCVP_geo_mode.eq(mode)]
        status=past["native_introduced_groups"][mode]
        require((t.n_shared_genera.gt(0).sum()==status["n_cell_has_shared_genus"]),
                "Old shared genus cells count not replicated "+mode)
        require((t.n_shared_genera.ge(5).sum()==status["n_cell_has_ge5_shared_genus"]),
                "Old >=5 shared genera support count not replicated "+mode)
        for threshold,n_expect in ((5,2516 if mode=="strict" else 2926),
                                   (10,2320 if mode=="strict" else 2697)):
            c=t[(t.n_shared_genera.ge(threshold))
                &t.shared_native_species_n.ge(10)
                &t.shared_intro_species_n.ge(10)]
            require(len(c)==n_expect,
                    f"Old common-genera source cohort at {mode} >= {threshold} disagrees")
    # Compare original species-weight source percent to common genus subcohort
    # percent at identical cell, not to unrelated region/time samples.
    require(d.observed_native_minus_introduced.notna().all(),
            "Signed original native-intro difference missing")
    good=d.n_shared_genera.ge(1)
    require(d.loc[good,"genus_equal_native_minus_intro"].between(-1,1).all(),
            "Within-cell genus-equal white share difference outside probability scale")
    require(d.loc[good,"matched_genera_support_weighted_difference"].between(-1,1).all(),
            "Within-cell support-weighted genus white share invalid")
    d["raw_minus_common_genus_equal_gap"]=(d.observed_native_minus_introduced-
                                            d.genus_equal_native_minus_intro)
    return d

def selected(d,mode,threshold):
    z=d[d.WCVP_geo_mode.eq(mode)]
    z=z[(z.n_shared_genera>=threshold)
       &(z.shared_native_species_n>=10)&(z.shared_intro_species_n>=10)].copy()
    require(len(z)>1700 and z.block_500km.nunique()>150,"Insufficient true shared-genera ecology sample")
    require(z[list(TARGETS.values())].notna().all().all(),
            "Missing common-genus response on selected matched cells")
    z["min_total_status_species"]=z[
        ["colour_species_n__native","colour_species_n__introduced"]].min(axis=1)
    z["shared_native_coverage"]=z.shared_native_species_n/z["colour_species_n__native"]
    z["shared_intro_coverage"]=z.shared_intro_species_n/z["colour_species_n__introduced"]
    require(z.shared_native_coverage.between(0,1).all() and
            z.shared_intro_coverage.between(0,1).all(),
            "Shared genus support > source native/intro species source")
    return z

def block_mean_bootstrap(d,weighting="equal"):
    w=(np.ones(len(d),float) if weighting=="equal" else
       d.min_total_status_species.to_numpy(float))
    g=d[["block_500km"]].copy()
    g["w"]=w
    for k,col in TARGETS.items():
        g[k]=d[col].to_numpy(float)*w
    b=g.groupby("block_500km",sort=False)[["w"]+list(TARGETS)].sum()
    arr=b.to_numpy(float)
    nb=arr.shape[0]
    require(nb>=80,"Only few spatially distinct source 500km clusters")
    rng=np.random.default_rng(SEED+(nb*17)+(0 if weighting=="equal" else 2222))
    sims=np.empty((N_BOOT,len(TARGETS)),float)
    for left in range(0,N_BOOT,250):
        take=min(250,N_BOOT-left)
        index=rng.integers(0,nb,size=(take,nb))
        total=arr[index].sum(axis=1)
        sims[left:left+take,:]=100*total[:,1:]/total[:,[0]]
    observed=100*np.average(d[list(TARGETS.values())].to_numpy(float),
                            axis=0,weights=w)
    ci=np.percentile(sims,[2.5,97.5],axis=0)
    require(abs(observed[0]-observed[1]-observed[-1])<1e-9,
            "Mean gap does not equal raw minus matched genus")
    return {name:{
         "difference_pp":float(observed[i]),
         "CI95_pp":[float(ci[0,i]),float(ci[1,i])],
         "n_spatial_blocks":int(nb),"n_cells":int(len(d))}
         for i,name in enumerate(TARGETS)}

def env_design(d):
    raw=pd.DataFrame({
      "log_current_precip":np.log1p(d.present_precip_mm),
      "NPP":d[NPP],
      "current_temp":d.present_temperature_C,
      "LGM_warming":d[LGM],
      "LGM_precip_change":d[RAIN],
      "NDVI":d.MODIS_NDVI_2001_2021_mean,
      "log_woody_AGB":np.log1p(d.ESA_AGB_woody_Mg_ha_2020),
      "abs_lat":d.lat.abs(),
      "log_native_species_n":np.log1p(d["colour_species_n__native"]),
      "log_intro_species_n":np.log1p(d["colour_species_n__introduced"]),
      "log_shared_genera_n":np.log1p(d.n_shared_genera),
      "shared_native_fraction":d.shared_native_coverage,
      "shared_intro_fraction":d.shared_intro_coverage,
    },index=d.index)
    require(np.isfinite(raw.to_numpy()).all(),"NaN or infinite ecological predictors")
    sd=raw.std(ddof=1)
    require(sd.gt(1e-9).all(),"Constant environmental exposure")
    return (raw-raw.mean())/sd

def fit_spatial(d,blockkm=500,weighting="equal"):
    x=env_design(d)
    b=d.block_500km.astype(str) if blockkm==500 else d.block_1000km.astype(str)
    w=(np.ones(len(d),float) if weighting=="equal" else
       d.min_total_status_species.to_numpy(float))
    # Correct weighted within-block transform, not unweighted centering + WLS.
    gsum=pd.Series(w,index=d.index).groupby(b).transform("sum")
    xavg=x.mul(w,axis=0).groupby(b).transform("sum").div(gsum,axis=0)
    xc=x-xavg
    valid=xc.abs().sum(axis=1)>1e-10
    xc=xc.loc[valid]
    bc=b.loc[valid]
    weights=w[valid.to_numpy()]
    require(len(xc)>800 and bc.nunique()>60,"insufficient fit support after spatial FE")
    fits={}
    for name,col in TARGETS.items():
        y=d[col].astype(float)
        meany=(y*w).groupby(b).transform("sum")/gsum
        yc=(y-meany).loc[valid]
        model=sm.WLS(yc,xc,weights=weights).fit(
            cov_type="cluster",cov_kwds={"groups":bc},use_t=False)
        effects={}
        for term in ("log_current_precip","NPP","LGM_warming","NDVI"):
            ci=model.conf_int().loc[term]
            effects[term]={
              "coef_pp_per_SD":float(100*model.params[term]),
              "CI95_pp":[float(100*ci.iloc[0]),float(100*ci.iloc[1])],
              "p_cluster":float(model.pvalues[term])}
        fits[name]={"effects":effects,"n_fit_cells":len(yc),
                    "n_fit_500km_clusters":int(bc.nunique()),
                    "r2_in_sample":float(model.rsquared),
                    "all_coefs":{k:float(v) for k,v in model.params.items()}}
    require(max(abs(fits["all_genera_native_minus_intro"]["all_coefs"][term]-
                    fits["shared_genera_equal_genus_native_minus_intro"]["all_coefs"][term]-
                    fits["raw_minus_shared_equal_gap"]["all_coefs"][term])
                for term in xc.columns)<1e-11,
            "Original-minus-matched regression must reconstruct exactly")
    # Pre-declared ecological primary covariates only: rain and NPP, two
    # outcomes (raw all genera; common genus equal), BH4 within strict.
    focus=("all_genera_native_minus_intro","shared_genera_equal_genus_native_minus_intro")
    tests=[(label,cov) for label in focus for cov in ("log_current_precip","NPP")]
    q=multipletests([fits[label]["effects"][cov]["p_cluster"]
                     for label,cov in tests],method="fdr_bh")[1]
    for (label,cov),qq in zip(tests,q):
        fits[label]["effects"][cov]["q_BH_four"]=float(qq)
    return {"blockkm":blockkm,"weighting":weighting,"responses":fits,
            "n_matched_cells":len(d),"n_fit_cells":len(xc),
            "n_fit_500km_clusters":int(bc.nunique())}

def run(common,paired,out):
    dst=Path(out)
    require(not dst.exists(),"No overwrite")
    old,current,past=load_signed(common,paired)
    full=assemble(old,current,past)
    results={}
    tab=[]
    outcells=[]
    for mode in ("strict","moderate"):
        results[mode]={}
        for nmin in (5,10):
            data=selected(full,mode,nmin)
            key=f"shared_genera_atleast_{nmin}_10species_each"
            result={"n_cells":len(data),
                    "n_500km_blocks":int(data.block_500km.nunique()),
                    "median_common_genera":float(data.n_shared_genera.median()),
                    "median_shared_native_species":float(data.shared_native_species_n.median()),
                    "median_shared_intro_species":float(data.shared_intro_species_n.median()),
                    "mean_shared_native_fraction":float(data.shared_native_coverage.mean()),
                    "mean_shared_intro_fraction":float(data.shared_intro_coverage.mean()),
                    "means_bootstrap":{},
                    "environmental_regressions":{}}
            for weighting in ("equal","min_status_species"):
                bt=block_mean_bootstrap(data,weighting)
                result["means_bootstrap"][weighting]=bt
                for response,entry in bt.items():
                    tab.append({"WCVP_mode":mode,"shared_genus_min":nmin,
                                "weighting":weighting,"target":response,**entry})
            if nmin==5:
                for weighting in ("equal","min_status_species"):
                    for km in (500,1000):
                        result["environmental_regressions"][f"{km}km_{weighting}"]=fit_spatial(data,km,weighting)
            results[mode][key]=result
            if nmin==5:
                k=["WCVP_geo_mode","grid_x","grid_y","block_500km",
                   "lat","n_shared_genera","shared_native_species_n",
                   "shared_intro_species_n","shared_native_coverage",
                   "shared_intro_coverage","min_total_status_species",
                   "observed_white__native","observed_white__introduced",
                   NPP,LGM,"present_precip_mm","present_temperature_C"]+list(TARGETS.values())
                outcells.append(data[k])
            summary={name:{
                "means":result["means_bootstrap"]["equal"][name],
                "support_weighted":result["means_bootstrap"]["min_status_species"][name]}
                for name in ("all_genera_native_minus_intro",
                             "shared_genera_equal_genus_native_minus_intro",
                             "shared_genera_min_support_native_minus_intro")}
            print("MOIS_SAME_CELL_SHARED_GENUS_ECOLOGY_FINISHED",json.dumps({
                "mode":mode,"min_genus":nmin,"cells":len(data),
                "median_common_genera":result["median_common_genera"],
                "means":summary},allow_nan=False),flush=True)
    dst.mkdir(parents=True)
    main=dst/"matched_common_genus_white_share_500km_block_bootstrap.csv"
    pd.DataFrame(tab).to_csv(main,index=False)
    ffull=dst/"same_cells_common_genera_5percell_env_and_photo.csv.gz"
    pd.concat(outcells,ignore_index=True).to_csv(ffull,index=False,compression="gzip")
    models=dst/"ecological_same_cell_common_genus_models.json"
    models.write_text(json.dumps(results,indent=2,allow_nan=False)+"\n")
    receipt={
      "schema":"MOIS_true_GBIF_WCVP_native_intro_photo_white_common_genus_20261118",
      "global_accepted_photo_hue_species":37104,
      "source_GBIF_species_cell_pairs":5493931,
      "real_WCVP_TDWG3_same_cell_source_all_cells":{"strict":3106,"moderate":3602},
      "primary_sample":"strict, >=5 shared genera, >=10 species from shared genera in each status",
      "primary_response":"equal-weight shared genera native-minus-introduced photo white species share",
      "n_bootstrap_500km":N_BOOT,
      "old_shared_genus_zip_SHA256":COMMON_SHA,
      "recent_paired_env_zip_SHA256":PAIRED_SHA,
      "current_GBIF_distribution_not_LGM_reconstruction":True,
      "photo_hue_not_verified_actual_flower_petals":True,
      "WCVP_TDWG3_native_introduced_not_individual_wild_status":True,
      "strict_and_moderate_not_independent_replicates":True,
      "genus_equal_vs_all_genera_gap_is_NOT_causal_genus_explanation":True,
      "findings":results,
      "files_SHA256":{p.name:sha(p) for p in dst.iterdir()}
    }
    (dst/"receipt.json").write_text(json.dumps(receipt,indent=2,allow_nan=False)+"\n")
    print("MOIS_SAME_CELL_COMMON_GENUS_ECOLOGY_SUCCESS",json.dumps({
       mode:{"n_cells":result["shared_genera_atleast_5_10species_each"]["n_cells"],
             "means":{k:v["difference_pp"] for k,v in
              result["shared_genera_atleast_5_10species_each"][
               "means_bootstrap"]["equal"].items()}}
       for mode,result in results.items()},allow_nan=False),flush=True)

def synthetic_test():
    rng=np.random.default_rng(11)
    n=1320
    d=pd.DataFrame({
      "block_500km":[str(i//8) for i in range(n)],
      "block_1000km":[str(i//16) for i in range(n)],
      "present_precip_mm":rng.uniform(200,3500,n),
      NPP:rng.uniform(.3,2,n),
      "present_temperature_C":rng.normal(13,6,n),
      LGM:rng.uniform(3,12,n),
      RAIN:rng.normal(0,.6,n),
      "MODIS_NDVI_2001_2021_mean":rng.uniform(.1,.9,n),
      "ESA_AGB_woody_Mg_ha_2020":rng.uniform(0,100,n),
      "lat":rng.uniform(-60,60,n),"n_shared_genera":rng.integers(5,35,n),
      "shared_native_species_n":rng.integers(10,80,n),
      "shared_intro_species_n":rng.integers(10,50,n),
      "colour_species_n__native":rng.integers(100,300,n),
      "colour_species_n__introduced":rng.integers(70,150,n)
    })
    d["shared_native_coverage"]=d.shared_native_species_n/d["colour_species_n__native"]
    d["shared_intro_coverage"]=d.shared_intro_species_n/d["colour_species_n__introduced"]
    d["min_total_status_species"]=d[["colour_species_n__native",
                                     "colour_species_n__introduced"]].min(axis=1)
    d["observed_native_minus_introduced"]=rng.normal(-.02,.06,n)
    d["genus_equal_native_minus_intro"]=rng.normal(0,.1,n)
    d["matched_genera_support_weighted_difference"]=rng.normal(0,.08,n)
    d["matched_genera_species_pooled_difference"]=rng.normal(0,.09,n)
    d["raw_minus_common_genus_equal_gap"]=(d.observed_native_minus_introduced-
                                            d.genus_equal_native_minus_intro)
    z=fit_spatial(d)
    require(z["n_fit_cells"]>=900,"Smoke test spatial FE collapsed")
    m=block_mean_bootstrap(d)
    require(len(m)==len(TARGETS),"Smoke test bootstrap responses incomplete")
    print("MOIS_SAMECELL_SHARED_GENUS_SYNTHETIC_PASS",flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--common");p.add_argument("--paired");p.add_argument("--out")
    p.add_argument("--selftest",action="store_true")
    a=p.parse_args()
    if a.selftest:synthetic_test()
    else:
        require(a.common and a.paired and a.out,"Need signed old/new sources and output")
        run(a.common,a.paired,a.out)
