#!/usr/bin/env python3
"""MOIS ecological introduction-filtering test from original GBIF/WCVP source.

Hypothesis (primary): within genera, species with a WHITE modal photographic
nine-colour label have different propensity to be *recorded as introduced in
>=1 WCVP TDWG3 botanical region* than nonwhite species, conditional on
native-range footprint, contemporary climate, vegetation, and sampling.
This is a species-level trait association, not evolutionary change following
introduction, introduction success, or true botanically verified petal colour.

Independent primary species (one row per COLXR accepted species);
WCVP matches unique-exact accepted binomials only. Same genus fixed-effect
linear probability model, genus cluster SE; sensitivity to genus balancing,
photo status and trait confidence; 3-fold genus-disjoint prediction increment.
"""
from __future__ import annotations
import argparse, hashlib, io, json, re
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, log_loss

ENV_SHA="d8b97528188d9f8d63cb802d5b3f6e8f7f251937fd30af42aeae47e43ebad1fc"
ENV_MEMBER="BIO19_weather_NDVI_ESACCIv7_AGB_37104_accepted_species.parquet"
NAME_FILE="mois350399_wcvp_source_taxon_exact_name_resolution.csv.gz"
DIST_FILE="mois350399_wcvp_accepted_native_distribution_TDWG3.csv.gz"
MANIFEST="mois350399_wcvp_TDWG3_source_only_taxon_native_receipt.json"
COLORS=("white","yellow","orange","red","pink","magenta","purple","blue","bronze")
VALID_STATES={"native_at_TDWG3","introduced_at_TDWG3","doubtful","extinct",
              "unknown_flag","contradictory_TDWG3_records"}
BINOMIAL=re.compile(r"^[A-Z][A-Za-z-]+ [a-z][A-Za-z-]+$")
SEED=20261118
COVARS=("native_region_count_log","BIO01","BIO12_log","vpd_kPa_approx",
        "MODIS_NDVI_2001_2021_mean","AGB_log","abs_lat",
        "GBIF_records_log","GBIF_cells_log","srad_monthly_mean",
        "wind_monthly_mean")

def require(cond,message):
    if not bool(cond):raise ValueError("MOIS_INTRO_FILTER_COLOUR: "+str(message))
def sha(v):
    if isinstance(v,(Path,str)):v=Path(v).read_bytes()
    return hashlib.sha256(v).hexdigest()

def signed_env(path):
    import pyarrow.parquet as pq
    require(sha(path)==ENV_SHA,"Genuine source environment ZIP SHA differs")
    with ZipFile(path) as z:
        original=z.read(ENV_MEMBER)
        meta=json.loads(z.read("receipt.json"))
    require(sha(original)==meta["output_SHA256"][ENV_MEMBER]
            and meta["source_global_accepted_species"]==37104,
            "Env Parquet SHA or COLXR accepted species denominator changed")
    d=pq.read_table(io.BytesIO(original)).to_pandas()
    d["key"]=d.key.astype(str)
    require(len(d)==37104 and d.key.is_unique and
            set(d.source_hue9)==set(COLORS),
            "Original accepted species/photo colour source changed")
    require(d.records_passed_spatial_QA.sum()==396142825,
            "Original GBIF valid occurrence total changed")
    return d

def signed_wcvp(folder):
    p=Path(folder)
    m=json.loads((p/MANIFEST).read_text())
    require(m["source_photo_attempts"]==350399
            and m["taxa_with_any_source_photos"]==74821
            and m["wcvp_distributions_resolution"]=="TDWG_level3_only",
            "Original WCVP imported release manifest no longer valid")
    for fname in (NAME_FILE,DIST_FILE):
        require(sha(p/fname)==m["data_files_sha256"][fname],
                "Kew WCVP original gzip SHA differs: "+fname)
    names=pd.read_csv(p/NAME_FILE,dtype=str,keep_default_na=False)
    dist=pd.read_csv(p/DIST_FILE,dtype=str,keep_default_na=False)
    require(len(names)==74821 and names.inat_taxon_id.is_unique,
            "WCVP exact original names count changed")
    require(len(dist)>600000 and not dist[
        ["wcvp_accepted_plant_name_id","tdwg_level3_code"]].duplicated().any(),
            "WCVP original TDWG3 observations non-unique")
    require(dist.native_introduced_state.isin(VALID_STATES).all(),
            "WCVP status changed or unexpected")
    return names,dist,m

def prepare(env,names,dist):
    unique=names.loc[
        names.wcvp_match_state.eq("unique_exact_name_to_accepted_id")
        &names.wcvp_accepted_plant_name_id.ne("")
        &names.source_species_exact_binomial.str.match(BINOMIAL),
        ["source_species_exact_binomial","wcvp_accepted_plant_name_id"]
    ].drop_duplicates()
    n=unique.groupby("source_species_exact_binomial").wcvp_accepted_plant_name_id.nunique()
    unique=unique.loc[unique.source_species_exact_binomial.isin(n.index[n.eq(1)])]
    unique=unique.drop_duplicates("source_species_exact_binomial")
    d=env.merge(unique,left_on="species",
                right_on="source_species_exact_binomial",how="left",
                validate="many_to_one")
    require(len(d)==37104 and d.key.is_unique,"WCVP exact binomial merge lost species")
    match_n=int(d.wcvp_accepted_plant_name_id.notna().sum())
    require(36000<=match_n<=37104,
            "Unexpected number of verified WCVP matchable accepted species")
    st=dist[dist.native_introduced_state.isin(
        ("native_at_TDWG3","introduced_at_TDWG3"))].copy()
    st["native"]=st.native_introduced_state.eq("native_at_TDWG3").astype(int)
    st["introduced"]=st.native_introduced_state.eq("introduced_at_TDWG3").astype(int)
    counts=st.groupby("wcvp_accepted_plant_name_id",sort=False).agg(
        native_region_count=("native","sum"),
        introduced_region_count=("introduced","sum"),
    ).reset_index()
    d=d.merge(counts,on="wcvp_accepted_plant_name_id",how="left",
              validate="many_to_one")
    d[["native_region_count","introduced_region_count"]]=d[
        ["native_region_count","introduced_region_count"]].fillna(0).astype(int)
    d["wcvp_matched"]=d.wcvp_accepted_plant_name_id.notna()
    d["has_native_WCVP"]=d.native_region_count.gt(0)
    d["recorded_introduced"]=d.introduced_region_count.gt(0).astype(int)
    d["photo_white"]=d.source_hue9.eq("white").astype(int)
    numeric=["BIO01","BIO12","vpd_kPa_approx",
      "MODIS_NDVI_2001_2021_mean","ESA_AGB_woody_Mg_ha_2020",
      "mean_abs_latitude_over_occupied_cells","records_passed_spatial_QA",
      "occupied_grid_cells","srad_monthly_mean","wind_monthly_mean"]
    needed=numeric+["genus","hue_status"]
    good=d.wcvp_matched & d.has_native_WCVP & d[needed].notna().all(axis=1)
    good &= d.BIO12.ge(0) & d.ESA_AGB_woody_Mg_ha_2020.ge(0)
    good &= d.occupied_grid_cells.gt(0)&d.records_passed_spatial_QA.gt(0)
    d=d.loc[good].copy().reset_index(drop=True)
    require(len(d)>22000 and d.genus.nunique()>4000,
            "Insufficient original native WCVP-referenced accepted species")
    require(d.key.is_unique and d.photo_white.nunique()==2
            and d.recorded_introduced.nunique()==2,
            "Introduced status or photographic white outcome unavailable")
    # WCVP accepted IDs are not always uniquely keyed by COLXR taxonomy:
    # source species still remain one independent accepted key per row.
    return d,{"accepted_source":37104,"unique_WCVP_matched":match_n,
              "accepted_with_native_and_full_current_features":len(d),
              "n_genera":int(d.genus.nunique()),
              "n_recorded_in_any_introduced_region":int(d.recorded_introduced.sum()),
              "n_photo_white":int(d.photo_white.sum())}

def design(d,photo_status=True):
    out=pd.DataFrame({
      "native_region_count_log":np.log1p(d.native_region_count),
      "BIO01":d.BIO01.astype(float),
      "BIO12_log":np.log1p(d.BIO12.astype(float)),
      "vpd_kPa_approx":d.vpd_kPa_approx.astype(float),
      "MODIS_NDVI_2001_2021_mean":d.MODIS_NDVI_2001_2021_mean.astype(float),
      "AGB_log":np.log1p(d.ESA_AGB_woody_Mg_ha_2020.astype(float)),
      "abs_lat":d.mean_abs_latitude_over_occupied_cells.astype(float),
      "GBIF_records_log":np.log1p(d.records_passed_spatial_QA.astype(float)),
      "GBIF_cells_log":np.log1p(d.occupied_grid_cells.astype(float)),
      "srad_monthly_mean":d.srad_monthly_mean.astype(float),
      "wind_monthly_mean":d.wind_monthly_mean.astype(float),
    },index=d.index)
    require(np.isfinite(out.to_numpy(float)).all(),
            "Original current covariates non-finite")
    sd=out.std(ddof=1)
    require(sd.gt(1e-10).all(),"No variation in source ecology")
    out=(out-out.mean())/sd
    if photo_status:
        out=pd.concat([out,pd.get_dummies(
            d.hue_status.astype(str),prefix="photo_status",drop_first=True,
            dtype=float)],axis=1)
    return out

def within_genus(d,white_col="photo_white",weighting="species",
                 photo_status=True,cohort="full"):
    y=d.recorded_introduced.astype(float)
    mixed=y.groupby(d.genus).transform("nunique").eq(2)
    d=d.loc[mixed].copy().reset_index(drop=True)
    require(len(d)>8000 and d.genus.nunique()>500,
            "Not enough mixed-introduction-status genera")
    x=design(d,photo_status=photo_status)
    x["colour_trait"]=d[white_col].astype(float)
    if white_col=="photo_white":
        require(x.colour_trait.nunique()==2,"No photographic white contrast")
    group=d.genus.astype(str)
    xx=x-x.groupby(group).transform("mean")
    yy=y.loc[mixed].reset_index(drop=True)
    yy=yy-yy.groupby(group).transform("mean")
    xx=xx.loc[:,xx.std()>1e-10]
    require("colour_trait" in xx.columns,"Colour trait has no within-genus variation")
    weights=(np.ones(len(d)) if weighting=="species" else
             1/d.genus.map(d.genus.value_counts()).to_numpy(float))
    fit=sm.WLS(yy,xx,weights=weights).fit(
        cov_type="cluster",cov_kwds={"groups":group},use_t=False)
    ci=fit.conf_int().loc["colour_trait"]
    return {"cohort":cohort,"photo_predictor":white_col,
            "weighting":weighting,"photo_status_adjusted":photo_status,
            "N_mixed_introduction_species":int(len(d)),
            "N_mixed_introduction_genera":int(group.nunique()),
            "white_photo_to_intro_status_risk_difference_pp":
              100*float(fit.params["colour_trait"]),
            "CI95_pp":[100*float(ci.iloc[0]),100*float(ci.iloc[1])],
            "p_genus_cluster":float(fit.pvalues["colour_trait"]),
            "n_covariates":len(xx.columns)}

def genus_folds_oof(d):
    df=design(d)
    y=d.recorded_introduced.to_numpy(int)
    c=d.source_hue9.to_numpy(str)
    genera=d.genus.to_numpy(str)
    extra=d.photo_white.to_numpy(float).reshape(-1,1)
    x=df.to_numpy(float)
    models={"ecology_only":x,"ecology_plus_white":np.column_stack((x,extra))}
    splits=list(StratifiedGroupKFold(
        n_splits=3,shuffle=True,random_state=SEED
    ).split(x,y,genera))
    results={k:np.full(len(d),np.nan) for k in models}
    for fold,(tr,te) in enumerate(splits):
        require(set(genera[tr]).isdisjoint(set(genera[te])),
                "Train/test genus overlap")
        for name,X in models.items():
            imputer=SimpleImputer(strategy="median").fit(X[tr])
            xtr=imputer.transform(X[tr])
            xte=imputer.transform(X[te])
            scaler=StandardScaler().fit(xtr)
            clf=LogisticRegression(C=.1,max_iter=750,solver="lbfgs")
            clf.fit(scaler.transform(xtr),y[tr])
            results[name][te]=clf.predict_proba(scaler.transform(xte))[:,1]
    require(all(np.isfinite(x).all() for x in results.values()),
            "Genus held-out predictions have missing species")
    metric={name:{"ROC_AUC":float(roc_auc_score(y,prob)),
                  "log_loss":float(log_loss(y,prob))}
            for name,prob in results.items()}
    p=results["ecology_only"];q=results["ecology_plus_white"]
    each=-(y*np.log(q)+(1-y)*np.log1p(-q))+(y*np.log(p)+(1-y)*np.log1p(-p))
    ids,g=np.unique(genera,return_inverse=True)
    sums=np.bincount(g,weights=each,minlength=len(ids))
    sizes=np.bincount(g,minlength=len(ids))
    rng=np.random.default_rng(SEED)
    bs=[]
    for _ in range(2000):
        wt=np.bincount(rng.integers(0,len(ids),len(ids)),minlength=len(ids))
        bs.append(float(np.dot(wt,sums)/np.dot(wt,sizes)))
    return {
        "n_species":len(d),"n_heldout_genera":len(ids),
        "folds":3,"models":metric,
        "delta_logloss_plus_white_minus_ecology_only":float(each.mean()),
        "delta_logloss_genus_bootstrap_CI95":
            [float(x) for x in np.percentile(bs,[2.5,97.5])],
        "introductions_not_directly_observed_as_field_establishment":True
    }

def run(env,wcvp,out):
    output=Path(out)
    require(not output.exists(),"Do not overwrite evidence")
    data=signed_env(env)
    nm,dist,manifest=signed_wcvp(wcvp)
    d,source=prepare(data,nm,dist)
    rates=d.groupby("photo_white").agg(
       species=("key","size"),
       recorded_introductions=("recorded_introduced","sum"),
       median_native_regions=("native_region_count","median"),
       mean_introduced_region_count=("introduced_region_count","mean"),
    ).reset_index()
    rates["intro_record_fraction"]=(rates.recorded_introductions/rates.species)
    estimates=[
      within_genus(d,weighting="species"),
      within_genus(d,weighting="genus"),
      within_genus(d,weighting="species",photo_status=False),
    ]
    # Strict three independent photographic observations is a highly
    # colour-selected cohort; report sensitivity and sample size, not truth.
    h=d.hue_status.eq("sampling_consensus_single")
    if int(h.sum())>850:
        d2=d.loc[h].copy().reset_index(drop=True)
        try:estimates.append(within_genus(
             d2,cohort="3plus_photo_consensus_selected"))
        except ValueError:
            estimates.append({"cohort":"3plus_photo_consensus_selected",
                              "status":"insufficient_mixed_genera"})
    # Exploratory 9-way phenotype contrasts: WCVP introduction propensity by hue.
    colour_est=[]
    for hue in COLORS:
        dc=d.copy()
        dc["one_vs_rest"]=dc.source_hue9.eq(hue).astype(int)
        r=within_genus(dc,white_col="one_vs_rest",
                       cohort="exploratory_9_colour_vs_rest")
        r["hue"]=hue
        colour_est.append(r)
    qs=multipletests([r["p_genus_cluster"] for r in colour_est],
                    method="fdr_bh")[1]
    for r,q in zip(colour_est,qs):r["BH_q_nine_hues"]=float(q)
    holdout=genus_folds_oof(d)
    output.mkdir(parents=True)
    table=output/"all_colour_intro_propensity_9class_OLS.csv"
    pd.DataFrame(colour_est).to_csv(table,index=False)
    rfile=output/"species_intro_rate_by_photo_white.csv"
    rates.to_csv(rfile,index=False)
    model=output/"source_verified_introduction_models.json"
    report={
        "study":"global WCVP TDWG3 introduced status x accepted-species photo hue",
        "source":source,"unadjusted_white_binary_table":rates.to_dict("records"),
        "primary_white_status_models":estimates,
        "exploratory_nine_hue_models":colour_est,
        "genus_disjoint_prediction":holdout,
        "source_env_zip_SHA256":ENV_SHA,
        "source_WCVP_names_gzip_SHA256":sha(Path(wcvp)/NAME_FILE),
        "source_WCVP_dist_gzip_SHA256":sha(Path(wcvp)/DIST_FILE),
        "native_region_count_adjusted":True,
        "source_has_native_region_requirement":True,
        "WCVP_no_introduced_region_record_not_true_absence_of_introduction":True,
        "WCVP_TDWG3_not_documented_field_establishment":True,
        "photo_hue_not_true_petal_colour_verified":True,
        "species_once_not_repeated_per_GBIF_cell":True,
        "shared_genus_controls_not_full_phylogeny_or_geography":True
    }
    model.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    (output/"receipt.json").write_text(json.dumps({
        **report,"numeric_output_sha256":{p.name:sha(p) for p in (table,rfile,model)}
    },indent=2,allow_nan=False)+"\n")
    print("MOIS_INTRO_FILTER_WITHIN_GENUS_REAL_SUCCESS",json.dumps({
       "source":source,
       "unadjusted":rates.to_dict("records"),
       "white":estimates,
       "colour_FDR":[{"hue":r["hue"],
                      "diff_pp":r["white_photo_to_intro_status_risk_difference_pp"],
                      "p":r["p_genus_cluster"],"q":r["BH_q_nine_hues"]}
                     for r in colour_est],
       "CV":holdout
    },allow_nan=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--env",required=True)
    p.add_argument("--wcvp",required=True)
    p.add_argument("--out",required=True)
    a=p.parse_args()
    run(a.env,a.wcvp,a.out)
