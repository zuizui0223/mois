#!/usr/bin/env python3
"""Descriptive woody biomass × photographed hue within genera across GBIF recorded ranges.

Tests if PHOTO white-label species live in recorded areas with different
SAR-inferred surrounding woody aboveground biomass; NOT that individual white
plants have more or less biomass. Paired genera, numeric climate covariates.
"""
from __future__ import annotations
import argparse,hashlib,io,json
from pathlib import Path
from zipfile import ZipFile
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import spearmanr

SOURCE_SHA="d8b97528188d9f8d63cb802d5b3f6e8f7f251937fd30af42aeae47e43ebad1fc"
FILE="BIO19_weather_NDVI_ESACCIv7_AGB_37104_accepted_species.parquet"

def require(ok,why):
    if not bool(ok):raise ValueError("MOIS_AGB_HUE_WITHIN_GENUS: "+why)

def read(path):
    import pyarrow.parquet as pq
    b=Path(path).read_bytes()
    require(hashlib.sha256(b).hexdigest()==SOURCE_SHA,"biomass photo hue ZIP digest changed")
    with ZipFile(io.BytesIO(b)) as z:
        receipt=json.loads(z.read("receipt.json"))
        p=z.read(FILE)
    require(hashlib.sha256(p).hexdigest()==receipt["output_SHA256"][FILE],
            "source biomass Parquet drift")
    d=pq.read_table(io.BytesIO(p)).to_pandas()
    require(len(d)==37104 and d.key.is_unique
        and d.records_passed_spatial_QA.sum()==396142825,
        "source accepted species roster incomplete")
    return d

def conditional_genus(d,resp,equal=False):
    c=pd.DataFrame(index=d.index)
    c["white_photo"]=d.source_hue9.eq("white").astype(float)
    c["log_GBIF_occ"]=np.log1p(d.records_passed_spatial_QA)
    c["abs_lat"]=d.mean_abs_latitude_over_occupied_cells
    c["temp_BIO01"]=d.BIO01
    c["precip_log_BIO12"]=np.log1p(d.BIO12.clip(lower=0))
    c["NDVI_mean"]=d.MODIS_NDVI_2001_2021_mean
    for s in sorted(d.hue_status.unique())[1:]:
        c["status_"+s]=d.hue_status.eq(s).astype(float)
    good=(c.notna().all(axis=1)&d[resp].notna())
    c=c.loc[good];sub=d.loc[good]
    outcome=sub[resp].astype(float)
    c=c-c.groupby(sub.genus).transform("mean")
    outcome=outcome-outcome.groupby(sub.genus).transform("mean")
    if equal:
        weights=(1/sub.groupby("genus").genus.transform("size")).to_numpy(float)
        f=sm.WLS(outcome,c,weights=weights).fit(
            cov_type="cluster",cov_kwds={"groups":sub.genus},use_t=False)
    else:
        f=sm.OLS(outcome,c).fit(
            cov_type="cluster",cov_kwds={"groups":sub.genus},use_t=False)
    ci=f.conf_int().loc["white_photo"]
    return {"beta_white_photo":float(f.params["white_photo"]),
       "95CI":[float(ci.iloc[0]),float(ci.iloc[1])],
       "p_genus_cluster":float(f.pvalues["white_photo"]),
       "n_species":len(sub),"n_genus":int(sub.genus.nunique()),
       "equal_genus_weight":equal}

def run(source,out):
    out=Path(out)
    require(not out.exists(),"refuse overwrite AGB source audit")
    d=read(source)
    d=d.loc[d.occupied_grid_cells.gt(0)
            &d.MODIS_NDVI_2001_2021_mean.notna()
            &d.ESA_AGB_woody_Mg_ha_2020.notna()
            &d.BIO01.notna()&d.BIO12.notna()
            &d.mean_abs_latitude_over_occupied_cells.notna()].copy()
    require(len(d)==36753,"matched 37104 accepted source AGB cohort drift")
    d["white_photo"]=d.source_hue9.eq("white")
    d["log_AGB_2020"]=np.log1p(d.ESA_AGB_woody_Mg_ha_2020)
    d["AGB_forest_cover"]=d.ESA_AGB_woody_positive_area_fraction_2020
    bywhite=d.groupby("white_photo").agg(
        species=("key","size"),
        median_recorded_woody_AGB_Mg_ha=("ESA_AGB_woody_Mg_ha_2020","median"),
        mean_recorded_woody_AGB_Mg_ha=("ESA_AGB_woody_Mg_ha_2020","mean"),
        median_recorded_forest_fraction=("AGB_forest_cover","median"),
        median_recorded_NDVI=("MODIS_NDVI_2001_2021_mean","median")).reset_index()
    byhue=d.groupby("source_hue9").agg(
        species=("key","size"),
        median_recorded_woody_AGB_Mg_ha=("ESA_AGB_woody_Mg_ha_2020","median"),
        median_recorded_NDVI=("MODIS_NDVI_2001_2021_mean","median")).reset_index()
    # Only genera that contain both white and nonwhite PHOTO categories
    mixed=d.groupby("genus").white_photo.nunique().eq(2)
    matched=d.loc[d.genus.map(mixed)].copy()
    require(len(matched)>27000 and matched.genus.nunique()>2000,
            "matching genera from global source unexpectedly low")
    tests={
       "log_2020_woody_AGB_all_genus_weighted":conditional_genus(matched,"log_AGB_2020"),
       "log_2020_woody_AGB_equal_genus_weight":conditional_genus(matched,"log_AGB_2020",equal=True),
       "fraction_woody_AGB_positive_equal_genus_weight":conditional_genus(matched,"AGB_forest_cover",equal=True)
    }
    spearman=spearmanr(d.ESA_AGB_woody_Mg_ha_2020,
                       d.MODIS_NDVI_2001_2021_mean)
    out.mkdir(parents=True)
    bywhite.to_csv(out/"AGB_2020_white_photo_vs_colour_descriptive.csv",index=False)
    byhue.to_csv(out/"AGB_2020_real_NDVI_nine_photo_hues_descriptive.csv",index=False)
    receipt={
      "schema":"mois_37104_accepted_photo9_woody_AGB_2020_genus_block_conditional_v1",
      "accepted_species_universe":37104,
      "species_with_all_real_AGB_NDVI_climate":len(d),
      "genera_white_and_nonwhite_comparison":int(matched.genus.nunique()),
      "species_within_mixed_genera":len(matched),
      "AGB_2020_vs_NDVI_Spearman_rho_species":float(spearman.statistic),
      "AGB_2020_vs_NDVI_Spearman_p_species":float(spearman.pvalue),
      "conditional_source_models":tests,
      "source_biomass_artifact_SHA256":SOURCE_SHA,
      "biomass_measures_surrounding_woody_vegetation_Mg_ha_not_flower_size_or_species_individual_mass":True,
      "flower_colour_label_not_expert_petal_verified":True,
      "GBIF_ranges_native_only":False,
      "outcome_is_not_flower_petal_biomass":True,
      "causal_interpretation_not_supported":True,
      "output_sha256":{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in out.iterdir()}
    }
    (out/"receipt.json").write_text(json.dumps(receipt,sort_keys=True,indent=2)+"\n")
    print("MOIS_ESACCI_WOODY_BIOMASS_37104_GLOBAL_HUE_WITHIN_GENUS_ASSOCIATION",
          json.dumps({"descriptive_white_vs_coloured":bywhite.to_dict("records"),
                      "by_hue":byhue.to_dict("records"),
                      "conditional_source_models":tests,
                      "NDVI_spearman":float(spearman.statistic)},
                     sort_keys=True),flush=True)
    return receipt

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--biomass",required=True);p.add_argument("--output",required=True)
    a=p.parse_args();run(a.biomass,a.output)
