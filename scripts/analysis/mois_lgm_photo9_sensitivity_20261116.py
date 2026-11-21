#!/usr/bin/env python3
"""MOIS LGM-3GCM signed output: photo-hue QA and predeclared sensitivity contrasts.

Primary question: is same-genus association of present-minus-LGM warming with
photographic white-label species robust to species weights, coverage, location
and independent observation hue consensus? NO claim of glacial refugia.
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
from statsmodels.stats.multitest import multipletests

from scripts.analysis.mois_lgm_climate_history_species_20261116 import genus_model

SOURCE_SHA = "db3b5b4fa2133194963be9233411bb6bb50bc86393d1dc50d996cc4cc78b25ab"
SOURCE = "GBIF37104_photo9_LGM_climate_species.parquet"
MEDIAN = ("threeGCM_delta_temperature_C_median",
          "threeGCM_delta_log_precip_median")
HUES = ("white","yellow","orange","red","pink","magenta","purple","blue","bronze")

def require(x,msg):
    if not bool(x):
        raise ValueError("MOIS_LGM_SENSITIVITY: "+msg)

def sha(b):
    if isinstance(b,(Path,str)):
        b=Path(b).read_bytes()
    return hashlib.sha256(b).hexdigest()

def read_signed(path):
    import pyarrow.parquet as pq
    require(sha(path)==SOURCE_SHA,"unaltered signed 3GCM analysis artifact expected")
    with ZipFile(path) as z:
        raw=z.read(SOURCE)
        audit=json.loads(z.read("receipt.json"))
    require(audit["gbif_accepted_species"]==37104
            and audit["species_x_grid_pairs"]==5493931
            and audit["source_valid_GBIF_occurrences"]==396142825
            and sha(raw)==audit["output_SHA256"][SOURCE],
            "source rows and signed parquet SHA drift")
    d=pq.read_table(io.BytesIO(raw)).to_pandas()
    require(len(d)==37104 and d.key.is_unique
            and d.source_hue9.isin(HUES).all(),
            "original accepted-species photo9 frame not preserved")
    return d,audit

def complete_frame(d,cover=.75):
    numeric = [*MEDIAN,"BIO01","BIO12","MODIS_NDVI_2001_2021_mean",
               "ESA_AGB_woody_Mg_ha_2020",
               "mean_abs_latitude_over_occupied_cells",
               "records_passed_spatial_QA","lgm_climate_coverage"]
    usable = (d[numeric].notna().all(axis=1)
              &(d.lgm_climate_coverage>=cover)
              &(d.BIO12>=0)
              &(d.ESA_AGB_woody_Mg_ha_2020>=0)
              &d.genus.notna()
              &d.hue_status.notna())
    return d.loc[usable].copy()

def same_genus_contrast(d, target="white"):
    has_both = d.assign(target=d.source_hue9.eq(target)).groupby(
        "genus", sort=False).target.transform("nunique").eq(2)
    return d.loc[has_both].copy()

def fit_contrast(d, target="white"):
    d=same_genus_contrast(d,target)
    if len(d)<200 or d.genus.nunique()<30:
        return {"status":"insufficient_genus_contrasts",
                "n_species":int(len(d)), "n_genera":int(d.genus.nunique())}
    # The source model tests white; re-code target-vs-others while preserving
    # original hue proportions in the provenance and summary reports.
    if target!="white":
        d=d.copy()
        d["source_hue9"]=np.where(d.source_hue9.eq(target),"white","yellow")
    result={}
    for label,equal in (("species_equal",False),("genus_equal",True)):
        try:
            r=genus_model(d,MEDIAN,equal_genus=equal)
        except (ValueError, np.linalg.LinAlgError) as ex:
            return {"status":"model_failed","reason":str(ex),
                    "n_species":int(len(d)),"n_genera":int(d.genus.nunique())}
        result[label]={
            "n_species":r["n_species"],"n_genera":r["n_genera"],
            "p_joint_historic":r["p_joint_historic"],
            "r_squared_within_genus":r["r_squared_within_genus"],
            "temperature":r["historic_effects"][MEDIAN[0]],
            "precipitation":r["historic_effects"][MEDIAN[1]]
        }
    return {"status":"fitted","model":result}

def sensitivity(d):
    base=complete_frame(d)
    original=base.copy()
    base["original_genus_n"] = base.genus.map(base.genus.value_counts())
    variants={
        "all_75pct_coverage":base,
        "atleast_90pct_LGM_grid_coverage":complete_frame(d,.90),
        "full_LGM_grid_coverage":complete_frame(d,1.0),
        "GBIF_at_least_100_records":base[base.records_passed_spatial_QA>=100],
        "genus_at_least_5_species":base[base.original_genus_n>=5],
        "genus_under_5_species":base[base.original_genus_n<5],
        "low_abs_lat_below_30":base[base.mean_abs_latitude_over_occupied_cells<30],
        "middle_high_abs_lat_atleast_30":base[base.mean_abs_latitude_over_occupied_cells>=30],
        "high_confidence_photo_consensus":base[
           base.hue_status.astype(str).eq("sampling_consensus_single")]
    }
    fits={}
    for label,s in variants.items():
        fits[label]={"sample_size_before_genus_contrast":int(len(s)),
                     "fit":fit_contrast(s)}
        z=fits[label]["fit"]
        print("MOIS_LGM_QA_SENSITIVITY",label,
              "n",z.get("model",{}).get("genus_equal",{}).get("n_species",z.get("n_species")),
              "beta",z.get("model",{}).get("genus_equal",{}).get("temperature",{}).get("beta_per_species_SD"),
              "p",z.get("model",{}).get("genus_equal",{}).get("temperature",{}).get("p"),flush=True)
    # Hues of the photographed proxy: exploratory nine one-vs-rest checks with
    # BH correction for individual temperature coefficients. No claims of
    # historical phenotype change or colour adaptation.
    nine={}
    for hue in HUES:
        nine[hue]=fit_contrast(original,target=hue)
    p=[nine[h]["model"]["genus_equal"]["temperature"]["p"]
       for h in HUES if nine[h]["status"]=="fitted"]
    q=multipletests(p,method="fdr_bh")[1].tolist() if p else []
    for h,qq in zip((h for h in HUES if nine[h]["status"]=="fitted"),q):
        nine[h]["model"]["genus_equal"]["temperature"]["q_BH_nine_hues"]=float(qq)
    return {"method":"one species=one row; within-genus linear-probability; genus-clustered",
            "all_species":len(d),"n_current_complete":len(base),
            "status_counts":d.hue_status.astype(str).value_counts().to_dict(),
            "hue_counts":d.source_hue9.astype(str).value_counts().to_dict(),
            "high_quality_label_is_true_petal_colour":False,
            "current_occupied_cells_are_glacial_refugia":False,
            "sensitivity_fits":fits,"exploratory_nine_hue_fits":nine}

def run(inp,out):
    out=Path(out)
    require(not out.exists(),"refuse overwrite")
    d,audit=read_signed(inp)
    results=sensitivity(d)
    out.mkdir(parents=True)
    (out/"LGM_colour_qa_sensitivity.json").write_text(
        json.dumps(results,indent=2,sort_keys=True,default=str)+"\n")
    summary={
        "source_run":"43868642184",
        "source_artifact_id":"14225733755",
        "source_zip_sha256":SOURCE_SHA,
        "source_run_success":True,
        "n_accepted_species":len(d),
        "n_analysed":results["sensitivity_fits"]["all_75pct_coverage"]["fit"]["model"]["genus_equal"]["n_species"],
        "summary":{k:{
           "status":v["fit"]["status"],
           "genus_equal_beta":v["fit"].get("model",{}).get("genus_equal",{}).get("temperature",{}).get("beta_per_species_SD"),
           "genus_equal_p":v["fit"].get("model",{}).get("genus_equal",{}).get("temperature",{}).get("p"),
           "species_equal_beta":v["fit"].get("model",{}).get("species_equal",{}).get("temperature",{}).get("beta_per_species_SD"),
           "species_equal_p":v["fit"].get("model",{}).get("species_equal",{}).get("temperature",{}).get("p")
        } for k,v in results["sensitivity_fits"].items()},
        "output_sha256":sha(out/"LGM_colour_qa_sensitivity.json")
    }
    (out/"receipt.json").write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n")
    print("MOIS_LGM_QA_COMPLETE",json.dumps(summary,sort_keys=True),flush=True)
    return summary

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--source",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    run(a.source,a.output)
