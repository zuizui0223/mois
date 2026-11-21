#!/usr/bin/env python3
"""Source-signed WorldClim coverage and photo-hue certainty sensitivities.

The fitted coef is for GBIF-recorded species climate proxies, NOT
inferred adaptation, native ranges or accurate flower petal color.
"""
from pathlib import Path
from zipfile import ZipFile
import argparse,hashlib,io,json
import numpy as np
import pandas as pd
import statsmodels.api as sm

SOURCE_SHA="82f06f69043a2080a961d5586ab4ecc2342b7d7d739c4ac819147971f4d8e069"
SOURCE_MEMBER="COLXR_37104_species_GBIF_grid_WorldClim_species_climate.csv.gz"

def require(cond,msg):
    if not bool(cond):raise ValueError("MOIS_WC_HUE_SENSITIVITY: "+msg)

def fit_genus_FE(d,outcome):
    d=d.copy()
    y=d[outcome].astype(float)
    X=pd.DataFrame(index=d.index)
    X["white_photo"]=d.source_hue9.eq("white").astype(float)
    X["log_GBIF_occ"]=np.log1p(d.records_passed_spatial_QA)
    X["mean_abs_latitude"]=d.mean_abs_latitude_over_occupied_cells.astype(float)
    if outcome=="mean_annual_temperature_C":
        X["log_BIO12_precip"]=np.log1p(d.annual_precipitation_mm)
    else:
        X["BIO1_temperature_C"]=d.mean_annual_temperature_C
    for s in sorted(d.hue_status.unique())[1:]:
        X["photo_status_"+s]=d.hue_status.eq(s).astype(float)
    g=d.genus
    y=y-y.groupby(g).transform("mean")
    X=X-X.groupby(g).transform("mean")
    m=sm.OLS(y,X).fit(cov_type="cluster",cov_kwds={"groups":g},use_t=False)
    CI=m.conf_int().loc["white_photo"]
    return {"beta":float(m.params["white_photo"]),
            "CI95_low":float(CI.iloc[0]),"CI95_high":float(CI.iloc[1]),
            "p_genus_clustered":float(m.pvalues["white_photo"])}

def run(source,out):
    out=Path(out);require(not out.exists(),"source-derived sensitivity output exists")
    raw=Path(source).read_bytes()
    require(hashlib.sha256(raw).hexdigest()==SOURCE_SHA,"original climate output ZIP SHA drift")
    with ZipFile(io.BytesIO(raw)) as z:
        csv=z.read(SOURCE_MEMBER)
        receipt=json.loads(z.read("receipt.json"))
    require(hashlib.sha256(csv).hexdigest()==receipt["output_SHA256"][SOURCE_MEMBER],
            "signed WorldClim photo9 climate source file changed")
    d=pd.read_csv(io.BytesIO(csv),compression="gzip",dtype={"key":str})
    require(len(d)==37104 and d.key.is_unique and
            d.records_passed_spatial_QA.sum()==396142825,"source source taxa or GBIF denominator changed")
    d["log_annual_precip"]=np.log1p(d.annual_precipitation_mm)
    rows=[]
    for minimum in (0.5,0.75,0.9,1.0):
        for omit_disagree in (False,True):
            s=d.loc[d.occupied_grid_cells.gt(0)
                 &d.climate_coverage_fraction.ge(minimum)
                 &d.mean_abs_latitude_over_occupied_cells.notna()
                 &d.mean_annual_temperature_C.notna()
                 &d.annual_precipitation_mm.notna()].copy()
            if omit_disagree:s=s.loc[~s.modal_mean9_disagree]
            for response in ("mean_annual_temperature_C","log_annual_precip"):
                rows.append({
                  "climate_coverage_threshold":minimum,
                  "omit_species_photo_modal_vs_mean9_conflicts":omit_disagree,
                  "response":response,
                  "taxa":len(s),"genera":s.genus.nunique(),
                  **fit_genus_FE(s,response)})
    result=pd.DataFrame(rows)
    require(len(result)==16 and result.taxa.min()>26000,"unexpected source quality subset")
    out.mkdir(parents=True)
    f=out/"BIO1_BIO12_white_photo_genusFE_coverage_quality_sensitivity.csv"
    result.to_csv(f,index=False)
    meta={"schema":"mois_37104_9photo_worldclim_2endpoint_coverage_sensitivity_v1",
      "signed_original_WorldClim_output_ZIP_SHA256":SOURCE_SHA,
      "source_accepted_COLXR_species":37104,
      "climate_coverage_thresholds":[0.5,0.75,0.9,1.0],
      "photo_modal_mean9_mismatch_excluded_as_sensitivity":True,
      "tests":16,"source_climate_years":"1970-2000",
      "native_status_verified":False,
      "petal_colour_at_GBIF_location_measured":False,
      "genus_cluster_covariance_not_HC3":True,
      "output_CSV_SHA256":hashlib.sha256(f.read_bytes()).hexdigest()}
    (out/"receipt.json").write_text(json.dumps(meta,indent=2,sort_keys=True)+"\n")
    print("MOIS_GBIF_WORLDCLIM_37104_PHOTO9_COVERAGE_HUE_QA_SENSITIVITY",
          result.to_json(orient="records"),flush=True)
    return result

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--source",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args();run(a.source,a.output)
