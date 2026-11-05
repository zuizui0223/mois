#!/usr/bin/env python3
"""Signed 2020 and 2024 spatial ALAN-white correspondence; NO change regression.

Population adjustment for 2024 means 2020 BASELINE density only.
MOIS white is continuous strict9 photographed taxon-equal ROI, not species loss.
"""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

SHA_2020="3562db1e5783e77d3bb408bd33817b4d7cccf31ad4d4ad42586469708fd407e3"
SHA_2024="1c5f03e5a7c355a0f7045133140c5cf2e39a5802c548f89d1c2f66c8a235db5f"
MODEL="white_pp ~ night + population + bio5 + log_rain + elev + photo_effort + C(block_id)"

def check(x,message):
    if not bool(x):raise ValueError("MOIS_VIIRS_SPATIAL: "+message)

def digest(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for part in iter(lambda:f.read(4<<20),b""):h.update(part)
    return h.hexdigest()

def source(earlier,later):
    check(digest(earlier)==SHA_2020 and digest(later)==SHA_2024,"Frozen source SHA drift")
    a=pd.read_csv(earlier,dtype={"cell_id":str,"block_id":str})
    b=pd.read_csv(later,dtype={"cell_id":str,"block_id":str})
    check(len(a)==1059 and a.source_complete.sum()==1044,"2020 original source changed")
    check(len(b)==302 and b.VIIRS_2024_n_valid_months.ge(3).all(),"2024 source changed")
    check(a.cell_id.is_unique and b.cell_id.is_unique,"Duplicate original geographic cell")
    a=a.loc[a.source_complete].copy()
    check(b.cell_id.isin(a.cell_id).all(),"2024 cell not found in 2020 geographic source")
    paired=b.merge(a[["cell_id","white_taxon_equal","bio5","bio12","elev"]],
                   on="cell_id",how="left",validate="one_to_one")
    paired["white2024"]=paired.white_taxon_equal+paired.delta_white_pp/100
    check(paired.white2024.between(0,1).all(),"Invalid reconstructed 2024 original white ROI")
    check(np.allclose(paired.strict2024/paired.attempts2024,
      b.n_strict9/b.n_attempted+b.delta_ROI_yield_pp/100,atol=1e-8),
      "2024 ROI success fraction inconsistent with source receipt")
    common=a.loc[a.cell_id.isin(b.cell_id)].copy()
    check(len(common)==302,"Common-cell support changed")
    return a,common,paired

def get_frame(x,response,light,attempts):
    p=x.copy()
    check(p.cell_id.is_unique and p.block_id.notna().all(),"Cell IDs/blocks invalid")
    check(p[response].between(0,1).all(),"ROI fraction outside [0,1]")
    check(p[attempts].ge(5).all(),"Photo opportunity cohort changed")
    p["white_pp"]=100*p[response]
    p["night"]=np.log1p(p[light].clip(lower=0))
    p["population"]=np.log1p(p.population_density.clip(lower=0))
    p["log_rain"]=np.log1p(p.bio12.clip(lower=0))
    p["photo_effort"]=np.log1p(p[attempts])
    check(np.isfinite(p[["white_pp","night","population","bio5",
                         "log_rain","elev","photo_effort"]].to_numpy(float)).all(),
          "Nonfinite predictor or outcome")
    return p

def fit(p,name,year,attempts):
    check(len(p)>=200 and p.block_id.nunique()>=40,"Insufficient geographic support")
    m=smf.ols(MODEL,data=p).fit(cov_type="cluster",
                               cov_kwds={"groups":p.block_id},use_t=True)
    ci=m.conf_int().loc["night"]
    return {"model":name,"photo_year":year,
            "baseline_population_year":2020,"n_cells":len(p),
            "n_geographic_blocks":int(p.block_id.nunique()),
            "original_photo_opportunities":int(p[attempts].sum()),
            "white_component_mean_percent":float(p.white_pp.mean()),
            "median_upward_VIIRS_radiance_nW_cm2_sr":float(np.expm1(p.night).median()),
            "night_white_coefficient_pp_per_log1p":float(m.params["night"]),
            "nominal_95ci_low_pp":float(ci.iloc[0]),
            "nominal_95ci_high_pp":float(ci.iloc[1]),
            "nominal_p":float(m.pvalues["night"]),
            "night_population_log_correlation":float(p[["night","population"]].corr().iloc[0,1]),
            "descriptive_model_R2":float(m.rsquared),
            "year2024_population_change_controlled":False}

def run(path2020,path2024,out):
    out=Path(out)
    check(not out.exists(),"Refuse to overwrite prior results")
    a,b,c=source(path2020,path2024)
    datasets=[
        ("2020_all",get_frame(a,"white_taxon_equal","viirs_quarterly_median","n_attempted"),2020,"n_attempted"),
        ("2020_common302",get_frame(b,"white_taxon_equal","viirs_quarterly_median","n_attempted"),2020,"n_attempted"),
        ("2024_common302",get_frame(c,"white2024","VIIRS_2024_quarter_median","attempts2024"),2024,"attempts2024"),
    ]
    fitted=[fit(x,name,year,effort) for name,x,year,effort in datasets]
    out.mkdir(parents=True)
    target=out/"VIIRS_white_spatial_correspondence_by_year.csv"
    pd.DataFrame(fitted).to_csv(target,index=False,float_format="%.15g")
    receipt={"schema":"mois_2020_2024_VIIRS_spatial_correspondence_NOT_change_v1",
             "2020_input_SHA256":SHA_2020,"2024_input_SHA256":SHA_2024,
             "analysis":"cross_sectional_level_not_temporal_change",
             "2020_population_density_used_for_2024_baseline_only":True,
             "population2024_measured":False,
             "2024_population_density_or_urbanization_change_measured":False,
             "photo_and_observation_from_frozen_source_only":True,
             "unknown_ROI9_not_white_nonwhite":True,
             "white_species_extinction_observed":False,
             "moth_pollination_trait_assigned":False,
             "no_new_GIS_requests_or_original_image_download":True,
             "MNAR_sensitivity_run":False,"formula":MODEL,
             "result_SHA256":digest(target),"models":len(fitted)}
    (out/"spatial_correspondence_receipt.json").write_text(
        json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print("MOIS_ALAN_SPATIAL_CORRESPONDENCE_EXECUTED",json.dumps(
        [{"n":r["n_cells"],"year":r["model"],"coefficient":r["night_white_coefficient_pp_per_log1p"],
          "interval":[r["nominal_95ci_low_pp"],r["nominal_95ci_high_pp"]]} for r in fitted]),flush=True)
    return fitted,receipt

if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source-2020",required=True,type=Path)
    ap.add_argument("--source-2024",required=True,type=Path)
    ap.add_argument("--out",required=True,type=Path)
    args=ap.parse_args()
    run(args.source_2020,args.source_2024,args.out)
