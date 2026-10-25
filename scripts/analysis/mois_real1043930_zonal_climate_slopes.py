#!/usr/bin/env python3
"""Genuine 390,541 source-photo ROI9 conditional climate effects by biogeography.

Univariate within acquisition source arm × named nominal genus × 5° site,
taxon-geocell-equal slope, bootstrapped by GENUS (nominal CI), of annual water
vapour, annual rain and DEM elevation on white fraction, white dominant
probability, mixed-colour entropy and bronze palette fraction. Do separately
for north/south, tropical/extratropical, polar and latitude subregions.
Heterogeneous source taxon/observer composition and ROI success can still
confound; no physiological anthocyanin or confirmed green/brown petals.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd

from scripts.analysis.mois_real1043930_global_zonal_colour import (
    main_colour_fraction_table,sha,check
)
SCHEMA="mois_real1043930_white_colour_hemisphere_latitude_regional_abiotic_slopes_v1"
FEATURES={
    "atmospheric_vapour":"clim_vapour",
    "log_annual_precipitation":"clim_rain",
    "elevation":"clim_elevation",
}
OUTCOMES={
    "white_fraction":"white_fraction",
    "major_white_probability":"dominant_white",
    "nine_palette_entropy":"entropy9",
    "bronze_fraction_NOT_brown_petals":"bronze_fraction",
}
PAIRS=[
    "north_all","south_all",
    "north_extratropics","south_extratropics",
    "north_tropical","south_tropical",
    "north_temperate","south_temperate",
    "north_subpolar","south_subpolar",
    "north_arctic","south_antarctic",
]
MIN_CELLS=400
MIN_GENUS=50
DRAWS=249
def sector(lat,name):
    north=lat>=0
    abs_lat=np.abs(lat)
    if name=="north_all":return north
    if name=="south_all":return lat<0
    if name=="north_extratropics":return north&(abs_lat>=23.5)
    if name=="south_extratropics":return (lat<0)&(abs_lat>=23.5)
    if name=="north_tropical":return north&(abs_lat<23.5)
    if name=="south_tropical":return (lat<0)&(abs_lat<23.5)
    if name=="north_temperate":return north&(abs_lat>=23.5)&(abs_lat<55)
    if name=="south_temperate":return (lat<0)&(abs_lat>=23.5)&(abs_lat<55)
    if name=="north_subpolar":return north&(abs_lat>=55)&(abs_lat<66.5622)
    if name=="south_subpolar":return (lat<0)&(abs_lat>=55)&(abs_lat<66.5622)
    if name=="north_arctic":return lat>=66.5622
    if name=="south_antarctic":return lat<=-66.5622
    raise ValueError("Undefined biogeographical region")

def gather(p):
    x=p.loc[p.strict & p.latitude.between(-90,90) &
            p.longitude.between(-180,180)].copy()
    x=x[x.taxon.str.fullmatch(r"\d+")&
        x.genus.str.fullmatch("[A-Z][A-Za-z-]+",na=False)]
    x["clim_rain"]=np.log1p(x.clim_rain.clip(lower=0))
    for col in FEATURES.values():
        x=x[np.isfinite(x[col])].copy()
    x["lat5"]=np.floor((x.latitude+90)/5).astype("int16")
    x["lon5"]=np.floor((x.longitude+180)/5).astype("int16")
    return x

def compute_region(x,zone,scales,draws=DRAWS):
    x=x[sector(x.latitude.to_numpy(dtype=float),zone)].copy()
    n_photo=len(x)
    support={"zone":zone,"strict_roi9_photo_count":n_photo,
             "distinct_taxa":int(x.taxon.nunique()),"genera":int(x.genus.nunique())}
    group=["source_cohort","genus","lat5","lon5","taxon"]
    if not n_photo:
        return [{**support,"outcome":o,"predictor":f,"slope_pp_per_global_source_SD":None,
                 "reason":"no_actual_original_photo_ROI_in_geographic_zone"}
                for o in OUTCOMES for f in FEATURES]
    agg={f:(col,"mean") for f,col in FEATURES.items()}
    agg.update({o:(col,"mean") for o,col in OUTCOMES.items()})
    t=x.groupby(group,observed=True).agg(**agg,
       n_photos=("photo_id","size")).reset_index()
    inner=group[:-1]
    num=t.groupby(inner,observed=True).taxon.transform("size")
    q=t.loc[num>=2].copy()
    support.update({"eligible_taxon_geocells":len(q),
                    "eligible_source_photos":int(q.n_photos.sum()),
                    "eligible_genus_clusters":int(q.genus.nunique()),
                    "eligible_nested_cohort_genus_5deg_groups":
                      q.groupby(inner,observed=True).ngroups})
    if len(q)<MIN_CELLS or q.genus.nunique()<MIN_GENUS:
        return [{**support,"outcome":o,"predictor":f,
                 "slope_pp_per_global_source_SD":None,
                 "reason":"too_little_repeat_taxon_within_5deg_genus_cohort_support"}
                for o in OUTCOMES for f in FEATURES]
    centered={}
    for c in [*OUTCOMES,*FEATURES]:
        centered[c]=(q[c]-q.groupby(inner,observed=True)[c].transform("mean")).to_numpy(float)
    out=[]
    for i,(outcome,feat) in enumerate(
        [(o,f) for o in OUTCOMES for f in FEATURES]):
        dx=centered[feat]/scales[feat]
        dy=centered[outcome]
        prod=pd.DataFrame({
            "genus":q.genus.to_numpy(dtype=str),
            "xy":dx*dy,
            "xx":dx*dx,
        }).groupby("genus",observed=True)[["xy","xx"]].sum()
        agg2=prod.to_numpy(dtype=float)
        if agg2[:,1].sum()<=1e-12:
            out.append({**support,"outcome":outcome,"predictor":feat,
                "slope_pp_per_global_source_SD":None,
                "reason":"predictor_constant_within_genus_cells"})
            continue
        rng=np.random.default_rng(20261025+sum(map(ord,zone))+i*137)
        idx=rng.integers(0,len(agg2),size=(draws,len(agg2)))
        sample=agg2[idx].sum(axis=1)
        ok=sample[:,1]>1e-12
        sim=100*sample[ok,0]/sample[ok,1]
        boot=[float(z) for z in np.quantile(sim,[.025,.975])] if len(sim)>10 else [None,None]
        b=100*float(agg2[:,0].sum()/agg2[:,1].sum())
        out.append({
            **support,"outcome":outcome,"predictor":feat,
            "slope_pp_per_global_source_SD":b,
            "nominal_genus_cluster_95_pp":boot,
            "predictor_global_original_ROI_photo_SD":scales[feat],
            "genus_cluster_bootstraps":int(len(sim)),
            "univariate_within_cohort_genus_5degree_NOT_independent_causality":True,
        })
    return out

def run(photo,receipt,climate,geo,geo_receipt,out):
    out=Path(out)
    check(not out.exists(),"Refuse overwrite of all-zones source climate colour analysis")
    p,pins=main_colour_fraction_table(photo,receipt,climate,geo,geo_receipt)
    x=gather(p)
    check(len(x)>=380000 and int(p.strict.sum())==390541,
          "Genuine nine palette ecological/photo-site overlap insufficient")
    scales={key:float(x[name].std(ddof=0)) for key,name in FEATURES.items()}
    check(all(np.isfinite(s) and s>0 for s in scales.values()),
          "Predictor source global standard deviation degenerates")
    results=[]
    for label in PAIRS:
        for v in compute_region(x,label,scales):
            results.append(v)
    check(len(results)==len(PAIRS)*len(OUTCOMES)*len(FEATURES),
          "Cannot silently omit hemisphere / Arctic / Antarctic or a colour outcome")
    report={
       "schema":SCHEMA,
       "original_pixel_source_attempts":1043930,
       "original_strict9_source_assays":390541,
       "strict_photo_ROI_eligible_for_all3climates_and_named_genus":len(x),
       "humidity_temperature_covariate_correlation_not_independently_resolved":True,
       "all_sectors_including_polar_even_when_insufficient":list(PAIRS),
       "outcomes_including_white_first_and_dominant_colour":list(OUTCOMES),
       "zonal_genus_geocell_cohort_source_slopes":results,
       "full_photo_source_SD":scales,
       "method":"source-arm×genus×5deg/nominal-taxon-equal univariate slope; 249 genus bootstrap; source site environmental SD globally",
       "inference_is_source_photo_ROI_composition_not_petal_plant_flora_or_abiotic_causation":True,
       "source9_SHA":pins["input_9roi_sha256"],
       "source_climate_SHA":pins["input_climate_sha256"],
       "source_GSHHG_SHA":pins["input_landmass_sha256"],
       "RGB_photo_downloads":0,
    }
    out.mkdir(parents=True)
    tab=out/"mois1043930_zonal_white_dominant_genus_local_slopes.csv"
    pd.DataFrame(results).to_csv(tab,index=False)
    (out/"mois1043930_zonal_white_dominant_climate_effects.json").write_text(
      json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+"\n")
    previews=[
      {k:v for k,v in row.items() if k in (
          "zone","outcome","predictor","slope_pp_per_global_source_SD",
          "nominal_genus_cluster_95_pp","eligible_taxon_geocells",
          "eligible_genus_clusters","strict_roi9_photo_count","reason")}
      for row in results if row["outcome"] in ("white_fraction","major_white_probability",
          "nine_palette_entropy") and
          row["predictor"] in ("atmospheric_vapour","elevation") and
          row["zone"] in ("north_extratropics","south_extratropics",
                         "north_arctic","south_antarctic","north_tropical","south_tropical")
    ]
    print("MOIS_LATEST1043930_TRUE_NORTH_SOUTH_POLAR_WET_DRY_ELEV_WITHIN_GENUS_SLOPES",
          json.dumps({"strict_source_climate":len(x),"major":previews,
                      "rows":len(results)},sort_keys=True,allow_nan=False),flush=True)
    return report

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("photo","receipt","climate","geo","geo-receipt","output"):
        p.add_argument("--"+key,required=True,type=Path)
    z=p.parse_args()
    run(z.photo,z.receipt,z.climate,z.geo,z.geo_receipt,z.output)

if __name__=="__main__":main()
