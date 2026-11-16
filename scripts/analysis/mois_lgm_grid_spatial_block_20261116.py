#!/usr/bin/env python3
"""Spatial-block independent-unit check of MOIS observed photo hue assemblages.

Inputs are signed sources: original GBIF 37k species x 100km cells, signed
three-GCM LGM climate cells/species, current WorldClim 1.4 10m, ESA CCI
woody above-ground biomass by GBIF cell. Species modal photo hues are imputed
onto their recorded distribution cells, NOT flower hues locally observed.
"""
from __future__ import annotations
import argparse,hashlib,io,json,tempfile
from pathlib import Path
from zipfile import ZipFile
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scripts.analysis.mois_lgm_climate_history_species_20261116 import (
    read_signed_sources, extract_two_bioclim, read_band, extract_at_points,
    coords_in_grid, sha, require, GBIF_SHA, ENV_SHA
)

LGM_ZIP_SHA="db3b5b4fa2133194963be9233411bb6bb50bc86393d1dc50d996cc4cc78b25ab"
CURRENT_ZIP_SHA="f289ff21c6c4bd24f0b7635dee9a57b0c93193d7194ae9e86b6ba8bcb442e6f1"
CELL_LGM="WorldClim1p4_LGM_threeGCM_18577_cells.parquet"
CELL_AGB="ESACCI_AGB_2010_2020_2024_18577_GBIF_cells.parquet"

def load_inputs(gbif_zip,lgm_zip,env_zip,current_zip):
    import pyarrow.parquet as pq
    # Signed first-run source contains true 5,493,931 species x cells.
    pairs,ignored=read_signed_sources(gbif_zip,env_zip)
    require(sha(lgm_zip)==LGM_ZIP_SHA,"original LGM output ZIP changed")
    with ZipFile(lgm_zip) as z:
        s=json.loads(z.read("receipt.json"))
        b=z.read(CELL_LGM)
    require(sha(b)==s["output_SHA256"][CELL_LGM]
            and s["n_current_occupied_global_cells"]==18577,
            "true LGM 100km cell source signature missing")
    climate=pq.read_table(io.BytesIO(b)).to_pandas()
    with ZipFile(env_zip) as z:
        s=json.loads(z.read("receipt.json"))
        b=z.read(CELL_AGB)
    require(sha(b)==s["output_SHA256"][CELL_AGB],"AGB recorded cell signature missing")
    agb=pq.read_table(io.BytesIO(b)).to_pandas()
    require(sha(current_zip)==CURRENT_ZIP_SHA, "WorldClim current baseline changed")
    require(len(agb)==18577 and len(climate)==18577
            and not agb[["grid_x","grid_y"]].duplicated().any()
            and not climate[["grid_x","grid_y"]].duplicated().any(),
            "GBIF spatial cell universe changed")
    return pairs,ignored,climate,agb

def present_climate_cells(climate,archive):
    x=climate[["grid_x","grid_y"]].copy()
    coords=coords_in_grid(x)
    with tempfile.TemporaryDirectory() as td:
        layers=extract_two_bioclim(archive,td,"CURRENT")
        for k,name in ((1,"present_temperature_C"),(12,"present_precip_mm")):
            vals,affine=read_band(layers[k][0],k,"CURRENT")
            pts=extract_at_points(vals,coords,affine)
            valid=np.isfinite(pts)
            cnt=valid.sum(axis=0)
            means=np.divide(np.nansum(pts,axis=0),cnt,
                            out=np.full(len(x),np.nan),where=cnt>0)
            x[name]=means
            x[f"present_BIO{k}_valid_samples"]=cnt
    return x

def photo_assemblage_cells(pairs,photo_species):
    hues=photo_species[["key","source_hue9"]].copy()
    hues["key"]=hues.key.astype(str)
    require(hues.key.is_unique,"duplicate source taxon ID")
    require(not pairs[["specieskey","grid_x","grid_y"]].duplicated().any(),
            "duplicate species-grid pair")
    q=pairs.merge(hues,left_on="specieskey",right_on="key",
                  how="left",validate="many_to_one")
    require(len(q)==5493931 and q.source_hue9.notna().all(),
            "all species geographic rows must have photo9 proxy")
    q["white"]=q.source_hue9.eq("white").astype(np.int8)
    s=q.groupby(["grid_x","grid_y"],sort=False).agg(
        photo_label_species_n=("white","size"),
        white_photo_label_species_n=("white","sum"),
        GBIF_occurrence_count=("occurrence_count","sum")).reset_index()
    s["observed_photo_label_white_fraction"]=(s.white_photo_label_species_n
                                              /s.photo_label_species_n)
    require(s.photo_label_species_n.sum()==5493931
            and s.GBIF_occurrence_count.sum()==396142825,
            "source occurrence loss during photo assemblage cell summary")
    return s

def select_cells(pairs,photo_species,climate,agb,present,min_photo_species=30):
    s=photo_assemblage_cells(pairs,photo_species)
    for tab in (climate,agb,present):
        s=s.merge(tab,on=["grid_x","grid_y"],how="left",validate="one_to_one")
    require(len(s)==18577,"invalid cell join")
    from pyproj import Transformer
    t=Transformer.from_crs("EPSG:6933","EPSG:4326",always_xy=True)
    lon,lat=t.transform((s.grid_x.to_numpy(float)+.5)*100000,
                        (s.grid_y.to_numpy(float)+.5)*100000)
    s["lon"]=lon
    s["lat"]=lat
    s["abs_lat"]=np.abs(lat)
    s["block_500km"]=s.grid_x.floordiv(5).astype(str)+"|"+s.grid_y.floordiv(5).astype(str)
    s["block_1000km"]=s.grid_x.floordiv(10).astype(str)+"|"+s.grid_y.floordiv(10).astype(str)
    columns=["threeGCM_delta_temperature_C_median",
             "threeGCM_delta_log_precip_median",
             "present_temperature_C","present_precip_mm",
             "ESA_AGB_woody_Mg_ha_2020","abs_lat","lon","lat",
             "observed_photo_label_white_fraction"]
    good=(s[columns].notna().all(axis=1)
          & s.threeGCM_valid_all_models
          & (s.photo_label_species_n>=min_photo_species)
          & (s.present_precip_mm>=0)
          & (s.ESA_AGB_woody_Mg_ha_2020>=0)
          & (s.photo_label_species_n>0))
    return s,s.loc[good].copy()

def features(d):
    z=pd.DataFrame(index=d.index)
    z["LGM_warming"]=d.threeGCM_delta_temperature_C_median.astype(float)
    z["LGM_precip_delta"]=d.threeGCM_delta_log_precip_median.astype(float)
    z["present_temp"]=d.present_temperature_C.astype(float)
    z["log_present_precip"]=np.log1p(d.present_precip_mm.astype(float))
    z["log_woody_AGB"]=np.log1p(d.ESA_AGB_woody_Mg_ha_2020.astype(float))
    z["abs_lat"]=d.abs_lat.astype(float)
    z["log_species_n"]=np.log1p(d.photo_label_species_n.astype(float))
    z["log_GBIF_occ"]=np.log1p(d.GBIF_occurrence_count.astype(float))
    # Standardize globally once for comparable slopes across fixed-effects
    # checks, and not after blocking. Report one original global SD.
    sd=z.std()
    require((sd>0).all() and np.isfinite(sd).all(),"nonvarying predictor")
    mean=z.mean()
    z=(z-mean)/sd
    return z,{k:float(sd[k]) for k in sd.index}

def fit_cells(d,block=None,where=None):
    t=d if where is None else d.loc[where].copy()
    if len(t)<350 or t.block_500km.nunique()<25:
        return {"status":"insufficient_cells","n_cells":len(t)}
    X,sd=features(t)
    y=t.observed_photo_label_white_fraction.astype(float)
    cluster=t.block_500km.astype(str)
    if block is not None:
        B=t[block].astype(str)
        X=X-X.groupby(B).transform("mean")
        y=y-y.groupby(B).transform("mean")
        keep=(X.abs().sum(axis=1)>0)
        X=X.loc[keep];y=y.loc[keep];cluster=cluster.loc[keep]
        if len(y)<200 or B.nunique()<20:
            return {"status":"underpowered_after_block_decomposition"}
        model=sm.OLS(y,X).fit(cov_type="cluster",cov_kwds={"groups":cluster},use_t=False)
    else:
        # Global level model: quadratic location controls, finite-sample
        # covariance grouped by 500km geographic blocks, not individual cells.
        X["lat2"]=(t.lat.astype(float)**2)
        X["lon2"]=(t.lon.astype(float)**2)
        for col in ("lat2","lon2"):
            X[col]=(X[col]-X[col].mean())/X[col].std()
        X=sm.add_constant(X,has_constant="add")
        model=sm.OLS(y,X).fit(cov_type="cluster",cov_kwds={"groups":cluster},use_t=False)
    k="LGM_warming"
    lo,hi=model.conf_int().loc[k].tolist()
    return {"status":"fitted","n_cells":int(len(y)),
            "blocks_500km":int(cluster.nunique()),
            "white_fraction_mean":float(y.mean()) if block is None else None,
            "beta_white_fraction_per_SD_LGM_warming":float(model.params[k]),
            "95CI":[float(lo),float(hi)],
            "p":float(model.pvalues[k]),
            "R2":float(model.rsquared),
            "LGM_warming_SD_C":sd[k],
            "block_500km_cluster_errors":True,
            "fixed_block":block,
            "outcome_is_incidence_of_photo_labeled_species_not_field_flower_color":True}

def run(gbif,lgm,env,current,out):
    out=Path(out)
    require(not out.exists(),"refuse overwrite")
    pairs,photo_species,lclim,agb=load_inputs(gbif,lgm,env,current)
    present=present_climate_cells(lclim,current)
    allcells,d=select_cells(pairs,photo_species,lclim,agb,present,30)
    result={"n_all_unique_global_GBIF_cells":len(allcells),
            "n_cells_with_enough_photo_label_taxa_and_all_environment":len(d),
            "n_GBIF_species_x_cell_pairs":len(pairs),
            "GBIF_occurrences":396142825,
            "present_source_period":"1960-1990",
            "LGM_source_epoch":"about 22ka",
            "colour_outcome":"species modal photographic hue imputed across current GBIF occupied cells",
            "climate_history_is_not_refugia":True}
    result["models"]={
      "global_quadratic_lat_lon_500km_cluster":fit_cells(d),
      "500km_within_block":fit_cells(d,"block_500km"),
      "1000km_within_block":fit_cells(d,"block_1000km"),
      "abs_lat_less_than_30":fit_cells(d,where=d.abs_lat<30),
      "abs_lat_atleast_30":fit_cells(d,where=d.abs_lat>=30)}
    out.mkdir(parents=True)
    col=["grid_x","grid_y","observed_photo_label_white_fraction",
         "white_photo_label_species_n","photo_label_species_n","block_500km",
         "block_1000km","present_temperature_C","present_precip_mm",
         "threeGCM_delta_temperature_C_median",
         "threeGCM_delta_log_precip_median","ESA_AGB_woody_Mg_ha_2020","lon","lat"]
    d[col].to_parquet(out/"LGM_white_photo9_species_assemblages_100km_cells.parquet",
                      index=False)
    (out/"spatial_block_models.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    r={
     "source_sha256":{"GBIF":GBIF_SHA,"LGM_3GCM":LGM_ZIP_SHA,
                      "current_WorldClim_v1p4":CURRENT_ZIP_SHA,
                      "ESA_biomass_climate":ENV_SHA},
     "source_lgm_run":"43868642184",
     "all_cell_rows":len(allcells),"retained_cell_rows":len(d),
     "results":result,
     "source_hues_are_local_petal_measurements":False,
     "source_species_presence_is_not_native_only":True,
     "source_current_climate_and_LGM_are_same_v1p4_baseline":True,
     "output_sha256":{x.name:sha(x) for x in out.iterdir()}}
    (out/"receipt.json").write_text(json.dumps(r,indent=2,sort_keys=True)+"\n")
    print("MOIS_LGM_TRUE_100KM_SPATIAL_BLOCK",json.dumps(result,sort_keys=True),flush=True)
    return r

if __name__=="__main__":
    p=argparse.ArgumentParser()
    for arg in ("gbif","lgm","env","current","output"):
        p.add_argument("--"+arg,required=True)
    a=p.parse_args()
    run(a.gbif,a.lgm,a.env,a.current,a.output)
