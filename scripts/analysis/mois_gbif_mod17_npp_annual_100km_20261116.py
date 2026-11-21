#!/usr/bin/env python3
"""Join source-verified real annual MODIS C6.1 NPP (2001-2023) to the frozen
18,577 EPSG:6933 100-km cells occupied by GBIF 37,104 accepted species.

Uses UHH/ICDC 0.5deg CMG quality-screened netCDF regridding of NASA
MOD17A3HGF V061, DOI: 10.25592/uhhfdm.14635, CC-BY-4.0.
This is yearly carbon flux (kg C/m²/year), not biomass or flower allocation.
It is NOT the original 500m spatial resolution and is NOT an LGM estimate.
No photographic hue appears in the NPP sampling or aggregation.
"""
from __future__ import annotations
import argparse, hashlib, io, json, re, tempfile
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pandas as pd

from scripts.analysis.mois_lgm_climate_history_species_20261116 import (
    coords_in_grid, sha, require
)

SOURCE_NPP_SHA = "be8e74e5cf3e4785c94a2799655efef903e3f0c5908af04f6fd114d28dab6c58"
SOURCE_NPP_MD5 = "d9109919b8ca9121c51e90a2ae8f12b6"
SOURCE_LGM_SHA = "db3b5b4fa2133194963be9233411bb6bb50bc86393d1dc50d996cc4cc78b25ab"
SOURCE_CELL_MEMBER = "WorldClim1p4_LGM_threeGCM_18577_cells.parquet"
YEARS = tuple(range(2001,2024))
NPP_DOI = "10.25592/uhhfdm.14635"
NASA_DOI = "10.5067/MODIS/MOD17A3HGF.061"
MIN_YEARS = 18

def read_global_GBIF_cells(lgm_zip):
    import pyarrow.parquet as pq
    require(sha(lgm_zip)==SOURCE_LGM_SHA,
            "true 18577 GBIF-cell signed climate history artifact changed")
    with ZipFile(lgm_zip) as z:
        b=z.read(SOURCE_CELL_MEMBER)
        receipt=json.loads(z.read("receipt.json"))
    require(sha(b)==receipt["output_SHA256"][SOURCE_CELL_MEMBER]
            and receipt["gbif_accepted_species"]==37104
            and receipt["species_x_grid_pairs"]==5493931,
            "LGM source no longer retains actual GBIF universe")
    p=pq.read_table(io.BytesIO(b),columns=["grid_x","grid_y"]).to_pandas()
    require(len(p)==18577 and not p.duplicated().any(),
            "frozen GBIF distinct EPSG:6933 100km cells altered")
    return p

def nearest_05degree_indices(ax, points):
    """Nearest 0.5 degree original CMG grid centers, ascending or descending."""
    axis=np.asarray(ax,float)
    require(axis.ndim==1 and len(axis)>100,"unexpected global CMG coordinate")
    stride=float(axis[1]-axis[0])
    require(abs(abs(stride)-.5)<1e-5 and
            np.allclose(np.diff(axis),stride,atol=1e-5),
            "source long/lat centers no longer regular 0.5 degrees")
    p=np.asarray(points,dtype=float)
    with np.errstate(invalid="ignore"):
        ix=np.rint((p-axis[0])/stride)
    valid=(np.isfinite(ix) & np.isfinite(p)
           &(ix>=0)&(ix<len(axis)))
    out=np.full(p.shape,-1,dtype=np.int32)
    out[valid]=ix[valid].astype(np.int32)
    return out

def nc_sample(npp,valid_500m_n,lat,lon,coords):
    require(npp.shape==valid_500m_n.shape==(360,720),
            "source npp/quality raster shape unexpected")
    require(np.asarray(lat).shape==(360,) and np.asarray(lon).shape==(720,),
            "coordinate axis mismatch")
    x,y=coords
    ix=nearest_05degree_indices(lon,x)
    iy=nearest_05degree_indices(lat,y)
    valid=(ix>=0)&(iy>=0)
    values=np.full(ix.shape,np.nan,dtype=float)
    count_500=np.full(ix.shape,np.nan,dtype=float)
    vix=np.where(valid)
    values[vix]=npp[iy[vix],ix[vix]]
    count_500[vix]=valid_500m_n[iy[vix],ix[vix]]
    # Quality-screened ICDC NPP: numberofvalidvalues>0 is the hard criterion.
    good=valid & np.isfinite(values) & np.isfinite(count_500) & (count_500>0)
    values[~good]=np.nan
    count_500[~good]=np.nan
    n=good.sum(axis=0).astype("int16")
    site_npp=np.divide(np.nansum(values,axis=0),n,
                       out=np.full(n.shape,np.nan),where=n>0)
    qc_count=np.divide(np.nansum(count_500,axis=0),n,
                       out=np.full(n.shape,np.nan),where=n>0)
    return site_npp,n,qc_count

def validate_source_nc(ds,year):
    require(set(("lon","lat","time","npp","numberofvalidvalues"))
            .issubset(ds.variables.keys()),
            f"{year}: original primary-production or QA variables missing")
    lat=np.ma.filled(ds.variables["lat"][:],np.nan)
    lon=np.ma.filled(ds.variables["lon"][:],np.nan)
    require(len(lat)==360 and len(lon)==720 and len(ds.dimensions["time"])==1,
            f"{year}: source netcdf dimensions drift")
    require(str(getattr(ds.variables["npp"],"units","")).replace(" ","")
            .lower()=="kgcm-2",f"{year}: unexpected NPP physical units")
    require(tuple(ds.variables["npp"].shape)==(1,360,720),
            "npp image not one annual global grid")
    return np.asarray(lat,float),np.asarray(lon,float)

def true_annual_NPP_by_GBIF_cell(nppzip,cells):
    import netCDF4
    raw=Path(nppzip)
    require(sha(raw)==SOURCE_NPP_SHA and
            hashlib.md5(raw.read_bytes()).hexdigest()==SOURCE_NPP_MD5,
            "official licensed open global MOD17 NPP source archive SHA/MD5 differs")
    coords=coords_in_grid(cells)
    n=len(cells)
    year_annual=np.full((len(YEARS),n),np.nan,dtype=np.float32)
    year_valid=np.zeros((len(YEARS),n),dtype=np.int16)
    year_qc_mean=np.full((len(YEARS),n),np.nan,dtype=np.float32)
    years_receipt={}
    with ZipFile(raw) as z:
        names=[name for name in z.namelist() if name.lower().endswith(".nc")]
        require(len(names)==23 and set(int(name.split("/")[0]) for name in names)
                ==set(YEARS),"missing or duplicate one of the 23 original annual files")
        for k,year in enumerate(YEARS):
            matches=[name for name in names if name.startswith(str(year)+"/")]
            require(len(matches)==1,f"{year}: ambiguous source member")
            member=matches[0]
            buf=z.read(member)
            require(len(buf)==12974212,f"{year}: original netCDF byte size changed")
            with tempfile.NamedTemporaryFile(suffix=".nc") as t:
                t.write(buf);t.flush()
                with netCDF4.Dataset(t.name,"r") as nc:
                    lat,lon=validate_source_nc(nc,year)
                    v=np.ma.filled(nc.variables["npp"][0,:,:],
                                   np.nan).astype(np.float32)
                    q=np.ma.filled(nc.variables["numberofvalidvalues"][0,:,:]
                                   .astype(np.float32),np.nan).astype(np.float32)
                # Missing_32767 is never considered biological 32767kg C/m².
                annual,valid,quality=nc_sample(v,q,lat,lon,coords)
            year_annual[k]=annual
            year_valid[k]=valid
            year_qc_mean[k]=quality
            years_receipt[str(year)]={
                "original_member":member,
                "netcdf_sha256":sha(buf),
                "n_GBIF_cells_valid":int(np.isfinite(annual).sum()),
                "n_median_land_samples_when_valid":float(
                    np.median(valid[valid>0])) if np.any(valid>0) else 0,
                "annual_NPP_median_kg_C_m2_year":float(np.nanmedian(annual))
            }
            print("MOIS_TRUE_MOD17_NPP_YEAR",year,
                  years_receipt[str(year)]["n_GBIF_cells_valid"],
                  years_receipt[str(year)]["annual_NPP_median_kg_C_m2_year"],
                  flush=True)
    count=np.sum(np.isfinite(year_annual),axis=0).astype(int)
    with np.errstate(invalid="ignore"):
        mean=np.divide(np.nansum(year_annual,axis=0),count,
                       out=np.full(n,np.nan),where=count>0)
        centered=np.where(np.isfinite(year_annual),
                          year_annual-mean[None,:],0)
        variance=np.divide(np.sum(centered**2,axis=0),
                           count-1,out=np.full(n,np.nan),where=count>1)
        sd=np.sqrt(variance)
        stdcv=np.divide(sd,np.abs(mean),
                        out=np.full(n,np.nan),where=np.abs(mean)>1e-6)
    cells=cells.copy()
    cells["NPP_valid_years_2001_2023"]=count
    cells["NPP_2001_2023_mean_kgC_m2_yr"]=mean
    cells["NPP_2001_2023_sd_kgC_m2_yr"]=sd
    cells["NPP_2001_2023_cv"]=stdcv
    cells["NPP_years_ge_18"]=count>=MIN_YEARS
    cells["NPP_2001_kgC_m2_yr"]=year_annual[0]
    cells["NPP_2010_kgC_m2_yr"]=year_annual[9]
    cells["NPP_2020_kgC_m2_yr"]=year_annual[19]
    cells["NPP_2023_kgC_m2_yr"]=year_annual[22]
    cells["NPP_2001_2023_quality_sample_avg"]=np.divide(
        np.nansum(year_valid,axis=0),count,
        out=np.full(n,np.nan),where=count>0)
    # Keep annual time series, not just pooled mean.
    yearly=pd.DataFrame({
        "grid_x":np.tile(cells.grid_x.to_numpy(),len(YEARS)),
        "grid_y":np.tile(cells.grid_y.to_numpy(),len(YEARS)),
        "year":np.repeat(YEARS,n),
        "annual_NPP_kgC_m2":year_annual.reshape(-1),
        "valid_25points":year_valid.reshape(-1),
        "valid_500m_QA_count_mean":year_qc_mean.reshape(-1)
    })
    require(len(yearly)==18577*23 and
            not yearly[["grid_x","grid_y","year"]].duplicated().any(),
            "annual true NPP source has duplicated/synthetic grid year")
    return cells,yearly,years_receipt

def run(lgm,npp,output):
    out=Path(output)
    require(not out.exists(),"refuse overwriting original NPP results")
    base=read_global_GBIF_cells(lgm)
    cell,annual,source_audit=true_annual_NPP_by_GBIF_cell(npp,base)
    good=int(cell.NPP_years_ge_18.sum())
    require(good>9000 and np.isfinite(
            cell.loc[cell.NPP_years_ge_18,
                "NPP_2001_2023_mean_kgC_m2_yr"]).all(),
            "too little real NPP year support on GBIF land")
    out.mkdir(parents=True)
    cellfile=out/"MOD17_C061_NPP_2001_2023_18577_true_GBIF100km_cells.parquet"
    yearfile=out/"MOD17_C061_NPP_23annual_18577_GBIF100km_cells.parquet"
    cell.to_parquet(cellfile,index=False)
    annual.to_parquet(yearfile,index=False)
    info={"schema":"mois_source_signed_mod17_c061_NPP_23years_18577_GBIF_100km_cells_v1",
          "source": {"data_doi":NPP_DOI,"derived_from_NASA_DOI":NASA_DOI,
                    "modis_original_resolution_m":500,
                    "downloaded_global_product_resolution_degrees":0.5,
                    "original_years":list(YEARS),
                    "npp_units":"kg carbon per square meter per year",
                    "quality_criteria":"source 500m confidence best/good, n_valid>0",
                    "source_NPP_zip_SHA256":SOURCE_NPP_SHA,
                    "source_GBIF_climate_zip_SHA256":SOURCE_LGM_SHA},
          "counts":{"source_accepted_GBIF_species":37104,
                    "source_GBIF_species_cell_pairs":5493931,
                    "source_100km_unique_cells":18577,
                    "source_years":len(YEARS),
                    "cell_years_output":len(annual),
                    "cells_with_18plus_valid_annual_NPP":good,
                    "cells_with_any_NPP":int((cell.NPP_valid_years_2001_2023>0).sum()),
                    "cells_all_23years":int((cell.NPP_valid_years_2001_2023==23).sum())},
          "annual_sources":source_audit,
          "no_photo_hue_in_NPP_measurement":True,
          "woody_AGB_not_equal_to_NPP":True,
          "NPP_is_not_Flower_Biomass":True,
          "NPP_is_not_LGM_or_refugia":True,
          "source_raw_NPP_zip_not_republished":True,
          "output_SHA256":{p.name:sha(p) for p in (cellfile,yearfile)}}
    (out/"receipt.json").write_text(json.dumps(info,indent=2,sort_keys=True)+"\n")
    print("MOIS_TRUE_GBIF_MOD17_23YEAR_NPP_COMPLETE",json.dumps({
        **info["counts"],
        "NPP_median_kgC_m2_year":float(cell.loc[cell.NPP_years_ge_18,
                                               "NPP_2001_2023_mean_kgC_m2_yr"].median()),
        "NPP_2001_2023_year_source_2001":source_audit["2001"]["netcdf_sha256"]
    },sort_keys=True),flush=True)
    return info

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--lgm",required=True)
    p.add_argument("--npp",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    run(a.lgm,a.npp,a.output)
