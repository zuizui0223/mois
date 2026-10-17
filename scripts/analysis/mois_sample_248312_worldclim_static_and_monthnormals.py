#!/usr/bin/env python3
"""MOIS: sample global WorldClim at ALL measured photo sites, separate date-varying normals.

This is geographic, static climate baseline and topography, not a weather
event, independently confirmed flowering, or FCP intraspecific polymorphism.
Original hues/data remain unchanged; missing originals remain missing.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio

ORIGINAL_COUNT = 248312
ORIGINAL_CLASSIFIED = 120807
ORIGINAL_SHA = "750c97390221d29770056c5b0c3b823d8126fc6629f5322c303a14979aadd26d"
# Note: 750c... is the SHA of the original uncompressed/encoded 248312
# photo-table entry within the verified Actions archive, not the outer ZIP.
CLASSIFIED = ("white", "yellow_orange", "red_pink", "blue_purple")
BIO_IDS = tuple(range(1, 20))
SOURCE = {
    "bio": "wc2.1_10m_bio.zip",
    "srad": "wc2.1_10m_srad.zip",
    "prec": "wc2.1_10m_prec.zip",
    "tavg": "wc2.1_10m_tavg.zip",
    "elev": "wc2.1_2.5m_elev.zip",
}
SOURCE_OFFICIAL = {
    name: "https://geodata.ucdavis.edu/climate/worldclim/2_1/base/" + file
    for name, file in SOURCE.items()
}
SOURCE_PREVIOUS_SHA = {
    "bio": "00513224583665ec0f2f955a4ec252730c4deb2004cce9e793492a3f26df4dcf",
    "srad": "c72ee7f4f9a0eb4b5f6cd7a003eddc05bda22a2e5666d968ed4e817fd36b9026",
}
STATIC_SOURCE = (
    "latitude", "longitude", "wc_elevation_m", "soil_pH", "soil_SOC",
    "soil_N", "soil_clay", "soil_available_water_proxy",
    "wc_bio1", "wc_bio5", "wc_bio12", "wc_bio15",
)


def file_sha(path: Path) -> str:
    d = hashlib.sha256()
    with Path(path).open("rb") as f:
        for piece in iter(lambda: f.read(1 << 20), b""):
            d.update(piece)
    return d.hexdigest()


def require(ok: bool, reason: str) -> None:
    if not ok:
        raise ValueError(reason)


def find_zip_tiff(root: Path, zkey: str, name: str) -> Path:
    folder = root / zkey
    matches = sorted(folder.rglob(name))
    require(len(matches) == 1, "WorldClim raster missing or duplicate: "+name)
    return matches[0]


def raster_positions(ds: rasterio.io.DatasetReader,
                     lon: np.ndarray, lat: np.ndarray) -> tuple[np.ndarray,np.ndarray,np.ndarray]:
    require(abs(ds.transform.b) < 1e-12 and abs(ds.transform.d) < 1e-12,
            "Rotated geotransforms require explicit geographic transformer")
    require(ds.crs is not None and ds.crs.to_epsg() == 4326,
            "WorldClim raster is not EPSG:4326")
    rows = np.zeros(len(lon), dtype=np.int64)
    cols = np.zeros(len(lon), dtype=np.int64)
    valid_geo = np.isfinite(lon) & np.isfinite(lat) & (np.abs(lon) <= 180) & (np.abs(lat) <= 90)
    col_float = np.zeros(len(lon))
    row_float = np.zeros(len(lon))
    col_float[valid_geo] = (lon[valid_geo]-ds.transform.c)/ds.transform.a
    row_float[valid_geo] = (lat[valid_geo]-ds.transform.f)/ds.transform.e
    cols[valid_geo] = np.floor(col_float[valid_geo]).astype(np.int64)
    rows[valid_geo] = np.floor(row_float[valid_geo]).astype(np.int64)
    inside = valid_geo & (rows >= 0) & (rows < ds.height) & (cols >= 0) & (cols < ds.width)
    return rows, cols, inside


def sample(arr: np.ndarray, ds: rasterio.io.DatasetReader,
           lon: np.ndarray, lat: np.ndarray,
           *, return_placements: bool = False):
    rows, cols, inside = raster_positions(ds, lon, lat)
    out = np.full(len(lon), np.nan, dtype=np.float32)
    idx = np.flatnonzero(inside)
    vals = np.asarray(arr[rows[idx],cols[idx]], dtype=np.float32)
    nodata = ds.nodata
    acceptable = np.isfinite(vals)
    if nodata is not None and np.isfinite(nodata):
        acceptable &= vals != nodata
    out[idx[acceptable]] = vals[acceptable]
    return (out, rows, cols, inside) if return_placements else out


def sample_tiff(tif: Path, lon: np.ndarray, lat: np.ndarray) -> np.ndarray:
    with rasterio.open(tif) as ds:
        arr = ds.read(1)
        return sample(arr, ds, lon, lat)


def raster_terrain(elev_tiff: Path, lon: np.ndarray, lat: np.ndarray) -> dict[str,np.ndarray]:
    """Explicit 2.5-arcminute DEM grid-scale (NOT field-scale) terrain metrics.

    Central finite differences and 3x3 elevation neighbors at ~5km spatial
    grain. Nothing is extrapolated off land or across nodata/Ocean cells.
    """
    with rasterio.open(elev_tiff) as ds:
        arr = ds.read(1)
        elev, rows, cols, inside = sample(arr, ds, lon, lat, return_placements=True)
        n = len(lon)
        stack = np.full((9,n), np.nan, dtype=np.float32)
        offsets = [(-1,-1),(-1,0),(-1,1),(0,-1),(0,0),(0,1),(1,-1),(1,0),(1,1)]
        for i, (dr, dc) in enumerate(offsets):
            rr, cc = rows+dr, cols+dc
            valid = inside & (rr >= 0) & (rr < ds.height) & (cc >= 0) & (cc < ds.width)
            ids = np.flatnonzero(valid)
            values = np.asarray(arr[rr[ids],cc[ids]],dtype=np.float32)
            ok = np.isfinite(values)
            if ds.nodata is not None and np.isfinite(ds.nodata):
                ok &= values != ds.nodata
            stack[i, ids[ok]] = values[ok]
        center = stack[4,:]
        finite = np.isfinite(stack).all(axis=0)
        # Near-pole resolution and DEM edge effects are explicitly not modelled.
        finite &= np.abs(lat) < 85
        dx = np.maximum(1.0,np.abs(ds.transform.a) * 111320.0 * np.cos(np.deg2rad(np.where(np.isfinite(lat),lat,0))))
        dy = np.abs(ds.transform.e) * 111132.0
        gx = (stack[5,:] - stack[3,:])/(2*dx)
        gy = (stack[1,:] - stack[7,:])/(2*dy) # north = row -1
        slope = np.rad2deg(np.arctan(np.hypot(gx,gy)))
        relief = np.max(np.where(np.isfinite(stack),stack,-np.inf),axis=0) - np.min(
            np.where(np.isfinite(stack),stack,np.inf),axis=0)
        mean_neigh = (np.nansum(stack,axis=0)-center)/8
        tpi = center-mean_neigh
        tri = np.sqrt(np.nanmean((stack-center[np.newaxis,:])**2,axis=0))
        # Downslope aspect north-clockwise. Flat pixels have undefined aspect.
        eastness = -gx/np.maximum(1e-20,np.hypot(gx,gy))
        northness = -gy/np.maximum(1e-20,np.hypot(gx,gy))
        flat = np.hypot(gx,gy) <= 1e-7
        out = {
            "wc2_elev_2p5m":elev,
            "terrain_slope_2p5m_degree":slope,
            "terrain_roughness_2p5m_m":relief,
            "terrain_TRI_2p5m_m":tri,
            "terrain_TPI_2p5m_m":tpi,
            "terrain_eastness_2p5m":eastness,
            "terrain_northness_2p5m":northness,
        }
        for key in out:
            a = np.asarray(out[key],dtype=np.float32)
            if key != "wc2_elev_2p5m":
                a[~finite] = np.nan
            if key.endswith("eastness_2p5m") or key.endswith("northness_2p5m"):
                a[flat] = np.nan
            out[key] = a
        return out


def run(photo_gzip: Path, source_receipt: Path, raster_dir: Path, outdir: Path) -> dict:
    require(not outdir.exists() or not any(outdir.iterdir()),
            "Do not overwrite already computed original source outputs")
    receipt=json.loads(source_receipt.read_text(encoding="utf-8"))
    require(receipt["final_photo_ids"]==ORIGINAL_COUNT and
            receipt["final_classified_hue_ids"]==ORIGINAL_CLASSIFIED,
            "Unexpected photograph/colour denominator")
    require(file_sha(photo_gzip)==receipt["result_source_photo_csv_sha256"],
            "Original 248312 measured FCP+MOIS photo file digest differs")
    d = pd.read_csv(photo_gzip, compression="gzip", low_memory=False,
                    dtype={"photo_id":str,"inat_taxon_id":str,
                           "observation_id":str})
    require(len(d)==ORIGINAL_COUNT and d.photo_id.nunique()==len(d),
            "Photo ID union incomplete or contains duplication")
    require(int(d.morph.isin(CLASSIFIED).sum())==ORIGINAL_CLASSIFIED,
            "Original four-colour classification labels have changed")
    require(all(raster_dir.joinpath(n).is_dir() for n in SOURCE),
            "Required physical WorldClim raster archive was not staged")
    lon = pd.to_numeric(d.longitude,errors="coerce").to_numpy(dtype=float)
    lat = pd.to_numeric(d.latitude,errors="coerce").to_numpy(dtype=float)
    n=len(d)
    date=pd.to_datetime(d.observed_on,errors="coerce")
    month=date.dt.month.to_numpy(dtype=float)
    # Official WorldClim covariates remain separate source columns.
    out=pd.DataFrame({
      "photo_id":d.photo_id,"inat_taxon_id":d.inat_taxon_id,
      "species":d.species,"morph":d.morph,
      "original_source":d.old_mois_or_new_fcp,
      "observed_on":d.observed_on,
      "latitude":lat,"longitude":lon,
    })
    for col in STATIC_SOURCE:
        if col not in ("latitude","longitude") and col in d:
            out["historical_original_"+col] = pd.to_numeric(d[col],errors="coerce").astype("float32")
    ground = raster_terrain(find_zip_tiff(raster_dir,"elev","wc2.1_2.5m_elev.tif"),lon,lat)
    for col,val in ground.items():
        out[col]=val
    for bio in BIO_IDS:
        path=find_zip_tiff(raster_dir,"bio",f"wc2.1_10m_bio_{bio}.tif")
        out[f"wc2_bio{bio:02d}_10m"]=sample_tiff(path,lon,lat)
    out["wc2_month_tavg_normal_10m_C"] = np.float32(np.nan)
    out["wc2_month_prec_normal_10m_mm"] = np.float32(np.nan)
    out["wc2_month_srad_normal_10m_Wm2"] = np.float32(np.nan)
    # Annual mean irradiance is a site-level long-term climate baseline for
    # ALL photos. Unlike photo-month irradiance it needs no observation date.
    # Keep it separate from genuinely observed-year/month solar radiation.
    radiation_sum = np.zeros(n,dtype=np.float64)
    radiation_count = np.zeros(n,dtype=np.int8)
    for m in range(1,13):
        take=np.flatnonzero(month==m)
        for key,field,unit in (
            ("tavg","wc2_month_tavg_normal_10m_C","tavg"),
            ("prec","wc2_month_prec_normal_10m_mm","prec"),
            ("srad","wc2_month_srad_normal_10m_Wm2","srad")
        ):
            name=f"wc2.1_10m_{unit}_{m:02d}.tif"
            tif=find_zip_tiff(raster_dir,key,name)
            if len(take):
                val=sample_tiff(tif,lon[take],lat[take])
                out.loc[take,field]=val
            if key=="srad":
                annual=sample_tiff(tif,lon,lat)
                valid=np.isfinite(annual)
                radiation_sum[valid]+=annual[valid]
                radiation_count[valid]+=1
    annual_irradiance=np.full(n,np.nan,dtype=np.float32)
    valid_annual=radiation_count==12
    annual_irradiance[valid_annual]=(radiation_sum[valid_annual]/12).astype(np.float32)
    out["wc2_annual_srad_normal_10m_Wm2"]=annual_irradiance
    # No season/time variable or monthly climate covariate is filled without an
    # original observed calendar date. Month-specific WorldClim is a NORMAL.
    require(out.loc[date.isna(),[
      "wc2_month_tavg_normal_10m_C","wc2_month_prec_normal_10m_mm",
      "wc2_month_srad_normal_10m_Wm2"]].notna().sum().sum()==0,
      "Photo without date was assigned an invented seasonal exposure")
    outdir.mkdir(parents=True,exist_ok=True)
    result=outdir/"mois_248312_original_photo_worldclim_static_monthly_normal_terrain.csv.gz"
    out.to_csv(result,index=False,compression={"method":"gzip","mtime":0},float_format="%.6f",lineterminator="\n")
    cover={k:int(out[k].notna().sum()) for k in (
        "wc2_elev_2p5m","terrain_slope_2p5m_degree",
        "terrain_roughness_2p5m_m","terrain_TRI_2p5m_m",
        "wc2_bio01_10m","wc2_bio05_10m","wc2_bio12_10m","wc2_bio15_10m",
        "wc2_month_tavg_normal_10m_C",
        "wc2_month_prec_normal_10m_mm",
        "wc2_month_srad_normal_10m_Wm2","wc2_annual_srad_normal_10m_Wm2")}
    previous={}
    for key,b in (
        ("wc_bio1","wc2_bio01_10m"),
        ("wc_bio5","wc2_bio05_10m"),
        ("wc_bio12","wc2_bio12_10m"),
        ("wc_bio15","wc2_bio15_10m"),
        ("wc_elevation_m","wc2_elev_2p5m")
    ):
        original="historical_original_"+key
        if original not in out:continue
        a=np.asarray(out[original],dtype=float)
        new=np.asarray(out[b],dtype=float)
        compared=np.isfinite(a)&np.isfinite(new)
        previous[key]={
            "n_both_sources":int(compared.sum()),
            "median_absolute_original_vs_worldclim_difference":
                float(np.median(np.abs(a[compared]-new[compared]))) if compared.any() else None,
            "mean_original_minus_worldclim":
                float(np.mean(a[compared]-new[compared])) if compared.any() else None,
        }
    white=out.morph.eq("white")
    classified=out.morph.isin(CLASSIFIED)
    diagnostics=[]
    for variable in ("wc2_elev_2p5m","terrain_slope_2p5m_degree",
                     "terrain_roughness_2p5m_m","wc2_bio01_10m","wc2_bio05_10m",
                     "wc2_bio12_10m","wc2_bio15_10m","wc2_annual_srad_normal_10m_Wm2"):
        for source in sorted(out.original_source.unique()):
            p=out.loc[(out.original_source==source)&classified&
                      out[variable].notna(),["inat_taxon_id",variable,"morph"]].copy()
            if len(p)<100:continue
            p["photo_white"]=p.morph.eq("white")
            try:p["quartile"]=pd.qcut(p[variable],4,duplicates="drop")
            except ValueError:continue
            for level,q in p.groupby("quartile",observed=True):
                sp=q.groupby("inat_taxon_id").photo_white.mean()
                diagnostics.append({
                    "source":source,"variable":variable,"quartile":str(level),
                    "classifiable_photos":int(len(q)),"named_taxa":int(len(sp)),
                    "photographic_white_fraction_photo_weighted":round(float(q.photo_white.mean()),6),
                    "photographic_white_fraction_species_equal_descriptive":round(float(sp.mean()),6),
                    "test_of_causality_or_frozen_global_coefficient":False,
                })
    quartile=outdir/"mois_248312_source_species_balanced_static_gradient_descriptives.csv"
    pd.DataFrame(diagnostics).to_csv(quartile,index=False,lineterminator="\n")
    source_checks={
        key:file_sha(raster_dir/key/SOURCE[key]) for key in SOURCE
        if (raster_dir/key/SOURCE[key]).exists()
    }
    for key,digest in SOURCE_PREVIOUS_SHA.items():
        require(source_checks.get(key)==digest,"Previously pinned WorldClim original source ZIP changed: "+key)
    result_receipt={
        "schema":"mois_248312_real_original_photo_worldclim_static_and_monthly_climate_normals_v1",
        "status":"REAL_248312_PHOTO_GEOGRAPHY_WORLDCLIM_NORMALS_NOT_OBSERVED_YEAR_WEATHER",
        "n_original_photo_ids":n,
        "n_four_hue_classified":int(classified.sum()),
        "n_original_dated_photo_ids":int(date.notna().sum()),
        "worldclim_2_1_bio_normal_period":"1970-2000",
        "monthly_normals_not_same_year_actual_weather":True,
        "terrain_elevation_dem_grain":"2.5 arcminutes ~4.6km at equator, central finite differences; not local-scale ecological terrain",
        "n_total_georeferenced":int((np.isfinite(lat)&np.isfinite(lon)).sum()),
        "coverage":cover,
        "existing_historical_source_vs_new_raster_overlap":previous,
        "original_source_sha256":file_sha(photo_gzip),
        "worldclim_source_zip_sha256":source_checks,
        "worldclim_official_download_endpoints":SOURCE_OFFICIAL,
        "full_photo_static_monthly_normal_csv_sha256":file_sha(result),
        "source_stratified_quartile_descriptive_sha256":file_sha(quartile),
        "FCP_within_species_polymorphism_statistics_imported":False,
        "no_genetic_fitness_selection_causality":True,
        "soil_not_newly_mapped_in_this_execution":True,
        "real_monthly_2025_or_previous_weather_not_yet_retrieved_for_full_global_248312":True,
        "time_of_day_requires_source_timestamp_timezone":True,
    }
    (outdir/"mois_248312_worldclim_static_month_normal_receipt.json").write_text(
        json.dumps(result_receipt,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("MOIS_248312_REAL_WORLDCLIM_STATIC_MONTHLY_NORMAL_JOIN",
          json.dumps({"coverage":cover,"known_original_colours":int(classified.sum()),
                      "dated_photos":int(date.notna().sum())},sort_keys=True),flush=True)
    return result_receipt


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--photo-gzip",type=Path,required=True)
    p.add_argument("--photo-source-receipt",type=Path,required=True)
    p.add_argument("--raster-dir",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    run(a.photo_gzip,a.photo_source_receipt,a.raster_dir,a.output)
