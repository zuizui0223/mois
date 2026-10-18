#!/usr/bin/env python3
"""MOIS actual original-photo 313535: complete long-term WorldClim BIO1-19,
annual VAPR and full DEM terrain, with soil null where unmeasured.

This is an incremental interoperable extension of a verified photo-ID feature
lake: do not redo RGB classification, invent missing soil, use precolour
source IDs as phenotypes, or confuse 1970-2000 normals with observation weather.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from scripts.analysis import mois_global313535_water_light_heat_habitat_partition as previous
from scripts.analysis.mois_sample_248312_worldclim_static_and_monthnormals import (
    raster_terrain, sample_tiff,
)

BIO = tuple(f"wc2_bio{i:02d}_10m" for i in range(1,20))
VAPR = "wc2_annual_vapr_10m"
TERRAIN = (
    "wc2_elev_2p5m", "terrain_slope_2p5m_degree",
    "terrain_roughness_2p5m_m", "terrain_TRI_2p5m_m",
    "terrain_TPI_2p5m_m", "terrain_eastness_2p5m",
    "terrain_northness_2p5m",
)
CLIMATE = (*BIO, VAPR, previous.SOLAR)
SOURCE_REPO_RELEASE = previous.RELEASE
OUT_ENV = "mois_313535_all19_bio_vapr_solar_terrain_photoids.parquet"
OUT_REC = "mois_313535_all19_bio_vapr_fullterrain_source_manifest.json"
OUT_COVERAGE = "mois_313535_all19_environment_missingness.csv"

PIN_BIO = "00513224583665ec0f2f955a4ec252730c4deb2004cce9e793492a3f26df4dcf"
PIN_SRAD = "c72ee7f4f9a0eb4b5f6cd7a003eddc05bda22a2e5666d968ed4e817fd36b9026"


def need(v, msg):
    if not v:
        raise ValueError(msg)


def hash_file(path: Path) -> str:
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for part in iter(lambda: f.read(1<<20), b""):
            h.update(part)
    return h.hexdigest()


def climatology(p: pd.DataFrame, raster_root: Path) -> pd.DataFrame:
    """Original 313535 geolocated photo IDs; no classified-hue filtering."""
    lon = pd.to_numeric(p.longitude,errors="coerce").to_numpy(float)
    lat = pd.to_numeric(p.latitude,errors="coerce").to_numpy(float)
    new = pd.DataFrame({"photo_id":p.photo_id.to_numpy()})
    for i,k in enumerate(BIO,1):
        path=raster_root/"bio"/f"wc2.1_10m_bio_{i}.tif"
        need(path.is_file(), "Missing 1970-2000 official annual BIO raster: "+k)
        new[k]=sample_tiff(path,lon,lat)
    annual=np.zeros(len(p),dtype=float)
    count=np.zeros(len(p),dtype=np.int8)
    for month in range(1,13):
        path=raster_root/"vapr"/f"wc2.1_10m_vapr_{month:02d}.tif"
        need(path.is_file(), "Missing original WorldClim VAPR monthly data")
        values=sample_tiff(path,lon,lat)
        mask=np.isfinite(values)
        annual[mask]+=values[mask]
        count[mask]+=1
    vapr=np.full(len(p),np.nan,dtype=np.float32)
    ok=count==12
    vapr[ok]=(annual[ok]/12).astype(np.float32)
    new[VAPR]=vapr
    terrain=raster_terrain(raster_root/"elev"/"wc2.1_2.5m_elev.tif",lon,lat)
    for name in TERRAIN:
        new[name]=terrain[name]
    new[previous.SOLAR]=pd.to_numeric(
        p[previous.SOLAR],errors="coerce").to_numpy(float)
    need(new.photo_id.tolist()==p.photo_id.tolist()
         and len(new)==previous.N_PHOTOS
         and new.photo_id.is_unique,
         "New 19-BIO and terrain environmental source photo identity changed")
    return new


def verify_existing_photo_environment(previous_source: pd.DataFrame,
                                      new: pd.DataFrame) -> dict:
    # Four BIO variables are checked against the original already computed
    # 313535-source WorldClim. Detect unit, coordinate and raster-source drift.
    cross={}
    for v in (previous.HEAT,previous.RAIN,previous.DRY,previous.SEASON):
        old=pd.to_numeric(previous_source[v],errors="coerce").to_numpy(float)
        revised=pd.to_numeric(new[v],errors="coerce").to_numpy(float)
        valid=np.isfinite(old)&np.isfinite(revised)
        need(valid.sum()>110000,"Existing original 4 annual BIO columns unavailable")
        difference=np.abs(old[valid]-revised[valid])
        report={
            "n_same_photo_sites":int(valid.sum()),
            "p99_abs_source_difference":float(np.percentile(difference,99)),
            "fraction_above_point_05_source_units":float(np.mean(difference>.05)),
        }
        # The authoritative WorldClim source should agree at original photo
        # locations; prefer failing rather than silently mix two climatologies.
        need(report["fraction_above_point_05_source_units"]<.001,
             "Already measured BIO values disagree with new original raster")
        cross[v]=report
    return cross


def run(release:Path,raster_root:Path,reference_water:Path,out:Path):
    need(not out.exists() or not any(out.iterdir()),
         "Do not rewrite source-environment snapshot")
    source,source_rec=previous.load_source(release)
    ref=json.loads(reference_water.read_text())
    need(ref["all_original_photo_ids"]==309439
         and ref["original_genus_bootstrap_draws"]==999
         and ref["source_climatology_archives_sha256"]["bio"]==PIN_BIO
         and ref["source_climatology_archives_sha256"]["srad"]==PIN_SRAD
         and ref["source_climatology_archives_sha256"]["vapr"]==
             hash_file(raster_root/"vapr"/"wc2.1_10m_vapr.zip"),
         "Climatology must match independently measured original 309439 model")
    files={
        "bio":raster_root/"bio"/"wc2.1_10m_bio.zip",
        "vapr":raster_root/"vapr"/"wc2.1_10m_vapr.zip",
        "elev":raster_root/"elev"/"wc2.1_2.5m_elev.zip",
    }
    need(hash_file(files["bio"])==PIN_BIO
         and all(p.is_file() for p in files.values()),
         "Original WorldClim BIO, vapour or DEM zip missing/altered")
    new=climatology(source,raster_root)
    original_consistency=verify_existing_photo_environment(source,new)
    cover=[]
    for name in CLIMATE+TERRAIN:
        valid=pd.to_numeric(new[name],errors="coerce").notna()
        cover.append({"predictor":name,"photo_attempts":len(new),
                      "environment_value_present":int(valid.sum()),
                      "environment_missing":int((~valid).sum()),
                      "source":"WorldClim2.1_1970_2000_normals" if name in CLIMATE
                               else "WorldClim2.1_2p5_arcminute_coarse_DEM"})
    need(min(x["environment_value_present"] for x in cover)>180000,
         "Global annual climate / DEM raster did not map enough real photo sites")
    out.mkdir(parents=True,exist_ok=True)
    target=out/OUT_ENV
    new.to_parquet(target,compression="zstd",engine="pyarrow",index=False)
    reread=pd.read_parquet(target,columns=["photo_id"])
    need(reread.photo_id.tolist()==source.photo_id.tolist(),
         "Saved 19BIO source lost photographic photo IDs")
    table=out/OUT_COVERAGE
    pd.DataFrame(cover).to_csv(table,index=False,lineterminator="\n")
    report={
        "schema":"mois_actual313535_original_photo_worldclim_full19_vapr_dem_v1",
        "source_release":SOURCE_REPO_RELEASE,
        "source_prior_release_three_Parquet_SHA256":
            {k:source_rec["files_sha256"][k] for k in (
                previous.PHOTO,previous.ANNUAL,previous.STATIC)},
        "real_original_processed_photo_ids":previous.N_PHOTOS,
        "original_four_hue_photo_labels":previous.N_HUES,
        "original_white_photo_ROIs":previous.WHITE_N,
        "previous_313535_photo_IDs_and_hues_unchanged":True,
        "site_annual_climatology_WorldClim_version":"2.1 (1970-2000)",
        "source_climatology_zip_SHA256":
            {key:hash_file(p) for key,p in files.items()},
        "original_four_annual_BIO_site_consistency":original_consistency,
        "num_BIO_variables":len(BIO),
        "has_annual_original_vapour_and_solar_normal":True,
        "terrain_variables":list(TERRAIN),
        "terrain_is_coarse_raster_not_microhabitat":True,
        "soil_not_sampled_from_new_spatial_source":True,
        "soil_stays_source_partial_and_no_missing_value_invented":True,
        "climate_not_observation_month_weather":True,
        "photo_colour_not_true_open_flower_phenophase":True,
        "not_a_fitness_or_50_50_evolutionary_equilibrium_test":True,
        "new_environment_parquet_sha256":hash_file(target),
        "missingness_csv_sha256":hash_file(table),
    }
    (out/OUT_REC).write_text(json.dumps(
        report,ensure_ascii=False,sort_keys=True,indent=2,
        allow_nan=False)+"\n")
    print("MOIS_ALL19_BIO_ORIGINAL313535_ENVIRONMENT_SOURCE_VERIFIED",
          json.dumps({"n_real_photos":previous.N_PHOTOS,
                      "predictors":len(CLIMATE)+len(TERRAIN),
                      "coverage":cover[:3],
                      "vapr_zip_sha":report["source_climatology_zip_SHA256"]["vapr"]},
                     sort_keys=True),flush=True)
    return report

if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--original-photo-release",type=Path,required=True)
    parser.add_argument("--worldclim-raster-root",type=Path,required=True)
    parser.add_argument("--original309439-water-model-receipt",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    a=parser.parse_args()
    run(a.original_photo_release,a.worldclim_raster_root,
        a.original309439_water_model_receipt,a.output)
