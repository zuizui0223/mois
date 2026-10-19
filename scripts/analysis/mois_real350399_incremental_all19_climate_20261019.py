#!/usr/bin/env python3
"""Immutable 317631 real-photo WorldClim lake + 32768 independently processed original photo sites.

Read actual original RGB/ROI measurement releases, never source-only IDs.
Point-sample WorldClim 2.1 1970–2000 climatologies/DEM ONLY for the newly
measured source photos; prior 317631 photo-site rasters must be byte-identical.
No soil imputation, flowering, climate-year anomaly, or selection/fitness claims.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from scripts.analysis.mois_actual_photo_cohort_season_audit_20261019 import (
    HUES, assemble, require, sha256,
)
from scripts.analysis import mois_original313535_full19_bio_vapr_terrain_lake as expanded
from scripts.analysis import mois_global313535_water_light_heat_habitat_partition as parent
from scripts.analysis.mois_sample_248312_worldclim_static_and_monthnormals import (
    raster_terrain, sample_tiff,
)

BASE=317631
NEW=32768
FULL=350399
OUT_ENV="mois350399_all19_annual_vapr_solar_terrain_source_photoids.parquet"
OUT_NEW="mois32768_new_only_real_photo_annual_climatology.parquet"
OUT_RECEIPT="mois350399_immutable_annual_climate_source_receipt.json"
OUT_COVERAGE="mois350399_original_photo_annual_environment_coverage.csv"
OUT_COLOUR="mois350399_original_photo_white_annual_rain_solar_descriptive.csv"
RASTER_ARCHIVES={
    "bio":"wc2.1_10m_bio.zip",
    "srad":"wc2.1_10m_srad.zip",
    "vapr":"wc2.1_10m_vapr.zip",
    "elev":"wc2.1_2.5m_elev.zip",
}
SOLAR=parent.SOLAR


def new_sites_climate(photos, raster, *, sampler=sample_tiff,
                      terrain_sampler=raster_terrain):
    photos=photos.reset_index(drop=True)
    n=len(photos)
    require(n>0 and photos.photo_id.is_unique,"Invalid new actual photo site cohort")
    lon=pd.to_numeric(photos.longitude,errors="coerce").to_numpy(float)
    lat=pd.to_numeric(photos.latitude,errors="coerce").to_numpy(float)
    out=pd.DataFrame({"photo_id":photos.photo_id.copy()})
    for i,name in enumerate(expanded.BIO, 1):
        tif=raster/"bio"/f"wc2.1_10m_bio_{i}.tif"
        require(tif.is_file(),"Missing WorldClim BIO original TIFF "+str(i))
        out[name]=sampler(tif,lon,lat)
    vap_sum=np.zeros(n, dtype=float); vap_n=np.zeros(n,dtype=np.int8)
    sun_sum=np.zeros(n, dtype=float); sun_n=np.zeros(n,dtype=np.int8)
    for month in range(1,13):
        a=raster/"vapr"/f"wc2.1_10m_vapr_{month:02d}.tif"
        b=raster/"srad"/f"wc2.1_10m_srad_{month:02d}.tif"
        require(a.is_file() and b.is_file(),"Official monthly VAPR/SRAD TIFF missing")
        va=sampler(a,lon,lat); sa=sampler(b,lon,lat)
        ok=np.isfinite(va); vap_sum[ok]+=va[ok]; vap_n[ok]+=1
        ok=np.isfinite(sa); sun_sum[ok]+=sa[ok]/86.4; sun_n[ok]+=1
    vap=np.full(n,np.nan,dtype=np.float32)
    sun=np.full(n,np.nan,dtype=np.float32)
    mask=vap_n==12; vap[mask]=(vap_sum[mask]/12).astype(np.float32)
    mask=sun_n==12; sun[mask]=(sun_sum[mask]/12).astype(np.float32)
    out[expanded.VAPR]=vap
    out[SOLAR]=sun
    tif=raster/"elev"/"wc2.1_2.5m_elev.tif"
    require(tif.is_file(),"Official original DEM source TIFF missing")
    terrain=terrain_sampler(tif,lon,lat)
    for name in expanded.TERRAIN:
        require(name in terrain,"Missing original terrain variable "+name)
        out[name]=terrain[name]
    require(out.photo_id.is_unique and len(out)==n
            and set(out.columns)==set(("photo_id",*expanded.BIO,
                                       expanded.VAPR,SOLAR,*expanded.TERRAIN)),
            "The new photo-climatology does not match original all19 schema")
    return out


def append_immutable(old_env, core_photos, new_photos, new_env):
    require(len(old_env)==BASE and len(core_photos)==BASE
            and len(new_env)==len(new_photos)==NEW,
            "Original previous climate lake and actual 32768 pixels differ")
    require(old_env.photo_id.is_unique and new_env.photo_id.is_unique
            and core_photos.photo_id.is_unique and new_photos.photo_id.is_unique,
            "Duplicate original photo ID in WorldClim source")
    require(old_env.photo_id.astype(str).tolist()==core_photos.photo_id.tolist()
            and new_env.photo_id.astype(str).tolist()==new_photos.photo_id.tolist(),
            "WorldClim row order is not aligned to source-verified original photos")
    require(not new_env.photo_id.isin(old_env.photo_id).any(),
            "New photo site already in old original annual climate cache")
    require(list(old_env.columns)==list(new_env.columns),
            "ABIOTIC variable units/columns do not match immutable prior source")
    result=pd.concat([old_env,new_env],ignore_index=True)
    require(len(result)==FULL and result.photo_id.is_unique
            and result.iloc[:BASE].reset_index(drop=True).equals(
                old_env.reset_index(drop=True)),
            "Previously measured WorldClim photo-site results were changed")
    return result


def descriptive(d, e):
    require(len(d)==len(e)==FULL and
            d.photo_id.tolist()==e.photo_id.astype(str).tolist(),
            "Photo hue source and actual annual source environmental IDs drifted")
    x=d[["photo_id","inat_taxon_id","species","arm","morph","classified","white"]].copy()
    for col in ("wc2_bio05_10m","wc2_bio12_10m",
                "wc2_bio14_10m",SOLAR,expanded.VAPR):
        x[col]=pd.to_numeric(e[col],errors="coerce").to_numpy(float)
    x["cohort"]=np.where(x.arm.eq("core"),"core","new_species_first")
    coverage=[]
    for (arm,source),v in x.groupby(["cohort","classified"],observed=True):
        coverage.append({
            "source_cohort":arm,"four_hue_classified":bool(source),
            "photo_attempts":len(v),
            "annual_rain_present":int(v.wc2_bio12_10m.notna().sum()),
            "annual_solar_present":int(v[SOLAR].notna().sum()),
            "annual_heat_present":int(v.wc2_bio05_10m.notna().sum()),
            "annual_vapour_present":int(v[expanded.VAPR].notna().sum())
        })
    classify=x[x.classified].copy()
    # Fixed across all pooled SOURCE PHOTOS BEFORE looking at white outcomes.
    # This still estimates *photo-label composition*, not a natural flower census.
    cuts={}
    for col,key in [("wc2_bio12_10m","rain"),(SOLAR,"solar")]:
        a=classify[col].dropna()
        require(len(a)>100 and a.nunique()>10,
                "No genuine environmental climate support for source photos")
        q=np.quantile(a,[1/3,2/3])
        require(q[0]<q[1],"Insufficient distinct abiotic environmental sample")
        classify[key]=pd.cut(classify[col],
                             [-np.inf,q[0],q[1],np.inf],
                             labels=("low","mid","high"))
        cuts[key]={"source_variable":col,
                   "low_middle_split":float(q[0]),
                   "middle_high_split":float(q[1]),
                   "all_classified_actual_source_photo_selection":True}
    rows=[]
    for (cohort,rain,solar),z in classify.groupby(
            ["cohort","rain","solar"],observed=True):
        taxa=z.groupby("inat_taxon_id",observed=True).white.mean()
        nominal_genus=z.species.str.split().str[0]
        p=z.assign(genus=nominal_genus).groupby(
            ["genus","inat_taxon_id"],observed=True).white.mean()
        equal_genus=p.groupby(level=0).mean()
        rows.append({
            "cohort":cohort,"rain_tertile":str(rain),
            "solar_tertile":str(solar),
            "four_hue_classified_real_photo_ROIs":len(z),
            "white_real_photo_ROIs":int(z.white.sum()),
            "white_photo_fraction":float(z.white.mean()),
            "taxon_equal_photo_white_fraction":float(taxa.mean()),
            "genus_equal_photo_white_fraction":float(equal_genus.mean()),
            "nominal_taxa":int(len(taxa)),
            "nominal_genera":int(len(equal_genus)),
            "median_annual_rain_mm":float(z.wc2_bio12_10m.median()),
            "median_long_term_solar_Wm2":float(z[SOLAR].median()),
            "median_warmest_month_temp_C":float(z.wc2_bio05_10m.median())
        })
    return pd.DataFrame(coverage),pd.DataFrame(rows),cuts


def run(core,broad,wave2,cached8192,tail16384,old_env,original_env_receipt,raster,out):
    out=Path(out)
    require(not out.exists() or not any(out.iterdir()),
            "Never overwrite previously signed original photographic environment source")
    source,arms=assemble(core,broad,wave2,cached8192,tail16384,strict=True)
    previous=pd.read_parquet(old_env)
    original=json.loads(Path(original_env_receipt).read_text())
    oldname=Path(old_env).name
    require(original["schema"]=="mois_real317631_fifth4096_incremental_RGB_annual_bio19_terrain_v1"
            and original["all_original_photo_attempts"]==BASE
            and sha256(old_env)==original["photo_ID_aligned_Zstd_files_sha256"][oldname],
            "Base 317631 immutable original photo all19 environment not verified")
    source_archives={}
    for name,file in RASTER_ARCHIVES.items():
        p=Path(raster)/name/file
        expected=original["source_original_raster_sha256"][name]
        require(p.is_file() and sha256(p)==expected,
                "Official raster source or physical solar units differ: "+name)
        source_archives[name]=expected
    require(source_archives["bio"]==expanded.PIN_BIO and
            source_archives["srad"]==expanded.PIN_SRAD,
            "Original WorldClim and solar source SHA unexpectedly changed")
    novel=pd.concat(arms[1:],ignore_index=True)
    require(len(novel)==NEW and novel.photo_id.is_unique
            and not novel.photo_id.isin(arms[0].photo_id).any(),
            "Species-first 32768 photographed pixels not disjoint")
    newenv=new_sites_climate(novel,Path(raster))
    newenv=newenv[previous.columns]
    full=append_immutable(previous,arms[0],novel,newenv)
    out.mkdir(parents=True,exist_ok=True)
    full.to_parquet(out/OUT_ENV,compression="zstd",index=False)
    newenv.to_parquet(out/OUT_NEW,compression="zstd",index=False)
    saved=pd.read_parquet(out/OUT_ENV,columns=["photo_id"])
    require(saved.photo_id.astype(str).tolist()==source.photo_id.tolist(),
            "Saved all350399 photo original IDs not aligned after Parquet output")
    cover,descript,cuts=descriptive(source,full)
    cover.to_csv(out/OUT_COVERAGE,index=False)
    descript.to_csv(out/OUT_COLOUR,index=False)
    receipt={
        "schema":"mois_full350399_source_photograph_original_worldclim_all19_annual_v1",
        "source_photo_attempts":FULL,
        "previous_unchanged_photo_climate":BASE,
        "new_only_original_image_photo_sites_sampled":NEW,
        "old_original_photo_environment_SHA256":sha256(old_env),
        "source_pixel_shas":{
            name:sha256(f) for name,f in [
                ("core",core),("broad_cumulative",broad),
                ("wave2",wave2),("cached8192",cached8192),
                ("tail16384",tail16384)]
        },
        "raster_source_original_ZIP_SHA256":source_archives,
        "source_annual_climate_period":"WorldClim 2.1 1970-2000 climatological normals",
        "solar_unit":"W/m2 after dividing official raw kJ/m2/day by 86.4 ONLY at newly measured photo sites",
        "source_soil_unmeasured_for_new_photos":True,
        "year_of_observation_weather_inferred_from_normals":False,
        "original_photo_ROI_four_hue_not_independently_botanically_reviewed":True,
        "sampling_design_preexisting_taxon_first_not_unbiased_world_floral_census":True,
        "original_photo_pixels_reprocessed":False,
        "additional_image_photo_pixels_downloaded":False,
        "original_317631_climate_data_immutable":True,
        "technical_unclassifiable_hue_not_called_colour":True,
        "annual_rain_solar_tertiles_from_classified_source":cuts,
        "source_environment_coverage":cover.to_dict("records"),
        "new_32768_photo_site_environment_variables":[c for c in newenv.columns if c!="photo_id"],
        "sha256":{k:sha256(out/k) for k in
            (OUT_ENV,OUT_NEW,OUT_COVERAGE,OUT_COLOUR)}
    }
    (out/OUT_RECEIPT).write_text(json.dumps(
        receipt,indent=2,sort_keys=True,ensure_ascii=False,allow_nan=False)+"\n")
    print("MOIS_REAL350399_ALL19_ANNUAL_WORLDCLIM_EXTENDED_ONLY_NEW32768",
          json.dumps({
            "photos":FULL,"old_unchanged":BASE,"new_sites":NEW,
            "classified_original_ROIs":int(source.classified.sum()),
            "env_coverage":receipt["source_environment_coverage"],
            "rain_solar_source_tertiles":len(descript),
            "written":receipt["sha256"]
          },sort_keys=True),flush=True)
    return receipt


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("core","broad","wave2","cached8192","tail16384",
                 "old-env","old-receipt","raster","output"):
        p.add_argument("--"+name,type=Path,required=True)
    a=p.parse_args()
    run(a.core,a.broad,a.wave2,a.cached8192,a.tail16384,
        a.old_env,a.old_receipt,a.raster,a.output)


if __name__=="__main__":
    main()
