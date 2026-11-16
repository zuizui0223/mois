#!/usr/bin/env python3
"""Real ESA CCI v7 10km woody above-ground biomass in MOIS 37104 species GBIF grid.

No conflation of biomass with NDVI or NPP. Biomass is woody AGB Mg/ha, and
zeros may include open nonforest landscapes, not missing flower colour.
Sample 25 points *on WorldClim land mask* in every 100km GBIF cell.
"""
from __future__ import annotations
import argparse,hashlib,io,json
from pathlib import Path
from zipfile import ZipFile
import numpy as np
import pandas as pd
from scripts.analysis.mois_global_bioclim_ndvi_rf_features_20261116 import (
    sample_lonlat,raster_probe,point_mean,sha
)
FEATURE_ZIP="e08e81d1c2afd1a2ee5312f3ebc4aeb5ab7678988c2453a2a9c63d28c82eb17a"
AGB_ZIP="4b90104549bc05eb3a9c2ca1cac84d56994a3f0a56a0fae1372d31662cf57407"
GLOBAL_ZIP="422f6d6489e61c6e57c5d51db614e47f613039781a901e3fd3481c13b2b085c2"
BIO_ZIP="00513224583665ec0f2f955a4ec252730c4deb2004cce9e793492a3f26df4dcf"
AGB_TIF="ESACCI-BIOMASS-L4-AGB-MERGED-10000m-fv7.0.tif"
AGB_SD_TIF="ESACCI-BIOMASS-L4-AGB_SD-MERGED-10000m-fv7.0.tif"
AGB_TIF_SHA="6b247c3dc3321dac747d65ca7a53f77753a7edfeeed5b02a697c48d52c9aee3e"
AGB_SD_TIF_SHA="b50d72eba78bf057576846c866369517c714251e424f00bd0b18de828328b9d9"
AGB_DOI="10.5285/6429d1aafe1e43b9b414e4a5a7f8b903"
PIXELS_PER_GBIF_CELL=25

def require(ok,why):
    if not bool(ok):raise ValueError("MOIS_GLOBAL_ESACCI_AGB: "+why)

def assert_source_file(path,expected):
    require(sha(path)==expected,"signed public source ZIP changed")
    return ZipFile(path)

def read_inputs(feature_zip,agb_zip,global_zip,bio_zip):
    import pyarrow.parquet as pq
    f=assert_source_file(feature_zip,FEATURE_ZIP)
    r=json.loads(f.read("receipt.json"))
    member="WorldClim19_weather_NDVI_18577_global_100km_cells.parquet"
    b=f.read(member)
    require(sha(b)==r["output_SHA256"][member],"original 18577 cell feature SHA changed")
    cell=pq.read_table(io.BytesIO(b)).to_pandas()
    require(len(cell)==18577 and not cell[["grid_x","grid_y"]].duplicated().any(),
            "unique GBIF equal-area cells changed")
    member="WorldClim19_weather_NDVI_37104_COLXR_photo9_species.parquet"
    b=f.read(member)
    require(sha(b)==r["output_SHA256"][member],"accepted 37k species climate SHA changed")
    species=pq.read_table(io.BytesIO(b)).to_pandas()
    require(len(species)==37104 and species.key.is_unique,"species taxonomy incomplete")

    z=assert_source_file(agb_zip,AGB_ZIP)
    receipt=json.loads(z.read("receipt.json"))
    require(receipt["no_auth_used"] and receipt["dataset_DOI"]==AGB_DOI,
            "ESA CCI true AGB v7 provenance not verified")
    a=z.read(AGB_TIF);sd=z.read(AGB_SD_TIF)
    require(sha(a)==AGB_TIF_SHA and sha(sd)==AGB_SD_TIF_SHA,
            "ESACCI woody biomass raw 10km raster SHA differs")

    z=assert_source_file(global_zip,GLOBAL_ZIP)
    raw=z.read("all_37104_species_x_EPSG6933_100km_cells.parquet")
    gbif_receipt=json.loads(z.read("receipt.json"))
    require(gbif_receipt["all_16_shards_source_verified"]
            and gbif_receipt["GBIF_occurrences_passed_all_filters"]==396142825,
            "global true GBIF source absent")
    pairs=pq.read_table(io.BytesIO(raw),
        columns=["specieskey","grid_x","grid_y","occurrence_count"]).to_pandas()
    require(len(pairs)==5493931 and pairs.occurrence_count.sum()==396142825,
            "global species-grid completeness changed")

    z=assert_source_file(bio_zip,BIO_ZIP)
    b=z.read("wc2.1_10m_bio_1.tif")
    return cell,species,pairs,{"AGB":a,"AGB_SD":sd,"land_mask_BIO1":b}

def get_band(pixels,year,member):
    from rasterio.io import MemoryFile
    with MemoryFile(pixels) as m:
        with m.open() as ds:
            require(ds.count==18 and str(ds.crs).upper()=="EPSG:4326"
                    and ds.width==3600 and ds.height==1800,
                    "ESA aggregated 18-band source geometry changed")
            bands={int(s.split("-")[-2]):i for i,s in enumerate(ds.descriptions,1)}
            require(set((2010,2020,2024)).issubset(bands),
                    "ESA CCI multi-year AGB bands missing")
            b=ds.read(bands[year])
            require(np.issubdtype(b.dtype,np.integer) and b.min()>=0 and b.max()<=1000,
                    "ESACCI woody biomass Mg/ha scale changed")
            return b,dict(transform=ds.transform,crs=str(ds.crs),year=year,
                          source_band_index=bands[year],source_band_name=ds.descriptions[bands[year]-1])

def point_samples_from_array(arr,coords):
    # ESA source and WorldClim are EPSG4326 geographic, but not same pixel grid.
    # Raster origin -180/90, resolution 0.1°, uint16 Mg/ha.
    n=len(coords[0][0]);out=np.full((len(coords),n),np.nan,dtype=float)
    for i,(lon,lat) in enumerate(coords):
        valid=np.isfinite(lon)&np.isfinite(lat)&(lon>=-180)&(lon<180)&(lat>=-90)&(lat<90)
        ix=np.flatnonzero(valid)
        xx=np.floor((lon[ix]+180)*10).astype(int)
        yy=np.floor((90-lat[ix])*10).astype(int)
        good=(xx>=0)&(xx<3600)&(yy>=0)&(yy<1800)
        out[i,ix[good]]=arr[yy[good],xx[good]]
    return out

def combine_sources(cell,rasters):
    coords=sample_lonlat(cell[["grid_x","grid_y"]])
    land=raster_probe(rasters["land_mask_BIO1"],{"id":"wc2_BIO1"},coords)
    land_mask=np.isfinite(land)
    require(land_mask.sum()>100000,"WorldClim land mask did not load")
    n=len(cell)
    extra=cell[["grid_x","grid_y"]].copy()
    source_info={}
    for year in (2010,2020,2024):
        arr,detail=get_band(rasters["AGB"],year,"AGB")
        values=point_samples_from_array(arr,coords)
        values[~land_mask]=np.nan
        mean,count=point_mean(values)
        extra[f"ESA_AGB_woody_Mg_ha_{year}"]=mean
        extra[f"ESA_AGB_woody_positive_area_fraction_{year}"]=np.divide(
            np.sum(values>0,axis=0),count,out=np.full(n,np.nan),where=count>0)
        if year==2020:
            extra["ESA_AGB_land_valid_25points"]=count
        source_info[str(year)]={"source_annual_band":detail,
            "cell_count_with_land_AGB":int(np.isfinite(mean).sum()),
            "zero_biomass_land_cell_count":int((mean==0).sum())}
    arr,detail=get_band(rasters["AGB_SD"],2020,"AGB_SD")
    sd=point_samples_from_array(arr,coords)
    sd[~land_mask]=np.nan
    extra["ESA_AGB_2020_uncertainty_Mg_ha"]=point_mean(sd)[0]
    extra["ESA_AGB_2020_minus_2010_Mg_ha"]=(
        extra.ESA_AGB_woody_Mg_ha_2020-extra.ESA_AGB_woody_Mg_ha_2010)
    extra["ESA_AGB_2024_minus_2020_Mg_ha"]=(
        extra.ESA_AGB_woody_Mg_ha_2024-extra.ESA_AGB_woody_Mg_ha_2020)
    source_info["2020_uncertainty"]={"source_band":detail}
    require(len(extra)==18577 and not extra[["grid_x","grid_y"]].duplicated().any(),
            "ESACCI original 100km cell keys altered")
    return extra,source_info

def attach_species(species,pairs,env):
    # Source features contain previously averaged weather+NDVI per species.
    g=pairs[["specieskey","grid_x","grid_y"]].merge(
        env,on=["grid_x","grid_y"],how="left",validate="many_to_one")
    require(len(g)==len(pairs),"AGB join duplicated GBIF species cell rows")
    x=g.groupby("specieskey",sort=False).agg({
        name:"mean" for name in env.columns if name not in ("grid_x","grid_y")
    })
    data=species.merge(x,left_on="key",right_index=True,
        how="left",validate="one_to_one")
    require(len(data)==37104 and data.key.is_unique
            and (data.records_passed_spatial_QA.sum()==396142825),
            "species or full source GBIF records lost during woody biomass join")
    return data

def run(feature_zip,agb_zip,global_zip,bio_zip,output):
    out=Path(output);require(not out.exists(),"refuse overwrite AGB provenance")
    cell,species,pairs,rasters=read_inputs(feature_zip,agb_zip,global_zip,bio_zip)
    agb,details=combine_sources(cell,rasters)
    combined=attach_species(species,pairs,agb)
    require((combined.ESA_AGB_woody_Mg_ha_2020.notna().sum()>30000)
            and (combined.MODIS_NDVI_2001_2021_mean.notna().sum()>30000),
            "biomass and NDVI joint source coverage too low")
    out.mkdir(parents=True)
    speciesfile=out/"BIO19_weather_NDVI_ESACCIv7_AGB_37104_accepted_species.parquet"
    cellfile=out/"ESACCI_AGB_2010_2020_2024_18577_GBIF_cells.parquet"
    combined.to_parquet(speciesfile,index=False)
    agb.to_parquet(cellfile,index=False)
    report={"schema":"mois_37104_global_9photo_hues_ESA_CCI_v7_biomass_real_10km_v1",
      "source_global_accepted_species":37104,
      "source_GBIF_global_snapshot":"2026-11-01",
      "source_GBIF_valid_occurrences":396142825,
      "source_GBIF_species_grid_pairs":5493931,
      "source_100km_cells":18577,
      "source_BIO19_NDVI_feature_artifact_SHA256":FEATURE_ZIP,
      "source_ESA_CCI_v7_AGB_artifact_SHA256":AGB_ZIP,
      "source_GBIF_full_grid_artifact_SHA256":GLOBAL_ZIP,
      "source_official_BIO19_ZIP_SHA256":BIO_ZIP,
      "source_ESA_CCI_DOI":AGB_DOI,
      "ESA_wood_biomass_units":"Mg dry woody mass per hectare",
      "ESA_AGB_years":[2010,2020,2024],
      "AGB_uncertainty_2020_provided":True,
      "ESA_AGB_v7_source_band_details":details,
      "species_with_2020_AGB_cell_samples":int(combined.ESA_AGB_woody_Mg_ha_2020.notna().sum()),
      "species_with_both_AGB_NDVI_BIO19":int((combined.ESA_AGB_woody_Mg_ha_2020.notna()
          &combined.MODIS_NDVI_2001_2021_mean.notna()&combined.BIO19.notna()).sum()),
      "global_AGB_woody_biomass_equals_all_grassland_biomass":False,
      "zero_forest_AGB_does_not_mean_no_vegetation":True,
      "AGB_biomass_is_not_NDVI":True,
      "source_photo_hue_is_not_local_flower_truth":True,
      "species_native_range_verification":False,
      "raw_raster_PNG_pixels_redistributed":False,
      "output_SHA256":{p.name:sha(p) for p in out.iterdir()}}
    (out/"receipt.json").write_text(json.dumps(report,indent=2,sort_keys=True,default=str)+"\n")
    print("MOIS_ACTUAL_ESA_CCI_WOODY_BIOMASS_JOINED_TO_37104_NDVI_WEATHER_SPECIES",
        json.dumps({k:report[k] for k in (
            "source_global_accepted_species","source_100km_cells",
            "species_with_2020_AGB_cell_samples","species_with_both_AGB_NDVI_BIO19"
        )},sort_keys=True),flush=True)
    return report

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--features",required=True);p.add_argument("--agb",required=True)
    p.add_argument("--gbif",required=True);p.add_argument("--worldclim-bio",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    run(a.features,a.agb,a.gbif,a.worldclim_bio,a.output)
