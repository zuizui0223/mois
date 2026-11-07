#!/usr/bin/env python3
"""Original WorldPop R2024B UC gridded population 2019,2021,2022,2023.

One year and ONE raster release per output, no CN/UC splicing.
Converts population counts per pixel to density using WGS84 geodesic area.
No original photo downloads, white imputation or nighttime moth inference.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from pathlib import Path
import numpy as np,pandas as pd,rasterio
from rasterio.windows import bounds,Window
from pyproj import Geod

SOURCES={
2019:("01f4e2d2e1c8bc9a1392f4143f693563bfa8ea1a0f65ba800e1caa450803831c",27020),
2021:("f0d7cb134c09c74cccfb69199396ac30bc5b663d39dc5369d6e296d8c3a1cdca",45966),
2022:("594fe8e729cf54aa44859a841b5d3e49e6a66fb85b217da8aff0160044e14df5",52075),
2023:("e2825386f45a77da6e3dacdf714fc940cf9ab9cc5174671672b4d57466d04fba",59862)}
def require(ok,msg):
    if not bool(ok):raise ValueError("MOIS_R2024B_POP: "+msg)
def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda:f.read(4<<20),b""):h.update(b)
    return h.hexdigest()
def pixel_area(ds,row,geod):
    west,south,east,north=bounds(Window(0,int(row),1,1),ds.transform)
    area,_=geod.polygon_area_perimeter([west,east,east,west],[south,south,north,north])
    km2=abs(area)/1e6
    require(math.isfinite(km2) and km2>0,"Nonsensical WGS84 pixel area")
    return km2
def sample_population(ds,source,year):
    require(year in SOURCES and len(source)==SOURCES[year][1] and
            source.target_id.is_unique and source.target_id.str.startswith(str(year)+"|").all(),
            "Wrong original source population/year target")
    require(ds.count==1 and ds.crs.to_epsg()==4326 and ds.width==43200 and ds.height==17280,
            "Not expected original 30 arcseconds global population raster")
    require(abs(ds.res[0]-1/120)<1e-6 and abs(ds.res[1]-1/120)<1e-6,
            "Wrong global raster resolution")
    locations=list(zip(source.query_longitude.astype(float),source.query_latitude.astype(float)))
    geod=Geod(ellps="WGS84");areas={}
    counts=[];status=[];km2=[]
    for (lon,lat),sample in zip(locations,ds.sample(locations,masked=True)):
        if not(ds.bounds.left<=lon<ds.bounds.right and ds.bounds.bottom<lat<=ds.bounds.top):
            counts.append(np.nan);status.append("OUTSIDE_RASTER");km2.append(np.nan);continue
        if np.ma.is_masked(sample[0]):
            counts.append(np.nan);status.append("NODATA");km2.append(np.nan);continue
        v=float(sample[0])
        if not math.isfinite(v) or (ds.nodata is not None and v==ds.nodata):
            counts.append(np.nan);status.append("NODATA");km2.append(np.nan);continue
        if v<0:
            counts.append(np.nan);status.append("NEGATIVE_INVALID");km2.append(np.nan);continue
        row,col=ds.index(lon,lat)
        require(0<=row<ds.height and 0<=col<ds.width,"Raster index outside supported source")
        if row not in areas:areas[row]=pixel_area(ds,row,geod)
        counts.append(v);status.append("OK");km2.append(areas[row])
    require(len(counts)==len(source),"Original photo-site observations lost")
    arr=np.asarray(counts,dtype=float)
    area=np.asarray(km2,dtype=float)
    densities=np.divide(arr,area,out=np.full(len(arr),np.nan),where=np.isfinite(arr)&np.isfinite(area))
    out=pd.DataFrame({"target_id":source.target_id,"land_class":source.land_class,
        "source_photo_year":year,"source_population_model":"WorldPop Global2 R2024B UC v1",
        "population_source_year":year,"population_status":status,
        "population_count_per_pixel":arr,"population_pixel_geodesic_km2":area,
        "population_density_per_km2":densities,
        "n_existing_photo_opportunities":source.n_attempted,
        "n_original_strict9_success":source.n_strict9})
    require(out.target_id.is_unique,"Population pixels replicated original target ID")
    return out,{"valid":int(out.population_status.eq("OK").sum()),
                "nodata":int(out.population_status.eq("NODATA").sum()),
                "outside":int(out.population_status.eq("OUTSIDE_RASTER").sum()),
                "negative_invalid":int(out.population_status.eq("NEGATIVE_INVALID").sum()),
                "population_real_zero":int((out.population_status.eq("OK")&
                                             out.population_count_per_pixel.eq(0)).sum()),
                "unique_latitude_pixel_areas":len(areas)}
def run(year,sites,raster,out):
    out=Path(out);require(not out.exists(),"Refuse overwriting a signed source")
    require(year in SOURCES,"Year not supported")
    src_sha,n=SOURCES[year]
    require(sha(sites)==src_sha,"Original source photo-cell SHA changed")
    tif=Path(raster)
    require(tif.name==f"global_pop_{year}_UC_1km_R2024B_v1.tif",
            "Year/raster name does not match requested official population release")
    pop_sha=sha(tif)
    d=pd.read_csv(sites,dtype={"target_id":str},low_memory=False)
    with rasterio.open(tif) as ds:
        result,quality=sample_population(ds,d,year)
        metadata={"crs":str(ds.crs),"width":ds.width,"height":ds.height,
                  "nodata":ds.nodata,"resolution":list(ds.res)}
    out.mkdir(parents=True)
    file=out/f"2024B_UC_{year}_worldpop_sameyear_photo_landmass_popdensity.csv"
    result.to_csv(file,index=False,float_format="%.12g")
    receipt={"schema":"MOIS_2019_2023_R2024B_UC_original_population_source_v1",
      "source_year":year,"original_photo_landmass_sites":n,
      "source_photo_site_manifest_SHA256":src_sha,
      "official_Global2_R2024B_UC_raster_name":tif.name,
      "official_original_raster_SHA256":pop_sha,
      "official_original_raster_bytes":tif.stat().st_size,
      "population_model_release":"WorldPop Global2 R2024B UC v1",
      "source_pixel_units":"estimated people per 30-arcsec grid pixel",
      "density_conversion":"WGS84 geodesic pixel area km2",
      "original_photo_year_equals_population_model_year":True,
      "only_one_population_release_in_this_file":True,
      "quality":quality,"raster_metadata":metadata,
      "new_photographs_or_ROI_inferred":0,
      "no_country_island_area_latitude_exclusion":True,
      "no_NODATA_imputed_as_zero":True,
      "VIIRS_will_be_joined_only_by_exact_year_and_target_ID":True,
      "population_values_SHA256":sha(file)}
    (out/"source_receipt.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print("MOIS_REAL_YEAR_R2024B_UC_POPULATION_RASTER",json.dumps({
      k:receipt[k] for k in ("source_year","original_photo_landmass_sites",
        "official_original_raster_bytes","official_original_raster_SHA256","quality",
        "population_values_SHA256")},sort_keys=True),flush=True)
    return receipt

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--year",type=int,required=True,choices=sorted(SOURCES))
    p.add_argument("--sites",type=Path,required=True)
    p.add_argument("--population-raster",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args();run(a.year,a.sites,a.population_raster,a.out)
