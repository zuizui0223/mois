#!/usr/bin/env python3
"""MOIS real 37,104-species GBIF + FULL 19 BIOCLIM + monthly weather + NDVI.

Each training unit is an accepted plant SPECIES, never a repeated species-grid row.
WorldClim v2.1 is 1970–2000 climate normals; NDVI is MODIS-derived 2001–2021.
Both are associated with GBIF *recorded habitat cells*, not species-specific petal
colour at those cells. NDVI is GRENESS, not biomass/NPP.
"""
from __future__ import annotations
import argparse,hashlib,io,json,math,zipfile,warnings
from pathlib import Path
import numpy as np
import pandas as pd

GLOBAL_ZIP="422f6d6489e61c6e57c5d51db614e47f613039781a901e3fd3481c13b2b085c2"
GEO_ZIP="63a260c0a9c68d5c7198c6f4a821bd391eae25954a31cf060e6b1fbac544e313"
NDVI_ZIP="eb6e44a22696f3098a0bdf9e5c28e346556b9b98435853812cd594190bc8f3d6"
NDVI_ARCHIVE_SHA="acc4646537e7c80e61d8d19f426eb2f556fe5aaa85b0a6763ae8558f3a35ce58"
BIO_ARCHIVE_SHA="00513224583665ec0f2f955a4ec252730c4deb2004cce9e793492a3f26df4dcf"
GRID_PARQUET_SHA="e00b0408055984abfe4ec0ff0a24833fc148888afd1ea0e0c9b270a077d08b11"
VARIABLES=("bio","prec","tavg","srad","vapr","wind")
OFFSETS=(-.4,-.2,0.,.2,.4)
HUES=("white","yellow","orange","red","pink","magenta","purple","blue","bronze")

def require(cond,msg):
    if not bool(cond):raise ValueError("MOIS_FULL_ENV_RF_FEATURES: "+str(msg))
def sha(src):
    if isinstance(src,(str,Path)):src=Path(src).read_bytes()
    return hashlib.sha256(src).hexdigest()

def get_original_sources(global_zip,geo_zip,ndvi_zip):
    import pyarrow.parquet as pq
    require(sha(global_zip)==GLOBAL_ZIP and sha(geo_zip)==GEO_ZIP and sha(ndvi_zip)==NDVI_ZIP,
        "signed GBIF, photo and NDVI artifact ZIPs changed")
    with zipfile.ZipFile(global_zip) as z:
        raw=z.read("all_37104_species_x_EPSG6933_100km_cells.parquet")
        receipt=json.loads(z.read("receipt.json"))
    require(sha(raw)==GRID_PARQUET_SHA and receipt["all_16_shards_source_verified"]
        and receipt["complete_public_S3_object_count"]==10346
        and receipt["GBIF_occurrences_passed_all_filters"]==396142825,
        "global GBIF source not authenticated by complete anonymous 16-shard audit")
    pairs=pq.read_table(io.BytesIO(raw),
        columns=["specieskey","grid_x","grid_y","occurrence_count"]).to_pandas()
    require(len(pairs)==5493931 and pairs.occurrence_count.sum()==396142825,
            "source species x grid count changed")
    with zipfile.ZipFile(geo_zip) as z:
        r=json.loads(z.read("receipt.json"))
        name="all_37104_accepted_species_global_GBIF_photo9_georange_latitudes.csv.gz"
        source=z.read(name)
    require(sha(source)==r["output_SHA256"][name],"photo9 species level geographic source changed")
    taxa=pd.read_csv(io.BytesIO(source),compression="gzip",dtype={"key":str,"inat_taxon_id":str},
        low_memory=False)
    require(len(taxa)==37104 and taxa.key.is_unique
      and taxa.records_passed_spatial_QA.sum()==396142825
      and set(taxa.source_hue9)==set(HUES),"source traits taxon universe not complete")
    with zipfile.ZipFile(ndvi_zip) as z:
        zreceipt=json.loads(z.read("receipt.json"))
        nested=z.read("ndvi_annual_climatology.zip")
    require(sha(nested)==NDVI_ARCHIVE_SHA and zreceipt["downloaded_NDVI_raster"] is True
      and zreceipt["target_DOI"]=="10.5281/zenodo.22967479",
      "public global MODIS NDVI original source not provenance-verified")
    return taxa,pairs,nested

def sample_lonlat(unique):
    from pyproj import Transformer
    proj=Transformer.from_crs("EPSG:6933","EPSG:4326",always_xy=True)
    x=unique.grid_x.to_numpy(float)
    y=unique.grid_y.to_numpy(float)
    out=[]
    for a in OFFSETS:
        for b in OFFSETS:
            lon,lat=proj.transform((x+.5+a)*100000,(y+.5+b)*100000)
            out.append((np.asarray(lon,float),np.asarray(lat,float)))
    require(len(out)==25,"25 sample points not created")
    return out

def raster_probe(buffer,geometry,coords):
    from rasterio.io import MemoryFile
    with MemoryFile(buffer) as m:
        with m.open() as ds:
            require(ds.count==1 and str(ds.crs).upper()=="EPSG:4326"
                and ds.transform.a>0 and ds.transform.e<0
                and ds.transform.b==0 and ds.transform.d==0,
                "not a north-up EPSG4326 raster")
            arr=ds.read(1)
            values=np.full((len(coords),len(coords[0][0])),np.nan,float)
            for j,(lon,lat) in enumerate(coords):
                valid=(np.isfinite(lon)&np.isfinite(lat)
                    &(lon>=ds.bounds.left)&(lon<ds.bounds.right)
                    &(lat>=ds.bounds.bottom)&(lat<ds.bounds.top))
                inds=np.flatnonzero(valid)
                if not len(inds):continue
                xx=np.floor((lon[inds]-ds.transform.c)/ds.transform.a).astype(int)
                yy=np.floor((lat[inds]-ds.transform.f)/ds.transform.e).astype(int)
                ok=(xx>=0)&(yy>=0)&(xx<arr.shape[1])&(yy<arr.shape[0])
                vv=arr[yy[ok],xx[ok]].astype(float)
                if ds.nodata is not None and np.isfinite(ds.nodata):
                    vv[vv==ds.nodata]=np.nan
                vv[~np.isfinite(vv)]=np.nan
                values[j,inds[ok]]=vv
            geometry[str(geometry.get("id","raster"))]={
                "height":ds.height,"width":ds.width,"dtype":str(ds.dtypes[0]),
                "bounds":tuple(ds.bounds),"crs":str(ds.crs),
                "nodata":str(ds.nodata)}
    return values

def point_mean(values):
    valid=np.isfinite(values)
    count=valid.sum(axis=0)
    out=np.divide(np.where(valid,values,0).sum(axis=0),count,
        out=np.full(values.shape[1],np.nan),where=count>0)
    return out,count

def feature_cells(unique,ndvi_nested,worldclim_files):
    from rasterio.io import MemoryFile
    coords=sample_lonlat(unique)
    data=unique.copy()
    raster_receipts={};metadata={}
    for var in VARIABLES:
        path=Path(worldclim_files[var])
        raw=path.read_bytes()
        if var=="bio":require(sha(raw)==BIO_ARCHIVE_SHA,"WorldClim19 source SHA drift")
        roster=[]
        with zipfile.ZipFile(io.BytesIO(raw)) as z:
            expect=range(1,20) if var=="bio" else range(1,13)
            all_means=[]
            for month in expect:
                filename=f"wc2.1_10m_{var}_{month:02d}.tif" if var!="bio" else f"wc2.1_10m_bio_{month}.tif"
                require(filename in z.namelist(),"WorldClim layer absent "+filename)
                tif=z.read(filename)
                temporary={"id":filename}
                values=raster_probe(tif,temporary,coords)
                mean,count=point_mean(values)
                if var=="bio":
                    # Data values as supplied by official 1970–2000 WorldClim.
                    data[f"BIO{month:02d}"]=mean
                    if month==1:data["BIO1_valid_subsamples"]=count
                else:
                    all_means.append(mean)
                roster.append({"tif":filename,"SHA256":sha(tif),
                    "valid_cells":int(np.isfinite(mean).sum())})
            if var!="bio":
                monthly=np.vstack(all_means)
                count=np.isfinite(monthly).sum(axis=0)
                sums=np.nansum(monthly,axis=0)
                av=np.divide(sums,count,out=np.full(len(count),np.nan),where=count>0)
                data[f"{var}_monthly_mean"]=av
                # Seasonal amplitude is observed max-min across 12 climate normal months.
                tmp=np.where(np.isfinite(monthly),monthly,np.nan)
                with warnings.catch_warnings():
                    warnings.filterwarnings("ignore",category=RuntimeWarning)
                    data[f"{var}_monthly_range"]=np.nanmax(tmp,axis=0)-np.nanmin(tmp,axis=0)
                    data[f"{var}_monthly_sd"]=np.nanstd(tmp,axis=0)
                if var=="prec":
                    data["prec_annual_sum_of_monthly_normal_mm"]=sums
                if var=="tavg":
                    data["tavg_warmest_month_C"]=np.nanmax(tmp,axis=0)
        raster_receipts[var]={"source_url":f"https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_10m_{var}.zip",
                             "zip_SHA256":sha(raw),"layers":roster}
    # Mean annual vapor pressure deficit: subtract measured water vapour pressure from
    # saturation pressure at WorldClim monthly-normal mean temperature.
    T=data.tavg_monthly_mean.to_numpy(float)
    saturated=0.6108*np.exp((17.27*T)/(T+237.3))
    data["vpd_kPa_approx"]=np.maximum(0.,saturated-data.vapr_monthly_mean.to_numpy(float))
    # Satellite NDVI from 2001–2021 (not WorldClim 1970–2000); source raster
    # encoded at 0.01 scale (01C in filenames). NDVI is not standing biomass.
    ndvi_members={
      "mean":"ANNUAL_NORMALS/NASA_MODIS_NDVI_MEAN_10K_01C_2001U2021.tif",
      "sd":"ANNUAL_NORMALS/NASA_MODIS_NDVI_STD_10K_01C_2001U2021.tif"}
    with zipfile.ZipFile(io.BytesIO(ndvi_nested)) as z:
        for stat,member in ndvi_members.items():
            require(member in z.namelist(),"NDVI mean/std source missing")
            pixels=z.read(member)
            values=raster_probe(pixels,{"id":"MODIS_NDVI_"+stat},coords)
            # source stores NDVI ×100 as float, with valid mean -1..1 (incl sparse water).
            values=values/100.
            values[(values<(-1.0001))|(values>1.0001)]=np.nan
            mean,n=point_mean(values)
            data["MODIS_NDVI_2001_2021_"+stat]=mean
            if stat=="mean":data["NDVI_valid_subsamples"]=n
            raster_receipts["MODIS_NDVI_"+stat]={"source_member":member,"tif_SHA256":sha(pixels),
              "valid_cells":int(np.isfinite(mean).sum()),
              "units":"dimensionless NDVI, source numeric scaled by 0.01"}
    require(len(data)==len(unique) and data[["grid_x","grid_y"]].duplicated().sum()==0,
            "global WorldClim and MODIS cell feature keys changed")
    return data,raster_receipts

def species_features(taxa,pairs,envcells):
    # Species-level aggregation: one row per accepted species, occupied cells equally
    # weighted. Never let abundance of S3 occurrence records multiply training rows.
    joined=pairs[["specieskey","grid_x","grid_y"]].merge(
        envcells,on=["grid_x","grid_y"],how="left",validate="many_to_one")
    require(len(joined)==len(pairs),"species-grid rows expanded")
    idx,keys=pd.factorize(joined.specieskey,sort=False)
    names=[c for c in envcells.columns if c not in ("grid_x","grid_y")]
    features={"specieskey":keys}
    denom=np.bincount(idx,minlength=len(keys))
    features["all_occupied_100km_cells"]=denom
    for name in names:
        vals=joined[name].to_numpy(float)
        mask=np.isfinite(vals)
        n=np.bincount(idx[mask],minlength=len(keys))
        sums=np.bincount(idx[mask],weights=vals[mask],minlength=len(keys))
        out=np.divide(sums,n,out=np.full(len(keys),np.nan),where=n>0)
        features[name]=out
        if name in ("BIO01","MODIS_NDVI_2001_2021_mean","vpd_kPa_approx"):
            features[name+"_cell_coverage"]=n/np.maximum(denom,1)
    data=taxa.merge(pd.DataFrame(features),how="left",
        left_on="key",right_on="specieskey",validate="one_to_one")
    require(len(data)==37104 and data.key.is_unique
      and data.all_occupied_100km_cells.fillna(0).sum()==len(pairs),
      "species climate source lost keys or duplicated locations")
    return data

def run(gbif,photo,ndvi,worldclim,out):
    out=Path(out);require(not out.exists(),"refuse overwriting earlier source")
    taxa,pairs,nested=get_original_sources(gbif,photo,ndvi)
    unique=pairs[["grid_x","grid_y"]].drop_duplicates().reset_index(drop=True)
    require(len(unique)==18577,"global GBIF cell count drift")
    cell,rasters=feature_cells(unique,nested,worldclim)
    species=species_features(taxa,pairs,cell)
    require("BIO19" in species and "MODIS_NDVI_2001_2021_mean" in species
      and "vpd_kPa_approx" in species,"requested climate+NDVI features missing")
    out.mkdir(parents=True)
    a=out/"WorldClim19_weather_NDVI_18577_global_100km_cells.parquet"
    b=out/"WorldClim19_weather_NDVI_37104_COLXR_photo9_species.parquet"
    cell.to_parquet(a,index=False)
    species.to_parquet(b,index=False)
    cols=[x for x in cell.columns if x not in ("grid_x","grid_y")]
    missing={x:round(float(species[x].isna().mean()),6) for x in cols}
    meta={
      "schema":"mois_37104_real_GBIF_9photo_hue_19BIO_monthly_weather_NDVI_features_v1",
      "global_accepted_species":37104,
      "global_S3_provenance_verified":True,
      "source_GBIF_snapshot_2026_11_01":True,
      "source_species_grid_pairs":len(pairs),
      "unique_100km_cells":len(cell),
      "GBIF_source_ZIP_SHA256":GLOBAL_ZIP,
      "species_photo9_source_ZIP_SHA256":GEO_ZIP,
      "MODIS_NDVI_2001_2021_ZIP_SHA256":NDVI_ZIP,
      "WorldClim_baseline":"1970-2000",
      "WorldClim_all_19_BIO":True,
      "WorldClim_monthly_variables":["prec","tavg","srad","vapr","wind"],
      "NDVI_is_biomass":False,
      "GPP_NPP_or_aboveground_biomass_features_in_this_stage":False,
      "raster_sources":rasters,
      "species_feature_missing_fraction":missing,
      "species_1plus_valid_NDVI_cell":int(species["MODIS_NDVI_2001_2021_mean"].notna().sum()),
      "species_1plus_valid_BIO1_cell":int(species["BIO01"].notna().sum()),
      "species_with_both_NDVI_and_BIO19":int(
          (species["MODIS_NDVI_2001_2021_mean"].notna()
           &species["BIO01"].notna()&species["BIO19"].notna()).sum()),
      "photo_colour_at_GBIF_location_measured":False,
      "native_status_verified":False,
      "GBIF_record_year_matches_WorldClim_year":False,
      "output_SHA256":{p.name:sha(p) for p in out.iterdir()}
    }
    (out/"receipt.json").write_text(json.dumps(meta,indent=2,sort_keys=True)+"\n")
    print("MOIS_FULL19_BIO_WEATHER_MODIS_NDVI_37104_SPECIES_REAL_FEATURES",
      json.dumps({k:meta[k] for k in (
       "global_accepted_species","unique_100km_cells",
       "species_1plus_valid_NDVI_cell","species_with_both_NDVI_and_BIO19",
       "WorldClim_monthly_variables")},sort_keys=True),flush=True)
    return meta

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--gbif",required=True);p.add_argument("--photo",required=True)
    p.add_argument("--ndvi",required=True);p.add_argument("--bio",required=True)
    p.add_argument("--prec",required=True);p.add_argument("--tavg",required=True)
    p.add_argument("--srad",required=True);p.add_argument("--vapr",required=True)
    p.add_argument("--wind",required=True);p.add_argument("--output",required=True)
    a=p.parse_args()
    run(a.gbif,a.photo,a.ndvi,{v:getattr(a,v) for v in VARIABLES},a.output)
