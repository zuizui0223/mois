#!/usr/bin/env python3
"""Join all actual GBIF 37,104 accepted species x 100km equal-area cells to WorldClim 2.1.

Data: 2026-11 GBIF S3 public full snapshot; photo nine-palette species traits.
Climate: WorldClim historical 1970-2000 10 arcminute BIO1 temperature (C)
and BIO12 precipitation (mm) sampled at 25 within-cell points per 100km
GBIF cell. This diagnoses climate of a recorded species' GBIF occupancy,
NOT local flower-petal colour, native-only range or climate adaptation.
"""
from __future__ import annotations
import argparse,hashlib,io,json,zipfile,math
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

SOURCE_GBIF_ZIP_SHA="422f6d6489e61c6e57c5d51db614e47f613039781a901e3fd3481c13b2b085c2"
SOURCE_GRID_PARQUET_SHA="e00b0408055984abfe4ec0ff0a24833fc148888afd1ea0e0c9b270a077d08b11"
SOURCE_GEO_ZIP_SHA="63a260c0a9c68d5c7198c6f4a821bd391eae25954a31cf060e6b1fbac544e313"
WORLDCLIM_URL="https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_10m_bio.zip"
WORLDCLIM_BASELINE="1970-2000"
HUES=("white","yellow","orange","red","pink","magenta","purple","blue","bronze")
OFFSETS=(-0.4,-0.2,0.0,0.2,0.4)

def require(ok,msg):
    if not bool(ok):raise ValueError("MOIS_WORLDCLIM_GLOBAL_PHOTO9: "+str(msg))
def sha(b):
    if isinstance(b,(str,Path)):b=Path(b).read_bytes()
    return hashlib.sha256(b).hexdigest()

def extract_tif_samples(values,lons,lats,x0,y0,dx,dy,nodata):
    a=np.asarray(values)
    lon=np.asarray(lons,float);lat=np.asarray(lats,float)
    require(lon.shape==lat.shape,"lat/lon source grid mismatch")
    out=np.full(lon.shape,np.nan,dtype=float)
    valid=(np.isfinite(lon)&np.isfinite(lat)&(lon>=-180)&(lon<180)&(lat>=-90)&(lat<=90))
    good=np.flatnonzero(valid)
    if not len(good):return out
    x=np.floor((lon[good]-x0)/dx).astype(np.int64)
    y=np.floor((lat[good]-y0)/dy).astype(np.int64)
    inr=(x>=0)&(x<a.shape[1])&(y>=0)&(y<a.shape[0])
    indexes=good[inr]
    got=a[y[inr],x[inr]].astype(float)
    if nodata is not None:got[np.isclose(got,nodata,rtol=0,atol=0)]=np.nan
    got[~np.isfinite(got)]=np.nan
    out[indexes]=got
    return out

def demean(frame,group):
    return frame-frame.groupby(group).transform("mean")

def climate_coverage_status(valid,total):
    if total==0:return "no_eligible_GBIF_grid"
    if valid==0:return "not_sampled_coastal_or_missing"
    return "complete" if valid==total else "partial"

def signed_sources(gbifzip,geojsonzip):
    import pyarrow.parquet as pq
    require(sha(gbifzip)==SOURCE_GBIF_ZIP_SHA,"full 2026 November GBIF output ZIP changed")
    require(sha(geojsonzip)==SOURCE_GEO_ZIP_SHA,"full recorded-range colour+latitude ZIP changed")
    with zipfile.ZipFile(gbifzip) as z:
        gridbytes=z.read("all_37104_species_x_EPSG6933_100km_cells.parquet")
        full=json.loads(z.read("receipt.json"))
    require(sha(gridbytes)==SOURCE_GRID_PARQUET_SHA and
            full["all_16_shards_source_verified"] and
            full["GBIF_accepted_species_with_usable_occurrences"]==36758,
            "original full GBIF source lacks signed completeness")
    with zipfile.ZipFile(geojsonzip) as z:
        metadata=json.loads(z.read("receipt.json"))
        name="all_37104_accepted_species_global_GBIF_photo9_georange_latitudes.csv.gz"
        b=z.read(name)
    require(sha(b)==metadata["output_SHA256"][name],
            "source photo-colour geography file SHA differs")
    taxa=pd.read_csv(io.BytesIO(b),compression="gzip",dtype={"key":str,"inat_taxon_id":str},
                     low_memory=False)
    require(len(taxa)==37104 and taxa.key.is_unique and set(taxa.source_hue9)==set(HUES)
            and taxa.records_passed_spatial_QA.sum()==396142825,
            "original photo9+global range roster drift")
    cells=pq.read_table(io.BytesIO(gridbytes),
                columns=["specieskey","grid_x","grid_y","occurrence_count"]).to_pandas()
    require(len(cells)==5493931 and cells.occurrence_count.sum()==396142825,
            "full GBIF spatial cell array missing or incomplete")
    return taxa,cells

def raster_worldclim(archive):
    import rasterio
    from rasterio.io import MemoryFile
    raw=Path(archive).read_bytes()
    require(1000000<len(raw)<100000000,"WorldClim 10 arcmin zip unexpected size")
    images={}
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        names=z.namelist()
        require(len([n for n in names if n.endswith(".tif")])==19,
                "expected exactly 19 WorldClim BIO GeoTIFF rasters")
        for index in (1,12):
            member=f"wc2.1_10m_bio_{index}.tif"
            require(member in names,"WorldClim expected BIO raster missing")
            blob=z.read(member)
            with MemoryFile(blob) as mf:
                with mf.open() as ds:
                    require(ds.count==1 and str(ds.crs).upper()=="EPSG:4326",
                            "WorldClim tif CRS or band changed")
                    require(ds.transform.a>0 and ds.transform.e<0 and
                            ds.transform.b==0 and ds.transform.d==0,
                            "WorldClim tif geometry not north-up")
                    images[index]={"data":ds.read(1),
                      "nodata":ds.nodata,"w":ds.width,"h":ds.height,
                      "x0":float(ds.transform.c),"y0":float(ds.transform.f),
                      "dx":float(ds.transform.a),"dy":float(ds.transform.e),
                      "original_SHA256":sha(blob),"raster_dtype":str(ds.dtypes[0]),
                      "raster_crs":str(ds.crs)}
    require(images[1]["data"].shape==images[12]["data"].shape,
            "WorldClim BIO rasters not aligned")
    return images,sha(raw),len(raw)

def sampled_climate_cells(grid,images):
    from pyproj import Transformer
    unique=grid[["grid_x","grid_y"]].drop_duplicates().reset_index(drop=True)
    require(unique.shape[0]>1000 and unique.grid_x.notna().all()
            and unique.grid_y.notna().all(),"too few GBIF global terrestrial cells")
    n=len(unique)
    x=unique.grid_x.to_numpy(float)
    y=unique.grid_y.to_numpy(float)
    sums_t=np.zeros(n,float);sums_p=np.zeros(n,float)
    counts=np.zeros(n,int)
    transform=Transformer.from_crs("EPSG:6933","EPSG:4326",always_xy=True)
    for u in OFFSETS:
        for v in OFFSETS:
            lon,lat=transform.transform((x+0.5+u)*100000,(y+0.5+v)*100000)
            t=extract_tif_samples(images[1]["data"],lon,lat,
                *[images[1][k] for k in ("x0","y0","dx","dy")],
                nodata=images[1]["nodata"])
            p=extract_tif_samples(images[12]["data"],lon,lat,
                *[images[12][k] for k in ("x0","y0","dx","dy")],
                nodata=images[12]["nodata"])
            good=(np.isfinite(t)&np.isfinite(p)&(t>=-100)&(t<=80)
                  &(p>=0)&(p<=30000))
            sums_t+=np.where(good,t,0.)
            sums_p+=np.where(good,p,0.)
            counts+=good.astype(int)
    unique["valid_25_point_samples"]=counts
    unique["cell_sampling_fraction"]=counts/25
    unique["mean_annual_temperature_C"]=np.divide(sums_t,counts,
         out=np.full(n,np.nan),where=counts>0)
    unique["annual_precipitation_mm"]=np.divide(sums_p,counts,
         out=np.full(n,np.nan),where=counts>0)
    unique["climate_valid_cell"]=counts>0
    unique["climate_source_is_GBIF_occurrence_year"]=False
    return unique

def robust_within(d,response,*,add_precip=False,add_temp=False,genus_equal=False):
    # Independent units are source species, never repeated species×grid rows.
    y=d[response].astype(float)
    x=pd.DataFrame(index=d.index)
    x["photo_white"]=d.source_hue9.eq("white").astype(float)
    x["log_GBIF_occ"]=np.log1p(d.records_passed_spatial_QA)
    x["mean_abs_latitude"]=d.mean_abs_latitude_over_occupied_cells.astype(float)
    if add_precip:x["log_BIO12_precip"]=np.log1p(d.annual_precipitation_mm)
    if add_temp:x["BIO1_temperature_C"]=d.mean_annual_temperature_C
    for status in sorted(d.hue_status.unique())[1:]:
        x["photo_status_"+status]=d.hue_status.eq(status).astype(float)
    grp=d.genus
    x=demean(x,grp)
    y=demean(y.to_frame("response"),grp)["response"]
    if genus_equal:
        w=1/d.groupby("genus")["genus"].transform("size").to_numpy(float)
        return sm.WLS(y,x,weights=w).fit(cov_type="cluster",cov_kwds={"groups":grp},use_t=False)
    return sm.OLS(y,x).fit(cov_type="cluster",cov_kwds={"groups":grp},use_t=False)

def beta(m,key="photo_white"):
    c=m.conf_int().loc[key]
    return {"beta":float(m.params[key]),"95CI":[float(c.iloc[0]),float(c.iloc[1])],
            "p_HC3_genus_cluster":float(m.pvalues[key])}

def run(gbifzip,geozip,worldclimzip,out):
    out=Path(out);require(not out.exists(),"refuse overwrite")
    taxa,grid=signed_sources(gbifzip,geozip)
    rasters,raster_sha,raster_size=raster_worldclim(worldclimzip)
    cells=sampled_climate_cells(grid,rasters)
    # Coastal/marine nulls are NOT zero temperature, zero precipitation, or
    # biological range absence. Missing species remain in the full output.
    climate=grid[["specieskey","grid_x","grid_y"]].merge(
        cells,on=["grid_x","grid_y"],how="left",validate="many_to_one")
    require(len(climate)==len(grid),"species grid duplicates after climate join")
    agg=climate.groupby("specieskey").agg(
        observed_grid_cells=("grid_y","size"),
        sampled_climate_cells=("climate_valid_cell","sum"),
        mean_annual_temperature_C=("mean_annual_temperature_C","mean"),
        annual_precipitation_mm=("annual_precipitation_mm","mean"),
        mean_cell_land_sample_fraction=("cell_sampling_fraction","mean"))
    agg["climate_coverage_fraction"]=agg.sampled_climate_cells/agg.observed_grid_cells
    d=taxa.merge(agg,left_on="key",right_index=True,how="left",validate="one_to_one")
    d["observed_grid_cells"]=d.observed_grid_cells.fillna(0).astype(int)
    d["sampled_climate_cells"]=d.sampled_climate_cells.fillna(0).astype(int)
    d["climate_coverage_fraction"]=d.climate_coverage_fraction.fillna(0)
    require((d.observed_grid_cells==d.occupied_grid_cells).all()
            and len(d)==37104 and d.key.is_unique,
            "GBIF species climate coverage no longer matches source")
    d["climate_coverage_status"]=[
        climate_coverage_status(v,n) for v,n in zip(
            d.sampled_climate_cells,d.observed_grid_cells)]
    primary=d[(d.observed_grid_cells.gt(0))
         &(d.climate_coverage_fraction>=0.75)
         &d.mean_abs_latitude_over_occupied_cells.notna()
         &np.isfinite(d.mean_annual_temperature_C)
         &np.isfinite(d.annual_precipitation_mm)].copy()
    require(len(primary)>10000 and primary.genus.nunique()>2000,
            "WorldClim shared eligible species too small for global inference")
    primary["log_annual_precip"]=np.log1p(primary.annual_precipitation_mm)
    mixed=primary.groupby("genus").source_hue9.apply(
        lambda x:x.eq("white").any() and x.ne("white").any())
    matched=primary[primary.genus.map(mixed)].copy()
    require(len(matched)>5000,"not enough within-genus white/non-white species")
    mtemp=robust_within(primary,"mean_annual_temperature_C",
                       add_precip=True)
    mprec=robust_within(primary,"log_annual_precip",add_temp=True)
    meqtemp=robust_within(matched,"mean_annual_temperature_C",
                         add_precip=True,genus_equal=True)
    meqprec=robust_within(matched,"log_annual_precip",
                         add_temp=True,genus_equal=True)
    p=[float(mtemp.pvalues["photo_white"]),float(mprec.pvalues["photo_white"])]
    q=multipletests(p,method="fdr_bh")[1]
    # Record full per-hue coverage and climate niche proxies, not petal adaptation.
    hue=primary.groupby("source_hue9").agg(
        species=("key","size"),
        median_temperature_C=("mean_annual_temperature_C","median"),
        median_annual_precip_mm=("annual_precipitation_mm","median"),
        median_GBIF_grid_cells=("occupied_grid_cells","median"),
        median_GBIF_climate_coverage=("climate_coverage_fraction","median")).reset_index()
    out.mkdir(parents=True)
    cells.to_csv(out/"WorldClim_2p1_10m_BIO1_BIO12_100km_cell_subsamples.csv.gz",
        index=False,compression={"method":"gzip","mtime":0})
    cols=["key","inat_taxon_id","species","genus","source_hue9","hue_status",
          "occupied_grid_cells","records_passed_spatial_QA",
          "mean_abs_latitude_over_occupied_cells","modal_mean9_disagree",
          "observed_grid_cells","sampled_climate_cells","climate_coverage_fraction",
          "climate_coverage_status","mean_annual_temperature_C","annual_precipitation_mm"]
    d[cols].to_csv(out/"COLXR_37104_species_GBIF_grid_WorldClim_species_climate.csv.gz",
        index=False,compression={"method":"gzip","mtime":0})
    hue.to_csv(out/"photo9_hue_WorldClim_BIO1_BIO12_global_descriptive.csv",index=False)
    receipt={"schema":"MOIS_37104_GBIF_100km_WorldClim21_climate_photo9_genus_FE_v1",
      "full_anonymous_public_GBIF_20261101_10346_files_verified":True,
      "accepted_COLXR_37104_species":37104,
      "source_global_GBIF_zip_SHA256":SOURCE_GBIF_ZIP_SHA,
      "source_photo9_geography_zip_SHA256":SOURCE_GEO_ZIP_SHA,
      "GBIF_valid_occurrence_count":396142825,
      "GBIF_species_100km_pairs":5493931,
      "distinct_GBIF_grid_cells":len(cells),
      "WorldClim_source_url":WORLDCLIM_URL,
      "WorldClim_source_ZIP_SHA256":raster_sha,
      "WorldClim_source_archive_bytes":raster_size,
      "WorldClim_BIO1_SHA256":rasters[1]["original_SHA256"],
      "WorldClim_BIO12_SHA256":rasters[12]["original_SHA256"],
      "WorldClim_resolution":"10 arcminutes",
      "WorldClim_historical_baseline":"1970-2000 (not source photo year)",
      "WorldClim_mean_temperature_units":"deg C (BIO1)",
      "WorldClim_annual_precipitation_units":"mm (BIO12)",
      "within_cell_climate_subsamples":25,
      "mean_climate_missing_coastal_grid_cells":int((~cells.climate_valid_cell).sum()),
      "species_with_any_valid_climate_cells":int(d.sampled_climate_cells.gt(0).sum()),
      "species_primary_climate_coverage_75percent":len(primary),
      "species_in_mixed_colored_genera_primary":len(matched),
      "primary_genus_count":int(primary.genus.nunique()),
      "matched_genus_count":int(matched.genus.nunique()),
      "models_primary":{
          "BIO1_within_genus_white_control_precip_lat_GBIFeffort_photoStatus":
             {**beta(mtemp),"q_BH_2_climate_endpoints":float(q[0])},
          "log_BIO12_within_genus_white_control_temp_lat_GBIFeffort_photoStatus":
             {**beta(mprec),"q_BH_2_climate_endpoints":float(q[1])},
          "BIO1_equal_weight_mixed_genus_sensitivity":beta(meqtemp),
          "log_BIO12_equal_weight_mixed_genus_sensitivity":beta(meqprec)},
      "methodological_boundaries":{
        "BIO_values_at_GBIF_occurrence_exact_coordinates":False,
        "BIO_values_at_all_100km_cell_land_sites":False,
        "source_photo9_is_true_local_petal_colour":False,
        "native_status_verified":False,
        "GBIF_records_are_observer_effort_adjusted_causes":False,
        "genus_FE_is_full_phylogenetic_correction":False,
        "climate_exposure_causes_flower_colour":False,
        "climate_observations_match_photo_date":False},
      "output_SHA256":{p.name:sha(p) for p in out.iterdir()}}
    (out/"receipt.json").write_text(json.dumps(receipt,sort_keys=True,indent=2)+"\n")
    print("MOIS_GBIF_FULL37104_WORLDCLIM_BIO1_BIO12_SPECIENICHE_AUDIT",
       json.dumps({k:receipt[k] for k in (
         "distinct_GBIF_grid_cells","mean_climate_missing_coastal_grid_cells",
         "species_primary_climate_coverage_75percent",
         "models_primary")},sort_keys=True),flush=True)
    return receipt

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--GBIF-global-zip",required=True)
    p.add_argument("--photo9-global-zip",required=True)
    p.add_argument("--WorldClim-bio10m-zip",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    run(a.GBIF_global_zip,a.photo9_global_zip,a.WorldClim_bio10m_zip,a.output)
