#!/usr/bin/env python3
"""MOIS: Island / elevation / Arctic / hemisphere gradients of photo9 assemblage.

This is a macroecological *geographical stress test* of the existing
37,104-photo-hue species × 5,493,931 GBIF occupied 100km grid pairs.
Outcome is the identical observed 9-colour Shannon minus the same strict
global 199-permutation genus × species occupancy-width null, not pixels.
No new photocolour label, local petal observation, or native classification.

Landmass is determined from original GSHHG 2.3.7 high-res L1 coastline
(verified source SHA) using 9 sampled points in the EPSG:6933 100km grid.
Only a dominant same sibling at >=7/9 subpoints is assigned (strict=9/9);
ocean/mixed cells stay unassigned. This only resolves large/majority islands,
not all small islands. Elevation is actual current WorldClim v2.1 10arcminute
global elevation sampled at the same points, not inferred from temperature.
Arctic circle is a purely astronomical boundary north of 66.562 degrees
N, not a full floristic Arctic ecoregion. Latitude and landmass class are
spatially confounded with many different biogeographic species pools.

Whole-catalogue genus×breadth labels are shuffled identically for every
spatial cell, 199 times exactly as the previous signed null. Group-level
permutation statistics use these WHOLE-catalogue permutations, not treating
12k neighbouring 100km cells as independent replicates.
"""
from __future__ import annotations
import argparse
import hashlib
import io
import json
import tempfile
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

from scripts.analysis.mois_global_hue9_community_assembly_genus_null_20261118 import (
    signed_original,baseline,counts_by_cell,empirical_entropy,SEED,GBIF_SHA,
    ENV_SHA,check,sha,
)
from scripts.analysis.mois_nonwhite8_assembly_white_fixed_null_20261118 import (
    load_strict_nine_null,NULL_SHA,
)

GSHHG_SHA="8dbbe7e071e77e9e75f2d639239099ebca8d5c16d6a07df8169729d49f15cf41"
GSHHG_VERSION="2.3.7"
GRID_METERS=100000
ARC_LAT=66.562
HUES=("white","yellow","orange","red","pink","magenta","purple","blue","bronze")
SEEDS={"africa":(20,0),"eurasia":(80,50),
       "north_america":(-100,40),"south_america":(-60,-15),
       "australia":(135,-25)}
SUBGRID=(.20,.50,.80)
N_PERM=199

def load_coast_and_continents(zippath,work):
    import geopandas as gpd
    from shapely.geometry import Point
    check(sha(zippath)==GSHHG_SHA,
          "Original GSHHG 2.3.7 coastline ZIP SHA differs from previous MOIS study")
    copied=[]
    with zipfile.ZipFile(zippath) as z:
        for item in z.namelist():
            filename=Path(item).name
            if (filename.startswith("GSHHS_h_L1.") and
                filename.split(".")[-1].lower() in ("shp","shx","dbf","prj","cpg")):
                dest=work/filename
                dest.write_bytes(z.read(item))
                copied.append(filename)
    check({"GSHHS_h_L1.shp","GSHHS_h_L1.shx","GSHHS_h_L1.dbf"}.issubset(copied),
          "High-resolution GSHHG source shapes absent")
    land=gpd.read_file(work/"GSHHS_h_L1.shp")
    check({"sibling_id","area","geometry"}.issubset(land.columns)
          and land.crs is not None,
          "Original GSHHG high-resolution landmass metadata drifted")
    check(len(land)>5000,"Unexpected sparse GSHHG landmass database")
    ids={}
    for name,(lon,lat) in SEEDS.items():
        pt=Point(lon,lat)
        mask=land.geometry.covers(pt)
        found=land.loc[mask,"sibling_id"].dropna().unique()
        check(len(found)==1,"GSHHG continental seed ambiguous "+name)
        # Geopandas spatial join coerces source integer sibling_id to float
        # because ocean points are NA. Canonical integer text is mandatory.
        ids[name]=str(int(found[0]))
    check(len(set(ids.values()))==5,
          "Five GSHHG continental siblings should be distinct")
    areas=land.groupby("sibling_id",sort=False).area.sum()
    return land,set(ids.values()),ids,areas.to_dict(),len(land)

def cell_coordinates(grid):
    from pyproj import Transformer
    coords=grid[["grid_x","grid_y"]].to_numpy(np.int64)
    n=len(coords)
    frac=np.array([(a,b) for a in SUBGRID for b in SUBGRID],dtype=float)
    shape=np.repeat(np.arange(n,dtype=np.int32),9)
    px=((coords[:,0,None]+frac[:,0])*GRID_METERS).reshape(-1)
    py=((coords[:,1,None]+frac[:,1])*GRID_METERS).reshape(-1)
    t=Transformer.from_crs("EPSG:6933","EPSG:4326",always_xy=True)
    lon,lat=t.transform(px,py)
    lon=np.asarray(lon);lat=np.asarray(lat)
    valid=np.isfinite(lon)&np.isfinite(lat)
    valid &= (lon>=-180)&(lon<=180)&(lat>=-90)&(lat<=90)
    return shape,lon,lat,valid

def georeference_cells(original,land,continents,areas,elevation_tif):
    import geopandas as gpd
    import rasterio
    n=len(original)
    cell_id,lon,lat,valid=cell_coordinates(original)
    # Read official 10' elevation in WGS84; masked ocean/NODATA stays missing.
    with rasterio.open(elevation_tif) as raster:
        check("4326" in raster.crs.to_string() or
              raster.crs.to_epsg()==4326,
              "WorldClim raster must be WGS84 geographical")
        check(raster.count==1 and raster.width>1000,
              "WorldClim v2.1 elevation raster not proper global band")
        vals=raster.sample(list(zip(lon[valid],lat[valid])),masked=True)
        z=np.full(len(cell_id),np.nan,dtype=float)
        observed=np.array([v[0] for v in vals],dtype=float)
        nodata=raster.nodata
        if nodata is not None:
            observed[np.isclose(observed,nodata)]=np.nan
        observed[~np.isfinite(observed)]=np.nan
        observed[(observed < -500)|(observed > 9000)]=np.nan
        z[valid]=observed
        profile={"width":raster.width,"height":raster.height,
                 "CRS":raster.crs.to_string(),"nodata":nodata}
    selected=np.flatnonzero(valid)
    pts=gpd.GeoDataFrame({"sample_id":selected},
         geometry=gpd.points_from_xy(lon[selected],lat[selected]),
         crs="EPSG:4326")
    # The old MOIS photo-island assignment uses this same GSHHG high L1 layer
    # and shapely "within" (not coastal snapping).
    hits=gpd.sjoin(pts,land[["sibling_id","geometry"]],
                   how="left",predicate="within")
    uniq=hits.groupby("sample_id").sibling_id.nunique()
    ambiguous=set(uniq.index[uniq.gt(1)])
    hits=hits.loc[~hits.sample_id.isin(ambiguous)]
    hits=hits.drop_duplicates("sample_id",keep="first")
    good=hits.sibling_id.notna()
    sample_land=np.full(len(cell_id),"",dtype=object)
    sample_land[hits.loc[good,"sample_id"].to_numpy(int)]=(
        pd.to_numeric(hits.loc[good,"sibling_id"],errors="raise")
          .astype(np.int64).astype(str).to_numpy())
    sample_land[np.fromiter(ambiguous,dtype=int,count=len(ambiguous))]=""
    frame=pd.DataFrame({
       "cell_index":cell_id,"sample_landmass":sample_land,
       "elevation_m":z})
    frame["is_land"]=frame.sample_landmass.ne("")
    count=frame.groupby("cell_index").is_land.sum().reindex(
        np.arange(n),fill_value=0).to_numpy(int)
    # Frequency for exact same GSHHG island/continental sibling.
    on_land=frame[frame.is_land]
    count_mass=(on_land.groupby(["cell_index","sample_landmass"],sort=False)
                    .size().rename("support").reset_index())
    best=(count_mass.sort_values(
            ["cell_index","support","sample_landmass"],
            ascending=[True,False,True])
          .drop_duplicates("cell_index",keep="first"))
    actual=pd.DataFrame({
       "cell_index":np.arange(n),
       "n_land_subpoints":count,
    }).merge(best,on="cell_index",how="left",validate="one_to_one")
    actual["n_dominant_sibling"]=actual.support.fillna(0).astype(int)
    actual["landmass_sibling_id"]=actual.sample_landmass.fillna("")
    actual["gshhg_land_class"]="unassigned"
    assigned=actual.n_dominant_sibling.ge(7)
    actual.loc[assigned,"gshhg_land_class"]=np.where(
        actual.loc[assigned,"landmass_sibling_id"].isin(continents),
        "mainland","island"
    )
    actual["strict_all9_same_landmass"]=actual.n_dominant_sibling.eq(9)
    # Entire continents must not become islands simply because sjoin
    # returned float 0.0, 1.0 etc alongside NaN ocean observations.
    check(actual.loc[actual.n_dominant_sibling.ge(7),
                     "gshhg_land_class"].eq("mainland").any(),
          "GSHHG continental sibling IDs did not match spatial join")
    check(actual.loc[actual.n_dominant_sibling.ge(7),
                     "gshhg_land_class"].eq("island").any(),
          "GSHHG returned no real islands despite global land sample")
    actual["landmass_area_km2"]=pd.to_numeric(
        actual.landmass_sibling_id.map(
            {str(int(k)):v for k,v in areas.items()}),errors="coerce")
    # Height is the source raster median of land subpoints within cell;
    # it is not the elevation of all individual source records/flowers.
    med=(frame[frame.is_land].groupby("cell_index",sort=False)
          .elevation_m.median())
    actual["sampled_land_median_elevation_m"]=actual.cell_index.map(med)
    actual["n_elevation_subpoints"]=(frame[
         frame.is_land&frame.elevation_m.notna()]
         .groupby("cell_index").size()
         .reindex(actual.cell_index,fill_value=0).to_numpy())
    # Geographic coordinate of 100km equal-area cell center; do not use
    # nearest recorded occurrence's point as a representative location.
    midpoint=4
    actual["center_longitude"]=lon.reshape(n,9)[:,midpoint]
    actual["center_latitude"]=lat.reshape(n,9)[:,midpoint]
    actual=actual.drop(columns=["sample_landmass","support"])
    source_info={
       "n_GSHHG_land_polygons":int(len(land)),
       "n_island_majority7":int(actual.gshhg_land_class.eq("island").sum()),
       "n_mainland_majority7":int(actual.gshhg_land_class.eq("mainland").sum()),
       "n_unassigned_mixed_or_ocean":int(actual.gshhg_land_class.eq("unassigned").sum()),
       "n_strict_island9":int((
            actual.gshhg_land_class.eq("island")&
            actual.strict_all9_same_landmass).sum()),
       "n_cells_with_median_land_elevation":int(
            actual.sampled_land_median_elevation_m.notna().sum()),
       "raster_profile":profile,
       "n_ambiguous_sample_points_multiple_GSHHG_siblings":len(ambiguous),
    }
    return actual,source_info

def region_masks(d):
    n=d.n_accepted_species.ge(30)
    lat=d.center_latitude
    h=d.sampled_land_median_elevation_m
    island=d.gshhg_land_class.eq("island")
    main=d.gshhg_land_class.eq("mainland")
    return {
       "global_recorded_richness30":n,
       "arctic_circle_N66p562plus":n&lat.ge(ARC_LAT),
       "boreal_N50_to_66p562":n&lat.between(50,ARC_LAT,inclusive="left"),
       "northern_midlatitude_N30_to_55":n&lat.ge(30)&lat.lt(55),
       "southern_midlatitude_S30_to_55":n&lat.le(-30)&lat.gt(-55),
       "northern_hemisphere":n&lat.ge(0),
       "southern_hemisphere":n&lat.lt(0),
       "north_tropical_0_to_23p44":n&lat.ge(0)&lat.lt(23.44),
       "south_tropical_minus23p44_to_0":n&lat.lt(0)&lat.ge(-23.44),
       "islands_GSHHG_majority7":n&island,
       "mainlands_GSHHG_majority7":n&main,
       "strict_100km_islands_all9":n&island&d.strict_all9_same_landmass,
       "medium_to_large_islands_1000_to_100000km2":n&island&
             d.landmass_area_km2.ge(1000)&d.landmass_area_km2.lt(100000),
       "large_islands_area100000km2plus":n&island&d.landmass_area_km2.ge(100000),
       "low_elevation_below500m":n&h.ge(-430)&h.lt(500),
       "intermediate_elevation_500to1500m":n&h.ge(500)&h.lt(1500),
       "high_elevation_above1500m":n&h.ge(1500),
       "very_high_elevation_above2500m":n&h.ge(2500),
    }

def evaluate_group_null(d,x,nperm=199):
    check(nperm==199,"Signed parent null has exactly 199 permutations")
    mask=region_masks(d)
    kept={k:np.flatnonzero(v.to_numpy()) for k,v in mask.items() if int(v.sum())>=25}
    check(kept["arctic_circle_N66p562plus"].size>=20 and
          kept["southern_hemisphere"].size>2000,
          "High-latitude or southern source support changed")
    obs=d.observed_hue9_Shannon.to_numpy(float)
    expected=d.expected_hue9_Shannon_given_genus.to_numpy(float)
    ordered={name:{"indices":ids,
                   "observed":float(obs[ids].mean()),
                   "expected_signed":float(expected[ids].mean()),
                   "null_draws":np.empty(nperm)}
             for name,ids in kept.items()}
    rng=np.random.default_rng(SEED)
    for i in range(nperm):
        perm=x["hue"].copy()
        for g in x["variable_groups"]:
            perm[g]=x["hue"][g][rng.permutation(len(g))]
        h=empirical_entropy(counts_by_cell(x,perm))
        for a in ordered.values():
            a["null_draws"][i]=float(h[a["indices"]].mean())
        if (i+1)%50==0:
            print("MOIS_GEO_SOURCE_NINECOLOUR_GLOBAL_PERM",i+1,flush=True)
    results={}
    for name,a in ordered.items():
        null=a["null_draws"]
        observed=a["observed"]
        expected=a["expected_signed"]
        check(abs(null.mean()-expected)<1e-10,
              "True original genus×breadth permuted expectation not reproduced "+name)
        delta=observed-expected
        p_low=(1+np.count_nonzero(null<=observed+1e-12))/(nperm+1)
        p_two=(1+np.count_nonzero(
            abs(null-null.mean())>=abs(delta)-1e-12))/(nperm+1)
        dat=d.iloc[a["indices"]]
        results[name]={
           "n_100km_cells":len(dat),
           "n_distinct_500km_blocks":int(dat.block_500km.nunique())
                if "block_500km" in dat else None,
           "observed_mean_Shannon":observed,
           "genus_breadth_null_mean_Shannon":expected,
           "observed_minus_null_Shannon":delta,
           "empirical_one_sided_more_clustered_p":float(p_low),
           "empirical_two_sided_p":float(p_two),
           "null_95pct_global_label_permutation_interval":
               [float(v) for v in np.percentile(null,[2.5,97.5])],
           "mean_white_photo_species_share":float(
               (dat.n_white_photo_species/dat.n_accepted_species).mean()),
           "median_recorded_species_per_cell":float(
               dat.n_accepted_species.median()),
           "mean_effective_number_of_photo_colours":float(
               np.exp(dat.observed_hue9_Shannon).mean()),
           "null_is_global_species_label_exchangeability_within_genus_and_occupancy_width":True,
        }
    comp=[("island_minus_mainland",
           "islands_GSHHG_majority7","mainlands_GSHHG_majority7"),
          ("high_elevation_minus_low",
           "high_elevation_above1500m","low_elevation_below500m"),
          ("arctic_minus_boreal","arctic_circle_N66p562plus","boreal_N50_to_66p562"),
          ("south_minus_north_hemisphere","southern_hemisphere","northern_hemisphere"),
          ("southern_mid_minus_northern_mid_30to55",
           "southern_midlatitude_S30_to_55","northern_midlatitude_N30_to_55")]
    differences={}
    for name,A,B in comp:
        if A not in ordered or B not in ordered:continue
        a,b=ordered[A],ordered[B]
        p=(a["observed"]-b["observed"])
        q=a["null_draws"]-b["null_draws"]
        diff=p-float(q.mean())
        differences[name]={
          "observed_difference_Shannon":p,
          "null_expected_geographic_difference":float(q.mean()),
          "difference_in_genus_breadth_anomalies":diff,
          "permutation_two_sided_p":float((1+np.count_nonzero(
               abs(q-q.mean())>=abs(diff)-1e-12))/(nperm+1)),
          "95pct_null_difference_distribution":[float(v) for v
               in np.percentile(q,[2.5,97.5])],
          "two_groups":{"A":A,"B":B},
          "unmatched_composition_not_causal_location_or_elevation":True
        }
    return results,differences

def run(args):
    out=Path(args.out)
    check(not out.exists(),"No overwrite of verified source analytical results")
    pairs,species=signed_original(args.GBIF,args.ENV)
    x=baseline(pairs,species,null_model="genus_breadth")
    null=load_strict_nine_null(args.NINE)
    coords=np.asarray(x["cell_key"].tolist(),dtype=np.int64)
    check(np.array_equal(coords,null[["grid_x","grid_y"]].to_numpy()),
          "Source original 18577 GBIF cells mismatch stricter 9-hue null")
    check(sha(args.GSHHG)==GSHHG_SHA,
          "This must use the same GSHHG coastline source as old MOIS island photo analysis")
    with tempfile.TemporaryDirectory() as scratch:
        land,continents,seeds,areas,npoly=load_coast_and_continents(
            args.GSHHG,Path(scratch))
        geo,geo_info=georeference_cells(null,land,continents,areas,args.ELEV_TIF)
    derived=null.copy()
    geo=geo.drop(columns="cell_index")
    check(len(geo)==len(derived),"GSHHG grid assignment changed source n")
    derived=pd.concat([derived,geo],axis=1)
    check(len(derived)==18577 and
          not derived[["grid_x","grid_y"]].duplicated().any(),
          "Duplicate derived GSHHG photo hue cells")
    # The original parent null was based on the same GBIF species once per cell.
    check(int(derived.n_accepted_species.ge(30).sum())==12728,
          "Original global photo colour richness threshold source drift")
    # Spatial 500km key from original EPSG6933 equal-area grid coordinates,
    # not geodetic 5° or photo observation coordinate.
    derived["block_500km"]=(derived.grid_x.floordiv(5).astype(str)+"|"+
                             derived.grid_y.floordiv(5).astype(str))
    groups,comparison=evaluate_group_null(derived,x,N_PERM)
    out.mkdir(parents=True)
    csv=out/"GBIF18577_real_photo9_island_alpine_arctic_hemisphere_cells.csv.gz"
    saved=[
      "grid_x","grid_y","n_accepted_species","n_white_photo_species",
      "observed_hue9_Shannon","expected_hue9_Shannon_given_genus",
      "observed_minus_genus_null_Shannon","SES_genus_conditioned",
      "n_unique_observed_photo_hues","block_500km",
      "center_latitude","center_longitude","gshhg_land_class",
      "landmass_sibling_id","n_land_subpoints","n_dominant_sibling",
      "strict_all9_same_landmass","landmass_area_km2",
      "sampled_land_median_elevation_m","n_elevation_subpoints"]
    derived[saved].to_csv(csv,compression="gzip",index=False)
    rfile=out/"geographic_group_photographic_hue_199globalpermutation.json"
    result={
        "analysis":"global 9-colour photos × GSHHG island/continent × WorldClim elevation × latitude in identical GBIF 100km cells",
        "source_accepted_species":37104,
        "source_species_grid_pairs":5493931,
        "n_global_100km_cells":18577,
        "n_expected_full_photo9_source_cells_richness30":12728,
        "n_source_genus_breadth_global_shuffles":199,
        "geographic_input_sha256":{
            "GBIF":GBIF_SHA,"ENV":ENV_SHA,"strict_genus_breadth_9hue_null":NULL_SHA,
            "GSHHG_coastline":GSHHG_SHA,
            "WorldClim2p1_10m_elev_TIF":sha(args.ELEV_TIF),
            "WorldClim2p1_10m_elev_archive":sha(args.ELEV_ZIP)},
        "elevation_source":"WorldClim v2.1 10arcminute topographic raster, nine 100km EPSG6933 subpoint median on GSHHG land",
        "arctic_operational":"latitude>=66.562 N, not floristic Arctic",
        "GSHHG_island_operational":"dominant same landmass sibling>=7/9 100km equal-area subpoints, no coastal snapping",
        "GSHHG_continental_seed_ids":seeds,
        "GSHHG_geometry":{"n_polygons":npoly,**geo_info},
        "geographic_group_null":groups,
        "geographic_comparisons":comparison,
        "these_are_species_community_photo_label_gradients_not_verified_petals":True,
        "GSHHG_large_island_bias_for_100km_cells":True,
        "WCVP_native_only_not_yet_filtered":True,
        "observed_grids_are_spatially_dependent":True,
        "different_temporal_baselines_are_not_causal_evolution":True
    }
    rfile.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    outputsha={p.name:sha(p) for p in (csv,rfile)}
    (out/"receipt.json").write_text(json.dumps({
         **result,"output_sha256":outputsha
    },indent=2,allow_nan=False)+"\n")
    print("MOIS_NINECOLOUR_REAL_GEO_GRADIENTS_SUCCESS",json.dumps({
       "GSHHG":geo_info,
       "region_n":{k:v["n_100km_cells"] for k,v in groups.items()},
       "geography_group_summary":{k:{
           "white_photo_share":v["mean_white_photo_species_share"],
           "entropy_anomaly":v["observed_minus_null_Shannon"],
           "p_two":v["empirical_two_sided_p"]}
           for k,v in groups.items()},
       "geographic_comparison_summary":{k:{
           "anomaly_diff":v["difference_in_genus_breadth_anomalies"],
           "p_two":v["permutation_two_sided_p"]}
           for k,v in comparison.items()},
    },allow_nan=False),flush=True)

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    for k in ("GBIF","ENV","NINE","GSHHG","ELEV_ZIP","ELEV_TIF","out"):
        parser.add_argument("--"+k,required=True)
    run(parser.parse_args())
