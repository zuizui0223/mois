#!/usr/bin/env python3
"""NEWEST 1,043,930 genuine original photo island + elevation analysis.

Reuses SHA-locked GSHHG 2.3.7 coastline classification for its exact 313,535
old original photo IDs; assigns island/continental sibling to every additional
latest genuine photo with exact same high-res land polygons, never country
proxies or coastal snapping. Re-estimates raw, spatial/climate conditioned,
optimal caliper-1 matching, and same-landmass nearby high/low elevations.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from scipy.optimize import linear_sum_assignment
from scipy.stats import t as t_dist
import statsmodels.api as sm
import statsmodels.formula.api as smf

from scripts.analysis.mois_global_geography_composition import (
    load_land,classify_land,community_table,fit_contrast,
)
PHOTO_SHA="9a4f3dae6fbcee967b71cb29576e71d52c48e1f1f790b7d3ea4750c4bcb0e6d4"
CLIM_SHA="ec39b49649acc2ca4a468c00c86d05800e33ea2db1f91a3b119076be38432735"
COAST_SHA="8dbbe7e071e77e9e75f2d639239099ebca8d5c16d6a07df8169729d49f15cf41"
HUES=("white","yellow_orange","red_pink","blue_purple")
SCHEMA="mois_latest1043930_source_photo_GSHHG_island_and_nearby_altitude_v1"

def require(ok,reason):
 if not ok:raise ValueError("MOIS_LATEST1043930_ISLAND: "+reason)
def sha(path):
 h=hashlib.sha256()
 with Path(path).open("rb") as f:
  for piece in iter(lambda:f.read(4<<20),b""):h.update(piece)
 return h.hexdigest()
def read_cached_land(cache):
 cache=Path(cache)
 original=json.loads((cache/"result.json").read_text())
 require(original["source_photo_ids"]==313535
         and original["gshhg_source_sha256"]==COAST_SHA,
         "Existing exact GSHHG source split changed")
 for name,digest in original["output_sha256"].items():
  require(sha(cache/name)==digest,"Original source coastline/photo result SHA drift "+name)
 d=pd.read_csv(cache/"photo_landmass_cache.csv.gz",
     dtype={"photo_id":str,"landmass":str},keep_default_na=False)
 require(len(d)==313535 and d.photo_id.is_unique,
         "Old 313535 classified-original-photo coastal cache absent")
 d.area_km2=pd.to_numeric(d.area_km2,errors="coerce")
 return d[["photo_id","land_type","landmass","area_km2"]],original

def all_geography(p,geo_cache,coast_zip,out,block_size=100000):
 prior,_=read_cached_land(geo_cache)
 known=p.photo_id.isin(prior.photo_id)
 require(int(known.sum())==313535,
    "Frozen GSHHG 313535 photo IDs must all exist within latest genuine 1043930 photo cohort")
 unused=p.loc[~known,["photo_id","latitude","longitude"]].copy()
 valid=unused.latitude.between(-90,90)&unused.longitude.between(-180,180)
 unique=unused.loc[valid,["latitude","longitude"]].drop_duplicates().reset_index(drop=True)
 require(sha(coast_zip)==COAST_SHA,"Not SHA-verified island GSHHG high-res shoreline ZIP")
 root=out/"shapefile";root.mkdir(parents=True,exist_ok=True)
 land,continents=load_land(coast_zip,root)
 pieces=[]
 for start in range(0,len(unique),block_size):
  frame=unique.iloc[start:start+block_size].copy()
  spatial=classify_land(frame,land,continents).reset_index(drop=True)
  piece=pd.concat([frame.reset_index(drop=True),spatial],axis=1)
  pieces.append(piece)
  if start%(block_size*3)==0:
   print("MOIS_LATEST1043930_GSHHG_NEW_EXACT_SITES",start,"of",len(unique),flush=True)
 coords=pd.concat(pieces,ignore_index=True) if pieces else pd.DataFrame(
  columns=["latitude","longitude","land_type","landmass","area_km2"])
 require(not coords[["latitude","longitude"]].duplicated().any(),
         "Distinct original photo coordinates repeated in landmass index")
 added=unused.merge(coords,on=["latitude","longitude"],how="left",
                    sort=False,validate="many_to_one")
 added["land_type"]=added.land_type.fillna("unassigned")
 added["landmass"]=added.landmass.fillna("")
 geo=pd.concat([prior,added[["photo_id","land_type","landmass","area_km2"]]],
               ignore_index=True)
 require(len(geo)==1043930 and geo.photo_id.is_unique,
         "Island category source photo-ID coverage lost or duplicated")
 geo=geo.set_index("photo_id").loc[p.photo_id.tolist()].reset_index()
 require(geo.photo_id.tolist()==p.photo_id.tolist(),
         "Lost original ordered source-image photo ID mapping")
 require(geo.land_type.isin(("island","mainland","unassigned")).all(),
         "No false unknown island label is acceptable")
 return geo,{"previous_coastline_assigned_photo_IDs_reused":313535,
    "new_coastline_photo_IDs_assigned_or_preserved_unknown":len(unused),
    "new_exact_lonlat_unique_sites_tested":len(unique),
    "original_seed_continental_sibling_ids":sorted(str(v) for v in continents),
    "all_photo_ids":len(geo)}

def coarse_matched(cells):
 c=cells.loc[cells.n_taxa.ge(5)&cells.n_photos.ge(10)].dropna(
     subset=["heat","rain","solar","elevation"]).copy()
 both=c.groupby("region").island.nunique()
 c=c[c.region.isin(both[both.eq(2)].index)]
 require(len(c)>60 and c.island.sum()>=10,
         "Missing island mainland shared coarse source geographic support")
 for col in ("heat","solar","elevation"):
  sd=c[col].std()
  require(sd>0,"Constant original predictor "+col)
  c["z_"+col]=(c[col]-c[col].mean())/sd
 z=np.log1p(c.rain.clip(lower=0));c["z_rain"]=(z-z.mean())/z.std()
 model=smf.ols(
    "p_white ~ island + z_heat + z_rain + z_solar + z_elevation + C(region)",c
 ).fit(cov_type="cluster",cov_kwds={"groups":c.block},use_t=True)
 ci=model.conf_int().loc["island"]
 predictors=["z_heat","z_rain","z_solar","z_elevation"]
 pairs=[]
 for region,g in c.groupby("region",sort=False):
  island=g[g.island.eq(1)]
  mainland=g[g.island.eq(0)]
  if island.empty or mainland.empty:continue
  distance=cdist(island[predictors].to_numpy(float),
                 mainland[predictors].to_numpy(float))
  cost=np.concatenate([np.where(distance<=1,distance,10000.0),
          np.full((len(island),len(island)),100.0)],axis=1)
  ii,jj=linear_sum_assignment(cost)
  for i,j in zip(ii,jj):
   if j<len(mainland) and distance[i,j]<=1:
    a=island.iloc[i];b=mainland.iloc[j]
    pairs.append({
      "region":str(region),"island_landmass":str(a.landmass),
      "island_white":float(a.p_white),
      "mainland_white":float(b.p_white),
      "difference":float(a.p_white-b.p_white),
      "matched_environment_distance":float(distance[i,j]),
    })
 paired=pd.DataFrame(pairs)
 require(len(paired)>=20 and paired.region.nunique()>=5,
         "No sufficient within-region ecologically matched island mainland photo cells")
 fit=sm.OLS(paired.difference.to_numpy(float),
            np.ones((len(paired),1))).fit(
              cov_type="cluster",cov_kwds={"groups":paired.region},use_t=True)
 return {
    "cells_in_geography_regression":int(len(c)),
    "island_cells":int(c.island.sum()),"mainland_cells":int((1-c.island).sum()),
    "coarse_regions":int(c.region.nunique()),"spatial_blocks":int(c.block.nunique()),
    "classified_photo_support_in_cells":int(c.n_photos.sum()),
    "spatial_climate_adjusted_island_minus_mainland_white_pp":
        100*float(model.params["island"]),
    "adjusted_nominal_spatial_cluster_95_pp":
        [100*float(v) for v in ci],
    "adjusted_nominal_p":float(model.pvalues["island"]),
    "optimal_match_pairs":int(len(paired)),
    "matched_island_landmasses":int(paired.island_landmass.nunique()),
    "matched_regions":int(paired.region.nunique()),
    "matched_island_minus_mainland_white_pp":100*float(fit.params[0]),
    "matched_nominal_region_cluster_95_pp":
        [100*float(v) for v in fit.conf_int()[0]],
    "matched_nominal_p":float(fit.pvalues[0]),
    "matched_island_white_percent":100*float(paired.island_white.mean()),
    "matched_mainland_white_percent":100*float(paired.mainland_white.mean()),
    "not_same_taxa_mandated_on_both_landmasses":True,
 },paired

def neighbourhood_elevation(cells,width=5):
 c=cells[cells.n_taxa.ge(5)&cells.n_photos.ge(10)].copy()
 c=c[np.isfinite(c.elevation)&np.isfinite(c.latitude)&np.isfinite(c.longitude)]
 c["latblock"]=np.floor((c.latitude+90)/width).astype(int)
 c["lonblock"]=np.floor((c.longitude+180)/width).astype(int)
 keys=["landmass","latblock","lonblock"]
 def summary(s):
  return s.groupby(keys,observed=True).agg(
    white=("p_white","mean"),n_cells=("p_white","size"),
    n_photo=("n_photos","sum")).reset_index()
 lo=summary(c[c.elevation.ge(0)&c.elevation.lt(500)])
 hi=summary(c[c.elevation.ge(1500)])
 pair=lo.merge(hi,on=keys,suffixes=("_low","_high"),validate="one_to_one")
 require(len(pair)>=15,"Nearby same landmass high-low elevation overlap missing")
 dif=(pair.white_high-pair.white_low).to_numpy(float)*100
 mean=float(np.mean(dif))
 se=float(np.std(dif,ddof=1)/np.sqrt(len(dif)))
 ci=t_dist.interval(.95,df=len(dif)-1,loc=mean,scale=se)
 return {
    "neighbourhood_width_deg":width,"paired_same_landmass_neighbourhoods":len(pair),
    "mean_low_0_500m_white_percent":100*float(pair.white_low.mean()),
    "mean_high_ge1500m_white_percent":100*float(pair.white_high.mean()),
    "high_minus_low_pp":mean,
    "nominal_neighbourhood_t95_pp":[float(x) for x in ci],
    "neighbourhoods_high_less_white":int((dif<0).sum()),
    "source_classifiable_photos_low":int(pair.n_photo_low.sum()),
    "source_classifiable_photos_high":int(pair.n_photo_high.sum()),
    "source_cells_low":int(pair.n_cells_low.sum()),"source_cells_high":int(pair.n_cells_high.sum()),
 },pair

def run(photo,climate,cache,coast,output):
 photo,climate,cache,coast,output=map(Path,(photo,climate,cache,coast,output))
 require(not output.exists(),"Cannot overwrite latest genuine island photographic result")
 require(sha(photo)==PHOTO_SHA and sha(climate)==CLIM_SHA,
         "Latest real photo ID source/WorldClim SHA changed")
 p=pd.read_parquet(photo);e=pd.read_parquet(climate)
 require(len(p)==len(e)==1043930 and p.photo_id.tolist()==e.photo_id.tolist(),
         "Newest image photo ID and environmental rows no longer align")
 output.mkdir(parents=True)
 geo,provenance=all_geography(p,cache,coast,output)
 geo_path=output/"latest1043930_all_original_photo_GSHHG_landmass.parquet"
 geo.to_parquet(geo_path,index=False,compression="zstd")
 d=p[["photo_id","inat_taxon_id","species","morph","latitude","longitude","arm"]].copy()
 d=d.merge(geo,on="photo_id",how="left",validate="one_to_one",sort=False)
 require(d.photo_id.tolist()==p.photo_id.tolist(),
         "Photo identity shifted by GSHHG coast join")
 for source,dest in [("wc2_bio05_10m","heat"),("wc2_bio12_10m","rain"),
       ("annual_solar_Wm2","solar"),("wc2_elev_2p5m","elevation")]:
  d[dest]=pd.to_numeric(e[source],errors="coerce").to_numpy(float)
 known=d.morph.isin(HUES)
 count={}
 for land in ("island","mainland","unassigned"):
  z=d[d.land_type.eq(land)]
  k=z[z.morph.isin(HUES)]
  count[land]={
    "source_original_photo_attempts":int(len(z)),
    "four_hue_classified":int(len(k)),
    "white_original_photo_ROI":int(k.morph.eq("white").sum()),
    "white_fraction_among_classified":
       float(k.morph.eq("white").mean()) if len(k) else None,
    "classifiability_photo_fraction":float(len(k)/len(z)) if len(z) else None,
    "hue_counts":{h:int(k.morph.eq(h).sum()) for h in HUES},
  }
 require(sum(x["source_original_photo_attempts"] for x in count.values())==1043930,
         "Source photo ID geography classification coverage does not sum to total")
 cells=community_table(d)
 adjusted,matched=coarse_matched(cells)
 e5,pair5=neighbourhood_elevation(cells,5)
 e10,pair10=neighbourhood_elevation(cells,10)
 result={
    "schema":SCHEMA,"real_original_source_photo_attempts":1043930,
    "latest_four_hue_classified":459235,
    "source_roi_success":int(known.sum()),
    "source_shas":{"colour":PHOTO_SHA,"climate":CLIM_SHA,
                   "old_GSHHG":sha(cache/"result.json"),"highres_GSHHG_zip":COAST_SHA},
    "island_mainland_photo_source_raw":count,
    "island_minus_mainland_raw_white_pp":
      100*(count["island"]["white_fraction_among_classified"]-
           count["mainland"]["white_fraction_among_classified"]),
    "island_mainland_environment_matched_and_spatial_regression":adjusted,
    "elevation_same_landmass_nearby_5deg":e5,
    "elevation_nearby_10deg_sensitivity":e10,
    "geographic_photo_ID_assignments":provenance,
    "output_geography_photo_ID_parquet_SHA256":sha(geo_path),
    "all_original_photo_ids_retained_once":True,
    "technical_image_4hue_unknowns_not_imputed":True,
    "wild_native_flora_and_pollinator_scarcity_identified":False,
    "landmass_area_not_island_nativity_verified":True,
    "all_image_downloads_zero":True,
 }
 for name,frame in [("latest1043930_island_environment_matched_pairs.csv",matched),
                    ("latest1043930_elevation5deg_same_landmass_pairs.csv",pair5),
                    ("latest1043930_elevation10deg_same_landmass_pairs.csv",pair10)]:
  frame.to_csv(output/name,index=False)
 (output/"mois_latest1043930_island_altitude_result.json").write_text(
    json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+"\n")
 print("MOIS_LATEST1043930_TRUE_ISLAND_AND_SAME_LANDMASS_ALTITUDE",json.dumps({
   "source":1043930,"mapped_original_photo_counts":{
      k:v["source_original_photo_attempts"] for k,v in count.items()},
   "raw_island_mainland_pp":result["island_minus_mainland_raw_white_pp"],
   "spatial_climate_adjusted_island_pp":
      adjusted["spatial_climate_adjusted_island_minus_mainland_white_pp"],
   "island_matched_pairs":adjusted["optimal_match_pairs"],
   "island_matched_pp":adjusted["matched_island_minus_mainland_white_pp"],
   "matched95":adjusted["matched_nominal_region_cluster_95_pp"],
   "high_low_elevation_5deg":e5
 },sort_keys=True),flush=True)
 return result

def main():
 a=argparse.ArgumentParser(description=__doc__)
 for k in ("photo","climate","cache","coast","output"):
  a.add_argument("--"+k,required=True,type=Path)
 x=a.parse_args()
 run(x.photo,x.climate,x.cache,x.coast,x.output)
if __name__=="__main__":main()
