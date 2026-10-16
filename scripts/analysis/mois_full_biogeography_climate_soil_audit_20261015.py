#!/usr/bin/env python3
"""Source audit: biogeographic sampling, climate, soil, altitude, heldout genus.

Runs ONLY on owner-held original FCP SHA256-verified ZIPs. Exports aggregate
numbers, not source photos, site coordinates, original individual taxon names
or per-genus results. Source elevation is NOT slope/aspect/TPI/TRI terrain.
This is retrospective observation, not causal evolutionary pigment fitness.
"""
from __future__ import annotations
import argparse, hashlib, json, zipfile
from pathlib import Path
import numpy as np
import pandas as pd

PHOTO_SHA="2d6097ba1b0434632e13bbe79dd78180608a287ff5dabebc770b5c0a2cfb52a0"
ENV_SHA="6e287294b3a05eaa383be3de644551825915bfecca9c4bde6950b4890f5e16ca"
DIR="fcp42111_real_site_abiotic/"
SITE="fcp42111_original_photo_all_climate_soil/all_original_100543_photo_sites_expanded_environment.csv.gz"
C=[*(f"wc_bio{i}" for i in range(1,20)),"wc_srad_annual_mean","wc_srad_monthly_cv","wc_wind_annual_mean",
   "wc_wind_monthly_cv","wc_vapr_annual_mean","wc_vapr_monthly_cv","wc_elevation_m"]
S=["soil_pH","soil_SOC","soil_N","soil_clay","soil_available_water_proxy",
   "soil_cec_0_30cm_source_raw","soil_sand_0_30cm_source_raw","soil_silt_0_30cm_source_raw",
   "soil_bdod_0_30cm_source_raw","soil_cfvo_0_30cm_source_raw"]
G=["abs_lat","lon_sin","lon_cos"];W=["wc_bio12","wc_vapr_annual_mean"]
TEMP=["wc_bio5","wc_srad_annual_mean"];ELEV=["wc_elevation_m"];SOIL5=S[:5]
MODELS={
 "geography":G,
 "geography_elevation":G+ELEV,
 "geography_moisture":G+W,
 "geography_moisture_elevation":G+W+ELEV,
 "geography_climate_elevation":G+W+TEMP+ELEV,
 "geography_climate_elevation_soil5":G+W+TEMP+ELEV+SOIL5,
 "geography_climate_elevation_soil10":G+W+TEMP+ELEV+S,
 "geography_elevation_soil10":G+ELEV+S
}
REPS=1999
SEED=202610151702

def sha(path):
 h=hashlib.sha256()
 with open(path,"rb") as f:
  for chunk in iter(lambda:f.read(1024*1024),b""):h.update(chunk)
 return h.hexdigest()

def read(path,member,cols):
 with zipfile.ZipFile(path) as z:
  with z.open(member) as stream:
   return pd.read_csv(stream,compression="gzip",usecols=cols,low_memory=False)

def original(photo_zip,env_zip):
 for path,required in ((photo_zip,PHOTO_SHA),(env_zip,ENV_SHA)):
  if sha(path)!=required:raise ValueError("Unexpected original FCP ZIP SHA256")
  with zipfile.ZipFile(path) as z:
   if z.testzip() is not None:raise ValueError("Invalid original archive member")
 cols=["inat_taxon_id","observation_id","photo_id","species","morph","measurement_status"]
 a=read(photo_zip,DIR+"breadth_flower_colour_climate_soil.csv.gz",cols)
 b=read(photo_zip,DIR+"taxon_cell_flower_colour_climate_soil.csv.gz",cols)
 p=pd.concat([a,b],ignore_index=True)
 if (len(a),len(b),int(p.photo_id.duplicated().sum()))!=(42111,85337,26905):raise ValueError("Original photo overlap drift")
 if p.groupby("photo_id").inat_taxon_id.nunique().max()!=1:raise ValueError("Identical photo used for differing taxa")
 p=p.drop_duplicates("photo_id")
 if (len(p),p.inat_taxon_id.nunique())!=(100543,42111):raise ValueError("Original source universe drift")
 env=read(env_zip,SITE,["inat_taxon_id","observation_id","photo_id","latitude","longitude",*C,*S])
 p=p.merge(env,on=["inat_taxon_id","observation_id","photo_id"],validate="one_to_one")
 ok=(p.measurement_status.eq("classified_four_state_morph")&
     p.morph.isin(["white","yellow_orange","red_pink","blue_purple"])&
     p.latitude.between(-90,90)&p.longitude.between(-180,180)&p[C].notna().all(axis=1))
 climate=p.loc[ok].copy()
 if len(climate)!=45836:raise AssertionError("Original climate26 label denominator drift")
 climate["genus"]=climate.species.str.split().str[0]
 climate["white"]=climate.morph.eq("white").astype(float)
 climate["abs_lat"]=climate.latitude.abs()
 climate["lon_sin"]=np.sin(np.deg2rad(climate.longitude))
 climate["lon_cos"]=np.cos(np.deg2rad(climate.longitude))
 soil=climate.loc[climate[S].notna().all(axis=1)].copy()
 if len(soil)!=31925:raise AssertionError("Original soil36 photo count drift")
 norm={k:float(soil[k].std(ddof=0)) for k in G+C+S}
 cols={k:"mean" for k in [*G,*C,*S,"latitude","longitude"]}
 cols.update({"genus":"first","white":"mean","photo_id":"size"})
 def species(rows):
  s=rows.groupby("inat_taxon_id",sort=True).agg(cols).rename(columns={"photo_id":"n_photos"}).reset_index()
  n=s.genus.value_counts()
  return s.loc[s.genus.isin(n.index[n>=3])].reset_index(drop=True)
 climsp=species(climate); soilsp=species(soil)
 if (len(climsp),climsp.genus.nunique(),int(climsp.n_photos.sum()))!=(20694,1934,38498):
  raise AssertionError("Climate congeneric species drift")
 if (len(soilsp),soilsp.genus.nunique(),int(soilsp.n_photos.sum()))!=(16801,1665,26163):
  raise AssertionError("Soil congeneric species drift")
 return climate,soil,climsp,soilsp,norm

def moments(s,norm,cols):
 genera,gi=np.unique(s.genus.to_numpy(),return_inverse=True)
 n=np.bincount(gi).astype(float)
 w=1/n[gi];X=np.stack([s[k].to_numpy(float)/norm[k] for k in cols],axis=1)
 y=s.white.to_numpy(float)
 mx=np.stack([np.bincount(gi,weights=X[:,j])/n for j in range(len(cols))],axis=1)
 my=np.bincount(gi,weights=y)/n
 x=X-mx[gi];y=y-my[gi]
 A=np.zeros((len(genera),len(cols),len(cols)));b=np.zeros((len(genera),len(cols)))
 np.add.at(A,gi,x[:,:,None]*x[:,None,:]*w[:,None,None])
 np.add.at(b,gi,x*y[:,None]*w[:,None])
 return A,b,float(np.sum(w*y*y))

def within_models(s,norm,genus_draws):
 records={}
 cases=[("geography_rain",G+W[:1]),
        ("climate_elevation",MODELS["geography_climate_elevation"]),
        ("with_soil5",MODELS["geography_climate_elevation_soil5"])]
 for label,cols in cases:
  Am,Bm,total=moments(s,norm,cols);A=Am.sum(axis=0);b=Bm.sum(axis=0)
  point=100*np.linalg.solve(A,b)
  if label=="geography_rain" and abs(point[-1]-2.256936929825337)>1e-8:
   raise AssertionError("Original soil36 genus-equal white/BIO12 slope not reproduced")
  mats=np.einsum("bg,gij->bij",genus_draws,Am,optimize=True)
  rhs=np.einsum("bg,gi->bi",genus_draws,Bm,optimize=True)
  boot=100*np.linalg.solve(mats,rhs[...,None])[:,:,0]
  keys=[key for key in [*W,*SOIL5] if key in cols]
  ix=np.array([cols.index(key) for key in keys])
  sd=np.std(boot[:,ix]-point[None,ix],axis=0,ddof=1)
  crit=float(np.quantile(np.max(abs((boot[:,ix]-point[None,ix])/sd[None,:]),axis=1),.95))
  result={"weighted_within_genus_in_sample_R2":float(1-(total-b@np.linalg.solve(A,b))/total),
          "normal_equation_condition_number":float(np.linalg.cond(A)),
          "selected_within_model_maxT_axes":keys,
          "within_model_maxT95_critical":crit,"coefficients":{}}
  for j,name in enumerate(cols):
   row={"point_pp_per_soil36_site_source_sd":float(point[j]),
        "ordinary_genus_bootstrap95":np.quantile(boot[:,j],[.025,.975]).tolist()}
   if name in keys:
    k=keys.index(name)
    row["selected_within_model_maxT95"]=[float(point[j]-crit*sd[k]),float(point[j]+crit*sd[k])]
   result["coefficients"][name]=row
  records[label]=result
 return records

def fold(g):
 return int.from_bytes(hashlib.sha256(str(g).encode()).digest()[:8],"big")%5

def predict_cv(s,nboot=REPS,seed=SEED):
 labels=["intercept",*MODELS]
 gids,gi=np.unique(s.genus.to_numpy(),return_inverse=True);n=np.bincount(gi)
 w=1/n[gi];y=s.white.to_numpy(float)
 F=sorted(set(x for v in MODELS.values() for x in v));X=s[F].to_numpy(float)
 split=s.genus.map(fold).to_numpy(int)
 errors=np.empty((len(s),len(labels)))
 for k in range(5):
  train=split!=k;test=~train;ww=w[train]/w[train].sum()
  intercept=float(np.sum(y[train]*ww))
  errors[test,0]=(y[test]-intercept)**2
  for j,name in enumerate(labels[1:],1):
   cols=MODELS[name];xx=X[:,[F.index(c) for c in cols]]
   mean=np.sum(xx[train]*ww[:,None],axis=0)
   sd=np.maximum(np.sqrt(np.sum(ww[:,None]*(xx[train]-mean)**2,axis=0)),1e-9)
   Z=(xx[train]-mean)/sd;T=(xx[test]-mean)/sd
   reg=np.linalg.solve((Z.T*ww)@Z+.01*np.eye(len(cols)),(Z.T*ww)@(y[train]-intercept))
   predicted=np.clip(intercept+T@reg,0,1)
   errors[test,j]=(y[test]-predicted)**2
 Eg=np.stack([np.bincount(gi,weights=errors[:,j])/n for j in range(len(labels))],axis=1)
 rng=np.random.default_rng(seed)
 draw=rng.integers(0,len(gids),size=(nboot,len(gids)),dtype=np.int32)
 boots=Eg[draw].mean(axis=1);points=Eg.mean(axis=0)
 geo=labels.index("geography");clim=labels.index("geography_climate_elevation")
 ret={"fold_species_counts":np.bincount(split,minlength=5).tolist(),"n_heldout_genera":len(gids),
      "n_bootstrap":nboot,"genus_fold_hash":"SHA256 first 8 bytes, modulo 5",
      "ridge_lambda":.01,"loss":"equal genus average squared error of classified photo-white species fraction, predictions clipped 0–1",
      "confidence_intervals":"1999 bootstrap resamples of observed genus-heldout losses WITHOUT nested refitting",
      "models":{}}
 for j,name in enumerate(labels):
  ret["models"][name]={
   "heldout_equal_genus_squared_error":float(points[j]),
   "gain_vs_geography":float(points[geo]-points[j]),
   "gain_vs_geography95":np.quantile(boots[:,geo]-boots[:,j],[.025,.975]).tolist(),
   "additional_gain_vs_climate_elevation":float(points[clim]-points[j]),
   "additional_gain_vs_climate_elevation95":np.quantile(boots[:,clim]-boots[:,j],[.025,.975]).tolist(),
  }
 return ret

def sector_transfer(s):
 sector=np.clip(np.floor((s.longitude.to_numpy(float)+180)/90),0,3).astype(int)
 genus_ct=s.genus.value_counts();w=s.genus.map(1/genus_ct).to_numpy(float)
 y=s.white.to_numpy(float)
 names=["geography","geography_moisture","geography_climate_elevation",
        "geography_climate_elevation_soil5","geography_climate_elevation_soil10"]
 out={}
 for j,sectorname in enumerate(("W180_W90","W90_E0","E0_E90","E90_E180")):
  test=sector==j;sampled=set(s.loc[test,"genus"])
  train=(~test)&(~s.genus.isin(sampled).to_numpy())
  ww=w[train]/w[train].sum();baseline=np.sum(y[train]*ww)
  errors={}
  for name in names:
   cols=MODELS[name];X=s[cols].to_numpy(float)
   mu=np.sum(X[train]*ww[:,None],axis=0)
   sd=np.maximum(np.sqrt(np.sum(ww[:,None]*(X[train]-mu)**2,axis=0)),1e-9)
   T=(X[train]-mu)/sd;Q=(X[test]-mu)/sd
   b=np.linalg.solve((T.T*ww)@T+.01*np.eye(len(cols)),(T.T*ww)@(y[train]-baseline))
   prediction=np.clip(baseline+Q@b,0,1)
   errors[name]=float(np.average((y[test]-prediction)**2,weights=w[test]))
  out[sectorname]={
   "n_test_species":int(test.sum()),"n_test_genera":len(sampled),
   "n_train_species_excluding_all_test_genera":int(train.sum()),
   "n_train_genera_excluding_all_test_genera":int(s.loc[train,"genus"].nunique()),
   "water_gain_vs_geography":errors["geography"]-errors["geography_moisture"],
   "soil5_increment_vs_climate_elevation":errors["geography_climate_elevation"]-errors["geography_climate_elevation_soil5"],
   "model_errors":errors,
  }
 return out

def source_support(climate,soil):
 sector=np.clip(np.floor((climate.longitude.to_numpy(float)+180)/90),0,3).astype(int)
 eligible=climate[S].notna().all(axis=1).to_numpy()
 r={}
 for j,label in enumerate(("W180_W90","W90_E0","E0_E90","E90_E180")):
  sel=sector==j
  r[label]={"n_classified_climate26":int(sel.sum()),
            "n_soil36_complete":int(np.sum(sel&eligible)),
            "soil36_complete_fraction":float(np.mean(eligible[sel])),
            "photographed_white_fraction":float(np.mean(climate.white.to_numpy(float)[sel]))}
 elevation=climate.wc_elevation_m.to_numpy(float)
 edges=np.quantile(elevation,[0,.25,.5,.75,1])
 qr=[]
 for i in range(4):
  mask=(elevation>=edges[i])&((elevation<=edges[i+1]) if i==3 else (elevation<edges[i+1]))
  qr.append({"quartile":i+1,"elevation_bounds_m":edges[i:i+2].tolist(),
             "classified_climate26_photos":int(mask.sum()),
             "complete_soil36_photos":int(np.sum(mask&eligible)),
             "soil36_completion":float(np.mean(eligible[mask]))})
 return {"soil36_completeness":len(soil)/len(climate),"longitude_sectors":r,"elevation_quartiles":qr}

def continent_summary(species,shapefile):
 if shapefile is None:return {"status":"NOT_TESTED_NO_110M_CONTINENT_POLYGON_SOURCE"}
 import geopandas as gpd
 world=gpd.read_file(shapefile)
 if "continent" not in world.columns or world.crs is None:raise ValueError("Natural Earth 110m polygon CRS and continent attributes needed")
 world=world.to_crs("EPSG:4326")
 world=world.loc[~world.continent.isin(["Antarctica","Seven seas (open ocean)"]),["continent","geometry"]]
 points=gpd.GeoDataFrame(species[["inat_taxon_id","latitude","longitude"]].copy(),
                       geometry=gpd.points_from_xy(species.longitude,species.latitude),crs="EPSG:4326")
 joined=gpd.sjoin(points,world,how="left",predicate="intersects")
 joined=joined.drop_duplicates("inat_taxon_id").set_index("inat_taxon_id")
 data=species.copy()
 data["continent"]=data.inat_taxon_id.map(joined.continent).fillna("unresolved_110m")
 groups={str(key):{"n_species":len(g),"mean_photo_white_species_equal":float(g.white.mean())}
         for key,g in data.groupby("continent")}
 return {"source_polygon":"Natural Earth 110m administered land country continent polygon, not a true biogeographic realm",
         "source_shapefile_sha256":sha(shapefile),"n_source_species":len(species),
         "unresolved_species":groups.get("unresolved_110m",{}).get("n_species",0),
         "interpretation":"Unresolved includes small islands/coasts and mean photo-site point mismatches. Global genus bias and photo selection uncorrected.",
         "groups":groups}

def full_report(photo,env,shapefile):
 climate,soil,climatespp,soilspp,norm=original(photo,env)
 rng=np.random.default_rng(202610151701)
 draws=rng.multinomial(1665,np.repeat(1/1665,1665),size=REPS)
 wls=within_models(soilspp,norm,draws)
 # Soil texture is a closed composition: 10x percent clay + raw sand + raw silt.
 texture=soilspp.soil_sand_0_30cm_source_raw+soilspp.soil_silt_0_30cm_source_raw+10*soilspp.soil_clay
 return {
  "schema":"mois_original_biogeography_climate_soil_elevation_full_audit_v1",
  "date_jst":"2026-10-15","source_photo_sha256":PHOTO_SHA,"source_environment_sha256":ENV_SHA,
  "sample":{"photo_opportunities":100543,"source_hue_classified":45953,"climate26_photos":len(climate),
            "soil36_photos":len(soil),"climate26_species":len(climatespp),"climate26_genera":climatespp.genus.nunique(),
            "soil36_species":len(soilspp),"soil36_genera":soilspp.genus.nunique(),
            "soil36_qualified_photos":int(soilspp.n_photos.sum())},
  "replayed_original_soil36_geography_adjusted_white_BIO12_pp_per_source_soil36_photo_sd":
      wls["geography_rain"]["coefficients"]["wc_bio12"]["point_pp_per_soil36_site_source_sd"],
  "soil36_genus_centered_wls":wls,
  "genus_heldout_fivefold":predict_cv(soilspp),
  "longitude_sector_heldout_excluding_target_genera":sector_transfer(soilspp),
  "soil_coverage_and_geographic_selection":source_support(climate,soil),
  "continent_110m_descriptive":continent_summary(climatespp,shapefile),
  "soil_texture_composition_diagnostic":{
   "closure_expression":"sand_raw+silt_raw+10*clay_percent approximately 1000",
   "closure_mean":float(texture.mean()),"closure_std":float(texture.std(ddof=0)),
   "closure_max_abs_distance_from_1000":float(np.max(abs(texture-1000)))},
  "true_terrain_metrics":{
   "source_terrain_point_elevation_m":True,"independent_projected_DEM_supplied":False,
   "slope":None,"aspect":None,"terrain_ruggedness_TRI":None,"topographic_position_TPI":None,
   "verdict":"NOT TESTED: point elevation cannot identify relief/slope/aspect or terrain geometry."},
  "limits":[
   "Source images, genus counts and photo colour classes were not randomly sampled; photo-label missingness and soil missingness are geographic and altitudinal.",
   "Within-source conditional water/soil associations are NOT pigment, pollinator or natural-selection causal effects.",
   "SOIL5 selected-coefficient maxT is conditional on this retrospective model; historical source 24-contrast inference is separate.",
   "Out-of-genus CV uses training-only centring, fixed ridge alpha=0.01 and genus-bootstrapped fixed heldout losses, not nested repeated model refitting.",
   "Five-variable soil increment and 10-soil increment refer to identical soil36 species; negative CV gain is NOT a proof that soil never affects flower colour.",
   "Longitude sectors are not independently verified biogeographical realms; 110m country polygons miss small islands and species mean photo-sites are not range maps.",
   "Original soil water availability is a historical numerical proxy, not an in situ measured field capacity.",
   "No original individual photo/taxon/site identities, geographic coordinates or attributed images exported."],
 }
def main():
 p=argparse.ArgumentParser()
 p.add_argument("--photo-zip",type=Path,required=True)
 p.add_argument("--env-zip",type=Path,required=True)
 p.add_argument("--continent-shp",type=Path)
 p.add_argument("--output",type=Path,required=True)
 mode=p.add_mutually_exclusive_group(required=True)
 mode.add_argument("--write",action="store_true")
 mode.add_argument("--check",action="store_true")
 args=p.parse_args()
 report=full_report(args.photo_zip,args.env_zip,args.continent_shp)
 if args.write:
  args.output.write_text(json.dumps(report,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
  print("WROTE original source-climate soil spatial audit",args.output)
 else:
  # Existing frozen safe summary is a compact audit receipt; source-level claims
  # are checked against its exact fixed scientific endpoints, not file equality.
  old=json.loads(args.output.read_text(encoding="utf8"))
  assert old["source_scope"]["soil36_species"]==report["sample"]["soil36_species"]==16801
  assert old["source_scope"]["soil36_genera"]==report["sample"]["soil36_genera"]==1665
  assert abs(old["source_scope"]["original_soil36_BIO12_genus_equal_coef_pp_per_photo_SD"]-
             report["replayed_original_soil36_geography_adjusted_white_BIO12_pp_per_source_soil36_photo_sd"])<1e-8
  x=report["genus_heldout_fivefold"]["models"]
  y=old["heldout_genus_cv"]["models"]
  for mine,historical in [
   ("geography_moisture","geography_moisture"),
   ("geography_climate_elevation","geography_climate_elevation"),
   ("geography_climate_elevation_soil5","geography_climate_elevation_soil5"),
   ("geography_climate_elevation_soil10","geography_climate_elevation_soil10")]:
   assert abs(x[mine]["heldout_equal_genus_squared_error"]-y[historical]["loss"])<1e-10
  assert abs(x["geography_climate_elevation_soil5"]["additional_gain_vs_climate_elevation"]-
             y["geography_climate_elevation_soil5"]["increment_vs_climate_elevation"])<1e-10
  aw=report["soil36_genus_centered_wls"]["with_soil5"]["coefficients"]["soil_available_water_proxy"]["point_pp_per_soil36_site_source_sd"]
  assert abs(aw-old["global_same_soil36_taxa_OLS_congeneric"]["adding_soil5"]["soil_available_water_proxy"]["point_pp_per_source_soil36_photo_site_SD"])<.001
  assert report["true_terrain_metrics"]["slope"] is None
  print("PASS source-only climate/soil/biogeographic matched receipt and real-terrain nonclaim")

if __name__=="__main__":main()
