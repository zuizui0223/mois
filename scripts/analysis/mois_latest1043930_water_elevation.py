#!/usr/bin/env python3
"""Genuine latest 1,043,930 original-photo MOIS moisture and elevation ecology.

458,174 automatically four-hue classified photos with complete re-sampled
WorldClim annual rainfall/BIO5/BIO14/VAPR/SRAD/DEM. Classifier nonresponse
remains missing. Photo-cohort, nominal genus, coarse region conditioned
taxon-geocell slopes (not within-individual causal pigment adaptation).
"""
from __future__ import annotations
import argparse,json,hashlib,re
from pathlib import Path
import numpy as np
import pandas as pd

PHOTO_SHA="9a4f3dae6fbcee967b71cb29576e71d52c48e1f1f790b7d3ea4750c4bcb0e6d4"
ENV_SHA="ec39b49649acc2ca4a468c00c86d05800e33ea2db1f91a3b119076be38432735"
HUES=("white","yellow_orange","red_pink","blue_purple")
PRED=("log1p_annual_rainfall","mean_atmospheric_vapour","annual_solar",
      "warmest_month_temperature","elevation")
SRC={
 "log1p_annual_rainfall":"wc2_bio12_10m",
 "mean_atmospheric_vapour":"wc2_annual_vapr_10m",
 "annual_solar":"annual_solar_Wm2",
 "warmest_month_temperature":"wc2_bio05_10m",
 "elevation":"wc2_elev_2p5m"}
SCHEMA="mois_original1043930_four_hue_latest_climate_source_conditional_water_altitude_v1"
def need(ok,why):
 if not ok:raise ValueError("MOIS_LATEST1043930_WATER_ELEVATION: "+why)
def sha(path):
 h=hashlib.sha256()
 with Path(path).open("rb") as s:
  for b in iter(lambda:s.read(4<<20),b""):h.update(b)
 return h.hexdigest()
def prepare(photo,clim):
 p=pd.read_parquet(photo);e=pd.read_parquet(clim)
 need(sha(photo)==PHOTO_SHA and sha(clim)==ENV_SHA,"Original image+environment SHAs changed")
 need(len(p)==len(e)==1043930 and p.photo_id.astype(str).tolist()==e.photo_id.astype(str).tolist(),
      "Not photo-ID-aligned latest original ROI-v4 + WorldClim")
 n=p.known_image_roi.eq(1)
 need(int(n.sum())==459235,"Latest genuine photo classifier terminal changed")
 x=p.loc[n,["photo_id","inat_taxon_id","species","latitude","longitude","arm","morph","white"]].copy()
 for dest,col in SRC.items():x[dest]=pd.to_numeric(e.loc[n,col],errors="coerce").to_numpy(float)
 for k in PRED:x=x[np.isfinite(x[k])].copy()
 x=x[x.latitude.between(-90,90)&x.longitude.between(-180,180)]
 x["inat_taxon_id"]=x.inat_taxon_id.fillna("").astype(str)
 x=x[x.inat_taxon_id.str.fullmatch(r"\d+")]
 x["genus"]=x.species.fillna("").astype(str).str.strip().str.split().str[0]
 x=x[x.genus.str.fullmatch("[A-Z][A-Za-z-]+",na=False)].copy()
 x["cohort"]=np.where(x.arm.str.startswith("FCP5_"),"FCP5",x.arm).astype(str)
 x["log1p_annual_rainfall"]=np.log1p(x["log1p_annual_rainfall"].clip(lower=0))
 x["white_float"]=x.white.astype(float)
 need(len(x)>=300000,"Latest classified-photo environmental and taxon support unexpectedly insufficient")
 return x

def aggregate(x,degree):
 d=x.copy()
 d["latblock"]=np.floor((d.latitude+90)/degree).astype(np.int16)
 d["lonblock"]=np.floor((d.longitude+180)/degree).astype(np.int16)
 grouping=["cohort","genus","latblock","lonblock","inat_taxon_id"]
 t=d.groupby(grouping,observed=True).agg(
    p_white=("white_float","mean"),n_photos=("photo_id","size"),
    **{z:(z,"mean") for z in PRED}).reset_index()
 need(t.n_photos.sum()==len(x),"Not all exact source photos accounted by named taxon/cell")
 return t

def simple_slopes(t,degree,min_taxa=2,draws=399,seed=20261025):
 group=["cohort","genus","latblock","lonblock"]
 size=t.groupby(group,observed=True).inat_taxon_id.transform("size")
 v=t[size>=min_taxa].copy()
 need(len(v)>=5000,"Too few genus/cohort support cells")
 scales={col:float(t[col].std(ddof=0)) for col in PRED}
 need(all(np.isfinite(z) and z>0 for z in scales.values()),"Degenerate WorldClim predictor")
 for col in ["p_white",*PRED]:
  v["c_"+col]=v[col]-v.groupby(group,observed=True)[col].transform("mean")
 dy=v.c_p_white.to_numpy(float)
 results=[]
 for i,col in enumerate(PRED):
  dx=v["c_"+col].to_numpy(float)/scales[col]
  xy=dx*dy;xx=dx*dx
  b=100*float(xy.sum()/xx.sum())
  g=pd.DataFrame({"genus":v.genus.to_numpy(),"xy":xy,"xx":xx})
  g=g.groupby("genus",observed=True)[["xy","xx"]].sum()
  rng=np.random.default_rng(seed+degree*91+i*53)
  z=g.to_numpy(float);idx=rng.integers(0,len(z),(draws,len(z)))
  agg=z[idx].sum(axis=1);bootstrap=100*agg[:,0]/agg[:,1]
  ci=np.quantile(bootstrap,[.025,.975]).tolist()
  results.append({
   "grid_degrees":degree,"min_named_taxa_per_genus_cell_arm":min_taxa,
   "predictor":col,"white_pp_per_original_taxon_cell_SD":b,
   "nominal_95_genus_cluster_bootstrap_pp":ci,
   "original_taxon_geocells":len(v),
   "source_cohort_genus_region_groups":v.groupby(group,observed=True).ngroups,
   "genera_for_bootstrap":len(g),"predictor_original_SD":scales[col],
   "selected_photo_ROIs":int(v.n_photos.sum()),
   "cohort_arm_within_stratum":True,"observational_univariate_not_causal":True,
  })
 # Environmental attributes are not independent: compute a simultaneous
 # within-source/cohort/geogenus partial model as a diagnostic, not direct selection.
 X=np.column_stack([v["c_"+col].to_numpy(float)/scales[col] for col in PRED])
 y=dy
 beta=np.linalg.lstsq(X,y,rcond=None)[0]*100
 terms=np.column_stack([np.einsum("ij,ik->ijk",X,X).reshape(len(v),len(PRED)**2),
                        X*y[:,None]])
 f=pd.DataFrame(terms)
 f["genus"]=v.genus.values
 z=f.groupby("genus",observed=True).sum().to_numpy(float)
 rng=np.random.default_rng(seed+degree*99+555)
 estimates=[]
 for indexes in rng.integers(0,len(z),(draws,len(z))):
  summed=z[indexes].sum(axis=0)
  A=summed[:len(PRED)**2].reshape(len(PRED),len(PRED))
  v2=summed[len(PRED)**2:]
  estimates.append(100*np.linalg.lstsq(A,v2,rcond=None)[0])
 ci=np.quantile(np.array(estimates),[.025,.975],axis=0)
 conditional=[{
   "grid_degrees":degree,"predictor":col,
   "cohort_genus_region_fully_partial_white_pp_per_SD":float(beta[i]),
   "cluster_genus_nominal95_partial_ci_pp":ci[:,i].tolist(),
   "sources_and_sample":int(len(v)),"model_rank":int(np.linalg.matrix_rank(X)),
   "not_a_causal_independent_climate_effect":True
  } for i,col in enumerate(PRED)]
 return results,conditional

def quintiles(x):
 t=aggregate(x,5);ans={}
 for name in ("log1p_annual_rainfall","mean_atmospheric_vapour","elevation"):
  q=t.copy()
  q["env_quintile"]=pd.qcut(q[name],5,labels=False,duplicates="drop")
  ans[name]=[
   {"low_to_high_quintile":int(v)+1,
    "taxon_geocells":int(len(g)),"white_percent":100*float(g.p_white.mean()),
    "median_annual_original_climate":float(g[name].median())}
   for v,g in q.groupby("env_quintile",observed=True)]
 return ans

def run(photo,climate,out):
 out=Path(out)
 need(not out.exists(),"Refuse overwrite of real newest photo gradients")
 x=prepare(photo,climate)
 outputs=[];conditional=[]
 for deg in (5,10):
  t=aggregate(x,deg)
  y,z=simple_slopes(t,deg,2 if deg==5 else 3)
  outputs.extend(y);conditional.extend(z)
 d={
  "schema":SCHEMA,"complete_original_actual_photo_attempts":1043930,
  "original_automatic_four_hue_ROI_success":459235,
  "actual_annual_climate_complete_named_genus_classifiable_photos":len(x),
  "source_four_hue_parquet_SHA256":PHOTO_SHA,"source_annual_climate_parquet_SHA256":ENV_SHA,
  "moisture_elevation_within_cohort_genus_geocell_univariate":outputs,
  "simultaneous_climate_correlation_diagnostic":conditional,
  "source_taxon_cell_environment_quintiles":quintiles(x),
  "source_unknowns_included_in_denominator_not_relabelled":1043930-459235,
  "open_flower_or_native_status_independently_verified":False,
  "causal_flower_colour_adaptation_identified":False,
  "old350399_results_are_not_new_data":True,
 }
 out.mkdir(parents=True)
 (out/"latest1043930_climate_macroecology.json").write_text(
  json.dumps(d,indent=2,sort_keys=True,allow_nan=False)+"\n")
 pd.DataFrame(outputs).to_csv(out/"latest1043930_water_elev_local_slopes.csv",index=False)
 pd.DataFrame(conditional).to_csv(out/"latest1043930_multivariable_climate_sensitivity.csv",index=False)
 print("MOIS_REAL_LATEST1043930_WATER_ELEVATION_TAXON_GENUS",json.dumps({
  "photo_attempts":1043930,"classifiable":459235,"climate_complete_named_genus":len(x),
  "5deg_within_cohort_genus":[v for v in outputs if v["grid_degrees"]==5
    and v["predictor"] in ("log1p_annual_rainfall","mean_atmospheric_vapour","elevation")],
  "5deg_joint_climate_diagnostic":[v for v in conditional
    if v["grid_degrees"]==5 and v["predictor"] in
     ("mean_atmospheric_vapour","elevation","log1p_annual_rainfall")]
 },sort_keys=True),flush=True)
 return d

def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument("--photo",type=Path,required=True)
 p.add_argument("--climate",type=Path,required=True)
 p.add_argument("--output",type=Path,required=True)
 a=p.parse_args();run(a.photo,a.climate,a.output)
if __name__=="__main__":main()
