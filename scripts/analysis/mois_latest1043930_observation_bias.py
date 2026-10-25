#!/usr/bin/env python3
"""Actual 1,043,930 original MOIS photograph colour observation bias by four axes.

Colour nonresponse is NOT no-pigment. Explicitly report both P(classified|source)
and P(WHITE|classified), plus source observed-white and upper/lower unknown
bounds for each climate quantile, altitude band, island/mainland and season.
No replacement of unknowns and no source ecological prevalence inference.
"""
from __future__ import annotations
import json,hashlib,argparse
from pathlib import Path
import numpy as np,pandas as pd
HUES=("white","yellow_orange","red_pink","blue_purple")
SCHEMA="mois_real1043930_five_axis_photo_colour_classifiability_bias_v1"
PHOTO_SHA="9a4f3dae6fbcee967b71cb29576e71d52c48e1f1f790b7d3ea4750c4bcb0e6d4"
CLIM_SHA="ec39b49649acc2ca4a468c00c86d05800e33ea2db1f91a3b119076be38432735"
def need(cond,msg):
 if not cond:raise ValueError("MOIS1043930_SELECTION: "+msg)
def sha(path):
 h=hashlib.sha256()
 with Path(path).open("rb") as f:
  for chunk in iter(lambda:f.read(4<<20),b""):h.update(chunk)
 return h.hexdigest()
def group_rows(d,key):
 ans=[]
 for value,x in d.groupby(key,dropna=False,observed=True):
  n=len(x); c=int(x.classified.sum()); white=int(x.white.sum())
  row={
   "stratum_variable":key,"stratum":str(value),
   "source_original_photo_attempts":n,
   "four_hue_classified":c,
   "automated_four_hue_unknown":n-c,
   "classified_fraction":float(c/n) if n else None,
   "observed_white_classified_photo_ROIs":white,
   "observed_white_per_source_photo":float(white/n) if n else None,
   "white_share_AMONG_classifiable":float(white/c) if c else None,
   "white_photo_fraction_lower_bound_all_unknown_nonwhite":float(white/n) if n else None,
   "white_photo_fraction_upper_bound_all_unknown_white":float((white+n-c)/n) if n else None,
   "taxa_with_original_four_hue":int(x.loc[x.classified,"inat_taxon_id"].nunique()),
  }
  ans.append(row)
 return ans
def run(photo,clim,geo,out):
 photo,clim,geo,out=map(Path,(photo,clim,geo,out))
 need(not out.exists(),"Cannot overwrite original source photo observation bias result")
 need(sha(photo)==PHOTO_SHA and sha(clim)==CLIM_SHA,
      "Original latest photo/climate SHA changed")
 p=pd.read_parquet(photo);c=pd.read_parquet(clim)
 g=pd.read_parquet(geo)
 need(len(p)==len(c)==len(g)==1043930,"Original row count not all 1.043m")
 need(p.photo_id.tolist()==c.photo_id.tolist()==g.photo_id.tolist(),
      "Original source photo ID alignment no longer identical")
 d=p[["photo_id","inat_taxon_id","date","latitude","longitude","morph","arm"]].copy()
 d["classified"]=d.morph.isin(HUES)
 d["white"]=d.morph.eq("white")
 need(int(d.classified.sum())==459235 and int(d.white.sum())==219600,
      "Original four-hue terminal unexpectedly changed")
 d["land_type"]=g.land_type.to_numpy()
 d["area_km2"]=g.area_km2.to_numpy()
 d["water"]=pd.to_numeric(c.wc2_annual_vapr_10m,errors="coerce").to_numpy(float)
 d["rain"]=pd.to_numeric(c.wc2_bio12_10m,errors="coerce").to_numpy(float)
 d["altitude"]=pd.to_numeric(c.wc2_elev_2p5m,errors="coerce").to_numpy(float)
 d["date"]=pd.to_datetime(d.date,errors="coerce")
 d["cohort"]=np.where(d.arm.str.startswith("FCP5_"),"FCP5",d.arm).astype(str)
 for name,col in (("vapour_quintile","water"),("rain_quintile","rain")):
  x=d[col].to_numpy(float)
  finite=np.isfinite(x)
  out_q=np.zeros(len(x),dtype=np.int8)
  out_q[finite]=np.asarray(pd.qcut(x[finite],5,labels=False,duplicates="drop"),dtype=np.int8)+1
  d[name]=np.where(out_q==0,"unknown",out_q.astype(str))
 # Elevation category must preserve invalid values distinctly.
 a=d.altitude.to_numpy(float)
 d["altitude_band"]=np.select(
  [(a>=0)&(a<500),(a>=500)&(a<1500),(a>=1500)&(a<2500),
   (a>=2500)&(a<10000)],
  ["0_500m","500_1500m","1500_2500m","2500m_plus"],
  default="unknown")
 m=d.date.dt.month
 north=d.latitude.ge(23.5);south=d.latitude.le(-23.5)
 summer=((north&m.isin((6,7,8)))|(south&m.isin((12,1,2))))
 winter=((north&m.isin((12,1,2)))|(south&m.isin((6,7,8))))
 d["season"]=np.select([north&summer,north&winter,south&summer,south&winter],
  ["north_summer","north_winter","south_summer","south_winter"],
  default="other_or_undated")
 outrows=[]
 for key in ("vapour_quintile","rain_quintile","altitude_band",
             "land_type","season","cohort"):
  outrows.extend(group_rows(d,key))
 overall={
  "original_true_photo_attempts":1043930,
  "four_hue_classified":459235,
  "unknown_four_hue":1043930-459235,
  "all_original_white_ROIs":219600,
  "observed_white_rate_among_classifiable":219600/459235,
  "automatic_four_colour_classification_rate":459235/1043930,
  "source_white_fraction_lower_bound_under_unknown_extreme":219600/1043930,
  "source_white_fraction_upper_bound_under_unknown_extreme":(219600+1043930-459235)/1043930,
  "no_unknown_imputed_as_white_or_chromatic":True
 }
 # Report differences in selection probability by geography/season; do not
 # renormalize unknown image labels into white or nonwhite.
 out.mkdir(parents=True)
 tab=pd.DataFrame(outrows)
 tab.to_csv(out/"latest1043930_classifiability_white_by_water_altitude_island_season.csv",index=False)
 r={"schema":SCHEMA,"source_photo_sha256":PHOTO_SHA,
    "source_climate_sha256":CLIM_SHA,"source_GSHHG_photo_parquet_SHA256":sha(geo),
    "overall":overall,"classifiability_and_observed_white_strata":outrows,
    "this_is_photo_sampling_and_ROI_success_not_botanical_organ_fitness":True}
 (out/"latest1043930_selection_bias_report.json").write_text(
  json.dumps(r,indent=2,sort_keys=True,allow_nan=False)+"\n")
 print("MOIS_LATEST1043930_CLASSIFICATION_SELECTION_BY_FOUR_GRADIENTS",json.dumps({
  "overall":overall,
  "strata":[x for x in outrows if (
    (x["stratum_variable"]=="land_type" and x["stratum"] in ("island","mainland")) or
    (x["stratum_variable"]=="season" and x["stratum"]!="other_or_undated") or
    (x["stratum_variable"]=="vapour_quintile" and x["stratum"] in ("1","5")) or
    (x["stratum_variable"]=="altitude_band" and x["stratum"] in ("0_500m","2500m_plus"))
  )]
 },sort_keys=True),flush=True)
 return r
def main():
 p=argparse.ArgumentParser(description=__doc__)
 for k in ("photo","clim","geo","output"):p.add_argument("--"+k,required=True,type=Path)
 x=p.parse_args()
 run(x.photo,x.clim,x.geo,x.output)
if __name__=="__main__":main()
