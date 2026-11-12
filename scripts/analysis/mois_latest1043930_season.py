#!/usr/bin/env python3
"""Latest 1,043,930 genuinely measured original MOIS photo seasonal composition.

All original photo IDs and four-hue labels signed, unknowns preserved.
The Kitagawa four terms are photographic taxon/observation composition,
not a trait plasticity or independent phenology causal effect.
"""
from pathlib import Path
import argparse, json, hashlib
import pandas as pd
from scripts.analysis.mois_real350399_taxon_season_decomposition_20261019 import estimate

HUES=("white","yellow_orange","red_pink","blue_purple")
FULL_SHA="9a4f3dae6fbcee967b71cb29576e71d52c48e1f1f790b7d3ea4750c4bcb0e6d4"
def sha(path):
 h=hashlib.sha256()
 with Path(path).open("rb") as f:
  for b in iter(lambda:f.read(4<<20),b""):h.update(b)
 return h.hexdigest()
def run(photo,receipt,out):
 photo,receipt,out=Path(photo),Path(receipt),Path(out)
 r=json.loads(receipt.read_text())
 assert r["source_pixel_attempts"]==1043930 and r["source_four_hue_classified"]==459235
 assert r["photo_parquet_sha256"]==FULL_SHA==sha(photo)
 d=pd.read_parquet(photo)
 assert len(d)==1043930 and d.photo_id.is_unique
 d["classified"]=d.known_image_roi.eq(1)
 d["white"]=d.morph.eq("white")
 d["date"]=pd.to_datetime(d.date,errors="coerce")
 assert int(d.classified.sum())==459235
 output=estimate(d,seed=20261025,draws=299)
 result={
  "schema":"mois_full1043930_real_photos_original_season_taxon_decomposition_v1",
  "source_sha256":FULL_SHA,"full_photo_attempts":len(d),
  "classified_four_hues":int(d.classified.sum()),
  "unknown_photo_ROIs":int((~d.classified).sum()),
  "classifiable_dated_extratropical_summer_winter":
   output["classified_extratropical_dated_summer_winter_photos"],
  "hemispheres":output["hemispheres"],
  "observation_equal_sensitivity":output["observation_equal_sensitivity"],
  "source_arm_seasonal_photo_counts":output["source_arm_seasonal_photo_counts"],
  "colour_plasticity_or_blooming_seasonality_proved":False,
  "unknowns_imputed_as_nonwhite":False,
  "all_original_image_bytes_newly_downloaded":0
 }
 out.parent.mkdir(parents=True,exist_ok=True)
 out.write_text(json.dumps(result,sort_keys=True,indent=2,allow_nan=False)+"\n")
 print("MOIS_LATEST1043930_REAL_SEASON_TAXON_DECOMPOSITION",json.dumps({
  "source":len(d),"classified":int(d.classified.sum()),
  "hemispheres":{k:{
      "summer_minus_winter_white_pp":v["raw_photo_white_summer_minus_winter_pp"],
      "photo_frequency_component_pp":v["changed_common_taxon_relative_photo_frequency_component_pp"],
      "shared_taxon_rate_component_pp":v["within_common_taxon_photographic_white_rate_component_pp"],
      "summer_photos":v["summer_classified_original_photos"],
      "winter_photos":v["winter_classified_original_photos"]
   } for k,v in output["hemispheres"].items()}
 },sort_keys=True),flush=True)

if __name__=="__main__":
 p=argparse.ArgumentParser()
 p.add_argument("--photo",type=Path,required=True)
 p.add_argument("--receipt",type=Path,required=True)
 p.add_argument("--output",type=Path,required=True)
 x=p.parse_args();run(x.photo,x.receipt,x.output)
