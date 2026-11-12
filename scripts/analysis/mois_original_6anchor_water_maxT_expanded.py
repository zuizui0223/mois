#!/usr/bin/env python3
"""MOIS global humidity/white-vs-three-other-hues: ORIGINAL 2026-10-13 estimator.

Reuse actual MOIS 248312 photo sites (120807 genuine automated four-hue labels),
identical original species-mean -> nominal-genus geographic and 999-genus
bootstrap maxT estimator, and original six fixed anchors. Not an FCP within-
species D/W statistic, independent dataset, or proof of causal adaptation.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
import rasterio
from scripts.analysis.fcp_abiotic_hue_direction_macro_20261013 import (
    CLASSES, GEO, ANCHOR_CLIMATE, BOOT_SEED, species_frame,
    cluster_suff, regress_boot, geo_cell
)
from scripts.analysis.mois_sample_248312_worldclim_static_and_monthnormals import (
    sample_tiff, find_zip_tiff
)

SOURCE_N=248312
HUE_N=120807
VAPR_ZIP="wc2.1_10m_vapr.zip"
CLIMATE={
  "wc_bio5":"wc2_bio05_10m",
  "wc_bio12":"wc2_bio12_10m",
  "wc_bio14":"wc2_bio14_10m",
  "wc_bio15":"wc2_bio15_10m",
  "wc_srad_annual_mean":"wc2_annual_srad_normal_10m_Wm2"
}

def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda:f.read(1<<20),b""):h.update(b)
    return h.hexdigest()

def require(flag,why):
    if not flag:raise ValueError(why)

def photo_and_moisture(site_gz:Path,receipt:Path,vapr_dir:Path):
    source=json.loads(receipt.read_text(encoding="utf-8"))
    require(sha(site_gz)==source["full_photo_static_monthly_normal_csv_sha256"],
            "Full authentic 248k original photographic site/raster source changed")
    require((source["n_original_photo_ids"],source["n_four_hue_classified"])==(SOURCE_N,HUE_N),
            "Original frozen 248312 photos and 120807 hues changed")
    columns=["photo_id","inat_taxon_id","species","morph","original_source",
             "latitude","longitude",*CLIMATE.values()]
    df=pd.read_csv(site_gz,compression="gzip",usecols=columns,
                   low_memory=False,dtype={"photo_id":str,"inat_taxon_id":str})
    require(len(df)==SOURCE_N and df.photo_id.nunique()==SOURCE_N,
            "Literal original photo IDs duplicated")
    require(int(df.morph.isin(CLASSES).sum())==HUE_N,"Classifiable hue original denominator changed")
    require(set(ANCHOR_CLIMATE)==set((*CLIMATE.keys(),"wc_vapr_annual_mean")),
            "Not the original six anchored ecology gradients")
    src=df.loc[df.morph.isin(CLASSES)].copy()
    coords_lon=pd.to_numeric(src.longitude,errors="coerce").to_numpy(float)
    coords_lat=pd.to_numeric(src.latitude,errors="coerce").to_numpy(float)
    n=len(src)
    vapr_sum=np.zeros(n,np.float64)
    vapr_n=np.zeros(n,np.int8)
    # Must use EVERY monthly normal (mean pressure across the year), not
    # observation-date weather. WorldClim vapr values are in kPa.
    for m in range(1,13):
        path=find_zip_tiff(vapr_dir,"vapr",f"wc2.1_10m_vapr_{m:02d}.tif")
        value=sample_tiff(path,coords_lon,coords_lat)
        good=np.isfinite(value)
        vapr_sum[good]+=value[good]
        vapr_n[good]+=1
    vapr=np.full(n,np.nan)
    complete=vapr_n==12
    vapr[complete]=vapr_sum[complete]/12.
    src["wc_vapr_annual_mean"]=vapr
    for official,new in CLIMATE.items():
        src[official]=pd.to_numeric(src[new],errors="coerce")
    # Historical cached value was labelled W/m² while original WorldClim
    # solar source is kJ/m²/day. 1 W/m² = 86.4 kJ/m²/day.
    src["wc_srad_annual_mean"]/=86.4
    src["genus"]=src.species.str.split().str[0]
    src["abs_lat"]=pd.to_numeric(src.latitude,errors="coerce").abs()
    src["lon_sin"]=np.sin(np.deg2rad(pd.to_numeric(src.longitude,errors="coerce")))
    src["lon_cos"]=np.cos(np.deg2rad(pd.to_numeric(src.longitude,errors="coerce")))
    ok_genus=src.genus.str.match(r"^[A-Z][a-zA-Z]+$",na=False)
    src=src.loc[ok_genus].copy()
    required=[*GEO,*ANCHOR_CLIMATE]
    src=src.dropna(subset=required)
    require(len(src)>90000,"Real global humidity-climate source coverage catastrophically reduced")
    return src,{
       "n_original_classified_four_hues":HUE_N,
       "n_real_248k_photo_sites_with_full_six_anchored_climate":len(src),
       "n_named_taxa_before_congeneric_filter":int(src.inat_taxon_id.nunique()),
       "n_source_photos_original_historical":
           int(src.original_source.eq("historical").sum()),
       "n_source_photos_added_from_FCP":
           int(src.original_source.ne("historical").sum()),
       "official_WorldClim_vapr_zip_sha256":sha(vapr_dir/"vapr"/VAPR_ZIP),
       "vapr_12_months_complete_photo_count":int(complete.sum()),
       "worldclim_bio12_annual_rainfall_units":"mm/year",
       "worldclim_vapr_annual_mean_units":"kPa",
       "worldclim_srad_annual_mean_units":"W/m2, converted from original kJ/m2/day by /86.4",
       "WorldClim1970_2000_monthly_normals_not_2025_weather":True,
    }

def one_scope(src,scope,nboot,extra_seed):
    means=src[[*GEO,*ANCHOR_CLIMATE]].mean()
    scales=src[[*GEO,*ANCHOR_CLIMATE]].std(ddof=0).clip(lower=1e-12)
    species=species_frame(src,ANCHOR_CLIMATE)
    if scope=="matched_coarse_cell":
        groups,meta=cluster_suff(species,ANCHOR_CLIMATE,means,scales,
                                 geos=False,cell_conditioned=True)
    else:
        groups,meta=cluster_suff(species,ANCHOR_CLIMATE,means,scales,geos=True)
    if groups is None:
        return {"scope":scope,"status":"INSUFFICIENT_CONGENERIC_REPLICATES",**meta}
    outcome=regress_boot(groups,ANCHOR_CLIMATE,nboot,
                         BOOT_SEED+extra_seed,
                         equal_genus=True,cell_conditioned=scope=="matched_coarse_cell")
    require(outcome["status"]=="OK","Real original six-anchor maxT failed: "+str(outcome))
    return {"scope":scope,
            "n_source_photo_ids":len(src),
            "n_source_taxa_before_genus_filter":int(src.inat_taxon_id.nunique()),
            "genus_equal":True,
            "same_source_999_bootstrap_across_all_4x6_anchors":True,
            **meta,**outcome}

def run(site:Path,site_receipt:Path,vapr:Path,out:Path,nboot:int):
    require(nboot==999,"Only 999 original-method bootstrap resamples admitted")
    require(not out.exists() or not any(out.iterdir()),"Frozen outcomes cannot be overwritten")
    src,support=photo_and_moisture(site,site_receipt,vapr)
    allworld=one_scope(src,"world_genus_geo",nboot,44)
    old=src.loc[src.original_source.eq("historical")]
    historic=one_scope(old,"old_historical_only_genus_geo",nboot,96)
    matched=one_scope(src,"matched_coarse_cell",nboot,944)
    report={
      "schema":"mois_expanded_original_six_anchor_four_hue_water_maxT_v1",
      "status":"ORIGINAL_GENUS_EQUAL_999_BOOTSTRAP_MAXT_EXECUTED_NEW_PHOTO_SOURCE",
      "source_original_photo_table_SHA256":sha(site),
      "source_raster_receipt_SHA256":sha(site_receipt),
      "source":{"main":support,"base_historical_frame_photo_count":100543,
                "expanded_combined_photo_source_count":SOURCE_N},
      "models":{"global_expanded":allworld,"historical_only":historic,"coarse_cell":matched},
      "original_null_and_estimator_import":"scripts/analysis/fcp_abiotic_hue_direction_macro_20261013.py::species_frame,cluster_suff,regress_boot",
      "exact_original_six_abiotic_anchors":list(ANCHOR_CLIMATE),
      "four_original_floral_photo_hues":list(CLASSES),
      "original_familywise_maxT":"24 four-hue x fixed-six-abiotic-anchor contrasts, joint genus multinomial draws, 999 iterations, within scope",
      "new_source_BIO12_white_result_is_not_independent_of_previous_MOIS":True,
      "FCP_within_species_D_W_local_photo_pair_statistic_not_reused":True,
      "strict_image_quality_source_columns_not_restored_in_current_expanded_photo_union":True,
      "causal_water_selection_not_established":True,
      "month_of_observation_temperature_or_VPD_NOT_fitted":True,
    }
    out.mkdir(parents=True,exist_ok=True)
    (out/"mois_248k_original_six_anchor_moisture_maxT.json").write_text(
         json.dumps(report,indent=2,ensure_ascii=False,sort_keys=True,allow_nan=False)+"\n")
    lines=[]
    for k,mod in report["models"].items():
        if mod["status"]!="OK":continue
        for e in mod["estimates"]:
            lines.append({"scope":k,
                "n_source_photo_ids":mod["n_source_photo_ids"],
                "n_qualified_genus":mod["n_named_genus"],
                "n_qualified_species":mod["n_species"],
                "bootstrap":mod["n_bootstrap"],
                **e})
    import csv
    from scripts.analysis.fcp_abiotic_hue_direction_macro_20261013 import flatten
    import pandas as pd
    table=pd.DataFrame(lines)
    table.to_csv(out/"mois_248k_original_maxT_water_four_hue_ledger.csv",index=False)
    for block in ("global_expanded","historical_only"):
        m=report["models"][block]
        p=next(z for z in m["estimates"] if z["anchor"]=="wc_bio12" and z["hue"]=="white")
        print("MOIS_ORIGINAL_999BOOT_6ANCHOR_RAINFALL_WHITE",block,
              "photos",m["n_source_photo_ids"],
              "species",m["n_species"],"genera",m["n_named_genus"],
              "pp",p["point_pp_per_original_photo_env_sd"],
              "maxT",p["scope_anchor_hue_maxT_ci95_pp"],flush=True)
    return report

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--real-photo-worldclim-gzip",type=Path,required=True)
    p.add_argument("--real-worldclim-receipt",type=Path,required=True)
    p.add_argument("--vapr-raster-dir",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--bootstrap",type=int,default=999)
    a=p.parse_args()
    run(a.real_photo_worldclim_gzip,a.real_worldclim_receipt,a.vapr_raster_dir,
        a.output,a.bootstrap)
