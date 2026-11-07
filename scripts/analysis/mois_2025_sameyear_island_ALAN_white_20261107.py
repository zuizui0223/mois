#!/usr/bin/env python3
"""2025 original photo white ROI against TRUE 2025 VIIRS and 2025 UC population.

Uses signed 2025 original-photo 0.05-degree GSHHG landmass source and 73,463
same-year VIIRS+R2024B UC population targets; no older-year population proxies.
No photo downloads, new ROI inference, botanical native-flora or moth claims.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
from scripts.analysis.mois_2019_2023_sameyear_island_ALAN_white_20261107 import (
 PHOTO_SHA,WCVP_SHA,sha,require,make_photo_cells,fit_if_supported,COVARS)

GEO_SHA="de2cd14690383f7e66cff68c945cf2a6812b33b0511f0189aaea9f6d62ea37fc"
SITES_SHA="45c5152e1edad997ff8749c645818996fb2b0bc03a2ab7c0753af446afaab143"
def run(paths,out):
    out=Path(out)
    require(not out.exists(),"Refuse replacing old source results")
    require(sha(paths["canonical"])==PHOTO_SHA and sha(paths["wcvp"])==WCVP_SHA,
            "Frozen original photo or WCVP accepted-taxa source changed")
    geo_proof=json.loads(Path(paths["geo_receipt"]).read_text())
    require(geo_proof.get("all_knownland_site_SHA256")==SITES_SHA and
            geo_proof.get("full_source_reference_SHA256")==GEO_SHA and
            geo_proof.get("full_source_SHA256")==sha(paths["geo"]) and
            geo_proof.get("known_mainland_targets")==63408 and
            geo_proof.get("known_island_targets")==10055,
            "2025 source reference or verified landmass/coordinate identity drift")
    # Exact known-land source photo ID / land class / coordinates were pinned.
    # Auxiliary GEOS area floating point bytes may vary; do not fake SHA identity.
    require(geo_proof.get("full_source_bytes_identical_to_original") or
            geo_proof.get("area_roundoff_exception_only_with_byte_identical_knownland_manifest"),
            "Unverified 2025 full geometry source discrepancy")
    source=json.loads(Path(paths["numeric_receipt"]).read_text())
    require(source["source_photo_year"]==source["VIIRS_year"]==
            source["population_raster_year"]==2025 and
            source["original_source_photo_grid_sha256"]==SITES_SHA,
            "Original 2025 image/ALAN/population year is not identical")
    require(source["original_target_count"]==73463 and
            sha(paths["numeric"])==source["output_values_SHA256"],
            "Original 2025 environmental values byte SHA drift")
    geo=pd.read_csv(paths["geo"],dtype={"target_id":str,"cell_id":str,"land_id":str})
    env=pd.read_csv(paths["numeric"],dtype={"target_id":str},low_memory=False)
    require(len(geo)==76630 and len(env)==73463 and
            geo.target_id.is_unique and env.target_id.is_unique,
            "Year-2025 GSHHG originals not full or unique")
    require(set(env.target_id)==set(geo.loc[
            geo.land_class.isin(("mainland","island")),"target_id"]),
            "2025 true population/light target was not same source landmass")
    require(env.target_id.str.startswith("2025|").all(),
            "Non2025 photographic ALAN source in 2025 year model")
    require(env.loc[env.population_source_status.ne("OK"),
         "population_density"].isna().all(),"Population NoData replaced with zero")
    require(env.loc[env.population_source_status.eq("OK"),
         "population_density"].ge(0).all(),"Invalid actual population density")
    env=env.rename(columns={
       "population_source_status":"population_status",
       "population_density":"population_density_per_km2"})
    photo=pd.read_csv(paths["canonical"],compression="gzip",low_memory=False,
        dtype={"photo_id":str,"observation_id":str,"inat_taxon_id":str},
        usecols=["photo_id","observation_id","inat_taxon_id","source_subarm",
                 "observed_on","latitude","longitude","strict9","strict9_white",
                 "wc2_bio05_10m","wc2_bio12_10m","wc2_elev_2p5m"])
    wcvp=pd.read_csv(paths["wcvp"],compression="gzip",dtype=str,
        usecols=["inat_taxon_id","wcvp_match_state","wcvp_accepted_plant_name_id"])
    photos,allcells,d,ambiguous=make_photo_cells(2025,photo,wcvp,geo,env)
    base="white_pp ~ log_light + log_pop + "+COVARS+" + C(region5)"
    form="white_pp ~ log_light * is_island + log_pop + "+COVARS+" + C(region5)"
    fits={
       "geographic_region_global":fit_if_supported(d,base,"region5"),
       "island_by_light_between_landmass_classes":fit_if_supported(
         d,form,"region5",interaction=True)}
    island=d.loc[d.is_island.eq(1)].copy()
    g=island.groupby("land_id").agg(cells=("target_id","nunique"),
                                    levels=("viirs_median","nunique"))
    ids=g.index[(g.cells>=2)&(g.levels>=2)]
    within=island.loc[island.land_id.isin(ids)].copy()
    fits["within_same_island"]=fit_if_supported(within,
        "white_pp ~ log_light + log_pop + "+COVARS+" + C(land_id)",
        "land_id",within=True)
    desc=d.groupby("land_class").agg(cells=("target_id","size"),
         mean_white_ROI_pp=("white_pp","mean"),
         median_VIIRS=("viirs_median","median")).reset_index()
    output={
       "schema":"MOIS_2025_original_2025_VIIRS_pop_photo_white_island_20261107",
       "year":2025,
       "population_source":"2025 WorldPop Global2 R2024B UC v1",
       "photo_year_VIIRS_year_population_year_all_2025":True,
       "sameyear_original_environment_known_land":len(env),
       "sameyear_original_environment_complete":int(env.complete.sum()),
       "historical_main_photo_attempts_in_unambiguous_landmass":len(photos),
       "historic_strict9_success_in_unambiguous_landmass":int(photos.strict9.sum()),
       "photographic_eligible_sameyear_alanspatial_population_cells":len(d),
       "island_eligible_photo_cells":len(island),
       "within_island_supported_cells":len(within),
       "independent_islands_with_exposure_variation":len(ids),
       "cross_year_2020_population_never_used":True,
       "minimum_island_area_exclusion":None,
       "ambiguous_multilandmass_cell_rows_not_falsely_allocated":ambiguous,
       "models":fits,
       "source_SHA256":{k:sha(v) for k,v in paths.items()},
       "same_photo_original_strict9_unchanged":True,
       "technical_unknown_relabels":0,
       "new_photo_downloads_or_ROI_inference":0,
       "no_causal_moth_or_native_species_extirpation_claim":True}
    out.mkdir(parents=True)
    desc.to_csv(out/"2025_sameyear_photo_island_mainland_descriptive.csv",index=False)
    output["descriptive_sha256"]=sha(out/"2025_sameyear_photo_island_mainland_descriptive.csv")
    (out/"source_receipt.json").write_text(json.dumps(output,indent=2,sort_keys=True)+"\n")
    print("MOIS_2025_REAL_SAMEYEAR_ISLAND_ALAN_WHITE",json.dumps({
        k:output[k] for k in ["year",
            "photographic_eligible_sameyear_alanspatial_population_cells",
            "island_eligible_photo_cells",
            "independent_islands_with_exposure_variation",
            "models"]},sort_keys=True),flush=True)
    return output
if __name__=="__main__":
    p=argparse.ArgumentParser()
    for key in ("canonical","wcvp","geo","geo_receipt","numeric","numeric_receipt"):
        p.add_argument("--"+key.replace("_","-"),dest=key,required=True,type=Path)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    run({k:getattr(a,k) for k in ("canonical","wcvp","geo","geo_receipt","numeric","numeric_receipt")},a.out)
