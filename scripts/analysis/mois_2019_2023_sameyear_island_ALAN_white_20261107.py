#!/usr/bin/env python3
"""Actual 2019/2021/2022/2023 MOIS island × ALAN photograph ROI refit.

Only source SHA-locked, same observation-year VIIRS and R2024B UC population.
Year-specific photo assemblages, WCVP accepted-species-equal original strict9
white ROI fractions and GSHHG landmass identities, no photo downloads.
2022 October VIIRS is provider missing, not zero or another year's exposure.
All associations observational; no botanical white species disappearance claim.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scripts.analysis.mois_2024_sameyear_global_island_photo_white_20261107 import (
    COVARS, cell_id, cluster_fit)
# Frozen identifiers from the prior GSHHG source manifest (no heavy GIS import).
YEAR={
  2019:{"all":28082,"known":27020,"known_sha":"01f4e2d2e1c8bc9a1392f4143f693563bfa8ea1a0f65ba800e1caa450803831c"},
  2021:{"all":47714,"known":45966,"known_sha":"f0d7cb134c09c74cccfb69199396ac30bc5b663d39dc5369d6e296d8c3a1cdca"},
  2022:{"all":54163,"known":52075,"known_sha":"594fe8e729cf54aa44859a841b5d3e49e6a66fb85b217da8aff0160044e14df5"},
  2023:{"all":62289,"known":59862,"known_sha":"e2825386f45a77da6e3dacdf714fc940cf9ab9cc5174671672b4d57466d04fba"}
}

PHOTO_SHA="f033661d76e07898a5c13199a55c5d0ce6e55353daa48b0058c456979e90644c"
WCVP_SHA="dc5bb0e3686fde0cfd327e4336245d0ceae3c7657b3bd12841bb814cc94d1eb2"
YEAR_SET=(2019,2021,2022,2023)
def require(condition,message):
    if not bool(condition):
        raise ValueError("MOIS_YEAR_ISLAND_ALAN: "+message)
def sha(file):
    h=hashlib.sha256()
    with Path(file).open("rb") as f:
        for data in iter(lambda:f.read(4<<20),b""):
            h.update(data)
    return h.hexdigest()

def verify_environment(year,geo,viirs,pop,source_receipt,viirs_receipt,pop_receipt):
    require(year in YEAR_SET,"outside original year-supported source")
    contract=YEAR[year]
    require(sha(geo)==source_receipt["full_source_SHA256"],
            "Source original photo landmass SHA drift")
    require(source_receipt["knownland_source_SHA256"]==contract["known_sha"],
            "Original known-land photo target cohort changed")
    require(source_receipt["year"]==year and
            source_receipt["original_known_land_targets"]==contract["known"],
            "Original year geographic source count changed")
    require(sha(viirs)==viirs_receipt["values_SHA256"] and
            viirs_receipt["source_year"]==year,
            "Original source-year satellite numeric SHA drift")
    require(sha(pop)==pop_receipt["population_values_SHA256"] and
            pop_receipt["source_year"]==year and
            pop_receipt["source_photo_site_manifest_SHA256"]==contract["known_sha"],
            "Population source version or observed year drift")
    land=pd.read_csv(geo,dtype={"target_id":str,"cell_id":str,"land_id":str})
    v=pd.read_csv(viirs,dtype={"target_id":str},low_memory=False)
    p=pd.read_csv(pop,dtype={"target_id":str},low_memory=False)
    require(len(land)==contract["all"] and len(v)==len(p)==contract["known"],
            "Old year photo landmass/source population count drift")
    require(land.target_id.is_unique and v.target_id.is_unique and p.target_id.is_unique,
            "Source target identity duplicated")
    require(set(v.target_id)==set(p.target_id)==set(land.loc[
        land.land_class.isin(("mainland","island")),"target_id"]),
        "Same-year photo-site/landmass VIIRS and population target mismatch")
    require(v.target_id.str.startswith(str(year)+"|").all(),
            "VIIRS requested for wrong photographic year")
    require(p.population_source_year.eq(year).all() and
            p.population_source_year.notna().all() and
            p.source_population_model.eq("WorldPop Global2 R2024B UC v1").all(),
            "Population other year or release substituted")
    pcol=p[["target_id","population_status","population_density_per_km2"]]
    env=v.merge(pcol,on="target_id",validate="one_to_one")
    require(env.target_id.is_unique,"Duplicate population match")
    require(env.loc[env.population_status.ne("OK"),
         "population_density_per_km2"].isna().all(),
         "NoData population must not be converted to zero")
    require(env.loc[env.population_status.eq("OK"),
         "population_density_per_km2"].ge(0).all(),
         "Population zero and negative pixel values mixed")
    return land,env

def make_photo_cells(year,canonical,wcvp,land,env):
    require(len(canonical)==1032882 and canonical.photo_id.is_unique
            and canonical.observation_id.is_unique,
            "Original frozen photo/observation ID population changed")
    require(wcvp.inat_taxon_id.is_unique,"WCVP accepted taxon matches not unique")
    # A 0.05-degree grid may contain more than one island/mainland landmass.
    # Do not arbitrarily attribute its photographs to one landmass.
    ambiguous=land.cell_id.duplicated(keep=False)
    safe=land.loc[~ambiguous &
                  land.land_class.isin(("island","mainland")),
                  ["cell_id","target_id","land_class","land_id"]].copy()
    require(safe.cell_id.is_unique,"Geography safety gate was bypassed")
    photos=canonical.loc[
        canonical.observed_on.astype(str).str[:4].eq(str(year))&
        canonical.source_subarm.eq("historical_main")&
        canonical.latitude.between(-90,90)&
        canonical.longitude.between(-180,180)].copy()
    photos["cell_id"]=cell_id(photos.longitude,photos.latitude)
    photos=photos.merge(safe,on="cell_id",how="inner",validate="many_to_one")
    require(photos.photo_id.is_unique and photos.observation_id.is_unique,
            "Photographic unit duplicated by geography")
    accept=wcvp.loc[wcvp.wcvp_match_state.eq("unique_exact_name_to_accepted_id"),
        ["inat_taxon_id","wcvp_accepted_plant_name_id"]]
    photos=photos.merge(accept,on="inat_taxon_id",how="left",validate="many_to_one")
    attempt=photos.groupby("target_id").agg(
        photo_attempts=("photo_id","size"),
        strict9_success=("strict9","sum"),
        bio5=("wc2_bio05_10m","mean"),
        bio12=("wc2_bio12_10m","mean"),
        elevation=("wc2_elev_2p5m","mean")).reset_index()
    good=photos.loc[
        photos.strict9 & photos.strict9_white.notna() &
        photos.wcvp_accepted_plant_name_id.notna(),
        ["target_id","wcvp_accepted_plant_name_id","strict9_white"]].copy()
    require(good.strict9_white.between(0,1).all(),"Invalid original ROI white mass")
    species=good.groupby(["target_id","wcvp_accepted_plant_name_id"],
                          observed=True).strict9_white.mean().reset_index()
    per_cell=species.groupby("target_id",observed=True).agg(
        white_equal=("strict9_white","mean"),
        accepted_species_n=("wcvp_accepted_plant_name_id","nunique")).reset_index()
    d=env.merge(attempt,on="target_id",how="left",validate="one_to_one")
    d=d.merge(per_cell,on="target_id",how="left",validate="one_to_one")
    # Numeric satellite tables retain target_id but not the GSHHG sibling
    # land_id: recover that ID from the verified original year geometry.
    d=d.merge(safe[["target_id","land_id"]],on="target_id",how="inner",
              validate="one_to_one")
    require(d.target_id.is_unique and d.land_id.notna().all(),
            "Selected light cells lack authentic original landmass identity")
    d["photo_eligible"]=(d.photo_attempts.ge(5)&d.strict9_success.ge(3)&
                          d.accepted_species_n.ge(3))
    d["light_eligible"]=(d.viirs_valid_months.ge(3)&d.viirs_median.notna())
    d["pop_eligible"]=(d.population_status.eq("OK")&
                        d.population_density_per_km2.notna())
    d["white_pp"]=100*d.white_equal
    d["log_light"]=np.log1p(d.viirs_median.clip(lower=0))
    d["log_pop"]=np.log1p(d.population_density_per_km2.clip(lower=0))
    d["log_rain"]=np.log1p(d.bio12.clip(lower=0))
    d["log_attempts"]=np.log1p(d.photo_attempts.clip(lower=0))
    d["region5"]=(np.floor((d.query_longitude+180)/5).astype("int64").astype(str)+
                  "_"+np.floor((d.query_latitude+90)/5).astype("int64").astype(str))
    selected=d.loc[d.photo_eligible&d.light_eligible&d.pop_eligible].dropna(
        subset=["white_pp","log_light","log_pop","log_rain","bio5",
                "elevation","log_attempts","region5"]).copy()
    selected["is_island"]=selected.land_class.eq("island").astype(int)
    return photos,d,selected,int(ambiguous.sum())

def fit_if_supported(df,form,cluster,*,interaction=False,within=False):
    clusters=int(df[cluster].nunique())
    if clusters<8 or len(df)<40:
        return {"estimable":False,"reason":"too_few_independent_blocks_or_cells",
                "cells":len(df),"clusters":clusters}
    if interaction:
        islands=df.loc[df.land_class.eq("island")]
        if len(islands)<30 or islands.land_id.nunique()<10:
            return {"estimable":False,"reason":"too_few_island_cells_or_islands",
                    "island_cells":len(islands),
                    "independent_islands":int(islands.land_id.nunique()),
                    "minimum_island_cells":30,"minimum_islands":10,
                    "cells":len(df),"clusters":clusters}
    if within and clusters<10:
        return {"estimable":False,"reason":"too_few_independent_within_islands",
                "cells":len(df),"clusters":clusters}
    result=cluster_fit(df,form,cluster)
    return {"estimable":True,"cells":len(df),"clusters":clusters,"terms":result}

def run(year,paths,out):
    require(year in YEAR_SET,"Only original 2019/21/22/23 source years")
    out=Path(out)
    require(not out.exists(),"Never overwrite existing source result")
    for key,wanted in [("canonical",PHOTO_SHA),("wcvp",WCVP_SHA)]:
        require(sha(paths[key])==wanted,"Frozen source hash drift "+key)
    source_receipt=json.loads(Path(paths["source_receipt"]).read_text())
    viirs_receipt=json.loads(Path(paths["viirs_receipt"]).read_text())
    pop_receipt=json.loads(Path(paths["pop_receipt"]).read_text())
    geo,env=verify_environment(year,paths["geo"],paths["viirs"],paths["population"],
                               source_receipt,viirs_receipt,pop_receipt)
    photo=pd.read_csv(paths["canonical"],compression="gzip",low_memory=False,
        dtype={"photo_id":str,"observation_id":str,"inat_taxon_id":str},
        usecols=["photo_id","observation_id","inat_taxon_id","source_subarm",
                 "observed_on","latitude","longitude","strict9","strict9_white",
                 "wc2_bio05_10m","wc2_bio12_10m","wc2_elev_2p5m"])
    wcvp=pd.read_csv(paths["wcvp"],compression="gzip",dtype=str,
        usecols=["inat_taxon_id","wcvp_match_state","wcvp_accepted_plant_name_id"])
    photos,allcells,d,ambiguous=make_photo_cells(year,photo,wcvp,geo,env)
    # Report 2022 provider 400 explicitly, rather than claim four-month parity.
    october_provider_missing=bool(year==2022)
    if october_provider_missing:
        require(int(env.viirs_10_status.eq("PROVIDER_RASTER_400").sum())>0,
                "Cannot silently treat 2022 October provider failure as real nightlight")
    else:
        require(not env.viirs_10_status.eq("PROVIDER_RASTER_400").any(),
                "Unannounced missing year-specific VIIRS quarter")
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
    descriptive=d.groupby("land_class").agg(cells=("target_id","size"),
        white_roi_mean_pp=("white_pp","mean"),
        median_viirs=("viirs_median","median"),
        accepted_species_cell_sum=("accepted_species_n","sum")).reset_index()
    output={
        "schema":"MOIS_source_photo_original_sameyear_ALAN_islands_R2024B_UC_v1",
        "year":year,"population_release":"WorldPop Global2 R2024B UC v1",
        "VIIRS_year":year,"photographic_original_arm":"historical_main",
        "original_source_lands":len(geo),
        "original_knownland_environment_sites":len(env),
        "source_VIIRS_3months_sites":int((env.viirs_valid_months>=3).sum()),
        "source_population_valid_sites":int(env.population_status.eq("OK").sum()),
        "source_mixed_landmass_cell_rows_omitted":ambiguous,
        "original_historical_photos_geolocated_and_matched":len(photos),
        "original_strict9_photo_success_matched":int(photos.strict9.sum()),
        "photographic_eligible_with_sameyear_light_and_pop_cells":len(d),
        "photographic_eligible_island_cells":len(island),
        "photographic_eligible_mainland_cells":len(d)-len(island),
        "within_island_supported_cells":len(within),
        "within_island_independent_islands":len(ids),
        "2022_Oct_VIIRS_missing_not_replaced":october_provider_missing,
        "four_month_comparability_to_other_years":not october_provider_missing,
        "year_specific_exploratory_not_independent_replication":True,
        "models":fits,
        "original_photo_white_ROI_not_botanical_white_species_fraction":True,
        "no_new_photos_or_ROI_inferences":True,
        "technical_unknown_colour_imputed":False,
        "island_minimum_area_or_latitude_filter":None,
        "no_night_pollination_or_causal_species_loss_claim":True,
        "source_SHA256":{k:sha(v) for k,v in paths.items()}}
    out.mkdir(parents=True)
    descriptive.to_csv(out/"sameyear_island_mainland_photo_white_descriptive.csv",index=False)
    output["output_descriptive_SHA256"]=sha(out/"sameyear_island_mainland_photo_white_descriptive.csv")
    (out/"source_receipt.json").write_text(json.dumps(output,indent=2,sort_keys=True)+"\n")
    print("MOIS_YEAR_REAL_PHOTO_ALAN_ISLAND_WHITE",json.dumps({
        key:output[key] for key in ["year","photographic_eligible_with_sameyear_light_and_pop_cells",
       "photographic_eligible_island_cells","within_island_independent_islands",
       "four_month_comparability_to_other_years","models"]},sort_keys=True),flush=True)
    return output

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--year",type=int,required=True)
    for key in ("canonical","wcvp","geo","source_receipt","viirs",
                "viirs_receipt","population","pop_receipt"):
        p.add_argument("--"+key.replace("_","-"),dest=key,type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    run(a.year,{k:getattr(a,k) for k in ("canonical","wcvp","geo","source_receipt",
        "viirs","viirs_receipt","population","pop_receipt")},a.out)
