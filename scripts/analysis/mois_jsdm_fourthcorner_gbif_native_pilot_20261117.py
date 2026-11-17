#!/usr/bin/env python3
"""MOIS Fourth-corner JSDM design matrix pilot for GBIF/photo9 flower colour.

Same-source true 5,493,931 distinct GBIF species×EPSG6933 100km cell pairs,
18,577 global cells, signed WCVP broad TDWG3 native-range statuses, true
WorldClim LGM (3 GCM), 2001–2023 annual MOD17 NPP, NDVI, woody AGB.

THIS GENERATES AN **APPARENT RECORDING** RESPONSE, NOT A SURVEY
OCCUPANCY OR TRUE PRESENCE/ABSENCE. A zero records 'eligible-to-occur by
WCVP TDWG3 but not recorded in this GBIF snapshot'; source coverage,
detection and any true absences remain unidentifiable. No claim of model
fit or true ancient refugium from this pilot.

Panels predeclared by latitude and longitude, never by flower colour outcome:
  North America 25°–55°N, -125°–-65°E
  Europe 35°–65°N, -15°–40°E
  East Asia 25°–55°N, 95°–145°E

Species chosen on native GBIF presence support and immutable stable hash,
never on hue, genus, historical response, or fitted p-value. No genus term.
"""
from __future__ import annotations
import argparse, hashlib, io, json
from pathlib import Path
from zipfile import ZipFile
import numpy as np
import pandas as pd

from scripts.analysis.mois_lgm_climate_history_species_20261116 import (
    read_signed_sources,require,sha
)
from scripts.analysis.mois_lgm_wcvp_native_regional_block_20261116 import (
    get_wcvp_source, assign_cells_broad_region, match_gbif_to_wcvp
)
from scripts.analysis.mois_lgm_grid_spatial_block_20261116 import (
    LGM_ZIP_SHA,CELL_LGM
)
from scripts.analysis.mois_gbif_family_excess_biomass_ndvi_20261116 import (
    signed_parquet,SPATIAL_SHA,SPATIAL_MEMBER,NDVI_SHA,NDVI_MEMBER
)
from scripts.analysis.mois_gbif_family_genus_npp_AGB_NDVI_20261116 import (
    NPP_SHA,NPP_MEMBER,NPP_FIELD
)

REGIONS={
  "North_America_temperate":(25,55,-125,-65),
  "Europe_temperate":(35,65,-15,40),
  "East_Asia_temperate":(25,55,95,145),
}
MAX_SITES=300
MAX_SPECIES=180
MIN_ELIGIBLE_SITES=30
MIN_RECORDINGS=5
MIN_SITES=65
SEED="mois-fourthcorner-20261117-predeclared-no-photo-selection"

def stable_rank(v):
    return hashlib.sha256((SEED+"|"+str(v)).encode("utf-8")).hexdigest()

def read_paired_environment(lgm,spatial,ndvi,npp):
    import pyarrow.parquet as pq
    require(sha(lgm)==LGM_ZIP_SHA,"LGM source signature lost")
    with ZipFile(lgm) as z:
        b=z.read(CELL_LGM)
        receipt=json.loads(z.read("receipt.json"))
    require(sha(b)==receipt["output_SHA256"][CELL_LGM],"LGM cells SHA drift")
    original=pq.read_table(io.BytesIO(b)).to_pandas()
    require(len(original)==18577,"true global climate cells changed")
    spatial=signed_parquet(spatial,SPATIAL_SHA,SPATIAL_MEMBER,"output_sha256")
    ndvi=signed_parquet(ndvi,NDVI_SHA,NDVI_MEMBER,"output_SHA256")
    npp=signed_parquet(npp,NPP_SHA,NPP_MEMBER,"output_SHA256")
    require(len(spatial)==12560 and len(ndvi)==len(npp)==18577,
            "environment site counts not exact")
    cols=["grid_x","grid_y","lat","lon",
          "block_500km","block_1000km",
          "present_temperature_C","present_precip_mm",
          "ESA_AGB_woody_Mg_ha_2020"]
    original=original.merge(spatial[cols],on=["grid_x","grid_y"],
                             how="inner",validate="one_to_one")
    original=original.merge(ndvi[[
        "grid_x","grid_y","MODIS_NDVI_2001_2021_mean",
        "MODIS_NDVI_2001_2021_sd"]],on=["grid_x","grid_y"],
        validate="one_to_one")
    original=original.merge(npp[[
        "grid_x","grid_y",NPP_FIELD,
        "NPP_valid_years_2001_2023"]],on=["grid_x","grid_y"],
        validate="one_to_one")
    required=[
        "threeGCM_delta_temperature_C_median",
        "threeGCM_delta_log_precip_median",
        "present_temperature_C","present_precip_mm",
        "ESA_AGB_woody_Mg_ha_2020",
        "MODIS_NDVI_2001_2021_mean",NPP_FIELD,"lat","lon"]
    good=(original[required].notna().all(axis=1)
        &original.NPP_valid_years_2001_2023.ge(18)
        &original.ESA_AGB_woody_Mg_ha_2020.ge(0)
        &original.present_precip_mm.ge(0))
    return original.loc[good].copy()

def source_species_native_records(pairs,species,names,dist,site):
    """Conservative native species-recording rows with distinct species/cell."""
    if site.empty:return pd.DataFrame()
    species_map=match_gbif_to_wcvp(species,names)
    site=site[["grid_x","grid_y","code"]].drop_duplicates()
    q=pairs.merge(site,on=["grid_x","grid_y"],how="inner",
                  validate="many_to_one")
    q=q.merge(species_map[["key","source_hue9",
             "hue_status","wcvp_accepted_plant_name_id"]],
            left_on="specieskey",right_on="key",validate="many_to_one")
    q=q.merge(dist[["wcvp_accepted_plant_name_id","tdwg_level3_code",
                    "native_introduced_state"]],
          left_on=["wcvp_accepted_plant_name_id","code"],
          right_on=["wcvp_accepted_plant_name_id","tdwg_level3_code"],
          how="left",validate="many_to_one")
    q=q[q.native_introduced_state.eq("native_at_TDWG3")].copy()
    require(not q[["grid_x","grid_y","specieskey"]].duplicated().any(),
            "GBIF native recording duplicated")
    return q[["grid_x","grid_y","specieskey","occurrence_count"]].copy()

def panel_sites(env,regions,box):
    latlo,lathi,lonlo,lonhi=box
    geo=regions[["grid_x","grid_y","code","region_strict"]].copy()
    geo=geo[geo.region_strict]
    x=env.merge(geo[["grid_x","grid_y","code"]],on=["grid_x","grid_y"],
                how="inner",validate="one_to_one")
    x=x[x.lat.between(latlo,lathi)&x.lon.between(lonlo,lonhi)]
    x=x.drop_duplicates(["grid_x","grid_y"])
    x["site_rank"]=[
        stable_rank(f"cell:{g}:{h}") for g,h in zip(x.grid_x,x.grid_y)
    ]
    x=x.sort_values("site_rank").head(MAX_SITES).reset_index(drop=True)
    return x.drop(columns=["site_rank"])

def build_panel(pairs,species,names,dist,sites,region_name):
    if len(sites)<MIN_SITES:
        return None,{"region":region_name,"status":"insufficient_native_climate_sites",
                     "n_eligible_sites":len(sites)}
    natives=source_species_native_records(pairs,species,names,dist,sites)
    if natives.empty:
        return None,{"region":region_name,"status":"no_native_recorded_species",
                     "n_sites":len(sites)}
    freq=natives.groupby("specieskey").agg(
        n_recorded_sites=("grid_x","size"),
        source_occurrences=("occurrence_count","sum")
    ).reset_index()
    # Species selection uses only native recorded range and a stable hash.
    freq=freq[(freq.n_recorded_sites>=MIN_RECORDINGS)
              &(freq.n_recorded_sites<=.9*len(sites))].copy()
    freq["priority"]=freq.specieskey.map(
         lambda k:stable_rank("species:"+str(k)))
    freq=freq.sort_values("priority").head(MAX_SPECIES)
    if len(freq)<60:
        return None,{"region":region_name,"status":"too_few_native_species",
                     "n_sites":len(sites),"n_source_species":len(freq)}
    species_map=match_gbif_to_wcvp(species,names)
    trait=species_map[species_map.key.isin(freq.specieskey)].copy()
    trait=trait[["key","species","source_hue9","hue_status",
                  "wcvp_accepted_plant_name_id"]].sort_values("key")
    require(trait.key.is_unique and len(trait)==len(freq),
            "species trait matching changed")
    # Eligibility means 'listed as native somewhere in the site's broad TDWG3',
    # not true ecological possibility or individual wild status.
    xx=pd.MultiIndex.from_product(
        [sites.index,trait.key.to_list()],
        names=["site_row","key"]).to_frame(index=False)
    xx=xx.merge(sites[["code"]].reset_index().rename(
          columns={"index":"site_row"}),on="site_row",validate="many_to_one")
    xx=xx.merge(trait[["key","wcvp_accepted_plant_name_id"]],
                on="key",validate="many_to_one")
    xx=xx.merge(dist[["wcvp_accepted_plant_name_id",
                       "tdwg_level3_code","native_introduced_state"]],
                left_on=["wcvp_accepted_plant_name_id","code"],
                right_on=["wcvp_accepted_plant_name_id","tdwg_level3_code"],
                how="left",validate="many_to_one")
    xx["eligible_native_range"]=xx.native_introduced_state.eq("native_at_TDWG3")
    sites_index=sites[["grid_x","grid_y"]].reset_index().rename(
        columns={"index":"site_row"})
    present=natives[natives.specieskey.isin(trait.key)].copy()
    present=present.merge(sites_index,on=["grid_x","grid_y"],
                          how="left",validate="many_to_one")
    present["recorded"]=1
    xx=xx.merge(present[["site_row","specieskey","recorded"]],how="left",
                left_on=["site_row","key"],right_on=["site_row","specieskey"],
                validate="one_to_one")
    xx["apparent_recorded"]=xx.recorded.fillna(0).astype("int8")
    require(not xx[(xx.apparent_recorded==1)&(~xx.eligible_native_range)].shape[0],
            "a true native GBIF recording was marked non-native eligible")
    count=xx.groupby("key").agg(
         n_eligible=("eligible_native_range","sum"),
         n_recorded=("apparent_recorded","sum"))
    keep=count.index[
        (count.n_eligible>=MIN_ELIGIBLE_SITES)
        &(count.n_recorded>=MIN_RECORDINGS)
        &(count.n_recorded/count.n_eligible<=.85)
    ]
    trait=trait[trait.key.isin(keep)].copy()
    xx=xx[xx.key.isin(keep)].copy()
    if len(trait)<60:return None,{
         "region":region_name,"status":"too_few_native_eligible_species",
         "n_species":len(trait),"n_sites":len(sites)}
    data=xx.pivot(index="site_row",columns="key",
                  values="apparent_recorded")
    elig=xx.pivot(index="site_row",columns="key",
                  values="eligible_native_range")
    data=data.reindex(index=sites.index,columns=trait.key.to_list())
    elig=elig.reindex(index=sites.index,columns=trait.key.to_list())
    require(data.notna().all().all() and elig.notna().all().all(),
            "matrix join missing native eligibility values")
    Y=data.to_numpy(np.int8)
    eligible=elig.to_numpy(bool)
    require(Y.max()==1 and Y.min()==0
            and not (Y.astype(bool)&~eligible).any(),
            "GBIF observed and WCVP native matrices inconsistent")
    outcome=Y.astype(np.float32)
    outcome[~eligible]=np.nan
    info={"region":region_name,
      "status":"source_authenticated_apparent_GBIF_recorded_matrix",
      "n_sites":len(sites),"n_species":len(trait),
      "n_500km_blocks":int(sites.block_500km.nunique()),
      "n_eligible_native_species_sites":int(eligible.sum()),
      "n_observed_native_species_sites":int(Y.sum()),
      "n_unrecorded_but_WCVP_eligible_species_sites":int(eligible.sum()-Y.sum()),
      "n_unavailable_region_species_sites":int((~eligible).sum()),
      "unrecorded_does_NOT_establish_true_absence":True,
      "site_selection_independent_of_flower_colour":True,
      "species_selection_independent_of_flower_colour":True,
      "native_TDWG3_does_NOT_establish_wild_individual":True,
      "source_photo9_does_NOT_establish_local_petal_colour":True,
      "GBIF_observation_effort_detection_unmodelled":True}
    return (sites,trait,outcome,eligible),info

def run(gbif,env,lgm,spatial,ndvi,npp,wcvp,tdwg,out):
    out=Path(out)
    require(not out.exists(),"refuse overwrite")
    pairs,photo=read_signed_sources(gbif,env)
    environmental=read_paired_environment(lgm,spatial,ndvi,npp)
    # Explicit detection/reporting effort controls are crucial for opportunistic
    # GBIF recording outcomes. These describe record intensity, not survey effort.
    effort=pairs.groupby(["grid_x","grid_y"],sort=False).agg(
        GBIF_all_selected_species_record_count=("occurrence_count","sum"),
        GBIF_all_selected_species_recorded_richness=("specieskey","size")
    ).reset_index()
    require(int(effort.GBIF_all_selected_species_record_count.sum())==396142825
            and int(effort.GBIF_all_selected_species_recorded_richness.sum())==5493931,
            "GBIF recording opportunity denominator drift")
    environmental=environmental.merge(effort,on=["grid_x","grid_y"],
                                      how="left",validate="one_to_one")
    require(environmental.GBIF_all_selected_species_record_count.notna().all(),
            "source effort missing from a GBIF selected cell")
    names,dist,_=get_wcvp_source(wcvp)
    regions=assign_cells_broad_region(
        # Read original complete climate 18577 table, not only the ~12k
        # model-eligible environmental sample.
        read_lgm_cells(lgm),tdwg)
    out.mkdir(parents=True)
    panels={}
    for name,box in REGIONS.items():
        site=panel_sites(environmental,regions,box)
        model,info=build_panel(pairs,photo,names,dist,site,name)
        panels[name]=info
        if model is not None:
            sites,traits,y,elig=model
            directory=out/name
            directory.mkdir(parents=True)
            x=sites[[
               "grid_x","grid_y","lat","lon","code",
               "block_500km","block_1000km",
               "present_temperature_C","present_precip_mm",
               "GBIF_all_selected_species_record_count",
               "GBIF_all_selected_species_recorded_richness",
               "threeGCM_delta_temperature_C_median",
               "threeGCM_delta_log_precip_median",NPP_FIELD,
               "MODIS_NDVI_2001_2021_mean",
               "ESA_AGB_woody_Mg_ha_2020"]].copy()
            x["site_id"]=np.arange(len(x))
            x.to_csv(directory/"X_true_GBIF_100km_environment.csv",index=False)
            tr=traits[["key","species","source_hue9","hue_status"]].copy()
            tr["white_proxy"]=tr.source_hue9.eq("white").astype(int)
            tr.to_csv(directory/"TR_photo9_species_traits.csv",index=False)
            pd.DataFrame(y,columns=traits.key.to_list()).to_csv(
                directory/"Y_apparent_recorded_not_absence.csv",
                index=False,na_rep="NA")
            pd.DataFrame(elig.astype(np.int8),columns=traits.key.to_list()).to_csv(
                directory/"Y_species_native_TDWG3_eligibility.csv",index=False)
        print("MOIS_JSDM_FOURTHCORNER_PANEL_PREFLIGHT",
              json.dumps(info,sort_keys=True),flush=True)
    require(sum(r.get("status","").startswith(
        "source_authenticated") for r in panels.values())>=1,
        "none of the three preregistered regions has a data-ready panel")
    receipt={
     "schema":"MOIS_GBIF_WCVP_photo9_fourthcorner_JSDM_preflight_v1",
     "source_accepted_species":37104,
     "source_GBIF_species_100km_pairs":5493931,
     "source_GBIF_valid_occurrences":396142825,
     "public_NPP_units":"kg C m-2 year-1",
     "native_status_scale":"WCVP broad TDWG3 only",
     "input_sha256":{"GBIF":sha(gbif),"ENV":sha(env),
             "LGM":sha(lgm),"spatial":sha(spatial),
             "NDVI":sha(ndvi),"NPP":sha(npp),"TDWG3":sha(tdwg)},
     "preselected_regions":REGIONS,
     "sampling_not_selected_for_phenotype_result":True,
     "no_genus_effect_or_genus_selection":True,
     "response_is_apparent_recorded_species_not_true_occupancy":True,
     "GBIF_reporting_effort_proxy_in_X":True,
     "GBIF_reporting_effort_proxy_is_not_survey_effort":True,
     "unrecorded_cells_must_not_be_called_true_absence":True,
     "not_fitted_JSDM_yet":True,
     "panels":panels}
    receipt["output_sha256"]={str(p.relative_to(out)):sha(p)
       for p in out.rglob("*.csv")}
    (out/"receipt.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    return receipt

def read_lgm_cells(lgm):
    import pyarrow.parquet as pq
    with ZipFile(lgm) as z:
        src=z.read(CELL_LGM)
        receipt=json.loads(z.read("receipt.json"))
    require(sha(src)==receipt["output_SHA256"][CELL_LGM],
            "source TDWG region lookup cell IDs changed")
    return pq.read_table(io.BytesIO(src)).to_pandas()

if __name__=="__main__":
    p=argparse.ArgumentParser()
    for name in ("gbif","env","lgm","spatial","ndvi","npp",
                 "wcvp","tdwg","out"):
        p.add_argument("--"+name,required=True)
    args=p.parse_args()
    run(args.gbif,args.env,args.lgm,args.spatial,args.ndvi,
        args.npp,args.wcvp,args.tdwg,args.out)
