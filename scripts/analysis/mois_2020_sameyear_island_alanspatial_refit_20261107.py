#!/usr/bin/env python3
"""Same-year 2020 MOIS photographic colour vs nightlight on new global land targets.

Use all five colour-blind 2020 acquisition increments: 36,070 known-land
target cells, WITHOUT reusing the earlier outcome-selected 1,014 legacy cells.
An original 2020 photographic point must be assigned to its true GSHHG
landmass, not inferred from its cell centre or an island size threshold.
The fitted outcome is species-equal automated ROI9 white mass, not white-
flowered species prevalence. No causal, native-flora or moth-loss inference.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import pandas as pd
import geopandas as gpd
from scripts.analysis.mois_global_2020_GSHHG_scope_20261107 import (
    prepare_land,assign_points)
from scripts.analysis.mois_2024_sameyear_global_island_photo_white_20261107 import (
    cell_id,cluster_fit,COVARS)

PINS={
    "canonical":"f033661d76e07898a5c13199a55c5d0ce6e55353daa48b0058c456979e90644c",
    "wcvp":"dc5bb0e3686fde0cfd327e4336245d0ceae3c7657b3bd12841bb814cc94d1eb2",
    "first256":"303424afe40b41e7bda2ec0ae456dd27820a50e17f7ca23204c81036adf21466",
    "next1024":"7aa0193178d664818bca777749074fc88fe0b9e594e83e5eb34d6f72110efce9",
    "next64":"a513250fa9911a35c7f40380684ae97f52d7cdfea7bf828cdcd7cbe1be34192c",
    "last30630":"174de5a28b62f36e456504062d099b8f1a55f3cc7a2d7fdf6fed003e511446d0"}
COUNTS={"first256":256,"next1024":1024,"next64":64,
        "first4096":4096,"last30630":30630}
COAST_SHA="8dbbe7e071e77e9e75f2d639239099ebca8d5c16d6a07df8169729d49f15cf41"
def check(ok,msg):
    if not bool(ok):raise ValueError("MOIS_2020_SAMEYEAR_ISLAND_ALAN: "+msg)
def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda:f.read(4<<20),b""):h.update(b)
    return h.hexdigest()
def parse_source_ids(df):
    check(df.target_id.is_unique,"duplicate original source target IDs")
    parts=df.target_id.astype(str).str.split("|",regex=False,expand=True)
    check(parts.shape[1]==3 and parts[0].eq("2020").all(),
          "Target ID must identify same-year 2020 original landmass")
    check(parts[1].str.fullmatch(r"\d+_\d+").all(),"invalid source grid identity")
    check(df.land_class.isin(["island","mainland"]).all(),
          "no unresolved landmass may be invented")
    df=df.copy()
    df["cell_id"]=parts[1].values
    df["land_id"]=parts[2].values
    check(df.land_id.str.fullmatch(r"\d+").all(),"not a GSHHG original landmass")
    return df
def load_environment(paths):
    blocks=[]
    for key,n in COUNTS.items():
        file=paths[key]
        if key in PINS:check(sha(file)==PINS[key],"Source file byte SHA drift: "+key)
        else:
            # Reuse the verified original workflow receipt for the 4,096
            # actual data values, whose source target manifest is separately pinned.
            receipt=Path(file).parent/"receipt.json"
            check(receipt.exists(),"Source 4,096 actual original numeric receipt missing")
            r=json.loads(receipt.read_text())
            expected=r.get("data_SHA256") or r.get("output_values_sha256")
            check(expected==sha(file),"Original 4,096 value receipt SHA mismatch")
            check(r.get("input_manifest_SHA256")==
                  "fa06170ce59e41b6d7e486c236acf7d004091a035f7cbddbb3632bcd788ffa65",
                  "Original 4,096 source manifest identity mismatch")
        rows=pd.read_csv(file,dtype={"target_id":str},low_memory=False)
        check(len(rows)==n,"Original increment target count drift: "+key)
        for col in ["target_id","land_class","viirs_valid_months","viirs_median",
                    "population_density","query_longitude","query_latitude"]:
            check(col in rows,"Source lacks required field "+key+" "+col)
        rows["source_increment"]=key
        blocks.append(rows)
    env=pd.concat(blocks,ignore_index=True)
    check(len(env)==36070,"Must use exact 36,070 real new 2020 targets")
    return parse_source_ids(env)
def make_frame(canonical,wcvp,env,shoreline):
    check(len(canonical)==1032882 and canonical.photo_id.is_unique
          and canonical.observation_id.is_unique,"Canonical original photo IDs drifted")
    original=canonical.loc[
        canonical.observed_on.astype(str).str[:4].eq("2020") &
        canonical.latitude.between(-90,90) &
        canonical.longitude.between(-180,180)].copy()
    original["unit_id"]=original.photo_id
    coast=Path(shoreline)
    proof=json.loads((coast/"source_receipt.json").read_text())
    check(proof["original_archive_sha256"]==COAST_SHA,"GSHHG original source changed")
    for filename,expected in proof["extracted_sha256"].items():
        check(sha(coast/filename)==expected,"GSHHG extracted part SHA changed "+filename)
    land=prepare_land(gpd.read_file(coast/"GSHHS_h_L1.shp"))
    lakes=gpd.read_file(coast/"GSHHS_h_L2.shp")
    geo=assign_points(original,land,lakes)
    original=original.merge(geo[["unit_id","land_id","land_class"]],
                            on="unit_id",how="left",validate="one_to_one")
    original=original.loc[original.source_subarm.eq("historical_main") &
                          original.land_class.isin(["island","mainland"])].copy()
    original["cell_id"]=cell_id(original.longitude,original.latitude)
    original["target_id"]="2020|"+original.cell_id+"|"+original.land_id.astype(str)
    selected=env[["target_id","land_class"]].rename(columns={
       "land_class":"environment_land_class"})
    original=original.merge(selected,on="target_id",how="inner",validate="many_to_one")
    check(original.land_class.eq(original.environment_land_class).all(),
          "Source original photos do not match environmental true landmass")
    check(original.photo_id.is_unique and original.observation_id.is_unique,
          "Photographic observations duplicated after coast/target join")
    accepted=wcvp.loc[wcvp.wcvp_match_state.eq(
       "unique_exact_name_to_accepted_id"),
       ["inat_taxon_id","wcvp_accepted_plant_name_id"]].copy()
    check(accepted.inat_taxon_id.is_unique,"WCVP one-to-many source error")
    original=original.merge(accepted,on="inat_taxon_id",how="left",
                            validate="many_to_one")
    per_target=original.groupby("target_id").agg(
       photo_attempts=("photo_id","size"),
       strict9_success=("strict9","sum"),
       bio5=("wc2_bio05_10m","mean"),
       bio12=("wc2_bio12_10m","mean"),
       elevation=("wc2_elev_2p5m","mean")).reset_index()
    good=original.loc[
        original.strict9 & original.strict9_white.notna() &
        original.wcvp_accepted_plant_name_id.notna(),
        ["target_id","wcvp_accepted_plant_name_id","strict9_white"]].copy()
    check(good.strict9_white.between(0,1).all(),"Invalid original colour mass")
    species=good.groupby(["target_id","wcvp_accepted_plant_name_id"],
                         observed=True).strict9_white.mean().reset_index()
    per_cell=species.groupby("target_id",observed=True).agg(
        white_equal=("strict9_white","mean"),
        accepted_species_n=("wcvp_accepted_plant_name_id","nunique")).reset_index()
    d=env.merge(per_target,on="target_id",how="left",validate="one_to_one")
    d=d.merge(per_cell,on="target_id",how="left",validate="one_to_one")
    d["photo_eligible"]=(d.photo_attempts.ge(5)&d.strict9_success.ge(3)&
                          d.accepted_species_n.ge(3))
    d["viirs_eligible"]=(d.viirs_valid_months.ge(3)&d.viirs_median.notna())
    d["pop_eligible"]=d.population_density.notna()
    d["white_pp"]=100*d.white_equal
    d["log_light"]=np.log1p(d.viirs_median.clip(lower=0))
    d["log_pop"]=np.log1p(d.population_density.clip(lower=0))
    d["log_rain"]=np.log1p(d.bio12.clip(lower=0))
    d["log_attempts"]=np.log1p(d.photo_attempts.clip(lower=0))
    d["region5"]=(np.floor((d.query_longitude+180)/5).astype(int).astype(str)+
                  "_"+np.floor((d.query_latitude+90)/5).astype(int).astype(str))
    with_light=d.loc[d.photo_eligible & d.viirs_eligible].dropna(
       subset=["white_pp","log_light","log_rain","bio5","elevation",
               "log_attempts","region5"]).copy()
    with_pop=with_light.loc[with_light.pop_eligible].dropna(
       subset=["log_pop"]).copy()
    with_light["is_island"]=with_light.land_class.eq("island").astype(int)
    with_pop["is_island"]=with_pop.land_class.eq("island").astype(int)
    return original,d,with_light,with_pop
def available_fit(d,form,cluster,*,require_island_support=False):
    if require_island_support:
        island=d.loc[d.land_class.eq("island")]
        n_cells=len(island)
        n_islands=int(island.land_id.nunique())
        # This is a minimum-support check, not a p-value threshold.
        # A tiny selected set of island cells can produce unstable
        # extrapolating interaction coefficients under regional fixed effects.
        if n_cells<30 or n_islands<10:
            return {"estimable":False,"reason":"island_interaction_below_prespecified_support",
                    "minimum_island_cells":30,"minimum_independent_islands":10,
                    "island_cells":n_cells,"independent_islands":n_islands,
                    "cells":len(d),"clusters":int(d[cluster].nunique())}
    if d[cluster].nunique()<5:
        return {"estimable":False,"reason":"fewer_than_five_independent_clusters",
                "clusters":int(d[cluster].nunique()),"cells":len(d)}
    vals=cluster_fit(d,form,cluster)
    return {"estimable":True,"cells":len(d),"clusters":int(d[cluster].nunique()),
            "terms":vals}
def run(paths,shoreline,out):
    out=Path(out);check(not out.exists(),"Do not overwrite prior ecological results")
    for name in ["canonical","wcvp"]:
        check(sha(paths[name])==PINS[name],"Original photo/trait SHA drift "+name)
    canonical=pd.read_csv(paths["canonical"],compression="gzip",low_memory=False,
       dtype={"photo_id":str,"observation_id":str,"inat_taxon_id":str},
       usecols=["photo_id","observation_id","inat_taxon_id","source_subarm",
           "observed_on","latitude","longitude","strict9","strict9_white",
           "wc2_bio05_10m","wc2_bio12_10m","wc2_elev_2p5m"])
    wcvp=pd.read_csv(paths["wcvp"],compression="gzip",dtype=str,
       usecols=["inat_taxon_id","wcvp_match_state","wcvp_accepted_plant_name_id"])
    env=load_environment(paths)
    photos,allcells,b,m=make_frame(canonical,wcvp,env,shoreline)
    global_form="white_pp ~ log_light + log_pop + "+COVARS+" + C(region5)"
    interaction_form="white_pp ~ log_light * is_island + log_pop + "+COVARS+" + C(region5)"
    results={
      "global_pop2020":available_fit(m,global_form,"region5"),
      "global_pop2020_island_light_interaction":available_fit(
            m,interaction_form,"region5",require_island_support=True)}
    island=m.loc[m.land_class.eq("island")].copy()
    info=island.groupby("land_id").agg(n=("target_id","nunique"),
                       light_n=("viirs_median","nunique"))
    ids=set(info.loc[info.n.ge(2)&info.light_n.ge(2)].index)
    island=island.loc[island.land_id.isin(ids)].copy()
    results["within_island_pop2020"]=available_fit(island,
       "white_pp ~ log_light + log_pop + "+COVARS+" + C(land_id)","land_id")
    desc=b.groupby("land_class").agg(
       eligible_cells=("target_id","nunique"),
       mean_photo_white_roi_pp=("white_pp","mean"),
       median_viirs=("viirs_median","median")).reset_index()
    receipt={
       "schema":"MOIS_2020_new36k_sameyear_whiteROI_island_ALAN_v1",
       "original_source_count_new_targets":len(env),
       "original_photo_2020_historical_main_linked":len(photos),
       "original_strict9_measured_linked":int(photos.strict9.sum()),
       "new2020_targets_with_photo_gate_and_VIIRS":len(b),
       "new2020_targets_with_photo_gate_VIIRS_population":len(m),
       "new2020_island_pop_eligible_cells":int(m.is_island.sum()),
       "new2020_within_island_cells":len(island),
       "new2020_within_island_independent_islands":len(ids),
       "models":results,
       "island_area_latitude_filters_used":False,
       "target_selection_colour_blind_original_2020_only":True,
       "source_population_product":"WorldPop Global1 requested for 2020; no response-side year echo",
       "sameyear_VIIRS_is_four_raw_month_medians_not_annual_composite":True,
       "automated_photo_white_ROI_fraction_not_white_species_fraction":True,
       "no_photo_RGB_downloads_or_ROI_recomputation":True,
       "unknown_colours_imputed":False,
       "no_ALAN_causal_native_species_extinction_or_moth_claim":True,
       "source_sha256":{k:sha(v) for k,v in paths.items()}}
    out.mkdir(parents=True)
    desc.to_csv(out/"2020_island_mainland_descriptive_source.csv",index=False)
    (out/"results.json").write_text(json.dumps(receipt,sort_keys=True,indent=2)+"\n")
    print("MOIS_2020_GLOBAL_NEW_SAMEYEAR_ISLAND_ALAN",json.dumps({
        k:receipt[k] for k in ("new2020_targets_with_photo_gate_and_VIIRS",
        "new2020_targets_with_photo_gate_VIIRS_population",
        "new2020_island_pop_eligible_cells",
        "new2020_within_island_independent_islands","models")},sort_keys=True),flush=True)
    return receipt
if __name__=="__main__":
    p=argparse.ArgumentParser()
    for name in ("canonical","wcvp",*COUNTS):
        p.add_argument("--"+name.replace("_","-"),dest=name,type=Path,required=True)
    p.add_argument("--shoreline",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    run({k:getattr(a,k) for k in ("canonical","wcvp",*COUNTS)},a.shoreline,a.out)
