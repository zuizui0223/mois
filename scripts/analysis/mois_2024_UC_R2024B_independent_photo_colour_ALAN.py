#!/usr/bin/env python3
"""2024 source-aligned MOIS white ROI–ALAN association with independent UC population.

Unconstrained WorldPop R2024B is a SEPARATE full-cohort source: no filling
missing R2025A constrained sites with estimates from a different model.
No source RGB or ROI9 reinference, no true nocturnal pollination labels.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import numpy as np,pandas as pd
from scripts.analysis import mois_2024_sameyear_global_island_photo_white_20261107 as old

INPUT_SHA={
 "canonical":"f033661d76e07898a5c13199a55c5d0ce6e55353daa48b0058c456979e90644c",
 "wcvp":"dc5bb0e3686fde0cfd327e4336245d0ceae3c7657b3bd12841bb814cc94d1eb2",
 "photo_landmass":"f7f65e6bc167089b3034555c655e777f528ecdd29952ab7ee9833d8722f39267",
 "CN_2024_environment":"6e0a077d496e90848d05d043bffae57e4772519f4572a457221c5ef4da0fb813"}
FORMULA_COVARS="bio5 + log_rain + elevation + log_attempts"

def require(ok,msg):
    if not bool(ok):raise ValueError("MOIS_UC_PHOTO_LIGHT: "+msg)

def digest(p):
    h=hashlib.sha256()
    with Path(p).open("rb") as f:
        for b in iter(lambda:f.read(4<<20),b""):h.update(b)
    return h.hexdigest()

def validate_sources(cn,uc,receipt):
    require(len(cn)==len(uc)==66456 and cn.target_id.is_unique and uc.target_id.is_unique,
            "photo-landmass source target IDs are duplicated")
    require(set(cn.target_id)==set(uc.target_id),"not same source 2024 target IDs")
    require(receipt.get("source_releases_NOT_mixed_in_one_covariate") is True,
            "population release cannot be mixed with CN input")
    require(receipt.get("original_2024_constrained_2024_values_SHA256")==INPUT_SHA["CN_2024_environment"],
            "unexpected source-year CN receipt")
    cols=["viirs_median","viirs_valid_months","land_class","n_attempted","n_strict9"]
    joined=cn[["target_id",*cols]].merge(uc[["target_id",*cols]],on="target_id",
                                      suffixes=("_CN","_UC"),validate="one_to_one")
    for col in cols:
        a,b=joined[col+"_CN"],joined[col+"_UC"]
        if col=="land_class":require(a.equals(b),"land status drift")
        else:require(np.allclose(pd.to_numeric(a),pd.to_numeric(b),equal_nan=True,rtol=1e-10,atol=1e-7),
                     "original same-year source changed: "+col)
    require(int(cn.complete.sum())==43595 and int(cn.population_source_status.eq("OK").sum())==44667,
            "original constrained population complete count changed")
    require((uc.UC_source_complete==(uc.viirs_valid_months.ge(3)&uc.population_UC_R2024B_status.eq("OK"))).all(),
            "UC source incorrectly promotes missing population to zero")
    require(uc.loc[uc.population_UC_R2024B_status.ne("OK"),
             "population_UC_R2024B_density_per_km2"].isna().all(),
            "Technical population NoData has a numerical value")
    require(uc.loc[uc.population_UC_R2024B_status.eq("OK"),
             "population_UC_R2024B_density_per_km2"].ge(0).all(),
            "Negative alternative population density")
    return True

def run(canonical,wcvp,photo_landmass,CN_environment,UC_environment,UC_receipt,out):
    out=Path(out);require(not out.exists(),"Do not overwrite confirmed ecology results")
    for label,path in (("canonical",canonical),("wcvp",wcvp),
                       ("photo_landmass",photo_landmass),("CN_2024_environment",CN_environment)):
        require(digest(path)==INPUT_SHA[label],"SHA-frozen original source changed: "+label)
    uc_receipt=json.loads(Path(UC_receipt).read_text())
    require(digest(UC_environment)==uc_receipt["output_SHA256"][Path(UC_environment).name],
            "UC original values do not match full GeoTIFF receipt")
    require(uc_receipt["UC_original_GeoTIFF"]["bytes"]==868021365 and
            uc_receipt["source_cohort"]==66456,
            "R2024B official source file size or cohort changed")
    photo=pd.read_csv(canonical,compression="gzip",low_memory=False,
       dtype={"photo_id":str,"observation_id":str,"inat_taxon_id":str},
       usecols=["photo_id","observation_id","inat_taxon_id","source_subarm",
                "observed_on","latitude","longitude","strict9","strict9_white",
                "wc2_bio05_10m","wc2_bio12_10m","wc2_elev_2p5m"])
    taxon=pd.read_csv(wcvp,compression="gzip",dtype=str,
       usecols=["inat_taxon_id","wcvp_match_state","wcvp_accepted_plant_name_id"])
    land=pd.read_csv(photo_landmass,low_memory=False,
       dtype={"target_id":str,"cell_id":str,"land_id":str})
    cn=pd.read_csv(CN_environment,dtype={"target_id":str},low_memory=False)
    uc=pd.read_csv(UC_environment,dtype={"target_id":str},low_memory=False)
    validate_sources(cn,uc,uc_receipt)
    original,_,b_old,m_old,_,_=old.make_frames(photo,taxon,land,cn)
    require(len(original)==105441 and len(b_old)==2091 and len(m_old)==1617,
            "confirmed original 2024 ROI9 photographic source model changed")
    # A separate full-cohort alternative population series, never a sitewise
    # hybrid of constrained CN and unconstrained UC:
    e=cn.merge(uc[["target_id","population_UC_R2024B_status",
       "population_UC_R2024B_density_per_km2","UC_source_complete"]],
       on="target_id",validate="one_to_one")
    e["population_source_status"]=e.pop("population_UC_R2024B_status")
    e["population_density"]=e.pop("population_UC_R2024B_density_per_km2")
    # original make_frames checks the CN source's old complete count for identity,
    # but eligibility is independently defined by the new UC pop fields.
    # The 'complete' column is ONLY for original source identity validation;
    # all new models consume data.pop_eligible created from UC.
    _,_,b,m,_,_=old.make_frames(photo,taxon,land,e)
    require(len(b)==len(b_old) and set(b.target_id)==set(b_old.target_id),
            "colour/VIIRS cohort changed while only population source should change")
    require(m.target_id.is_unique and m.pop_eligible.all(),"duplicate/invalid UC photographic populations")
    models=[]
    specs=[
      ("UC_R2024B_global_region_adjusted",m,
       "white_pp ~ log_light + log_pop + "+FORMULA_COVARS+" + C(region5)","region5"),
      ("UC_R2024B_island_mainland_light_interaction",m,
       "white_pp ~ log_light * is_island + log_pop + "+FORMULA_COVARS+" + C(region5)","region5")]
    for name,frame,formula,cluster in specs:
        for term,result in old.cluster_fit(frame,formula,cluster).items():
            models.append({"model":name,"term":term,"cells":len(frame),
                   "clusters":int(frame[cluster].nunique()),**result})
    island=m.loc[m.land_class.eq("island")].copy()
    island_group=island.groupby("land_id").agg(cells=("target_id","size"),
                                              lights=("viirs_median","nunique"))
    ids=set(island_group.loc[island_group.cells.ge(2)&island_group.lights.ge(2)].index)
    same=island.loc[island.land_id.isin(ids)].copy()
    require(len(ids)>=5 and len(same)>=30,"too few independent supported islands")
    f="white_pp ~ log_light + log_pop + "+FORMULA_COVARS+" + C(land_id)"
    stats=old.cluster_fit(same,f,"land_id")["log_light"]
    models.append({"model":"UC_R2024B_within_same_island_FE","term":"log_light",
           "cells":len(same),"clusters":len(ids),**stats})
    # Keep original and new source in separate named rows, not a mixed model.
    diag=pd.DataFrame([
      {"source":"CN_R2025A_2024_original","all_photo_VIIRS_cells":len(b_old),
       "population_eligible_photo_cells":len(m_old),
       "within_island_supported_cells":142,"within_island_supported_islands":17},
      {"source":"UC_R2024B_2024_alternative","all_photo_VIIRS_cells":len(b),
       "population_eligible_photo_cells":len(m),
       "within_island_supported_cells":len(same),
       "within_island_supported_islands":len(ids)}])
    out.mkdir(parents=True)
    model_csv=out/"2024_UC_R2024B_only_photo_white_nightlight_models.csv"
    support_csv=out/"CN_R2025A_versus_UC_R2024B_photo_source_support.csv"
    pd.DataFrame(models).to_csv(model_csv,index=False)
    diag.to_csv(support_csv,index=False)
    receipt={"schema":"MOIS_2024_original_ROI9_white_VIIRS_WorldPop_UC_separate_source_v1",
       "input_SHA256":INPUT_SHA,"new_UC_pop_source_values_SHA256":digest(UC_environment),
       "UC_GeoTIFF_SHA256":uc_receipt["UC_original_GeoTIFF"]["file_sha256"],
       "source_release_full_UC_R2024B_no_population_splicing":True,
       "original_2024_photo_VIIRS_eligible_cells":len(b),
       "previous_CN_photo_population_eligible_cells":len(m_old),
       "new_UC_photo_population_eligible_cells":len(m),
       "new_UC_supported_island_cells":len(same),
       "new_UC_independent_islands":len(ids),
       "new_UC_island_beta_pp_per_loglight":stats,
       "associations_exploratory_not_causal":True,
       "photographic_white_not_biologically_verified_species_white":True,
       "moth_effectiveness_not_measured":True,
       "new_original_photo_downloads_or_ROI_inferences":0,
       "no_MNAR_tipping_scenarios_run":True,
       "output_SHA256":{p.name:digest(p) for p in (model_csv,support_csv)}}
    (out/"2024_UC_photo_white_source_receipt.json").write_text(
       json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print("MOIS_2024_UC_ALAN_WHITE_COLOUR_SOURCE_ONLY",
       json.dumps({k:receipt[k] for k in ("original_2024_photo_VIIRS_eligible_cells",
       "previous_CN_photo_population_eligible_cells","new_UC_photo_population_eligible_cells",
       "new_UC_supported_island_cells","new_UC_independent_islands",
       "new_UC_island_beta_pp_per_loglight")},sort_keys=True),flush=True)
    return receipt

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for arg in ("canonical","wcvp","photo_landmass","CN_environment",
                "UC_environment","UC_receipt","out"):
        p.add_argument("--"+arg.replace("_","-"),dest=arg,required=True,type=Path)
    a=p.parse_args()
    run(a.canonical,a.wcvp,a.photo_landmass,a.CN_environment,a.UC_environment,
        a.UC_receipt,a.out)
