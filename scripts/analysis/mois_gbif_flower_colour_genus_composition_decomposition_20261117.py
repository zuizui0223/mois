#!/usr/bin/env python3
"""Do regional NPP/historic climate links to WHITE photo-hue species composition
arise from which genera occupy current GBIF grids, or species turnover WITHIN
genera? A true current-GBIF-cell, WCVP regional native/introduced decomposition.

For genus g, colour propensity p_g = fraction white modal-photo9 species in
signed global accepted-colour catalog (one species once, not photo count).
For region/cell assemblage C with n_g species of genus g:
  observed_white_C = sum white_i / n_C
  expected_due_to_genus_pool_C = sum n_g p_g / n_C
  residual_within_genus_species_selection_C = observed - expected
The last term describes current GBIF composition, NOT genetic colour change
or individual phenotype within species. Current occupied site ≠ paleo refugium.
All three responses are fitted on identical cells with identical covariates,
allowing exact slope reconstruction and spatial 500km clustering.
"""
from __future__ import annotations
import argparse,io,json
from pathlib import Path
from zipfile import ZipFile
import numpy as np
import pandas as pd

from scripts.analysis.mois_lgm_climate_history_species_20261116 import (
    read_signed_sources,require,sha,HUES
)
from scripts.analysis.mois_lgm_wcvp_native_regional_block_20261116 import (
    get_wcvp_source,assign_cells_broad_region
)
from scripts.analysis.mois_lgm_native_intro_common_genus_20261116 import (
    regional_species_cell_rows
)
from scripts.analysis.mois_gbif_family_excess_biomass_ndvi_20261116 import (
    signed_parquet, NATIVE_SHA, SPATIAL_SHA, NDVI_SHA,
    NATIVE_MEMBER, SPATIAL_MEMBER, NDVI_MEMBER
)
from scripts.analysis.mois_gbif_family_genus_npp_AGB_NDVI_20261116 import (
    NPP_SHA,NPP_MEMBER,NPP_FIELD
)
from scripts.analysis.mois_flower_colour_LGM_NPP_GBIF_20261117 import (
    cell_white_model,cell_cohort
)

LGM_SHA="db3b5b4fa2133194963be9233411bb6bb50bc86393d1dc50d996cc4cc78b25ab"
LGM_CELL="WorldClim1p4_LGM_threeGCM_18577_cells.parquet"
TARGET=["native","introduced"]

def signed_inputs(gbif,env,lgm,native,spatial,ndvi,npp):
    import pyarrow.parquet as pq
    pairs,species=read_signed_sources(gbif,env)
    require(sha(lgm)==LGM_SHA,"LGM GBIF 100km source ZIP drift")
    with ZipFile(lgm) as z:
        rec=json.loads(z.read("receipt.json"))
        b=z.read(LGM_CELL)
    require(sha(b)==rec["output_SHA256"][LGM_CELL]
            and rec["species_x_grid_pairs"]==5493931,
            "LGM source cell columns changed")
    climate=pq.read_table(io.BytesIO(b)).to_pandas()
    n=signed_parquet(native,NATIVE_SHA,NATIVE_MEMBER,"output_sha256")
    s=signed_parquet(spatial,SPATIAL_SHA,SPATIAL_MEMBER,"output_sha256")
    v=signed_parquet(ndvi,NDVI_SHA,NDVI_MEMBER,"output_SHA256")
    p=signed_parquet(npp,NPP_SHA,NPP_MEMBER,"output_SHA256")
    require(len(climate)==18577 and len(s)==12560
            and len(p)==len(v)==18577,"source present day cells changed")
    keep=["grid_x","grid_y","block_500km","block_1000km",
          "lat","lon","present_temperature_C","present_precip_mm",
          "threeGCM_delta_temperature_C_median",
          "threeGCM_delta_log_precip_median",
          "ESA_AGB_woody_Mg_ha_2020"]
    region=n.merge(s[keep],on=["grid_x","grid_y"],how="inner",
                   validate="many_to_one")
    region=region.merge(v[["grid_x","grid_y","MODIS_NDVI_2001_2021_mean"]],
                        on=["grid_x","grid_y"],validate="many_to_one")
    region=region.merge(p[["grid_x","grid_y",NPP_FIELD,
            "NPP_valid_years_2001_2023"]],
                 on=["grid_x","grid_y"],validate="many_to_one")
    require(not region[["grid_x","grid_y","tdwg3_assignment_mode"]].duplicated().any(),
            "region source not unique")
    return pairs,species,climate,region

def genus_colour_reference(species):
    require(len(species)==37104 and species.key.is_unique,
            "missing global photographic colour species catalogue")
    require(species.source_hue9.isin(HUES).all(),"source hue not in known nine")
    d=species[["key","genus","source_hue9"]].copy()
    d["white"]=d.source_hue9.eq("white").astype(int)
    g=d.groupby("genus",sort=False).agg(
        source_catalogue_species_n=("key","size"),
        source_catalogue_white_n=("white","sum")).reset_index()
    g["genus_global_white_rate"]=(g.source_catalogue_white_n
                                   /g.source_catalogue_species_n)
    require(g.genus.is_unique and np.isfinite(g.genus_global_white_rate).all(),
            "signed original 37k photo-genera not unique")
    return g

def local_genus_expected(source,ref):
    require(not source[["specieskey","grid_x","grid_y"]].duplicated().any(),
            "repeat species within one GBIF cell")
    s=source[["grid_x","grid_y","genus","status","white"]].merge(
        ref[["genus","genus_global_white_rate","source_catalogue_species_n"]],
        on="genus",how="left",validate="many_to_one")
    require(s.genus_global_white_rate.notna().all(),
            "not all GBIF source genus names in photographic colour census")
    s["white_source"]=s.white.astype(int)
    local=s.groupby(["grid_x","grid_y","status"],sort=False).agg(
        n_status_source_species=("white_source","size"),
        white_source_species=("white_source","sum"),
        genus_predicted_white=("genus_global_white_rate","sum"),
        n_mono_colour_genus_species=(
            "source_catalogue_species_n",lambda x:int((x==1).sum()))
    ).reset_index()
    local["observed_white"]=local.white_source_species/local.n_status_source_species
    local["expected_white_from_genus_composition"]=(
        local.genus_predicted_white/local.n_status_source_species)
    local["within_genus_colour_residual"]=(
        local.observed_white-local.expected_white_from_genus_composition)
    require(np.allclose(local.observed_white,
        local.expected_white_from_genus_composition
         +local.within_genus_colour_residual,atol=1e-12,rtol=0),
        "true GBIF within-genus colour residual fails decomposition")
    return local

def add_decomposition_to_signed_regional(region,local,mode,kind,minimum=30):
    status=("native_at_TDWG3" if kind=="native" else "introduced_at_TDWG3")
    z=region[region.tdwg3_assignment_mode.eq(mode)].copy()
    n,w=(("native_species_n","native_white_n") if kind=="native"
         else ("introduced_species_n","introduced_white_n"))
    z=z.merge(local[local.status.eq(status)].drop(columns=["status"]),
              on=["grid_x","grid_y"],how="left",validate="one_to_one")
    # Matching group totals MUST match the original independently created
    # signed WCVP native/introduced photo-hue GBIF cell aggregate.
    available=z.n_status_source_species.notna()
    require((z.loc[available,n].astype(int).to_numpy()
             ==z.loc[available,"n_status_source_species"].astype(int).to_numpy()).all(),
            "original 37k GBIF WCVP group species counts changed")
    require((z.loc[available,w].astype(int).to_numpy()
             ==z.loc[available,"white_source_species"].astype(int).to_numpy()).all(),
            "source white classification or WCVP native/intro drift")
    z["colour_species_n"]=z[n]
    fields=["observed_white","within_genus_colour_residual",
            "expected_white_from_genus_composition",
            NPP_FIELD,"threeGCM_delta_temperature_C_median",
            "threeGCM_delta_log_precip_median","present_temperature_C",
            "present_precip_mm","ESA_AGB_woody_Mg_ha_2020",
            "MODIS_NDVI_2001_2021_mean","GBIF_occurrence_count","lat"]
    good=z[fields].notna().all(axis=1)
    good &= z.colour_species_n.ge(minimum)
    good &= z.NPP_valid_years_2001_2023.ge(18)
    good &= z.ESA_AGB_woody_Mg_ha_2020.ge(0) & z.present_precip_mm.ge(0)
    z=z.loc[good].copy()
    require(len(z)>2500 if kind=="native" else len(z)>1200,
            "too few geographic sample cells matched by flower colour genus")
    return z

def matched_spatial_slopes(z,block="500",weighted=False):
    out={}
    for label,field in (
        ("raw_observed_white","observed_white"),
        ("expected_from_genus_pool","expected_white_from_genus_composition"),
        ("residual_within_genus","within_genus_colour_residual")
    ):
        d=z.copy()
        d["white_photo_species_fraction"]=d[field].astype(float)
        out[label]=cell_white_model(d,block=block,weight=weighted,npp=True)
    add={}
    for key in ("NPP","LGM_warming","LGM_rain",
                "log_woody_AGB","NDVI"):
        a=out["raw_observed_white"]["effects"][key]["beta_fraction_per_SD"]
        b=out["expected_from_genus_pool"]["effects"][key]["beta_fraction_per_SD"]
        c=out["residual_within_genus"]["effects"][key]["beta_fraction_per_SD"]
        err=a-b-c
        require(abs(err)<1e-10,
                f"{block} {key}: genus decomposition slope did not reconstruct")
        add[key]=float(err)
    return {"models":out,"additivity_max_abs":max(map(abs,add.values())),
            "additivity_error":add}

def analyze(pairs,species,climate,region,names,dist,tdwg_path):
    ref=genus_colour_reference(species)
    from scripts.analysis.mois_lgm_wcvp_native_regional_block_20261116 import (
        match_gbif_to_wcvp
    )
    region_keys=region[["grid_x","grid_y"]].drop_duplicates()
    geo=assign_cells_broad_region(climate,tdwg_path)
    source=regional_species_cell_rows(
        pairs,species,names,dist,geo,region_keys)
    local=local_genus_expected(source,ref)
    results={};output=[]
    for mode in ("strict","moderate"):
        results[mode]={}
        for kind in TARGET:
            d=add_decomposition_to_signed_regional(region,local,mode,kind)
            results[mode][kind]={
                "eligible_cells":len(d),
                "mean_observed_photo_white_fraction":float(d.observed_white.mean()),
                "mean_expected_due_to_genus":float(
                    d.expected_white_from_genus_composition.mean()),
                "mean_within_genus_species_composition_residual":float(
                    d.within_genus_colour_residual.mean()),
                "500km_equal_cell":matched_spatial_slopes(d),
                "1000km_equal_cell":matched_spatial_slopes(d,"1000"),
                "500km_species_support_weighted":matched_spatial_slopes(
                    d,weighted=True)
            }
            t=d[["grid_x","grid_y","observed_white",
                "expected_white_from_genus_composition",
                "within_genus_colour_residual",NPP_FIELD,
                "threeGCM_delta_temperature_C_median","block_500km",
                "n_status_source_species","colour_species_n"]].copy()
            t["WCVP_geo_mode"]=mode
            t["WCVP_status"]=kind
            output.append(t)
            models=results[mode][kind]["500km_equal_cell"]["models"]
            print("MOIS_GBIF_WHITE_COLOUR_GENUS_DECOMPOSITION",json.dumps({
                "mode":mode,"group":kind,"cells":len(d),
                "raw_NPP_beta":models["raw_observed_white"]["effects"][
                    "NPP"]["beta_fraction_per_SD"],
                "expected_NPP_beta":models["expected_from_genus_pool"]["effects"][
                    "NPP"]["beta_fraction_per_SD"],
                "within_genus_NPP_beta":models["residual_within_genus"]["effects"][
                    "NPP"]["beta_fraction_per_SD"],
                "within_genus_NPP_p":models["residual_within_genus"]["effects"][
                    "NPP"]["p"],
                "within_genus_LGM_beta":models["residual_within_genus"]["effects"][
                    "LGM_warming"]["beta_fraction_per_SD"],
                "within_genus_LGM_p":models["residual_within_genus"]["effects"][
                    "LGM_warming"]["p"],
                "additivity":results[mode][kind]["500km_equal_cell"][
                    "additivity_max_abs"]},sort_keys=True),flush=True)
    return results,pd.concat(output,ignore_index=True)

def run(gbif,env,lgm,native,spatial,ndvi,npp,wcvp,tdwg,out):
    out=Path(out)
    require(not out.exists(),"refuse overwrite biological colour results")
    pairs,species,climate,region=signed_inputs(
        gbif,env,lgm,native,spatial,ndvi,npp)
    names,dist,_=get_wcvp_source(wcvp)
    r,table=analyze(pairs,species,climate,region,names,dist,tdwg)
    out.mkdir(parents=True)
    f=out/"GBIF_native_introduced_genus_expected_photo_white_vs_NPP_LGM.parquet"
    table.to_parquet(f,index=False)
    mod=out/"models.json"
    mod.write_text(json.dumps(r,indent=2,sort_keys=True)+"\n")
    receipt={"study":"MOIS_photo9_species_white_current_GBIF_WCVP_native_vs_introduced",
        "true_GBIF_species":37104,"true_GBIF_species_cell_pairs":5493931,
        "flower_colour_outcome_photo9_species_modal_proxy":True,
        "genus_expected_uses_species_once_not_photo_counts":True,
        "LGM_paleoclimate_does_not_observe_refugium":True,
        "native_status_broad_region_not_local_wild":True,
        "NPP_is_annual_flux_not_floral_biomass":True,
        "source_zip_sha256":{"GBIF":sha(gbif),"ENV":sha(env),
                "LGM":sha(lgm),"NATIVE":sha(native),
                "SPATIAL":sha(spatial),"NDVI":sha(ndvi),"NPP":sha(npp),
                "TDWG3_geo":sha(tdwg)},
        "result":r,
        "output_SHA256":{f.name:sha(f),mod.name:sha(mod)}}
    (out/"receipt.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    return receipt

if __name__=="__main__":
    p=argparse.ArgumentParser()
    for k in ("gbif","env","lgm","native","spatial","ndvi","npp",
              "wcvp","tdwg","out"):
        p.add_argument("--"+k,required=True)
    a=p.parse_args()
    run(a.gbif,a.env,a.lgm,a.native,a.spatial,a.ndvi,a.npp,a.wcvp,a.tdwg,a.out)
