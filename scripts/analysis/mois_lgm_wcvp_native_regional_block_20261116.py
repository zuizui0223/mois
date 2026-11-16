#!/usr/bin/env python3
"""MOIS signed global GBIF+LGM floral colour: WCVP TDWG3-regional native sensitivity.

'Native' means accepted species has native range status in a broad TDWG3
botanical country; it does NOT establish a wild plant, local individual
nativeness, or island-specific nativeness. We do not call current occupied
GBIF cells an inferred ice-age refugium.

Uses original species×100km GBIF index without treating repeated species cells
as independent species; the cell-level model uses 500km spatial clustering and
block fixed effects. Species photographic modal hue is a fallible label.
"""
from __future__ import annotations
import argparse, hashlib, io, json, re, zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import geopandas as gpd
from pyproj import Transformer

from scripts.analysis.mois_lgm_grid_spatial_block_20261116 import (
    load_inputs, present_climate_cells, fit_cells
)
from scripts.analysis.mois_lgm_climate_history_species_20261116 import (
    coords_in_grid, require, sha
)
from scripts.analysis.mois_real350399_wcvp_geographic_native_status_20261019 import (
    polygons_read
)

NAME_FILE="mois350399_wcvp_source_taxon_exact_name_resolution.csv.gz"
DIST_FILE="mois350399_wcvp_accepted_native_distribution_TDWG3.csv.gz"
MANIFEST="mois350399_wcvp_TDWG3_source_only_taxon_native_receipt.json"
TDWG_GEO_SHA="c172bcf6aba20e19477adc60aebf0023068f0175dca8480f0760b090dcf64840"
BINOMIAL=re.compile(r"^[A-Z][A-Za-z-]+ [a-z][A-Za-z-]+$")
HUES={"white","yellow","orange","red","pink","magenta","purple","blue","bronze"}
VALID_STATES={"native_at_TDWG3","introduced_at_TDWG3","doubtful",
              "extinct","unknown_flag","contradictory_TDWG3_records"}

def get_wcvp_source(folder):
    folder=Path(folder)
    receipt=json.loads((folder/MANIFEST).read_text())
    require(receipt["source_photo_attempts"]==350399
            and receipt["taxa_with_any_source_photos"]==74821
            and receipt["wcvp_distributions_resolution"]=="TDWG_level3_only"
            and receipt["H3_native_island_hypothesis_tested"] is False,
            "original Kew regional/native source receipt failed")
    for fn in (NAME_FILE,DIST_FILE):
        require(sha(folder/fn)==receipt["data_files_sha256"][fn],
                "WCVP 350399-photo source data SHA differs")
    nm=pd.read_csv(folder/NAME_FILE,dtype=str,keep_default_na=False)
    di=pd.read_csv(folder/DIST_FILE,dtype=str,keep_default_na=False)
    require(len(nm)==74821 and nm.inat_taxon_id.is_unique
            and len(di)>=600000 and not di[[
                "wcvp_accepted_plant_name_id","tdwg_level3_code"
            ]].duplicated().any() and di.native_introduced_state.isin(VALID_STATES).all(),
            "WCVP regional accepted-name/distro rows drift")
    return nm,di,receipt

def match_gbif_to_wcvp(photo_species,names):
    """Conservative exact species binomials only, never assume taxonomic IDs map."""
    need=["source_species_exact_binomial","wcvp_match_state",
          "wcvp_accepted_plant_name_id"]
    require(set(need).issubset(names),"WCVP name source columns absent")
    match=names.loc[
        names.wcvp_match_state.eq("unique_exact_name_to_accepted_id")
        &names.wcvp_accepted_plant_name_id.ne("")
        &names.source_species_exact_binomial.str.match(BINOMIAL),need].drop_duplicates()
    # Reject any exact name mapping to multiple distinct WCVP IDs.
    cand=match.groupby("source_species_exact_binomial").wcvp_accepted_plant_name_id.nunique()
    unique_names=cand.index[cand.eq(1)]
    match=match[match.source_species_exact_binomial.isin(unique_names)]
    match=match.drop_duplicates("source_species_exact_binomial")
    species=photo_species[["key","species","source_hue9","hue_status"]].copy()
    species["key"]=species.key.astype(str)
    require(len(species)==37104 and species.key.is_unique,"signed COLXR accepted source changed")
    species=species.merge(match,how="left",left_on="species",
                          right_on="source_species_exact_binomial",
                          validate="many_to_one")
    species["wcvp_match_exact_binomial"]=species.wcvp_accepted_plant_name_id.notna()
    return species

def assign_cells_broad_region(lgm_cells,tdwg_path):
    """Find 25 subpoints/cell; select botanical-country codes only if unique.

    'strict' = >=20/25 samples uniquely in ONE region, none in other regions.
    'moderate' = >=10/25 with same rule, coastal holes remaining unresolved.
    No snapping GBIF coordinates to a convenient botanical country.
    """
    require(sha(tdwg_path)==TDWG_GEO_SHA,
            "Original attributed TDWG level-3 polygon SHA differs")
    polys=polygons_read(tdwg_path)
    cells=lgm_cells[["grid_x","grid_y"]].copy().reset_index(drop=True)
    require(len(cells)==18577 and not cells.duplicated().any(),"global cell count drift")
    lon,lat=coords_in_grid(cells)
    n=lon.shape[1]
    require(lon.shape==(25,n) and lat.shape==(25,n),"expected 25 exact subpoints")
    # A point's matching botanical country MUST be unique. A cell is assigned
    # only if every matched subpoint belongs to the same country.
    spots=pd.DataFrame({
        "subpoint_id":np.arange(n*25,dtype=np.int32),
        "cell_index":np.tile(np.arange(n,dtype=np.int32),25),
        "lon":lon.ravel(),
        "lat":lat.ravel()
    })
    valid=spots.lon.between(-180,180)&spots.lat.between(-90,90)
    geo=gpd.GeoDataFrame(
        spots.loc[valid,["subpoint_id","cell_index"]],
        geometry=gpd.points_from_xy(spots.loc[valid,"lon"],spots.loc[valid,"lat"]),
        crs="EPSG:4326")
    hits=gpd.sjoin(geo,polys[["LEVEL3_COD","geometry"]],
                   how="inner",predicate="within")
    hitcodes=hits[["subpoint_id","cell_index","LEVEL3_COD"]].drop_duplicates()
    perpoint=hitcodes.groupby("subpoint_id").LEVEL3_COD.nunique()
    unique_pt=perpoint.index[perpoint.eq(1)]
    hitcodes=hitcodes[hitcodes.subpoint_id.isin(unique_pt)]
    percell=hitcodes.groupby("cell_index").agg(
        num_points=("subpoint_id","nunique"),
        num_codes=("LEVEL3_COD","nunique"),
        code=("LEVEL3_COD","first"))
    result=cells.copy()
    result["code"]=pd.Series(pd.NA,index=result.index,dtype="string")
    result["tdwg3_points"]=0
    result["tdwg3_distinct_codes"]=0
    inds=percell.index.astype(int)
    result.loc[inds,"tdwg3_points"]=percell.num_points.to_numpy(int)
    result.loc[inds,"tdwg3_distinct_codes"]=percell.num_codes.to_numpy(int)
    ok=(percell.num_codes==1)
    result.loc[percell.index[ok],"code"]=percell.loc[ok,"code"].astype(str).to_numpy()
    result["region_moderate"]=result.code.notna() & result.tdwg3_points.ge(10)
    result["region_strict"]=result.code.notna() & result.tdwg3_points.ge(20)
    print("MOIS_WCVP_TDWG3_CELL_ASSIGNMENT",json.dumps({
        "total_cells":len(result),
        "moderate":int(result.region_moderate.sum()),
        "strict":int(result.region_strict.sum()),
        "ambiguous":int((result.tdwg3_distinct_codes>1).sum()),
        "unique_region_codes":int(result.code.nunique(dropna=True)),
        "source_invalid_polygons_repaired":polys.attrs.get("source_geometry_repair_audit",{})
    },default=str),flush=True)
    return result

def summarize_photo_species_by_tdwg(pairs,species,regions,distribution,mode="strict"):
    require(mode in ("strict","moderate"),"unknown region stringency")
    flag="region_"+mode
    r=regions.loc[regions[flag],["grid_x","grid_y","code"]].copy()
    require(r[["grid_x","grid_y"]].is_unique.all() if hasattr(r[["grid_x","grid_y"]],"is_unique") else True,
            "region cells duplicated")
    p=pairs.merge(species[[
        "key","source_hue9","wcvp_accepted_plant_name_id","wcvp_match_exact_binomial"]],
        left_on="specieskey",right_on="key",how="left",validate="many_to_one")
    p=p.merge(r,on=["grid_x","grid_y"],how="inner",validate="many_to_one")
    require(p.source_hue9.isin(HUES).all() and not p[[
        "specieskey","grid_x","grid_y"]].duplicated().any(),"source hue or grid duplicates")
    dist=distribution[["wcvp_accepted_plant_name_id","tdwg_level3_code",
                       "native_introduced_state"]]
    p=p.merge(dist,how="left",left_on=["wcvp_accepted_plant_name_id","code"],
              right_on=["wcvp_accepted_plant_name_id","tdwg_level3_code"],
              validate="many_to_one")
    native=p.native_introduced_state.eq("native_at_TDWG3")
    introduced=p.native_introduced_state.eq("introduced_at_TDWG3")
    p["photo_white"]=p.source_hue9.eq("white")
    p["native"]=native
    p["introduced"]=introduced
    p["native_white"]=native & p.photo_white
    p["introduced_white"]=introduced & p.photo_white
    summary=p.groupby(["grid_x","grid_y"],sort=False).agg(
        photo_label_species_n=("photo_white","size"),
        white_photo_label_species_n=("photo_white","sum"),
        native_species_n=("native","sum"),
        native_white_n=("native_white","sum"),
        introduced_species_n=("introduced","sum"),
        introduced_white_n=("introduced_white","sum"),
        GBIF_occurrence_count=("occurrence_count","sum")).reset_index()
    summary["unresolved_or_nonmodern_status_n"]=(summary.photo_label_species_n
                    -summary.native_species_n-summary.introduced_species_n)
    audit={"mode":mode,"species_cells_retained":len(p),
           "matched_WCVP_unique_exact_species":int(species.wcvp_match_exact_binomial.sum()),
           "native_at_TDWG3_species_cell_pairs":int(native.sum()),
           "introduced_at_TDWG3_species_cell_pairs":int(introduced.sum()),
           "uncertain_pairs":int((~(native|introduced)).sum()),
           "with_any_native":int(summary.native_species_n.gt(0).sum()),
           "with_atleast30_native":int(summary.native_species_n.ge(30).sum())}
    return summary,audit

def block_model_data(summary,lgm_cells,agb,present):
    d=summary.merge(lgm_cells,on=["grid_x","grid_y"],validate="one_to_one")
    d=d.merge(agb,on=["grid_x","grid_y"],validate="one_to_one")
    d=d.merge(present,on=["grid_x","grid_y"],validate="one_to_one")
    t=Transformer.from_crs("EPSG:6933","EPSG:4326",always_xy=True)
    lon,lat=t.transform((d.grid_x.to_numpy(float)+.5)*100000,
                        (d.grid_y.to_numpy(float)+.5)*100000)
    d["lon"]=lon;d["lat"]=lat;d["abs_lat"]=np.abs(lat)
    d["block_500km"]=d.grid_x.floordiv(5).astype(str)+"|"+d.grid_y.floordiv(5).astype(str)
    d["block_1000km"]=d.grid_x.floordiv(10).astype(str)+"|"+d.grid_y.floordiv(10).astype(str)
    return d

def fit_subset(d,response="native",min_species=30):
    if response=="native":
        chosen=d.copy()
        chosen["photo_label_species_n"]=d.native_species_n
        chosen["white_photo_label_species_n"]=d.native_white_n
    elif response=="introduced":
        chosen=d.copy()
        chosen["photo_label_species_n"]=d.introduced_species_n
        chosen["white_photo_label_species_n"]=d.introduced_white_n
    elif response=="all_same_region_cells":
        chosen=d.copy()
    else:
        raise ValueError("unknown proxy")
    chosen["observed_photo_label_white_fraction"]=(chosen.white_photo_label_species_n
                                        /chosen.photo_label_species_n.replace(0,np.nan))
    fields=["threeGCM_delta_temperature_C_median",
            "threeGCM_delta_log_precip_median",
            "present_temperature_C","present_precip_mm",
            "ESA_AGB_woody_Mg_ha_2020","abs_lat",
            "observed_photo_label_white_fraction"]
    good=(chosen.photo_label_species_n.ge(min_species)
          &chosen[fields].notna().all(axis=1)
          &chosen.threeGCM_valid_all_models
          &chosen.present_precip_mm.ge(0)
          &chosen.ESA_AGB_woody_Mg_ha_2020.ge(0))
    chosen=chosen.loc[good].copy()
    fits={"n_valid_cells":len(chosen)}
    if len(chosen)>=350 and chosen.block_500km.nunique()>=25:
        fits["global_geographic"]=fit_cells(chosen)
        fits["within_1000km"]=fit_cells(chosen,"block_1000km")
        fits["within_500km"]=fit_cells(chosen,"block_500km")
    else:
        fits["status"]="insufficient_cells"
    return fits,chosen

def run(gbif,env,lgm,current,wcvp_dir,tdwg,out):
    out=Path(out)
    require(not out.exists(),"refuse overwrite")
    pairs,photo_species,lgm_cells,agb=load_inputs(gbif,lgm,env,current)
    names,distribution,receipt=get_wcvp_source(wcvp_dir)
    x=match_gbif_to_wcvp(photo_species,names)
    regions=assign_cells_broad_region(lgm_cells,tdwg)
    present=present_climate_cells(lgm_cells,current)
    results={}
    tables=[]
    for mode in ("strict","moderate"):
        s,audit=summarize_photo_species_by_tdwg(
            pairs,x,regions,distribution,mode)
        d=block_model_data(s,lgm_cells,agb,present)
        results[mode]={"region_and_match_audit":audit,"models":{}}
        fits_data={}
        for label in ("native","all_same_region_cells","introduced"):
            fit,used=fit_subset(d,label,30)
            results[mode]["models"][label]=fit
            fits_data[label]=used
        # Paired comparison: global all-hues model on native-support cells,
        # rather than confounding a changed region sample with changing taxa.
        native_cells=fits_data["native"][["grid_x","grid_y"]]
        matched=d.merge(native_cells,on=["grid_x","grid_y"],validate="one_to_one")
        matched_fit,_=fit_subset(matched,"all_same_region_cells",30)
        results[mode]["models"]["all_exactly_native_support_cells"]=matched_fit
        s["tdwg3_assignment_mode"]=mode
        tables.append(s)
        print("MOIS_WCVP_REGIONAL_NATIVE_SPATIAL",mode,json.dumps(
            {k:{"n_cells":v.get("n_valid_cells"),"500km_beta":v.get(
                "within_500km",{}).get("beta_white_fraction_per_SD_LGM_warming"),
                "500km_p":v.get("within_500km",{}).get("p")}
             for k,v in results[mode]["models"].items()},sort_keys=True),flush=True)
    out.mkdir(parents=True)
    pd.concat(tables,ignore_index=True).to_parquet(
        out/"WCVP_TDWG3_native_introduced_photo9_LGM_100km_cell_summary.parquet",
        index=False)
    (out/"native_broad_region_spatial_test.json").write_text(
        json.dumps(results,indent=2,sort_keys=True,default=str)+"\n")
    summary={"schema":"mois_true_GBIF37104_LGM3GCM_WCVP_TDWG3_broad_region_native_v1",
      "source_wcvp_Kew_region_only":True,"individual_is_wild_native_inferred":False,
      "island_native_inferred":False,"LGM_refugium_inferred":False,
      "species_photo_hue_is_verified_petal_colour":False,
      "source_taxa":37104,"source_GBIF_species_cells":5493931,
      "source_GBIF_occurrences":396142825,
      "results":results,
      "source_sha256":{"Kew_names_csv":sha(Path(wcvp_dir)/NAME_FILE),
                       "Kew_distributions_csv":sha(Path(wcvp_dir)/DIST_FILE),
                       "TDWG_GeoJSON":sha(tdwg)},
      "output_sha256":{p.name:sha(p) for p in out.iterdir()}}
    (out/"receipt.json").write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n")
    return summary

if __name__=="__main__":
    p=argparse.ArgumentParser()
    for n in ("gbif","env","lgm","current","wcvp-dir","tdwg","output"):
        p.add_argument("--"+n,required=True)
    a=p.parse_args()
    result=run(a.gbif,a.env,a.lgm,a.current,a.wcvp_dir,a.tdwg,a.output)
    print("MOIS_LGM_WCVP_NATIVE_DATA_SUCCESS",json.dumps({
        "n_source_taxa":result["source_taxa"],
        "strict_native":result["results"]["strict"]["models"]["native"].get("n_valid_cells"),
        "moderate_native":result["results"]["moderate"]["models"]["native"].get("n_valid_cells"),
        "native_on_island_proven":False,
    },sort_keys=True),flush=True)
