#!/usr/bin/env python3
"""Source-audited GBIF native/introduced WHITE PHOTO-HUE at same cell, same genus.

Does broad TDWG3 native/introduced composition alone account for the paired
LGM warming association? All outcomes use species-level photographic modal
palette labels, NOT GBIF field-measured flowers. Current presence != LGM range,
and regional native != individual wild or island-native.
"""
from __future__ import annotations
import argparse,io,json
from pathlib import Path
from zipfile import ZipFile
import numpy as np
import pandas as pd

from scripts.analysis.mois_lgm_climate_history_species_20261116 import (
    read_signed_sources,require,sha
)
from scripts.analysis.mois_lgm_wcvp_native_regional_block_20261116 import (
    get_wcvp_source,match_gbif_to_wcvp,assign_cells_broad_region
)
from scripts.analysis.mois_lgm_native_introduced_matched_cells_20261116 import (
    signed_tables,same_cell_cohort,fit
)
LGM_SHA="db3b5b4fa2133194963be9233411bb6bb50bc86393d1dc50d996cc4cc78b25ab"
LGM_CELL="WorldClim1p4_LGM_threeGCM_18577_cells.parquet"

def regional_species_cell_rows(pairs,photo_species,names,distribution,regions,matched):
    """Only GBIF species presences in cells with both source-defined guilds."""
    ids=matched[["grid_x","grid_y"]].drop_duplicates()
    ids=ids.merge(regions[["grid_x","grid_y","code","region_strict","region_moderate"]],
                  on=["grid_x","grid_y"],validate="one_to_one")
    require(ids.region_moderate.all(),"matched cell lacked moderate TDWG3 assignment")
    s=match_gbif_to_wcvp(photo_species,names)
    s["genus"]=s.species.str.split().str[0]
    p=pairs.merge(ids,on=["grid_x","grid_y"],how="inner",validate="many_to_one")
    p=p.merge(s[["key","genus","source_hue9","wcvp_accepted_plant_name_id"]],
              left_on="specieskey",right_on="key",how="left",validate="many_to_one")
    dist=distribution[["wcvp_accepted_plant_name_id","tdwg_level3_code",
                       "native_introduced_state"]]
    p=p.merge(dist,left_on=["wcvp_accepted_plant_name_id","code"],
              right_on=["wcvp_accepted_plant_name_id","tdwg_level3_code"],
              how="left",validate="many_to_one")
    require(not p[["specieskey","grid_x","grid_y"]].duplicated().any(),
            "a plant species repeated in a single 100km GBIF cell")
    p["white"]=p.source_hue9.eq("white").astype("int8")
    p["status"]=p.native_introduced_state
    p=p.loc[p.status.isin(["native_at_TDWG3","introduced_at_TDWG3"]),
             ["grid_x","grid_y","region_strict","region_moderate",
              "specieskey","genus","white","status"]]
    require(not p.empty,"native/introduced species source rows absent")
    return p

def shared_genus_cells(source):
    """One record per geographic cell after within-cell genus matching.

    Per genus, both statuses must be present. Equal genus means each common
    genus contributes 1 comparison, regardless of size. Shared-species ratio
    compares only species in genera present on BOTH sides (but still can differ
    in their genus weights). All-genus source counts remain separate.
    """
    group=source.groupby(["grid_x","grid_y","genus","status"],observed=True,
                         sort=False).agg(n=("white","size"),
                                         white=("white","sum")).reset_index()
    n=group.loc[group.status.eq("native_at_TDWG3"),
                ["grid_x","grid_y","genus","n","white"]].rename(
        columns={"n":"native_n","white":"native_white"})
    i=group.loc[group.status.eq("introduced_at_TDWG3"),
                ["grid_x","grid_y","genus","n","white"]].rename(
        columns={"n":"intro_n","white":"intro_white"})
    paired=n.merge(i,on=["grid_x","grid_y","genus"],how="inner",
                   validate="one_to_one")
    paired["genus_diff"]=(paired.native_white/paired.native_n
                          -paired.intro_white/paired.intro_n)
    paired["min_support"]=paired[["native_n","intro_n"]].min(axis=1)
    paired["weighted_genus_diff"]=paired.genus_diff*paired.min_support
    cell=paired.groupby(["grid_x","grid_y"],sort=False).agg(
        n_shared_genera=("genus","size"),
        shared_native_species_n=("native_n","sum"),
        shared_intro_species_n=("intro_n","sum"),
        shared_native_white_n=("native_white","sum"),
        shared_intro_white_n=("intro_white","sum"),
        shared_genus_equal_sum=("genus_diff","sum"),
        weighted_genus_difference_sum=("weighted_genus_diff","sum"),
        weighted_genus_support_sum=("min_support","sum"),
    ).reset_index()
    cell["genus_equal_native_minus_intro"]=(cell.shared_genus_equal_sum
                                            /cell.n_shared_genera)
    cell["matched_genera_support_weighted_difference"]=(
        cell.weighted_genus_difference_sum/cell.weighted_genus_support_sum)
    cell["matched_genera_species_pooled_difference"]=(
        cell.shared_native_white_n/cell.shared_native_species_n
        -cell.shared_intro_white_n/cell.shared_intro_species_n)
    return cell,paired

def exactly_same_cell_reference(regional,signed_native,mode):
    """Groupby counts must reproduce previously frozen native/introduced totals."""
    q=regional[regional["region_"+mode]].copy()
    v=q.groupby(["grid_x","grid_y","status"],observed=True).agg(
        n=("white","size"),white=("white","sum")).reset_index()
    tab=v.pivot(index=["grid_x","grid_y"],columns="status",
                values=["n","white"]).fillna(0)
    tab.columns=["_".join(x) for x in tab.columns]
    tab=tab.reset_index()
    current=signed_native[
        signed_native.tdwg3_assignment_mode.eq(mode)].copy()
    z=current.merge(tab,on=["grid_x","grid_y"],how="inner",
                    validate="one_to_one")
    require(len(z)>500,"too few source native cells to audit")
    for key,new in (
        ("native_species_n","n_native_at_TDWG3"),
        ("introduced_species_n","n_introduced_at_TDWG3"),
        ("native_white_n","white_native_at_TDWG3"),
        ("introduced_white_n","white_introduced_at_TDWG3"),
    ):
        require(np.array_equal(z[key].to_numpy(int),z[new].to_numpy(int)),
                f"species-cell source crosswalk mismatch: {key}")
    return {"n_signed_region_cells_with_both_min30":int(len(z)),
            "exact_species_native_introduced_colour_counts_verified":True}

def fit_one_cohort(data,label):
    if len(data)<300:
        return {"status":"insufficient_common_genus_cell_support","n_cells":len(data)}
    require(data.n_shared_genera.ge(1).all(),"empty shared genera")
    results={"n_cell_support":len(data),
             "n_500km_blocks":int(data.block_500km.nunique()),
             "median_shared_genera":float(data.n_shared_genera.median()),
             "median_shared_native_species":float(data.shared_native_species_n.median()),
             "median_shared_intro_species":float(data.shared_intro_species_n.median())}
    for target,col in (
        ("original_all_species","native_minus_intro_white"),
        ("common_genus_equal","genus_equal_native_minus_intro"),
        ("common_genus_min_species_weighted",
            "matched_genera_support_weighted_difference"),
        ("common_genus_species_pooled",
            "matched_genera_species_pooled_difference")
    ):
        d=data.copy()
        d["native_minus_intro_white"]=d[col]
        # Controls reflect original group species support across all genera,
        # retaining identically specified same-cell environmental regression.
        results[target]={}
        for block in ("500","1000"):
            results[target][block+"_km_fixed"]=fit(d,kind="difference",block=block,
                                                  species_weight=False)
        results[target]["500km_cell_species_support_weighted"]=fit(
            d,kind="difference",block="500",species_weight=True)
    return results

def run(gbif,env,lgm,native,spatial,wcvp_dir,tdwg,out):
    out=Path(out)
    require(not out.exists(),"refuse overwrite")
    require(sha(lgm)==LGM_SHA,"signed LGM original ZIP modified")
    import pyarrow.parquet as pq
    with ZipFile(lgm) as z:
        receipt=json.loads(z.read("receipt.json"))
        b=z.read(LGM_CELL)
    require(sha(b)==receipt["output_SHA256"][LGM_CELL],
            "real LGM cell table not original verified file")
    lclim=pq.read_table(io.BytesIO(b)).to_pandas()
    region=assign_cells_broad_region(lclim,tdwg)
    wcvp_names,wcvp_dist,wcvp_receipt=get_wcvp_source(wcvp_dir)
    pairs,photo_species=read_signed_sources(gbif,env)
    signed_native,spatial_grid=signed_tables(native,spatial)

    # Signed same-cell sample fixes both original science and denominator.
    matched={}
    all_ids=[]
    for mode in ("strict","moderate"):
        d=same_cell_cohort(signed_native,spatial_grid,mode,min_each=30)
        matched[mode]=d
        all_ids.append(d[["grid_x","grid_y"]])
    ids=pd.concat(all_ids,ignore_index=True).drop_duplicates()
    source=regional_species_cell_rows(
        pairs,photo_species,wcvp_names,wcvp_dist,region,ids)
    summary,pairs_g=shared_genus_cells(source)
    result={"source_GBIF_distinct_species":37104,
            "source_GBIF_species_cell_pairs":5493931,
            "source_LGM_3GCM_zip_sha256":LGM_SHA,
            "comparison_unit":"current recorded GBIF 100km cell; within-genus native introduced pairs",
            "colour_is_source_photo_species_modal_not_GBIF_local_petal":True,
            "native_status_is_TDWG3_broad_region_not_wild_or_island":True,
            "LGM_warming_is_current_cell_paleoclimate_difference_not_refugia":True,
            "native_introduced_groups":{}}
    for mode in ("strict","moderate"):
        audit=exactly_same_cell_reference(source,signed_native,mode)
        d=matched[mode].merge(summary,on=["grid_x","grid_y"],
                              how="left",validate="one_to_one")
        d["n_shared_genera"]=d.n_shared_genera.fillna(0).astype(int)
        result["native_introduced_groups"][mode]={"audit":audit,
            "all_matched_cells_before_genus_filter":len(matched[mode]),
            "n_cell_has_shared_genus":int(d.n_shared_genera.gt(0).sum()),
            "n_cell_has_ge5_shared_genus":int(d.n_shared_genera.ge(5).sum()),
            "analyses":{}}
        for k in (1,5,10):
            x=d[(d.n_shared_genera>=k)
                &(d.shared_native_species_n>=10)
                &(d.shared_intro_species_n>=10)].copy()
            result["native_introduced_groups"][mode]["analyses"][
                f"atleast_{k}_common_genera_and_10_each_species"]=fit_one_cohort(x,k)
            m=result["native_introduced_groups"][mode]["analyses"][
                f"atleast_{k}_common_genera_and_10_each_species"]
            print("MOIS_LGM_SAMECELL_SHARED_GENUS",mode,k,json.dumps({
                "n_cells":len(x),"n_genus":int(x.n_shared_genera.sum()),
                "original_500_beta":m.get("original_all_species",{}).get(
                    "500_km_fixed",{}).get("effect_white_fraction_per_SD_LGM_warming"),
                "genus_equal_500_beta":m.get("common_genus_equal",{}).get(
                    "500_km_fixed",{}).get("effect_white_fraction_per_SD_LGM_warming"),
                "genus_equal_500_p":m.get("common_genus_equal",{}).get(
                    "500_km_fixed",{}).get("p")
            },sort_keys=True),flush=True)
    out.mkdir(parents=True)
    summary.to_parquet(out/"samecell_native_intro_shared_genera_100km_summary.parquet",
                       index=False)
    (out/"within_genus_native_intro_500km_models.json").write_text(
        json.dumps(result,indent=2,sort_keys=True,default=str)+"\n")
    review={"source_sha256":{
         "GBIF_full_cells":sha(gbif),"climate_AGB_species":sha(env),
         "LGM":LGM_SHA,
         "native_regions":sha(native),"spatial_base":sha(spatial),
         "WCVP_names":sha(Path(wcvp_dir)/"mois350399_wcvp_source_taxon_exact_name_resolution.csv.gz"),
         "WCVP_distribution":sha(Path(wcvp_dir)/"mois350399_wcvp_accepted_native_distribution_TDWG3.csv.gz"),
         "TDWG3_geometry":sha(tdwg)},
         "observational_only":True,
         "native_at_TDWG3_is_not_local_individual_status":True,
         "results":result,
         "output_sha256":{p.name:sha(p) for p in out.iterdir()}}
    (out/"receipt.json").write_text(json.dumps(review,indent=2,sort_keys=True)+"\n")
    print("MOIS_LGM_SHARED_GENUS_FINISHED",json.dumps({
        mode:{name:{
            "n_cells":v.get("n_cell_support"),
            "genus_equal_500_p":v.get("common_genus_equal",{}).get(
                "500_km_fixed",{}).get("p")}
            for name,v in result["native_introduced_groups"][mode]["analyses"].items()}
        for mode in ("strict","moderate")},sort_keys=True),flush=True)
    return review

if __name__=="__main__":
    p=argparse.ArgumentParser()
    for key in ("gbif","env","lgm","native","spatial","wcvp-dir","tdwg","output"):
        p.add_argument("--"+key,required=True)
    a=p.parse_args()
    run(a.gbif,a.env,a.lgm,a.native,a.spatial,a.wcvp_dir,a.tdwg,a.output)
