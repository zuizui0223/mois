#!/usr/bin/env python3
"""MOIS native-status restricted global nine-colour community assembly and
geographic gradients, on exact accepted GBIF species×100km cells.

Restore the same original GBIF 37,104 accepted species + 5,493,931 species-cell
pairs, same source iNaturalist photo-hue class, same 199 genus×range-breadth
permutation, and *conservatively* classify current species status in each
100-km cell with original Kew WCVP accepted unique-exact binomial × WGSRPD TDWG3
native/introduced source, rejecting ambiguous coastal botanical regions.

Report the genus×range-breadth null for recorded native species only in
island vs mainland, >1500m, Arctic circle, south and north hemisphere and
north/south midlatitudes. Strict and moderate TDWG3 region assignments are
overlapping sensitivity populations, not independent repetitions.

Native here means a species has TDWG3 regional native status, not proof of a
wild individual, island-specific endemicity or a biological flower observation.
Original photo-hue class is unverified as petal pigmentation.
"""
from __future__ import annotations
import argparse
import hashlib
import io
import json
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pandas as pd

from scripts.analysis.mois_global_hue9_community_assembly_genus_null_20261118 import (
    signed_original,baseline,empirical_entropy,counts_by_cell,
    SEED,GBIF_SHA,ENV_SHA,check,sha,
)
from scripts.analysis.mois_nonwhite8_assembly_white_fixed_null_20261118 import (
    load_strict_nine_null,NULL_SHA,
)
from scripts.analysis.mois_lgm_wcvp_native_regional_block_20261116 import (
    assign_cells_broad_region,get_wcvp_source,match_gbif_to_wcvp,TDWG_GEO_SHA,
)
from scripts.analysis.mois_lgm_native_intro_common_genus_20261116 import (
    regional_species_cell_rows,LGM_SHA,LGM_CELL,
)

GEO_SHA="8677b4e7bf6fa85f482272f86afcc87b153e4c6ec8f8d8030bc29f89fedfdbef"
GEO_CELL="GBIF18577_real_photo9_island_alpine_arctic_hemisphere_cells.csv.gz"
COHORT_NAMES=("strict_native","moderate_native",
              "strict_introduced","moderate_introduced")
N_PERM=199
ARCTIC_LAT=66.562

def load_geog(path):
    check(sha(path)==GEO_SHA,"Original GSHHG-WC21 island/elevation geography SHA changed")
    with ZipFile(path) as z:
        rec=json.loads(z.read("receipt.json"))
        for name,digest in rec["output_sha256"].items():
            check(sha(z.read(name))==digest,
                  "Original GSHHG geography internal CSV/JSON SHA mismatch")
        geom=pd.read_csv(io.BytesIO(z.read(GEO_CELL)),compression="gzip",
                         low_memory=False,dtype={"landmass_sibling_id":"string"})
    check(len(geom)==18577 and
          not geom[["grid_x","grid_y"]].duplicated().any(),
          "Signed full global geography cell source count changed")
    check(rec["GSHHG_geometry"]["n_island_majority7"]==516 and
          rec["GSHHG_geometry"]["n_mainland_majority7"]==11567,
          "Historical GIS mainland-island ID bug may have recurred")
    check(rec["n_source_genus_breadth_global_shuffles"]==199
          and rec["source_accepted_species"]==37104,
          "Previous signed colour null input no longer matches")
    return geom

def load_palaeoclimate_grid(path):
    import pyarrow.parquet as pq
    check(sha(path)==LGM_SHA,"WCVP region source true LGM ZIP SHA changed")
    with ZipFile(path) as z:
        raw=z.read(LGM_CELL)
        rec=json.loads(z.read("receipt.json"))
    check(sha(raw)==rec["output_SHA256"][LGM_CELL],
          "Original 18,577 WorldClim1p4 LGM present-grid parquet SHA wrong")
    cells=pq.read_table(io.BytesIO(raw)).to_pandas()
    check(len(cells)==18577 and not cells[
        ["grid_x","grid_y"]].duplicated().any(),
          "LGM 18,577 real geographic cells changed")
    return cells

def source_guild_rows(pairs,species,lgm,wcvp,tdwg,original_grid):
    check(sha(tdwg)==TDWG_GEO_SHA,"TDWG3 official polygons not source-matched")
    regions=assign_cells_broad_region(lgm,tdwg)
    check(len(regions)==18577 and not regions[
        ["grid_x","grid_y"]].duplicated().any(),
          "Source WCVP region assignment wrong denominator")
    check(int(regions.region_strict.sum())==9288
          and int(regions.region_moderate.sum())==10391,
          "True original strict/moderate TDWG3 geographic assignment drift")
    names,distribution,manifest=get_wcvp_source(wcvp)
    resolved=match_gbif_to_wcvp(species,names)
    check(int(resolved.wcvp_match_exact_binomial.sum())==36762,
          "Accepted COLXR to WCVP unique exact species mapping changed")
    selected=regions.loc[
        regions.region_moderate,["grid_x","grid_y"]].drop_duplicates()
    src=regional_species_cell_rows(pairs,species,names,distribution,
                                   regions,selected)
    src=src.merge(species[["key","source_hue9"]],left_on="specieskey",
                  right_on="key",how="left",validate="many_to_one")
    check(src.source_hue9.notna().all()
          and (src.white==src.source_hue9.eq("white")).all(),
          "Photo-hue species-key provenance failed within native status groups")
    src=src.drop(columns=["key"])
    original_index=pd.MultiIndex.from_frame(
        original_grid[["grid_x","grid_y"]])
    recorded_index=pd.MultiIndex.from_frame(src[["grid_x","grid_y"]])
    ii=original_index.get_indexer(recorded_index)
    spp_index=pd.Index(species.key.astype(str))
    jj=spp_index.get_indexer(src.specieskey.astype(str))
    check((ii>=0).all() and (jj>=0).all(),
          "Original signed WCVP native source GBIF cell/species not in null")
    check(np.array_equal(
        species.source_hue9.iloc[jj].to_numpy(str),
        src.source_hue9.to_numpy(str)),
        "Source original species hue order mapping does not reproduce")
    # Ensure crosswalk rows describe one original GBIF species in one cell.
    check(not src[["grid_x","grid_y","specieskey"]].duplicated().any(),
          "Repeated accepted species in same WCVP-native source cell")
    cohorts={}
    for name in COHORT_NAMES:
        moderate=name.startswith("moderate")
        native=name.endswith("_native")
        condition=(src.region_moderate if moderate else src.region_strict).copy()
        condition &= src.status.eq("native_at_TDWG3" if native else
                                   "introduced_at_TDWG3")
        keep=condition.to_numpy(bool)
        cohorts[name]={
            "cell_index":ii[keep].astype(np.int32),
            "species_index":jj[keep].astype(np.int32),
            "n_source_species_grid_pairs":int(keep.sum()),
        }
    # The strict historical native/introduced denominators are signed
    # by the independent original 2026-11-16 WCVP source audit.
    check(cohorts["strict_native"]["n_source_species_grid_pairs"]==2271852,
          "Original strict WCVP native species-cell source did not reproduce")
    check(cohorts["strict_introduced"]["n_source_species_grid_pairs"]==436422,
          "Original strict WCVP introduced species-cell source did not reproduce")
    return cohorts,{
      "n_WCVP_accepted_unique_exact_matched_species":36762,
      "n_regions_strict":9288,"n_regions_moderate":10391,
      "n_source_resolved_native_introduced_species_grid_pairs":len(src),
      "cohorts_source_species_grid_pairs":{
          k:v["n_source_species_grid_pairs"] for k,v in cohorts.items()},
      "source_WCVP_release_receipt_SHA256":sha(
          Path(wcvp)/"mois350399_wcvp_TDWG3_source_only_taxon_native_receipt.json"),
      "source_TDWG3_geojson_SHA256":sha(tdwg),
    }

def assemble_geographic_masks(geo,n_species):
    lat=geo.center_latitude
    elev=geo.sampled_land_median_elevation_m
    island=geo.gshhg_land_class.eq("island")
    mainland=geo.gshhg_land_class.eq("mainland")
    mask=n_species>=30
    return {
       "all_WCVP_status_species_richness_ge30":mask,
       "richness_ge100":n_species>=100,
       "islands_GSHHG_majority7":mask&island,
       "mainlands_GSHHG_majority7":mask&mainland,
       "islands_all9_same_sibling":mask&island&geo.strict_all9_same_landmass,
       "high_elevation_land_median1500plus":mask&elev.ge(1500),
       "low_elevation_land_median_below500":mask&elev.ge(-430)&elev.lt(500),
       "very_high_elevation_land_median2500plus":mask&elev.ge(2500),
       "Arctic_circle_lat66p562plus":mask&lat.ge(ARCTIC_LAT),
       "north_boreal_lat50to66p562":mask&lat.ge(50)&lat.lt(ARCTIC_LAT),
       "northern_hemisphere":mask&lat.ge(0),
       "southern_hemisphere":mask&lat.lt(0),
       "northern_midlatitude_30to55":mask&lat.ge(30)&lat.lt(55),
       "southern_midlatitude_minus55to_minus30":mask&lat.le(-30)&lat.gt(-55),
       "north_tropical_0to23p44":mask&lat.ge(0)&lat.lt(23.44),
       "south_tropical_minus23p44to0":mask&lat.lt(0)&lat.ge(-23.44),
    }

def start_cohorts(cohorts,n_global_cells,labels,geo):
    for name,c in cohorts.items():
        ident=c["cell_index"].astype(np.int64)*9+labels[c["species_index"]]
        counts=np.bincount(ident,minlength=n_global_cells*9).reshape((-1,9))
        n=counts.sum(axis=1)
        check(int(n.sum())==c["n_source_species_grid_pairs"],
              "Original WCVP native/introduced counts no longer match 9 colours")
        c["counts"]=counts
        c["cell_species_n"]=n
        c["original_Shannon"]=empirical_entropy(counts)
        c["mask"]=assemble_geographic_masks(geo,n)
        c["mask"]={k:np.flatnonzero(np.asarray(a,bool))
                   for k,a in c["mask"].items()
                   if int(np.sum(a))>=25}
        c["expected_Shannon_sum"]=np.zeros(n_global_cells)
        c["null_group_draws"]={k:np.empty(N_PERM)
                               for k in c["mask"]}
    return cohorts

def run_permutations(cohorts,original,source_species,strict_null,geo):
    x=original
    start_cohorts(cohorts,x["ngrids"],x["hue"],geo)
    rng=np.random.default_rng(SEED)
    expected_control=np.zeros(x["ngrids"])
    full_groups={}
    for i in range(N_PERM):
        perm=x["hue"].copy()
        for group in x["variable_groups"]:
            perm[group]=x["hue"][group][rng.permutation(len(group))]
        # Recompute source global expected Shannon to verify exact matching
        # to the previously signed 199-global-permutation null.
        h=empirical_entropy(counts_by_cell(x,perm))
        expected_control+=h
        for name,c in cohorts.items():
            encoded=(c["cell_index"].astype(np.int64)*9+
                     perm[c["species_index"]])
            ncounts=np.bincount(encoded,minlength=x["ngrids"]*9).reshape((-1,9))
            check(np.array_equal(ncounts.sum(axis=1),c["cell_species_n"]),
                  "Whole-catalog hue permutation changed WCVP cohort membership")
            H=empirical_entropy(ncounts)
            c["expected_Shannon_sum"]+=H
            for k,ids in c["mask"].items():
                c["null_group_draws"][k][i]=float(H[ids].mean())
        if (i+1)%50==0:
            print("MOIS_WCVP_NATIVE_NINECOLOUR_199_GLOBAL_PERMS",i+1,flush=True)
    expected_control/=N_PERM
    check(np.allclose(
        expected_control,
        strict_null.expected_hue9_Shannon_given_genus.to_numpy(float),
        atol=1e-10,rtol=0),
        "Original signed genus×breadth 199-permutation null was not reproduced")
    for c in cohorts.values():
        c["expected_Shannon_sum"]/=N_PERM
    return cohorts

def assess(cohorts,geo):
    results={}
    for mode,c in cohorts.items():
        reported={}
        for name,ids in c["mask"].items():
            sampled=geo.iloc[ids]
            y=c["original_Shannon"][ids]
            nulls=c["null_group_draws"][name]
            observed=float(y.mean())
            exp=float(nulls.mean())
            require_repro=float(c["expected_Shannon_sum"][ids].mean())
            check(abs(exp-require_repro)<1e-10,
                  "Native WCVP geographic group perm expected value mismatch")
            delta=observed-exp
            p_low=float((1+np.count_nonzero(nulls<=observed+1e-12))/(N_PERM+1))
            p_two=float((1+np.count_nonzero(abs(nulls-exp)>=
                                            abs(delta)-1e-12))/(N_PERM+1))
            counts=c["counts"][ids]
            n=counts.sum(axis=1)
            reported[name]={
              "n_real_100km_cells":len(ids),
              "n_photo_colour_species_median":float(np.median(n)),
              "observed_mean_Shannon":observed,
              "permuted_genus_breadth_expected_mean_Shannon":exp,
              "observed_minus_null_Shannon":delta,
              "empirical_p_one_sided_more_clustering":p_low,
              "empirical_p_two_sided":p_two,
              "null_95pct_interval":[float(v) for v in
                   np.percentile(nulls,[2.5,97.5])],
              "cell_equal_mean_white_photo_species_share":float(
                   np.mean(counts[:,0]/np.maximum(n,1))),
              "n_distinct_original_500km_blocks":int(
                   (sampled.grid_x.floordiv(5).astype(str)+"|"+
                    sampled.grid_y.floordiv(5).astype(str)).nunique())
            }
        contrasts={}
        groups=(("island_minus_mainland","islands_GSHHG_majority7",
                 "mainlands_GSHHG_majority7"),
                ("alpine_minus_lowland","high_elevation_land_median1500plus",
                 "low_elevation_land_median_below500"),
                ("Arctic_minus_boreal","Arctic_circle_lat66p562plus",
                 "north_boreal_lat50to66p562"),
                ("south_minus_north_hemisphere","southern_hemisphere",
                 "northern_hemisphere"),
                ("southern_mid_minus_northern_mid_30to55",
                 "southern_midlatitude_minus55to_minus30",
                 "northern_midlatitude_30to55"))
        for name,A,B in groups:
            if A not in c["mask"] or B not in c["mask"]:continue
            obs=(reported[A]["observed_mean_Shannon"]-
                 reported[B]["observed_mean_Shannon"])
            diff=c["null_group_draws"][A]-c["null_group_draws"][B]
            delta=obs-float(np.mean(diff))
            contrasts[name]={
                "n_cells_A":reported[A]["n_real_100km_cells"],
                "n_cells_B":reported[B]["n_real_100km_cells"],
                "difference_in_WCVP_genus_breadth_entropy_anomalies":delta,
                "empirical_p_two_sided":float((1+np.count_nonzero(
                    abs(diff-np.mean(diff))>=abs(delta)-1e-12))/(N_PERM+1)),
                "source_A":A,"source_B":B
            }
        results[mode]={"geographic_groups":reported,
                       "geographic_contrasts":contrasts}
    return results

def run(args):
    out=Path(args.out)
    check(not out.exists(),"Do not overwrite signed native WCVP ecology analysis")
    pairs,species=signed_original(args.GBIF,args.ENV)
    x=baseline(pairs,species,null_model="genus_breadth")
    fullnull=load_strict_nine_null(args.NINE)
    geo=load_geog(args.GEO)
    coords=np.asarray(x["cell_key"].tolist(),dtype=np.int64)
    check(np.array_equal(coords,geo[["grid_x","grid_y"]].to_numpy()),
          "Geographic GSHHG cell order not same original signed GBIF source")
    check(np.array_equal(coords,fullnull[["grid_x","grid_y"]].to_numpy()),
          "Genus breadth nine-hue photo null and geography grid differ")
    lgm=load_palaeoclimate_grid(args.LGM)
    src,audit=source_guild_rows(pairs,species,lgm,args.WCVP,args.TDWG,
                               geo[["grid_x","grid_y"]])
    print("MOIS_WCVP_NATIVE_PHOTO9_REGIONAL_SOURCE_VERIFIED",json.dumps(audit),flush=True)
    src=run_permutations(src,x,species,fullnull,geo)
    results=assess(src,geo)
    out.mkdir(parents=True)
    grid=geo[["grid_x","grid_y","center_latitude","center_longitude",
              "gshhg_land_class","strict_all9_same_landmass",
              "sampled_land_median_elevation_m"]].copy()
    for name,c in src.items():
        label=name+"_"
        grid[label+"n_photo9_species"]=c["cell_species_n"]
        grid[label+"white_species_fraction"]=np.divide(
            c["counts"][:,0],c["cell_species_n"],
            out=np.full(len(grid),np.nan),where=c["cell_species_n"]>0)
        grid[label+"Shannon_observed"]=c["original_Shannon"]
        grid[label+"Shannon_null"]=c["expected_Shannon_sum"]
        grid[label+"Shannon_anomaly"]=(c["original_Shannon"]-
                                      c["expected_Shannon_sum"])
    fp=out/"global_WCVP_TDWG3_native_only_photo9_island_alpine_polar_south_cells.csv.gz"
    grid.to_csv(fp,index=False,compression="gzip")
    result={
      "study":"WCVP native/introduced TDWG3 recorded assemblage 9-hue genus-by-range-breadth null across global islands altitudes latitude gradients",
      "source_37104_accepted_photo_hue_taxa":37104,
      "source_full_GBIF_species_cell_pairs":5493931,
      "source_null_199_same_global_genus_by_width_colour_label_permutations":True,
      "source_WCVP_exact_match_and_TDWG3_status":audit,
      "results":results,
      "source_sha256":{"GBIF":GBIF_SHA,"ENV":ENV_SHA,"LGM":LGM_SHA,
                       "TDWG":TDWG_GEO_SHA,"GSHHG_WC21_geography":GEO_SHA,
                       "original_strict_colour_null":NULL_SHA},
      "source_photo_hue_is_not_an_independently_verified_petal_colour":True,
      "WCVP_native_means_broad_TDWG3_not_island_endemic_or_wild_individual":True,
      "small_islands_unresolvable_at_100km":True,
      "strict_and_moderate_not_independent_samples":True,
      "community_assembling_organism_or_pollinator_causality_not_identified":True
    }
    rfile=out/"WCVP_native_only_geographical_gradients_and_nulls.json"
    rfile.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    result["output_sha256"]={p.name:sha(p) for p in (fp,rfile)}
    (out/"receipt.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    for name,mode in results.items():
        print("MOIS_NATIVE_WCVP_REAL_HUE9_GEO_SUCCESS",json.dumps({
            "status":name,
            "n_native_status_source_pairs":audit["cohorts_source_species_grid_pairs"][name],
            "group_results":{k:{
                "n":v["n_real_100km_cells"],
                "white":v["cell_equal_mean_white_photo_species_share"],
                "H9_anomaly":v["observed_minus_null_Shannon"],
                "permutation_p":v["empirical_p_two_sided"]}
                for k,v in mode["geographic_groups"].items()},
            "contrasts":mode["geographic_contrasts"]
        },allow_nan=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    for name in ("GBIF","ENV","NINE","GEO","LGM","WCVP","TDWG","out"):
        p.add_argument("--"+name,required=True)
    run(p.parse_args())
