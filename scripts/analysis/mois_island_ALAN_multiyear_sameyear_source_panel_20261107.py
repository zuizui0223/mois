#!/usr/bin/env python3
"""Original-photo 5-year, same-year population/nightlight global island panel.

One pre-specified association of strict9 photographed white ROI mass vs
log(1 + original VIIRS) * island, adjusting original geographic blocks,
calendar year, year-matched WorldPop R2024B UC, BIO5/BIO12, elevation,
photographic effort. Data from 2019, 2021, 2023, 2024, 2025.

2022's October source returned provider raster 400 and is NOT quietly treated
as a comparable four-quarter median. 2020's Global1 population is NOT mixed.
The region5 clusters persist across years: do not claim year samples independent.
No raw photographs or restricted source trait rows are distributed.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import pandas as pd
import numpy as np
from scripts.analysis.mois_2019_2023_sameyear_island_ALAN_white_20261107 import (
    PHOTO_SHA,WCVP_SHA,sha,require,verify_environment,make_photo_cells)
from scripts.analysis.mois_2024_sameyear_global_island_photo_white_20261107 import (
    COVARS,make_frames,cluster_fit)
from scripts.analysis.mois_2025_sameyear_island_ALAN_white_20261107 import (
    GEO_SHA,SITES_SHA)
YEAR_COUNTS={2019:547,2021:1223,2023:1788,2024:2085,2025:2339}
PHOTO_YEARS=(2019,2021,2023,2024,2025)
CN_SHA="6e0a077d496e90848d05d043bffae57e4772519f4572a457221c5ef4da0fb813"
UC_SHA="a07af611836d215fb7997c7f574729df10787f5a7cc9bfa1b67bf5a831d613dc"
GEOGRAPHY_2024_SHA="f7f65e6bc167089b3034555c655e777f528ecdd29952ab7ee9833d8722f39267"

def load_frozen_photos(paths):
    require(sha(paths["canonical"])==PHOTO_SHA and
            sha(paths["wcvp"])==WCVP_SHA,
            "Original immutable photo/accepted taxon source SHA drift")
    photo=pd.read_csv(paths["canonical"],compression="gzip",low_memory=False,
        dtype={"photo_id":str,"observation_id":str,"inat_taxon_id":str},
        usecols=["photo_id","observation_id","inat_taxon_id",
          "source_subarm","observed_on","latitude","longitude","strict9",
          "strict9_white","wc2_bio05_10m","wc2_bio12_10m","wc2_elev_2p5m"])
    wcvp=pd.read_csv(paths["wcvp"],compression="gzip",dtype=str,
        usecols=["inat_taxon_id","wcvp_match_state","wcvp_accepted_plant_name_id"])
    return photo,wcvp

def year_historic(year, paths, photos,wcvp, *, return_all=False):
    geo_file=paths[f"geo{year}"]
    sat_file=paths[f"viirs{year}"]
    pop_file=paths[f"pop{year}"]
    from scripts.analysis.mois_2019_2023_sameyear_island_ALAN_white_20261107 import YEAR
    base=Path(geo_file).parent
    source=json.loads((base/"source_receipt.json").read_text())
    satrec=json.loads((Path(sat_file).parent/"receipt.json").read_text())
    poprec=json.loads((Path(pop_file).parent/"source_receipt.json").read_text())
    geo,env=verify_environment(year,geo_file,sat_file,pop_file,
                               source,satrec,poprec)
    _,full,d,_=make_photo_cells(year,photos,wcvp,geo,env)
    return (full,d) if return_all else d

def year2024(paths,photos,wcvp, *, return_all=False):
    require(sha(paths["geo2024"])==GEOGRAPHY_2024_SHA and
            sha(paths["cn2024"])==CN_SHA and sha(paths["uc2024"])==UC_SHA,
            "Original 2024 matched CN/UC/landmass source byte drift")
    from scripts.analysis.mois_2024_UC_R2024B_independent_photo_colour_ALAN import validate_sources
    receipt=json.loads(Path(paths["uc2024receipt"]).read_text())
    geo=pd.read_csv(paths["geo2024"],dtype={"target_id":str,"cell_id":str,"land_id":str})
    cn=pd.read_csv(paths["cn2024"],dtype={"target_id":str},low_memory=False)
    uc=pd.read_csv(paths["uc2024"],dtype={"target_id":str},low_memory=False)
    validate_sources(cn,uc,receipt)
    data=cn.merge(uc[["target_id","population_UC_R2024B_status",
                      "population_UC_R2024B_density_per_km2","UC_source_complete"]],
                  on="target_id",validate="one_to_one")
    data["population_source_status"]=data.pop("population_UC_R2024B_status")
    data["population_density"]=data.pop("population_UC_R2024B_density_per_km2")
    _,full,_,m,_,_=make_frames(photos,wcvp,geo,data)
    require(len(m)==2085 and m.target_id.is_unique,
            "2024 original UC photographic model source changed")
    return (full,m) if return_all else m

def year2025(paths,photos,wcvp, *, return_all=False):
    rec=json.loads(Path(paths["geo2025receipt"]).read_text())
    require(rec.get("full_source_reference_SHA256")==GEO_SHA and
            rec.get("full_source_SHA256")==sha(paths["geo2025"]) and
            rec.get("all_knownland_site_SHA256")==SITES_SHA,
            "2025 GSHHG knownland ID/coast source not verified")
    receipt=json.loads(Path(paths["viirs2025"]).parent.joinpath("receipt.json").read_text())
    require(receipt["original_source_photo_grid_sha256"]==SITES_SHA and
            receipt["source_photo_year"]==receipt["VIIRS_year"]==
            receipt["population_raster_year"]==2025 and
            receipt["original_target_count"]==73463 and
            sha(paths["viirs2025"])==receipt["output_values_SHA256"],
            "2025 original population/nightlight source year or data SHA drift")
    geo=pd.read_csv(paths["geo2025"],dtype={"target_id":str,"cell_id":str,"land_id":str})
    env=pd.read_csv(paths["viirs2025"],dtype={"target_id":str},low_memory=False)
    require(len(env)==73463 and
            set(env.target_id)==set(geo.loc[geo.land_class.isin(
                ("mainland","island")),"target_id"]),
            "2025 original numeric/coast cohort drift")
    env=env.rename(columns={"population_source_status":"population_status",
                            "population_density":"population_density_per_km2"})
    require(env.loc[env.population_status.ne("OK"),
         "population_density_per_km2"].isna().all(),
         "2025 population no-data not retained")
    _,full,d,_=make_photo_cells(2025,photos,wcvp,geo,env)
    return (full,d) if return_all else d

def year_safeguards(year,d):
    require(len(d)==YEAR_COUNTS[year] and d.target_id.is_unique,
            f"Original {year} outcome/effort/raster photo cells changed")
    require(set(d.land_class)=={"island","mainland"},
            f"Original {year} island/mainland selection unsupported")
    require(d.log_light.notna().all() and d.log_pop.notna().all()
            and d.white_pp.between(0,100).all(),
            f"Missing source-colour or current year exposure in {year}")
    result=d[["target_id","land_id","land_class","is_island",
       "region5","white_pp","log_light","log_pop","bio5","log_rain",
       "elevation","log_attempts","photo_attempts",
       "accepted_species_n","strict9_success"]].copy()
    result["year"]=year
    # A source target is year|0.05°cell|landmass, but clustering must
    # treat repeated geographical blocks as dependent across years.
    coord=result.target_id.str.split("|",regex=False,expand=True)
    require(coord.shape[1]==3 and coord[0].eq(str(year)).all(),
            "Original photo-year geography target drift")
    result["cell_id_without_year"]=coord[1].values
    return result

def common_geographic_support_blocks(cells):
    """Only observed island+mainland comparisons within the SAME YEAR and 5° block."""
    by_year_block=cells.groupby(["year","region5"],observed=True).is_island.agg(["min","max"])
    return set(by_year_block.index[
        (by_year_block["min"]==0)&(by_year_block["max"]==1)])

def test_generalizability(cells,formula):
    """Two empirical questions beyond the pooled coefficient.

    (1) Does the island × light slope retain its sign when any one year is
        excluded? These five fits share much source support (NOT five
        independent replications).
    (2) Is that slope identifiable within 5-degree geographic blocks
        containing BOTH island and mainland photographic opportunities?
        Do not extrapolate island-only blocks into mainland comparators.

    Also audit how many island cell-years are concentrated on the largest
    original GSHHG landmass. No photo/coordinate rows are exported.
    """
    require(cells.year.nunique()==5 and len(cells)==sum(YEAR_COUNTS.values()),
            "Generalization diagnostics must use exact original 5-year population")
    leave_year={}
    for omitted in PHOTO_YEARS:
        kept=cells.loc[cells.year.ne(omitted)]
        stats=cluster_fit(kept,formula,"region5")["log_light:is_island"]
        leave_year[str(omitted)]={
            "cell_years":int(len(kept)),
            "island_cell_years":int(kept.is_island.sum()),
            "spatial_blocks":int(kept.region5.nunique()),
            **stats}
    both=common_geographic_support_blocks(cells)
    present=pd.MultiIndex.from_frame(cells[["year","region5"]])
    shared=cells.loc[present.isin(both)].copy()
    shared_island=int(shared.is_island.sum())
    shared_blocks=int(shared.region5.nunique())
    co_support={
      "year_and_block_groups_with_both_land_classes":len(both),
      "distinct_blocks_with_sameyear_island_mainland":shared_blocks,
      "total_spatial_blocks":int(cells.region5.nunique()),
      "co_supported_cell_years":len(shared),
      "co_supported_island_cell_years":shared_island,
      "only_sameyear_island_mainland_comparisons":True,
      "does_not_impute_unobserved_mainland_or_island":True,
    }
    if shared_blocks>=10 and shared_island>=40 and len(shared)>=100:
        co_support["estimable"]=True
        co_support["interaction"]=cluster_fit(
            shared,formula,"region5")["log_light:is_island"]
    else:
        co_support.update({
            "estimable":False,
            "reason":"insufficient_within_5degree_block_island_mainland_common_support"})
    island=cells.loc[cells.is_island.eq(1)]
    ranked=island.groupby("land_id").size().sort_values(ascending=False)
    require(len(ranked)>0,"No real GSHHG islands in original source")
    largest=str(ranked.index[0])
    omitted_island=cells.loc[~(cells.is_island.eq(1)&cells.land_id.eq(largest))]
    dominant={
        "independent_landmasses_with_any_island_photo_cell":int(len(ranked)),
        "largest_island_photo_cell_years":int(ranked.iloc[0]),
        "largest_island_share_of_island_cell_years":float(ranked.iloc[0]/len(island)),
        "five_largest_islands_share_of_island_cell_years":float(ranked.iloc[:5].sum()/len(island)),
        "largest_island_exclusion_interaction":cluster_fit(
            omitted_island,formula,"region5")["log_light:is_island"],
        "largest_island_excluded_cell_years":int(len(omitted_island)),
    }
    return {
        "leave_one_year_out":leave_year,
        "leave_one_year_out_all_negative":all(
            d["beta"]<0 for d in leave_year.values()),
        "geographically_co_supported_island_mainland_comparison":co_support,
        "original_island_sampling_concentration":dominant,
        "diagnostics_are_source_support_checks_not_causes_or_preregistered_confirmations":True,
    }

def run(paths,out):
    out=Path(out);require(not out.exists(),"Do not overwrite prior results")
    original,wcvp=load_frozen_photos(paths)
    frames={}
    for y in (2019,2021,2023):
        frames[y]=year_safeguards(y,year_historic(y,paths,original,wcvp))
    frames[2024]=year_safeguards(2024,year2024(paths,original,wcvp))
    frames[2025]=year_safeguards(2025,year2025(paths,original,wcvp))
    allcells=pd.concat([frames[y] for y in PHOTO_YEARS],ignore_index=True)
    require(allcells.target_id.is_unique and
            len(allcells)==sum(YEAR_COUNTS.values()),
            "Five-year original photographed eligible cell denominator drift")
    require(allcells.year.nunique()==5 and
            allcells.groupby("year").size().to_dict()==YEAR_COUNTS,
            "Five real photo-year sources are not all represented")
    reg="white_pp ~ log_light * is_island + log_pop + "+COVARS+" + C(region5) + C(year)"
    between=cluster_fit(allcells,reg,"region5")
    require("log_light:is_island" in between and
            all(np.isfinite(v) for v in between["log_light:is_island"].values()),
            "Multiyear main island/nightlight interaction not identified")
    island=allcells.loc[allcells.is_island.eq(1)].copy()
    # For a true within-island spatial contrast, at least two different
    # site cells on the same island, not simply one cell across multiple years.
    groups=island.groupby("land_id").agg(
      distinct_sites=("cell_id_without_year","nunique"),
      distinct_light=("log_light","nunique"),
      n_years=("year","nunique"))
    ids=set(groups.index[(groups.distinct_sites>=2)&(groups.distinct_light>=2)])
    within=island.loc[island.land_id.isin(ids)].copy()
    require(len(ids)>=10 and len(within)>=60,
            "Too few independent island units for repeated-year fixed effects")
    within_form="white_pp ~ log_light + log_pop + "+COVARS+" + C(land_id) + C(year)"
    island_fit=cluster_fit(within,within_form,"land_id")
    geographic_temporal_audit=test_generalizability(allcells,reg)
    desc=allcells.groupby(["year","land_class"],observed=True).agg(
        matched_photo_cells=("target_id","nunique"),
        mean_photo_white_ROI_pp=("white_pp","mean"),
        median_source_loglight=("log_light","median")).reset_index()
    result={
      "schema":"MOIS_REAL_5YEAR_ISLANDxALAN_sameyear_population_and_flowerROI_v1",
      "included_years":list(PHOTO_YEARS),
      "excluded_from_main_2022_due_original_Oct_source_unavailable":True,
      "excluded_from_main_2020_due_different_WorldPop_Global1_source":True,
      "current_year_R2024B_unconstrained_population_without_release_mixing":True,
      "original_5_year_cell_counts":YEAR_COUNTS,
      "total_photo_colour_eligible_cells":len(allcells),
      "total_island_cells":int(allcells.is_island.sum()),
      "geographic_clusters_persist_across_years":int(allcells.region5.nunique()),
      "independent_islands_in_supported_within_fit":len(ids),
      "within_island_cells":len(within),
      "original_white_ROI_photographic_not_biologically_verified_species_white":True,
      "multiyear_model":reg,
      "between_island_mainland_light_interaction_geoblock_cluster":between["log_light:is_island"],
      "between_model_full_light_and_island_terms":between,
      "within_island_FE_year_adjusted_light":island_fit["log_light"],
      "geographic_temporal_generalisability":geographic_temporal_audit,
      "year_repeated_same_regions_not_independent_experiments":True,
      "no_causal_ALAN_or_native_plant_extirpation_or_moth_service_claim":True,
      "no_ROI_unknown_imputed":True,
      "original_photographs_downloaded_or_reinferred":0,
      "original_input_SHA256":{k:sha(v) for k,v in paths.items()}}
    out.mkdir(parents=True)
    desc.to_csv(out/"five_year_island_mainland_source_support.csv",index=False)
    result["source_support_SHA256"]=sha(out/"five_year_island_mainland_source_support.csv")
    (out/"five_year_models_and_source_receipt.json").write_text(
        json.dumps(result,sort_keys=True,indent=2)+"\n")
    print("MOIS_REAL_FIVE_YEAR_ISLAND_ALAN_WHITE_ASSOCIATION",
        json.dumps({k:result[k] for k in [
          "total_photo_colour_eligible_cells","total_island_cells",
          "geographic_clusters_persist_across_years",
          "independent_islands_in_supported_within_fit",
          "between_island_mainland_light_interaction_geoblock_cluster",
          "within_island_FE_year_adjusted_light",
          "geographic_temporal_generalisability"]},sort_keys=True),flush=True)
    return result

if __name__=="__main__":
    p=argparse.ArgumentParser()
    source_keys=("canonical","wcvp",
       "geo2019","viirs2019","pop2019",
       "geo2021","viirs2021","pop2021",
       "geo2023","viirs2023","pop2023",
       "geo2024","cn2024","uc2024","uc2024receipt",
       "geo2025","geo2025receipt","viirs2025")
    for key in source_keys:
        p.add_argument("--"+key.replace("_","-"),dest=key,type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    run({k:getattr(a,k) for k in source_keys},a.out)
