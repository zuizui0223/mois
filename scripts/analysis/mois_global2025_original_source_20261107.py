#!/usr/bin/env python3
"""MOIS 2025 original photograph to world mainland/all-island source.

Do NOT require the photograph to pass ROI or any native trait filter.
Reconstruct exact 2025 photo/landmass/0.05-degree cells using the same
source GSHHG 2.3.7 and strict one-photo-per-observation canonical frame
as the completed 2020 cohort, without island-area or latitude exclusions.
"""
from __future__ import annotations
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np,pandas as pd,geopandas as gpd
from scripts.analysis.mois_global_2020_GSHHG_scope_20261107 import (
    assign_points,aggregate_targets,prepare_land)

CANONICAL_SHA="f033661d76e07898a5c13199a55c5d0ce6e55353daa48b0058c456979e90644c"
SOURCE_2025_FULL_MANIFEST_SHA="de2cd14690383f7e66cff68c945cf2a6812b33b0511f0189aaea9f6d62ea37fc"
SELECTED_4096_SHA="fda3a319dad406753bfeaa0bc379795b1a73d1e3991c215550bbbf0681f87a6d"
REMAINING_69367_SHA="68990f90ecc5157ffd21f0fb2269802dbd79be8d86de2e56ab80cfe8b7c72a65"
COAST_SHA="8dbbe7e071e77e9e75f2d639239099ebca8d5c16d6a07df8169729d49f15cf41"

def require(ok,msg):
    if not bool(ok):raise ValueError("MOIS_2025_GLOBAL_SOURCE: "+msg)
def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for part in iter(lambda:f.read(4<<20),b""):h.update(part)
    return h.hexdigest()

def run(canonical,shoreline,out,allow_knownland_verified_geo_roundoff=False):
    require(sha(canonical)==CANONICAL_SHA,"Original canonical image-ID SHA drift")
    source=Path(shoreline)
    rec=json.loads((source/"source_receipt.json").read_text())
    require(rec["original_archive_sha256"]==COAST_SHA,"Coastline release SHA drift")
    for file,digest in rec["extracted_sha256"].items():
        require(sha(source/file)==digest,"Modified GSHHG source file "+file)
    out=Path(out)
    require(not out.exists(),"Refuse overwriting existing 2025 target source")
    start=time.monotonic()
    df=pd.read_csv(canonical,compression="gzip",low_memory=False,
        dtype={"photo_id":str,"observation_id":str},
        usecols=["photo_id","observation_id","source_subarm","strict9",
                 "observed_on","latitude","longitude"])
    require(len(df)==1032882 and df.photo_id.is_unique and df.observation_id.is_unique,
            "Original photo identity changed")
    d=df.loc[df.observed_on.astype(str).str[:4].eq("2025")].copy()
    require(len(d)==172090,"2025 original photograph cohort changed")
    d["year"]=2025
    d["unit_id"]=d.photo_id
    land=prepare_land(gpd.read_file(source/"GSHHS_h_L1.shp"))
    lakes=gpd.read_file(source/"GSHHS_h_L2.shp")
    assignment=assign_points(d,land,lakes)
    d=d.merge(assignment,on="unit_id",validate="one_to_one")
    q=aggregate_targets(d)
    require(len(q)==76630,"2025 GSHHG 0.05-degree source cohort changed")
    ij=q.cell_id.str.split("_",expand=True).astype(int)
    q["longitude"]=(ij[0]+.5)*.05-180
    q["latitude"]=(ij[1]+.5)*.05-90
    centres=q[["cell_id","longitude","latitude"]].drop_duplicates("cell_id").rename(
        columns={"cell_id":"unit_id"})
    ctr=assign_points(centres,land,lakes).rename(columns={
         "unit_id":"cell_id","land_id":"centre_land_id",
         "land_class":"centre_land_class","area_km2":"centre_area_km2"})
    q=q.merge(ctr,on="cell_id",validate="many_to_one")
    same=q.land_id.notna()&q.land_id.eq(q.centre_land_id)&q.land_class.isin(
        ["island","mainland"])
    require(int(same.sum())==69387,"same-land 2025 center count changed")
    q["exposure_site_kind"]=np.where(same,"same_land_grid_centre","original_photo_anchor")
    q["query_longitude"]=np.where(same,q.longitude,q.anchor_longitude)
    q["query_latitude"]=np.where(same,q.latitude,q.anchor_latitude)
    q["land_area_km2"]=q.land_id.map(land.set_index("land_id").area_km2)
    q=q.sort_values("target_id").reset_index(drop=True)
    require(q.target_id.is_unique and int(q.n_attempted.sum())==172062,
            "2025 original IDs duplicate or lost")
    require(int(q.n_strict9.sum())==65546,
            "2025 technical ROI statuses changed")
    known=q.loc[q.land_class.isin(["mainland","island"])].copy()
    require(len(known)==73463 and int(known.land_class.eq("mainland").sum())==63408 and
            int(known.land_class.eq("island").sum())==10055,
            "2025 true mainland/island denominator changed")
    known["blind_hash"]=known.target_id.map(
        lambda s:hashlib.sha256(("MOIS-global-2025-v1|"+s).encode()).hexdigest())
    ordered=known.sort_values("blind_hash")
    selected=ordered.head(4096).copy()
    remaining=ordered.iloc[4096:].copy()
    require(len(remaining)==69367 and remaining.target_id.is_unique and
            not set(selected.target_id)&set(remaining.target_id),
            "Duplicate/missing original 2025 source sites")
    fields=["target_id","land_class","query_longitude","query_latitude",
            "n_attempted","n_strict9","source_subarms"]
    out.mkdir(parents=True)
    all_knownland=out/"mois_2025_all73463_knownland_source.csv"
    known[fields].to_csv(all_knownland,index=False,float_format="%.12g",lineterminator="\n")
    require(sha(all_knownland)=="45c5152e1edad997ff8749c645818996fb2b0bc03a2ab7c0753af446afaab143",
            "2025 full knownland site manifest differs from independent original extraction")
    full=out/"mois_2025_full_photo_landmass_source.csv"
    q.to_csv(full,index=False,float_format="%.12g",lineterminator="\n")
    full_byte_identical=(sha(full)==SOURCE_2025_FULL_MANIFEST_SHA)
    require(full_byte_identical or allow_knownland_verified_geo_roundoff,
            "Independent original 2025 known-land cell source SHA differs")
    if not full_byte_identical:
        # Frozen known-land CSV above is byte-identical: ID, land class,
        # observation count and query coordinate all checked. Some GEOS
        # versions differ only in the full auxiliary geometry-area values.
        # Do NOT report full-source byte identity in this case.
        require(sha(all_knownland)=="45c5152e1edad997ff8749c645818996fb2b0bc03a2ab7c0753af446afaab143",
                "Source IDs/landmass/coords differ, cannot accept area roundoff")
    target=out/"mois_2025_next4096_source.csv"
    selected[fields].to_csv(target,index=False,float_format="%.12g",lineterminator="\n")
    require(sha(target)==SELECTED_4096_SHA,
            "Independent 2025 4096 selected source SHA differs")
    rest=out/"mois_2025_remaining69367_source.csv"
    remaining[fields].to_csv(rest,index=False,float_format="%.12g",lineterminator="\n")
    require(sha(rest)==REMAINING_69367_SHA,
            "Independent 2025 remaining69367 source SHA differs")
    report={
      "schema":"MOIS_2025_GSHHG_original_global_land_all_islands_v1",
      "canonical_source_SHA256":CANONICAL_SHA,
      "GSHHG_original_source_SHA256":COAST_SHA,
      "all_canonical_original_photos":1032882,
      "original2025_photo_count":172090,
      "valid_coordinate_source_photo_count":172062,
      "strict9_preexisting_success":65546,
      "source_original2025_photo_cell_landmass_targets":len(q),
      "known_mainland_targets":63408,
      "known_island_targets":10055,
      "unresolved_source_geography":len(q)-len(known),
      "same_landmass_grid_centres":int(same.sum()),
      "selected_4096_mainland":int(selected.land_class.eq("mainland").sum()),
      "selected_4096_island":int(selected.land_class.eq("island").sum()),
      "selected_4096_island_below_5km2":int((selected.land_class.eq("island")&
                   selected.land_area_km2.lt(5)).sum()),
      "selected_4096_outside_plusminus60":int(selected.query_latitude.abs().gt(60).sum()),
      "selected_4096_attempted_original_photos":int(selected.n_attempted.sum()),
      "selected_4096_strict9_success":int(selected.n_strict9.sum()),
      "full_source_SHA256":sha(full),
      "full_source_reference_SHA256":SOURCE_2025_FULL_MANIFEST_SHA,
      "full_source_bytes_identical_to_original":full_byte_identical,
      "area_roundoff_exception_only_with_byte_identical_knownland_manifest":bool(
          not full_byte_identical and allow_knownland_verified_geo_roundoff),
      "selected_4096_SHA256":sha(target),
      "all_knownland_site_SHA256":sha(all_knownland),
      "remaining_69367_SHA256":sha(rest),
      "remaining_69367_mainland":int(remaining.land_class.eq("mainland").sum()),
      "remaining_69367_island":int(remaining.land_class.eq("island").sum()),
      "selection_independent_of_flower_colour_or_ROI":True,
      "no_minimum_island_area_km2":None,
      "latitude_exclusion":None,
      "original_photos_downloaded_or_reinferred":0,
      "2025_population_not_reused_from_2020":True,
      "elapsed_source_seconds":round(time.monotonic()-start,1)}
    (out/"2025_source_receipt.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print("MOIS_2025_GLOBAL_ORIGINAL_GSHHG_SOURCE_SHA_VERIFIED",
          json.dumps(report,sort_keys=True),flush=True)
    return report

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--canonical",type=Path,required=True)
    p.add_argument("--shoreline",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--allow-knownland-verified-geo-roundoff",action="store_true",
                  help="Only when frozen 73,463 knownland IDs, classes and coordinates SHA match")
    a=p.parse_args()
    run(a.canonical,a.shoreline,a.out,a.allow_knownland_verified_geo_roundoff)
