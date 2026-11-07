#!/usr/bin/env python3
"""MOIS 2024 original photograph to world mainland/all-island source.

Do NOT require the photograph to pass ROI or any native trait filter.
Reconstruct exact 2024 photo/landmass/0.05-degree cells using the same
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
SOURCE_2024_FULL_MANIFEST_SHA="f7f65e6bc167089b3034555c655e777f528ecdd29952ab7ee9833d8722f39267"
SELECTED_4096_SHA="541b6134090b61a285ffa9c45ee3df8a382e66ae5628c394ca0c4d47af9207a0"
REMAINING_62360_SHA="a792602a05735c5049c960a80fc6a75f4a1241c0a3abaa9769f452330559bf2d"
COAST_SHA="8dbbe7e071e77e9e75f2d639239099ebca8d5c16d6a07df8169729d49f15cf41"

def require(ok,msg):
    if not bool(ok):raise ValueError("MOIS_2024_GLOBAL_SOURCE: "+msg)
def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for part in iter(lambda:f.read(4<<20),b""):h.update(part)
    return h.hexdigest()

def run(canonical,shoreline,out):
    require(sha(canonical)==CANONICAL_SHA,"Original canonical image-ID SHA drift")
    source=Path(shoreline)
    rec=json.loads((source/"source_receipt.json").read_text())
    require(rec["original_archive_sha256"]==COAST_SHA,"Coastline release SHA drift")
    for file,digest in rec["extracted_sha256"].items():
        require(sha(source/file)==digest,"Modified GSHHG source file "+file)
    out=Path(out)
    require(not out.exists(),"Refuse overwriting existing 2024 target source")
    start=time.monotonic()
    df=pd.read_csv(canonical,compression="gzip",low_memory=False,
        dtype={"photo_id":str,"observation_id":str},
        usecols=["photo_id","observation_id","source_subarm","strict9",
                 "observed_on","latitude","longitude"])
    require(len(df)==1032882 and df.photo_id.is_unique and df.observation_id.is_unique,
            "Original photo identity changed")
    d=df.loc[df.observed_on.astype(str).str[:4].eq("2024")].copy()
    require(len(d)==151355,"2024 original photograph cohort changed")
    d["year"]=2024
    d["unit_id"]=d.photo_id
    land=prepare_land(gpd.read_file(source/"GSHHS_h_L1.shp"))
    lakes=gpd.read_file(source/"GSHHS_h_L2.shp")
    assignment=assign_points(d,land,lakes)
    d=d.merge(assignment,on="unit_id",validate="one_to_one")
    q=aggregate_targets(d)
    require(len(q)==69196,"2024 GSHHG 0.05-degree source cohort changed")
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
    require(int(same.sum())==62778,"same-land 2024 center count changed")
    q["exposure_site_kind"]=np.where(same,"same_land_grid_centre","original_photo_anchor")
    q["query_longitude"]=np.where(same,q.longitude,q.anchor_longitude)
    q["query_latitude"]=np.where(same,q.latitude,q.anchor_latitude)
    q["land_area_km2"]=q.land_id.map(land.set_index("land_id").area_km2)
    q=q.sort_values("target_id").reset_index(drop=True)
    require(q.target_id.is_unique and int(q.n_attempted.sum())==151347,
            "2024 original IDs duplicate or lost")
    require(int(q.n_strict9.sum())==59588,
            "2024 technical ROI statuses changed")
    known=q.loc[q.land_class.isin(["mainland","island"])].copy()
    require(len(known)==66456 and int(known.land_class.eq("mainland").sum())==57548 and
            int(known.land_class.eq("island").sum())==8908,
            "2024 true mainland/island denominator changed")
    known["blind_hash"]=known.target_id.map(
        lambda s:hashlib.sha256(("MOIS-global-2024-v1|"+s).encode()).hexdigest())
    ordered=known.sort_values("blind_hash")
    selected=ordered.head(4096).copy()
    remaining=ordered.iloc[4096:].copy()
    require(len(remaining)==62360 and remaining.target_id.is_unique and
            not set(selected.target_id)&set(remaining.target_id),
            "Duplicate/missing original 2024 source sites")
    fields=["target_id","land_class","query_longitude","query_latitude",
            "n_attempted","n_strict9","source_subarms"]
    out.mkdir(parents=True)
    full=out/"mois_2024_full_photo_landmass_source.csv"
    q.to_csv(full,index=False,float_format="%.12g",lineterminator="\n")
    require(sha(full)==SOURCE_2024_FULL_MANIFEST_SHA,
            "Independent original 2024 known-land cell source SHA differs")
    target=out/"mois_2024_next4096_source.csv"
    selected[fields].to_csv(target,index=False,float_format="%.12g",lineterminator="\n")
    require(sha(target)==SELECTED_4096_SHA,
            "Independent 2024 4096 selected source SHA differs")
    rest=out/"mois_2024_remaining62360_source.csv"
    remaining[fields].to_csv(rest,index=False,float_format="%.12g",lineterminator="\n")
    require(sha(rest)==REMAINING_62360_SHA,
            "Independent 2024 remaining62360 source SHA differs")
    report={
      "schema":"MOIS_2024_GSHHG_original_global_land_all_islands_v1",
      "canonical_source_SHA256":CANONICAL_SHA,
      "GSHHG_original_source_SHA256":COAST_SHA,
      "all_canonical_original_photos":1032882,
      "original2024_photo_count":151355,
      "valid_coordinate_source_photo_count":151347,
      "strict9_preexisting_success":59588,
      "source_original2024_photo_cell_landmass_targets":len(q),
      "known_mainland_targets":57548,
      "known_island_targets":8908,
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
      "selected_4096_SHA256":sha(target),
      "remaining_62360_SHA256":sha(rest),
      "remaining_62360_mainland":int(remaining.land_class.eq("mainland").sum()),
      "remaining_62360_island":int(remaining.land_class.eq("island").sum()),
      "selection_independent_of_flower_colour_or_ROI":True,
      "no_minimum_island_area_km2":None,
      "latitude_exclusion":None,
      "original_photos_downloaded_or_reinferred":0,
      "2024_population_not_reused_from_2020":True,
      "elapsed_source_seconds":round(time.monotonic()-start,1)}
    (out/"2024_source_receipt.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print("MOIS_2024_GLOBAL_ORIGINAL_GSHHG_SOURCE_SHA_VERIFIED",
          json.dumps(report,sort_keys=True),flush=True)
    return report

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--canonical",type=Path,required=True)
    p.add_argument("--shoreline",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    run(a.canonical,a.shoreline,a.out)
