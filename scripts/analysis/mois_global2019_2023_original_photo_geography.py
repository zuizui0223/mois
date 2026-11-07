#!/usr/bin/env python3
"""Source-frozen global 2019,2021,2022,2023 original MOIS photo landmass frames.

Exact original photo IDs, full GSHHG L1+L2, no area or latitude threshold,
no raw RGB downloads or inferred floral colour statuses.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import numpy as np,pandas as pd,geopandas as gpd
from scripts.analysis.mois_global_2020_GSHHG_scope_20261107 import (
    prepare_land,assign_points,aggregate_targets)
PHOTO_SHA="f033661d76e07898a5c13199a55c5d0ce6e55353daa48b0058c456979e90644c"
COAST_SHA="8dbbe7e071e77e9e75f2d639239099ebca8d5c16d6a07df8169729d49f15cf41"
YEAR={
  2019:{"photos":51311,"colour":20461,"all":28082,"known":27020,
        "island":3843,"mainland":23177,"known_sha":"01f4e2d2e1c8bc9a1392f4143f693563bfa8ea1a0f65ba800e1caa450803831c",
        "full_sha":"bb8cc60c8d8c3a4a4a51425cb09208124e9f89c1d2ba0639eb4ffa6a95dbb484"},
  2021:{"photos":95042,"colour":39113,"all":47714,"known":45966,
        "island":5131,"mainland":40835,"known_sha":"f0d7cb134c09c74cccfb69199396ac30bc5b663d39dc5369d6e296d8c3a1cdca",
        "full_sha":"4ffa96082987acad966e42849522edaae1047279408dda100abee2d605cd5b86"},
  2022:{"photos":110267,"colour":44140,"all":54163,"known":52075,
        "island":6813,"mainland":45262,"known_sha":"594fe8e729cf54aa44859a841b5d3e49e6a66fb85b217da8aff0160044e14df5",
        "full_sha":"bf2a7163d9c71c93150d9a53d2c8d62edc0fe92266c1b87e52497c81d627b480"},
  2023:{"photos":132793,"colour":52594,"all":62289,"known":59862,
        "island":7826,"mainland":52036,"known_sha":"e2825386f45a77da6e3dacdf714fc940cf9ab9cc5174671672b4d57466d04fba",
        "full_sha":"d644dd19cb01fe0787e63785072e3608f989becb475bacf7b0a4dff7b4ca05e8"}
}
def require(ok,msg):
    if not bool(ok):raise ValueError("MOIS_MULTIPLE_YEAR_SOURCE: "+msg)
def sha(p):
    h=hashlib.sha256()
    with Path(p).open("rb") as f:
        for b in iter(lambda:f.read(4<<20),b""):h.update(b)
    return h.hexdigest()
def run(photo_path,shoreline,out,years):
    out=Path(out);require(not out.exists(),"Source output exists; refuse overwrite")
    require(sha(photo_path)==PHOTO_SHA,"Original one-photo ID checksum changed")
    shore=Path(shoreline)
    rec=json.loads((shore/"source_receipt.json").read_text())
    require(rec["original_archive_sha256"]==COAST_SHA,"GSHHG release changed")
    for filename,dig in rec["extracted_sha256"].items():
        require(sha(shore/filename)==dig,"GSHHG source component changed "+filename)
    original=pd.read_csv(photo_path,compression="gzip",low_memory=False,
        dtype={"photo_id":str,"observation_id":str},
        usecols=["photo_id","observation_id","source_subarm","strict9","observed_on","latitude","longitude"])
    require(len(original)==1032882 and original.photo_id.is_unique and original.observation_id.is_unique,
            "Original unique photo/observation population changed")
    original["year"]=pd.to_numeric(original.observed_on.astype(str).str[:4],errors="coerce")
    land=prepare_land(gpd.read_file(shore/"GSHHS_h_L1.shp"))
    lakes=gpd.read_file(shore/"GSHHS_h_L2.shp")
    out.mkdir(parents=True)
    summaries=[]
    for year in years:
        contract=YEAR[year]
        d=original.loc[original.year.eq(year)].copy()
        require(len(d)==contract["photos"] and int(d.strict9.sum())==contract["colour"],
                "Original year source photo counts changed")
        d["unit_id"]=d.photo_id
        d=d.merge(assign_points(d,land,lakes),on="unit_id",validate="one_to_one")
        q=aggregate_targets(d)
        require(len(q)==contract["all"],"Original year landmass grid counts changed")
        xy=q.cell_id.str.split("_",expand=True).astype(int)
        q["longitude"]=(xy[0]+.5)*.05-180
        q["latitude"]=(xy[1]+.5)*.05-90
        centers=q[["cell_id","longitude","latitude"]].drop_duplicates("cell_id").rename(columns={"cell_id":"unit_id"})
        matched=assign_points(centers,land,lakes).rename(columns={
            "unit_id":"cell_id","land_id":"centre_land_id",
            "land_class":"centre_land_class","area_km2":"centre_area_km2"})
        q=q.merge(matched,on="cell_id",validate="many_to_one")
        same=q.land_id.notna()&q.land_id.eq(q.centre_land_id)&q.land_class.isin(("mainland","island"))
        q["exposure_site_kind"]=np.where(same,"same_land_grid_centre","original_photo_anchor")
        q["query_longitude"]=np.where(same,q.longitude,q.anchor_longitude)
        q["query_latitude"]=np.where(same,q.latitude,q.anchor_latitude)
        q["land_area_km2"]=q.land_id.map(land.set_index("land_id").area_km2)
        q=q.sort_values("target_id").reset_index(drop=True)
        good=q.loc[q.land_class.isin(("mainland","island"))].copy()
        require(len(good)==contract["known"] and
                int(good.land_class.eq("mainland").sum())==contract["mainland"] and
                int(good.land_class.eq("island").sum())==contract["island"],
                "Original 2020-compatible GSHHG landmass counts changed")
        folder=out/str(year);folder.mkdir()
        full=folder/f"mois_{year}_full_photo_landmass_source.csv"
        q.to_csv(full,index=False,float_format="%.12g",lineterminator="\n")
        fields=["target_id","land_class","query_longitude","query_latitude",
                "n_attempted","n_strict9","source_subarms"]
        sites=folder/f"mois_{year}_knownland_sites.csv"
        good[fields].to_csv(sites,index=False,float_format="%.12g",lineterminator="\n")
        require(sha(sites)==contract["known_sha"],
                "Original year target ID, landmass, precise coordinates and photo counts SHA mismatch")
        # Geographic area geometry may differ at the last decimal by GEOS,
        # but the signed source target geometry and photo IDs MUST match exactly.
        local_full_bytes_same=sha(full)==contract["full_sha"]
        summary={"year":year,"original_photo_records":len(d),
            "original_source_photo_strict9_success":int(d.strict9.sum()),
            "original_all_photo_landmass_targets":len(q),
            "original_known_land_targets":len(good),
            "original_island_targets":contract["island"],
            "original_mainland_targets":contract["mainland"],
            "unresolved_landmass_targets":len(q)-len(good),
            "high_latitude_targets":int(q.query_latitude.abs().gt(60).sum()),
            "small_island_lt5km2_targets":int((q.land_class.eq("island")&q.land_area_km2.lt(5)).sum()),
            "full_source_SHA256":sha(full),"knownland_source_SHA256":sha(sites),
            "independent_local_full_area_bytes_identical":local_full_bytes_same,
            "local_full_area_GEOS_roundoff_reference_SHA256":contract["full_sha"],
            "target_cell_ID_coordinate_original_photo_fields_exact_SHA":True,
            "no_new_original_photographs_or_ROI_inference":True,"island_area_cutoff":None,
            "latitude_cutoff":None}
        (folder/"source_receipt.json").write_text(json.dumps(summary,sort_keys=True,indent=2)+"\n")
        summaries.append(summary)
        print("MOIS_YEAR_GSHHG_SOURCE_ORIGINAL_VERIFIED",json.dumps(summary,sort_keys=True),flush=True)
    (out/"all_years_original_source_receipt.json").write_text(
        json.dumps({"schema":"MOIS_original_2019_2023_global_landmasses_source_v1",
                    "canonical_SHA256":PHOTO_SHA,"GSHHG_SHA256":COAST_SHA,
                    "year_reports":summaries},sort_keys=True,indent=2)+"\n")
    return summaries
if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--canonical",required=True,type=Path)
    p.add_argument("--shoreline",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    p.add_argument("--years",nargs="+",type=int,default=[2019,2021,2022,2023])
    a=p.parse_args();run(a.canonical,a.shoreline,a.out,a.years)
