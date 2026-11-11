#!/usr/bin/env python3
"""Small authentic GBIF occurrence-to-equal-area-grid pilot for MOIS photo species.

Only previously source-verified exact accepted COL XR species keys enter this
query. Enumerate complete filtered occurrence pages for small-record taxa,
skip taxa above cap (never mistake a truncated page for the full distribution).
The hue remains a SPECIES-LEVEL observed-photo trait, not a site hue.
"""
from __future__ import annotations
import argparse
import json
import math
import time
from collections import Counter
from pathlib import Path
from urllib.error import HTTPError,URLError
from urllib.parse import urlencode
from urllib.request import Request,urlopen

import pandas as pd
from pyproj import Transformer
from scripts.analysis.mois_gbif_colxr_nine_palette_name_pilot import COL_XR,require

API="https://api.gbif.org/v1/occurrence/search"
SCHEMA="mois_gbif_colxr_exact_species_100km_grid_pilot_v1"
GRID_METRES=100000


def api_page(taxon_key,offset=0,limit=300,*,attempts=3):
    assert 0<=limit<=300 and offset>=0
    args={"checklistKey":COL_XR,"taxonKey":str(taxon_key),
          "hasCoordinate":"true","hasGeospatialIssue":"false",
          "occurrenceStatus":"PRESENT","limit":str(limit),"offset":str(offset)}
    req=Request(API+"?"+urlencode(args),headers={
        "User-Agent":"MOIS-global-ecology-GBIF-grid-pilot/1.0",
        "Accept":"application/json"})
    last=None
    for i in range(attempts):
        try:
            with urlopen(req,timeout=30) as f:
                return json.loads(f.read(2500000))
        except HTTPError as e:
            last=e
            if e.code not in (429,500,502,503,504):break
        except (URLError,TimeoutError,ValueError) as e:
            last=e
        if i+1<attempts:time.sleep(2**i)
    raise RuntimeError("GBIF occurrence source request failed for "+
                       str(taxon_key)+": "+str(last))


def choose_species(name_matches,*,max_per_hue=1):
    require(max_per_hue>=1,"At least 1 species/hue required")
    required={"status","key","inat_taxon_id","source_species","dominant_hue9"}
    require(required<=set(name_matches.columns),"GBIF-name source crosswalk incomplete")
    a=name_matches.loc[name_matches.status.eq("colxr_exact_accepted_species") &
                       name_matches.key.fillna("").ne("")].copy()
    require(a.key.is_unique,"COL XR speciesKey must be distinct for pilot")
    return a.sort_values(["dominant_hue9","inat_taxon_id"],kind="stable")


def allowed_coordinates(record):
    try:
        lat=float(record["decimalLatitude"]);lon=float(record["decimalLongitude"])
        unc=record.get("coordinateUncertaintyInMeters")
        unc=float(unc) if unc not in (None,"") else None
        return (math.isfinite(lat) and math.isfinite(lon) and -90<=lat<=90 and
                -180<=lon<=180 and (lat!=0 or lon!=0) and
                (unc is None or (math.isfinite(unc) and unc<=50000)))
    except (KeyError,TypeError,ValueError):
        return False


def grid_100km(lon,lat,proj):
    x,y=proj.transform(float(lon),float(lat))
    require(math.isfinite(x) and math.isfinite(y),"Invalid projected occurrence point")
    return math.floor(x/GRID_METRES),math.floor(y/GRID_METRES)


def run(name_match_csv,output,*,max_per_hue=1,max_occurrences=1200):
    output=Path(output)
    require(not output.exists() and 100<=max_occurrences<=3000,
            "Refuse overwrite or unsupported completeness cap")
    names=pd.read_csv(name_match_csv,dtype=str,keep_default_na=False)
    accepted=choose_species(names,max_per_hue=max_per_hue)
    require(len(accepted)>0,"No actually matched COL XR exact accepted species")
    output.mkdir(parents=True)
    projection=Transformer.from_crs("EPSG:4326","EPSG:6933",always_xy=True)
    selected=Counter();source=[];gridrows=[]
    for row in accepted.itertuples(index=False):
        if selected[row.dominant_hue9]>=max_per_hue:continue
        selected[row.dominant_hue9]+=1  # at most one target per hue in pilot
        record={"inat_taxon_id":row.inat_taxon_id,
                "source_species":row.source_species,
                "GBIF_colxr_accepted_species_key":row.key,
                "photo_palette9":row.dominant_hue9,
                "colour_at_gbif_location_measured":False}
        try:
            first=api_page(row.key,0,1)
            count=int(first.get("count",0))
            record["gbif_filter_occurrence_count"]=count
            if count>max_occurrences:
                record["status"]="above_enumeration_cap_no_grid"
                source.append(record);continue
            if count==0:
                record["status"]="no_filtered_occurrences";source.append(record);continue
            seen=set();cells=Counter()
            for offset in range(0,count,300):
                batch=api_page(row.key,offset,min(300,count-offset))
                for hit in batch.get("results",[]):
                    gbifid=str(hit.get("key",""))
                    require(gbifid and gbifid not in seen,"Repeated/missing GBIF occurrence ID across pages")
                    seen.add(gbifid)
                    if not allowed_coordinates(hit):continue
                    x,y=grid_100km(hit["decimalLongitude"],hit["decimalLatitude"],projection)
                    cells[(x,y)]+=1
                time.sleep(0.3)
            require(len(seen)==count,"GBIF paginated occurrence count changed during pilot")
            record.update(status="all_filtered_occurrences_enumerated",
                          original_georeferenced_records=len(seen),
                          quality_retained_records=sum(cells.values()),
                          occupied_100km_cells=len(cells))
            for (x,y),n in sorted(cells.items()):
                gridrows.append({**record,
                                 "grid_epsg6933_x_index":x,
                                 "grid_epsg6933_y_index":y,
                                 "gbif_records_in_cell":n})
        except Exception as e:
            record.update(status="source_error_no_grid",reason=str(e)[:200])
        source.append(record)
        time.sleep(0.4)
    summary=pd.DataFrame(source)
    grids=pd.DataFrame(gridrows,columns=[
        "inat_taxon_id","source_species","GBIF_colxr_accepted_species_key",
        "photo_palette9","colour_at_gbif_location_measured","status",
        "gbif_filter_occurrence_count","original_georeferenced_records",
        "quality_retained_records","occupied_100km_cells",
        "grid_epsg6933_x_index","grid_epsg6933_y_index","gbif_records_in_cell"])
    summary.to_csv(output/"gbif_species_distributions_pilot_receipt.csv",index=False)
    grids.to_csv(output/"gbif_species_by_100km_grid_photo_hue_overlay_pilot.csv",index=False)
    report={"schema":SCHEMA,"source_GBIF_taxonomy":"COL XR",
            "source_checklistKey":COL_XR,"api_source":API,
            "sample_species_queried":len(summary),
            "hue_categories_queried":len(selected),
            "status_counts":{str(k):int(v) for k,v in summary.status.value_counts().items()},
            "successfully_enumerated_species":int(summary.status.eq("all_filtered_occurrences_enumerated").sum()),
            "species_grid_rows":len(grids),
            "raw_occurrences_not_republished":True,
            "grid_projection":"EPSG:6933",
            "grid_cell_length_m":GRID_METRES,
            "zero_grid_for_truncated_species":True,
            "photo_hue_at_GBIF_locations_measured":False,
            "species_sampling_biased_and_not_global_inference":True,
            "GBIF_download_DOI_created":False}
    (output/"receipt.json").write_text(json.dumps(report,sort_keys=True,indent=2)+"\n")
    print("MOIS_ACTUAL_GBIF_COLXR_100KM_GRID_PILOT",json.dumps(report,sort_keys=True),flush=True)
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--matches",required=True);p.add_argument("--output",required=True)
    p.add_argument("--max-per-hue",type=int,default=1)
    p.add_argument("--max-occurrences",type=int,default=1200)
    a=p.parse_args()
    run(a.matches,a.output,max_per_hue=a.max_per_hue,max_occurrences=a.max_occurrences)


if __name__=="__main__":
    main()
