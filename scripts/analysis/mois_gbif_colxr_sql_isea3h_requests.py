#!/usr/bin/env python3
"""Prepare (and optionally validate) GBIF SQL occurrence cube request JSONs.

GBIF SQL download is preferable to enumerating 28k species via 100k-window
occurrence search. Uses the current COL XR speciesKey and global equal-area
ISEA3H resolution 8 grid (~7774 km2 cells). Public validation only; actual
bulk download requires the user's own GBIF authentication.
"""
from __future__ import annotations
import argparse,hashlib,json,re
from pathlib import Path
from urllib.request import Request,urlopen
from urllib.error import HTTPError
import pandas as pd

SCHEMA="mois_species_first_COLXR_isea3h_r8_sql_cube_request_v1"
VALIDATE_URL="https://api.gbif.org/v1/occurrence/download/request/validate"
KEY=re.compile(r"^[A-Z0-9]{2,16}$")
GRID_RES=8


def require(condition,reason):
    if not bool(condition):raise ValueError("MOIS_GBIF_SQL_CUBE: "+reason)


def sql_for_species_keys(keys,*,resolution=GRID_RES):
    require(1<=len(keys)<=1000 and 1<=resolution<=22,
            "Invalid user source keys or global hex resolution")
    require(len(keys)==len(set(keys)) and all(KEY.fullmatch(str(x)) for x in keys),
            "COL XR taxon identifiers not unique or not valid alpha-numeric keys")
    values=", ".join("'"+k+"'" for k in sorted(keys))
    cell=f"GBIF_ISEA3HCode({resolution}, decimalLatitude, decimalLongitude, 0)"
    return (
      "SELECT specieskey, "+cell+" AS equal_area_hex, COUNT(*) AS observation_count "
      "FROM occurrence WHERE specieskey IN ("+values+") "
      "AND hascoordinate = true AND hasgeospatialissues = false "
      "AND occurrenceStatus = 'PRESENT' "
      "AND decimalLatitude IS NOT NULL AND decimalLongitude IS NOT NULL "
      "GROUP BY specieskey, "+cell)


def make_requests(crosswalk,*,chunk_size=250):
    require(1<=chunk_size<=1000,"GBIF SQL source key batch exceeds supported size")
    required={"status","key"}
    require(required<=set(crosswalk.columns),"Actual GBIF COL XR source absent")
    d=crosswalk.loc[crosswalk.status.eq("colxr_exact_accepted_species") &
                    crosswalk.key.fillna("").astype(str).ne("")].copy()
    keys=sorted(d.key.astype(str).unique())
    require(keys and len(keys)==len(d),
            "COL XR accepted species list empty or unresolved duplicates")
    result=[]
    for start in range(0,len(keys),chunk_size):
        sub=keys[start:start+chunk_size]
        payload={"sendNotification":True,"format":"SQL_TSV_ZIP",
                 "sql":sql_for_species_keys(sub)}
        result.append({"source_species_keys":sub,"GBIF_download_request":payload})
    return result


def validate_public(request):
    payload=json.dumps(request,separators=(",",":")).encode()
    r=Request(VALIDATE_URL,data=payload,method="POST",
              headers={"Accept":"application/json","Content-Type":"application/json",
                       "User-Agent":"MOIS-GBIF-SQL-cube-validation/1.0"})
    try:
        with urlopen(r,timeout=50) as response:
            return {"status_code":response.status,
                    "body":response.read(150000).decode("utf-8",errors="replace")}
    except HTTPError as e:
        return {"status_code":e.code,
                "body":e.read(150000).decode("utf-8",errors="replace")}


def run(source,output,*,chunk_size=250,do_validate=False):
    source=Path(source);output=Path(output)
    require(not output.exists(),"Refuse overwrite frozen GBIF request package")
    frame=pd.read_csv(source,dtype=str,keep_default_na=False,
                      compression="infer")
    requests=make_requests(frame,chunk_size=chunk_size)
    output.mkdir(parents=True)
    files={}
    for i,req in enumerate(requests):
        path=output/f"colxr_isea3h_r8_sql_request_{i:03d}.json"
        path.write_text(json.dumps(req["GBIF_download_request"],indent=2)+"\n")
        files[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
    tested=None
    if do_validate:
        # One small source-proven true-GBIF-key query validates SQL grammar.
        subset=requests[0]["source_species_keys"][:3]
        tested=validate_public({"sendNotification":False,"format":"SQL_TSV_ZIP",
                                "sql":sql_for_species_keys(subset)})
        (output/"public_sql_validation_response.json").write_text(
            json.dumps(tested,indent=2,sort_keys=True)+"\n")
    receipt={"schema":SCHEMA,"source_colxr_species_key_count":
             sum(len(r["source_species_keys"]) for r in requests),
             "source_actual_colxr_crosswalk_sha256":hashlib.sha256(source.read_bytes()).hexdigest(),
             "sql_request_count":len(requests),"sql_input_key_chunk_size":chunk_size,
             "spatial_grid":"ISEA3H global DGGS resolution 8",
             "approx_global_hex_area_km2":7774,
             "grid_coordinate_uncertainty_randomization_m":0,
             "source_taxonomy":"COL XR",
             "user_authenticated_GBIF_download_submitted":False,
             "GBIF_download_DOI":None,
             "public_validator_response":tested,
             "request_files_sha256":files}
    (output/"receipt.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print("MOIS_COLXR_GBIF_ISEA3H_SQL_REQUEST_PACKAGE",
          json.dumps({k:receipt[k] for k in ("source_colxr_species_key_count",
                  "sql_request_count","public_validator_response",
                  "user_authenticated_GBIF_download_submitted")},sort_keys=True),flush=True)
    return receipt


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source",required=True);p.add_argument("--output",required=True)
    p.add_argument("--chunk-size",type=int,default=250)
    p.add_argument("--validate",action="store_true")
    a=p.parse_args()
    run(a.source,a.output,chunk_size=a.chunk_size,do_validate=a.validate)


if __name__=="__main__":
    main()
