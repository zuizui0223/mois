#!/usr/bin/env python3
"""Test GBIF's actual public SQL validator with all 37,104 COL XR accepted species in ONE SQL.

No authenticated GBIF download submission, DOI issuance or bulk occurrence data transfer.
This is syntax/size validation, not a claim of feasible GBIF server execution.
"""
from __future__ import annotations
import argparse,hashlib,io,json,re,time
from pathlib import Path
from zipfile import ZipFile
from urllib.request import Request,urlopen
from urllib.error import HTTPError,URLError
import pandas as pd

SOURCE_ZIP_SHA="abcd3aee766c85a92b86a8f93f8ed696bbe49d0f6cf2f64e5f75b457abf78e6a"
SOURCE_CSV_SHA="5a77f3c968cf221287cd98a855559b3c3e8a8cf1322c08add28878adeb293df8"
COL_XR="7ddf754f-d193-4cc9-b351-99906754a03b"
API="https://api.gbif.org/v1/occurrence/download/request/validate"
COUNTS=(37104,30000,10000,3000,1000,250)

def require(ok,why):
    if not bool(ok):raise ValueError("MOIS_37104_SINGLE_SQL: "+why)
def sha(raw):return hashlib.sha256(raw).hexdigest()

def source_keys(path):
    raw=Path(path).read_bytes()
    require(sha(raw)==SOURCE_ZIP_SHA,"signed COL XR source ZIP SHA changed")
    with ZipFile(io.BytesIO(raw)) as z:
        csv=z.read("current_exact_accepted_COLXR_species.csv.gz")
    require(sha(csv)==SOURCE_CSV_SHA,"accepted roster CSV SHA changed")
    d=pd.read_csv(io.BytesIO(csv),compression="gzip",dtype=str,keep_default_na=False)
    require(len(d)==37104 and d.key.is_unique and d.inat_taxon_id.is_unique,
            "full 37104 species lost")
    require(d.status.eq("colxr_exact_accepted_species").all()
            and d.GBIF_kingdom.eq("Plantae").all()
            and d.COL_XR_checklistKey.eq(COL_XR).all()
            and d.candidate_rank.eq("SPECIES").all(),"accepted COL XR taxonomy invalid")
    keys=sorted(d.key.tolist())
    require(all(re.fullmatch(r"[A-Z0-9]{4,5}",k) for k in keys),"unexpected species key")
    return keys

def cube_sql(keys):
    require(len(keys)==len(set(keys)) and len(keys)>0
            and all(re.fullmatch(r"[A-Z0-9]{4,5}",k) for k in keys),
            "unsafe or duplicated accepted keys")
    cell="GBIF_ISEA3HCode(8, decimalLatitude, decimalLongitude, 0)"
    return ("SELECT specieskey, "+cell+" AS equal_area_hex, "
       "COUNT(*) AS observation_count FROM occurrence WHERE specieskey IN ("+
       ", ".join("'"+key+"'" for key in sorted(keys))+") "
       "AND hascoordinate = true AND hasgeospatialissues = false "
       "AND occurrenceStatus = 'PRESENT' "
       "AND decimalLatitude IS NOT NULL AND decimalLongitude IS NOT NULL "
       "GROUP BY specieskey, "+cell)

def payload(keys):
    return {"sendNotification":True,"format":"SQL_TSV_ZIP","sql":cube_sql(keys)}

def validate_public(req,expected):
    body=json.dumps(req).encode()
    request=Request(API,data=body,method="POST",headers={
        "Content-Type":"application/json","Accept":"application/json",
        "User-Agent":"MOIS-single-37104-COLXR-GBIF-syntax-validation/20261115"})
    try:
        with urlopen(request,timeout=110) as response:
            code=response.status; data=response.read(4_000_000)
    except HTTPError as e:
        code=e.code; data=e.read(400_000)
    except (TimeoutError,OSError,URLError) as e:
        return {"accepted":False,"status":"network_error_indeterminate",
                "error_type":type(e).__name__}
    try:
        obj=json.loads(data.decode("utf-8"))
    except (ValueError,UnicodeDecodeError):
        obj={"body_excerpt":data[:400].decode("utf-8","replace")}
    if code in (200,201):
        sql=obj.get("sql","")
        echoed=set(re.findall(r"'([A-Z0-9]{4,5})'",sql))
        require(echoed==set(expected),"validator changed or dropped species keys")
        return {"accepted":True,"status":"syntax_valid_not_executed","HTTP":code,
                "echoed_unique_species":len(echoed),"response_SHA256":sha(data)}
    if code in (408,425,429,500,502,503,504):
        return {"accepted":False,"status":"transient_or_server_error_inconclusive",
                "HTTP":code,"error_excerpt":str(obj)[:800]}
    return {"accepted":False,"status":"rejected_by_validator","HTTP":code,
            "error_excerpt":str(obj)[:900]}

def run(source,output,validate=True):
    output=Path(output)
    require(not output.exists(),"refuse to overwrite source receipt")
    keys=source_keys(source)
    allreq=payload(keys)
    output.mkdir(parents=True)
    filename=output/"COLXR_37104_one_ISEA3H8_SQL_request.json"
    filename.write_text(json.dumps(allreq,indent=2)+"\n")
    trials=[]
    if validate:
        for n in COUNTS:
            req=payload(keys[:n])
            outcome=validate_public(req,keys[:n])
            trial={"accepted_species":n,"sql_bytes":len(req["sql"].encode()),
                   "request_JSON_bytes":len(json.dumps(req).encode()),**outcome}
            trials.append(trial)
            print("MOIS_ONE_SQL_PUBLIC_VALIDATION",json.dumps(trial,sort_keys=True),flush=True)
            if outcome["accepted"] or outcome["status"]=="transient_or_server_error_inconclusive" or outcome["status"]=="network_error_indeterminate":
                break
            time.sleep(1.5)
    receipt={"schema":"mois_COLXR_37104_one_bulk_SQL_public_validator_v1",
       "signed_source_artifact_SHA256":SOURCE_ZIP_SHA,
       "signed_accepted_roster_SHA256":SOURCE_CSV_SHA,
       "COL_XR_accepted_species":37104,"single_SQL_species":len(keys),
       "single_SQL_bytes":len(allreq["sql"].encode()),
       "single_JSON_bytes":len(json.dumps(allreq).encode()),
       "GBIF_equal_area_isea3h_resolution":8,
       "public_SQL_validator_trials":trials,
       "all_37104_single_SQL_publicly_validated":bool(trials and trials[0]["accepted"]),
       "authenticated_SQL_download_submitted":False,"GBIF_download_key":None,
       "GBIF_download_DOI":None,"GBIF_species_grid_data_retrieved":False,
       "validation_does_not_guarantee_async_SQL_execution":True,
       "full_SQL_request_sha256":sha(filename.read_bytes())}
    (output/"receipt.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print("MOIS_37104_ONE_SQL_PUBLIC_VALIDATOR_COMPLETE",
          json.dumps(receipt,sort_keys=True),flush=True)
    return receipt

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--source",required=True)
    p.add_argument("--output",required=True)
    p.add_argument("--skip-public-validate",action="store_true")
    a=p.parse_args()
    run(a.source,a.output,not a.skip_public_validate)
