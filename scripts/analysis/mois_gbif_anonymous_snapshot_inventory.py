#!/usr/bin/env python3
"""Inventory EVERY Parquet object of a GBIF public snapshot via anonymous AWS S3.

This stage downloads no occurrence rows and needs no AWS/GBIF credentials.
The exact, auditable file roster is used for disjoint partitioned parquet scans.
"""
from __future__ import annotations
import argparse,hashlib,json,re
from pathlib import Path
import pandas as pd

BUCKET="gbif-open-data-us-east-1"
REGION="us-east-1"
DATE="2026-11-01"
CITATION_SHA="d9770dd7ef38ee8fb0e461f460e5c26abed60ed4a13319747b7f726f8c20fc0"

def require(x,msg):
    if not bool(x):raise ValueError("MOIS_GBIF_ANON_INVENTORY: "+msg)
def file_sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def distribute(sizes,groups=16):
    require(groups>=1,"partition count")
    bins=[0]*groups
    assigned=[]
    for key,size in sorted(sizes,key=lambda x:(-x[1],x[0])):
        g=min(range(groups),key=lambda x:(bins[x],x))
        assigned.append((key,int(size),g))
        bins[g]+=int(size)
    return assigned,bins

def run(out,groups=16):
    import boto3
    from botocore import UNSIGNED
    from botocore.config import Config
    out=Path(out);require(not out.exists(),"refuse overwritten inventory")
    s3=boto3.client("s3",region_name=REGION,config=Config(signature_version=UNSIGNED,
        connect_timeout=15,read_timeout=90,retries={"max_attempts":4}))
    prefix=f"occurrence/{DATE}/occurrence.parquet/"
    raw=[];n_pages=0
    for pg in s3.get_paginator("list_objects_v2").paginate(
        Bucket=BUCKET,Prefix=prefix,PaginationConfig={"PageSize":1000}):
        n_pages+=1
        for obj in pg.get("Contents",[]):
            key=obj["Key"]
            if int(obj["Size"])==0 and key.endswith("/"):continue
            require(re.fullmatch(re.escape(prefix)+r"\d{6,}",key),
                    "unexpected non-parquet file path "+key)
            raw.append({"key":key,"bytes":int(obj["Size"]),
                        "etag":str(obj.get("ETag","")),
                        "last_modified":obj["LastModified"].isoformat()})
    require(raw and len({x["key"] for x in raw})==len(raw),"no Parquet source or duplicate keys")
    # Disjoint, approximately equal total-compressed-byte work shares.
    assigned,totals=distribute([(x["key"],x["bytes"]) for x in raw],groups)
    bykey={x["key"]:x for x in raw}
    a=pd.DataFrame([{**bykey[key],"partition":part} for key,_,part in assigned]).sort_values("key")
    require(a.key.is_unique and len(a)==len(raw),"inventory lost files")
    require(a.groupby("partition").bytes.sum().sum()==a.bytes.sum(),"incomplete partition accounting")
    c=s3.get_object(Bucket=BUCKET,Key=f"occurrence/{DATE}/citation.txt")["Body"].read(100000)
    # Citation raw content is a stronger source reference than presumed snapshot DOI.
    out.mkdir(parents=True)
    a.to_csv(out/"GBIF_2026-11-01_all_anonymous_parquet_objects.csv.gz",index=False,
             compression={"method":"gzip","mtime":0})
    (out/"citation.txt").write_bytes(c)
    summary={"schema":"mois_anon_GBIF_S3_20261101_complete_file_inventory_v1",
       "snapshot_bucket":BUCKET,"snapshot_date":DATE,
       "snapshot_prefix":prefix,"record_files_inventory_complete":True,
       "S3_list_pages":n_pages,"Parquet_file_count":len(a),
       "total_Parquet_compressed_bytes":int(a.bytes.sum()),
       "total_Parquet_compressed_GB_decimal":round(float(a.bytes.sum()/1e9),3),
       "min_Parquet_file_bytes":int(a.bytes.min()),
       "max_Parquet_file_bytes":int(a.bytes.max()),
       "partitions":groups,
       "partition_bytes":{str(i):int(b) for i,b in enumerate(totals)},
       "citation_SHA256":hashlib.sha256(c).hexdigest(),
       "no_GBIF_account_used":True,"no_AWS_account_used":True,
       "occurrence_rows_downloaded_by_inventory":0,
       "global_species_grid_computed":False,
       "files_SHA256":{p.name:file_sha(p) for p in out.iterdir()}}
    (out/"receipt.json").write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n")
    print("MOIS_GBIF_FULL_PUBLIC_SNAPSHOT_FILE_INVENTORY_SUCCESS",
          json.dumps(summary,sort_keys=True),flush=True)
    return summary

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--output",required=True);p.add_argument("--partitions",type=int,default=16)
    a=p.parse_args();run(a.output,a.partitions)
