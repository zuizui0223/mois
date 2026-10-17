#!/usr/bin/env python3
"""Probe official iNaturalist monthly public S3 photo/observation metadata.

No images or metadata bulk files are downloaded here. Validate exact object
size/ETag/LastModified and gzip magic with a 64-byte HTTP Range request only.
The Open Dataset uses photo_id -> observation_uuid -> observed_on, NOT historic
numeric observation_id -> observed_on.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from datetime import datetime,timezone
from pathlib import Path
import requests

BUCKET="https://inaturalist-open-data.s3.amazonaws.com"
PATHS={"photo_table":"photos.csv.gz","observation_table":"observations.csv.gz"}
OFFICIAL_SCHEMA="https://github.com/inaturalist/inaturalist-open-data/blob/main/Metadata/README.md"


def probe(url:str, client:requests.Session) -> dict:
    r=client.head(url,allow_redirects=True,timeout=30)
    r.raise_for_status()
    expected_length=int(r.headers.get("Content-Length","0"))
    etag=r.headers.get("ETag")
    last_modified=r.headers.get("Last-Modified")
    if expected_length<=1024 or not etag or not last_modified:
        raise ValueError("Official monthly metadata object cannot be source-locked")
    with client.get(url,headers={"Range":"bytes=0-63","If-Match":etag},
                    stream=True,timeout=45) as tiny:
        tiny.raise_for_status()
        if tiny.status_code!=206:
            raise ValueError("Official S3 object ignored Range; refuse large source transfer")
        prefix=tiny.raw.read(64,decode_content=False)
        returned=tiny.headers.get("ETag")
        if returned and returned!=etag:
            raise ValueError("Official source changed between HEAD and ranged GET")
    if len(prefix)!=64 or prefix[:2]!=b"\x1f\x8b":
        raise ValueError("Official metadata is not a gzipped input")
    return {"url":url,"content_length_bytes":expected_length,
            "content_length_GiB":round(expected_length/(1024**3),5),
            "etag":etag,"last_modified":last_modified,
            "http_range_bytes":64,"gzip_magic_verified":True,
            "first_64_bytes_SHA256":hashlib.sha256(prefix).hexdigest()}


def run(output:Path) -> dict:
    if output.exists() and any(output.iterdir()):
        raise ValueError("Refuse to overwrite prior provider snapshot source verification")
    output.mkdir(parents=True,exist_ok=True)
    with requests.Session() as client:
        source={k:probe(BUCKET+"/"+v,client) for k,v in PATHS.items()}
    result={
        "schema":"mois_official_monthly_inat_open_data_metadata_source_probe_v1",
        "checked_at_UTC":datetime.now(timezone.utc).isoformat(),
        "official_metadata_schema_documentation":OFFICIAL_SCHEMA,
        "photo_id_source_join":"photos.csv.gz/photo_id -> observation_uuid",
        "current_observation_date_join":"observations.csv.gz/observation_uuid -> observed_on",
        "original_numeric_observation_id_in_open_csv":False,
        "need_full_original_photo_id_rejoin_with_citeable_monthly_snapshot":True,
        "not_all_restricted_or_nonopen_licence_photos_in_snapshot":True,
        "no_photo_pixels_downloaded":True,
        "no_full_metadata_files_downloaded":True,
        "source":source,
        "estimated_input_transfer_GiB":round(sum(v["content_length_bytes"]
               for v in source.values())/1024**3,5),
        "current_observed_on_is_editable_observation_metadata":True,
        "no_flowering_anthesis_inference":True,
    }
    target=output/"mois_official_inaturalist_monthly_open_s3_metadata_source_probe.json"
    target.write_text(json.dumps(result,sort_keys=True,ensure_ascii=False,indent=2)+"\n")
    print("MOIS_OFFICIAL_INAT_PUBLIC_OPEN_DATA_METADATA_SOURCE",json.dumps({
        k:{"GiB":v["content_length_GiB"],"ETag":v["etag"],"modified":v["last_modified"]}
        for k,v in source.items()},sort_keys=True),flush=True)
    return result


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,required=True)
    a=parser.parse_args()
    run(a.output)
