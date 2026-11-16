#!/usr/bin/env python3
"""Find and verify real public GLOBAL NDVI climatology without any service credentials.

Source is a third-party Zenodo raster compilation of MODIS NDVI, not biomass.
Fail closed if license/temporal provenance/units or downloadable size unknown.
"""
import argparse,hashlib,json,urllib.request,urllib.error,zipfile
from pathlib import Path

DOI="10.5281/zenodo.22967479"
API="https://zenodo.org/api/records/22967479"
TARGET="NDVI-10K_ANNUAL_NORMALS.zip"

def request(url):
    req=urllib.request.Request(url,headers={
        "User-Agent":"MOIS-GBIF-global-NDVI-remote-source-verification/20261116",
        "Accept":"application/json"})
    with urllib.request.urlopen(req,timeout=80) as r:return r.read(8_000_000)

def run(output):
    out=Path(output)
    if out.exists():raise RuntimeError("Refuse overwriting prior NDVI source audit")
    out.mkdir(parents=True)
    evidence={"schema":"mois_global_NDVI_Zenodo_public_source_provenance_v1",
       "target_DOI":DOI,"catalog_API":API,"target_filename":TARGET,
       "downloaded_NDVI_raster":False,"used_for_models":False,
       "NDVI_is_biomass":False,"GBIF_species_NDVI_joined":False,
       "credential_required":False,"source_snapshot_citation_validated":False}
    try:
        raw=request(API)
        source=json.loads(raw)
        (out/"public_dataset_metadata.json").write_bytes(raw)
        evidence["catalog_response_sha256"]=hashlib.sha256(raw).hexdigest()
        evidence["repository_record_id"]=source.get("id")
        evidence["publication_date"]=source.get("metadata",{}).get("publication_date")
        evidence["title"]=source.get("metadata",{}).get("title")
        evidence["license"]=source.get("metadata",{}).get("license",{})
        files=[]
        for f in source.get("files",[]):
            files.append({"key":f.get("key"),"size":f.get("size"),
                "checksum":f.get("checksum"),"self":f.get("links",{}).get("self"),
                "download":f.get("links",{}).get("content")})
        evidence["available_file_count"]=len(files)
        evidence["target_asset"]=next((f for f in files if f["key"]==TARGET),None)
        evidence["source_snapshot_citation_validated"]=bool(
            evidence["repository_record_id"]==22967479 and evidence["target_asset"] is not None)
        (out/"ndvi_candidate_source_files.json").write_text(json.dumps(files,indent=2)+"\n")
        if evidence["target_asset"]:
            f=evidence["target_asset"]
            link=f["download"] or f["self"]
            if not link or f["size"]>150_000_000:
                evidence["raster_status"]="available_metadata_only_unbounded_or_large"
            else:
                path=out/"ndvi_annual_climatology.zip"
                req=urllib.request.Request(link,headers={"User-Agent":"MOIS-NDVI-signed-source/20261116"})
                count=0; digest=hashlib.sha256()
                with urllib.request.urlopen(req,timeout=150) as incoming, path.open("wb") as target:
                    while True:
                        b=incoming.read(1<<20)
                        if not b:break
                        count+=len(b)
                        if count>150_000_000:raise RuntimeError("NDVI download exceeds fixed source cap")
                        target.write(b);digest.update(b)
                evidence["download_sha256"]=digest.hexdigest()
                evidence["download_bytes"]=count
                if f["size"] and count!=f["size"]:raise RuntimeError("NDVI size mismatch")
                with zipfile.ZipFile(path) as z:
                    members=[{"name":x.filename,"uncompressed":x.file_size} for x in z.infolist()]
                evidence["archive_members"]=members
                evidence["downloaded_NDVI_raster"]=True
                evidence["raster_status"]="actual_public_archive_fetched"
    except Exception as e:
        evidence["raster_status"]="source_access_failed_not_imputed"
        evidence["error_type"]=type(e).__name__
        evidence["error_detail_excerpt"]=str(e)[:300]
    (out/"receipt.json").write_text(json.dumps(evidence,sort_keys=True,indent=2,default=str)+"\n")
    print("MOIS_PUBLIC_NDVI_SOURCE_PROBE",json.dumps({k:v for k,v in evidence.items()
        if k not in ("archive_members","license")},sort_keys=True,default=str),flush=True)
    return evidence

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--output",required=True)
    args=p.parse_args();run(args.output)
