#!/usr/bin/env python3
"""NO GBIF ACCOUNT OR AWS ACCOUNT: inspect one authentic public GBIF monthly snapshot
from unsigned S3 and compare its name/key columns to signed MOIS COL XR 37,104.
This is a source/identifier compatibility audit, not all-snapshot extraction.
"""
from __future__ import annotations
import argparse,hashlib,io,json,re,datetime
from pathlib import Path
from zipfile import ZipFile
import pandas as pd
SOURCE_SHA="abcd3aee766c85a92b86a8f93f8ed696bbe49d0f6cf2f64e5f75b457abf78e6a"
CSV_SHA="5a77f3c968cf221287cd98a855559b3c3e8a8cf1322c08add28878adeb293df8"
BUCKET="gbif-open-data-us-east-1"
REGION="us-east-1"
SCHEMA="mois_anon_GBIF_S3_COLXR_taxonomy_1_parquet_provenance_v1"

def require(ok,msg):
    if not bool(ok):raise ValueError("MOIS_GBIF_ANON_SNAPSHOT: "+msg)

def sha(v):return hashlib.sha256(v).hexdigest()

def snapshot_name(prefix):
    m=re.fullmatch(r"occurrence/(\d{4}-\d{2}-\d{2})/",str(prefix))
    if not m:return None
    try:datetime.date.fromisoformat(m.group(1))
    except ValueError:return None
    return m.group(1)

def choose_recent_snapshot(prefixes,asof="2026-11-15"):
    a=sorted({date for p in prefixes if (date:=snapshot_name(p)) is not None and date<=asof})
    require(a,"no source snapshot as of date")
    return a[-1]

def key_overlap(occ_keys,col_keys):
    return len({str(x) for x in occ_keys if str(x).strip()}.intersection(
        str(x) for x in col_keys))

def canonical_pair(name):
    t=str(name).split()
    if len(t)<2 or not re.fullmatch(r"[A-Z][A-Za-z-]+",t[0]) or not re.fullmatch(r"[a-z][a-z-]+",t[1]):
        return None
    if t[1] in ("sp","spp","aff","cf"):return None
    return t[0]+" "+t[1]

def classify_join(key_matches,name_matches):
    if key_matches:return "actual_colxr_key_overlap"
    return "taxonomy_mismatch_key_bridge_required" if name_matches else "no_matched_species_in_sample"

def source_traits(path):
    raw=Path(path).read_bytes()
    require(sha(raw)==SOURCE_SHA,"MOIS COL XR ZIP changed")
    with ZipFile(io.BytesIO(raw)) as z:b=z.read("current_exact_accepted_COLXR_species.csv.gz")
    require(sha(b)==CSV_SHA,"MOIS COL XR CSV changed")
    d=pd.read_csv(io.BytesIO(b),compression="gzip",dtype=str,keep_default_na=False)
    require(len(d)==37104 and d.key.is_unique and d.inat_taxon_id.is_unique,
            "source accepted 37104 not intact")
    require(d.status.eq("colxr_exact_accepted_species").all()
            and d.GBIF_kingdom.eq("Plantae").all(),"source taxonomic status drift")
    d["canonical"]=d.species.map(canonical_pair)
    require(d.canonical.notna().all(),"some MOIS accepted species not binomial")
    return d

def run(source,out,asof="2026-11-15"):
    # boto3 UNSIGNED: NO AWS credentials and NO GBIF credentials.
    import boto3
    from botocore import UNSIGNED
    from botocore.config import Config
    import pyarrow as pa
    import pyarrow.parquet as pq
    from pyarrow import fs
    out=Path(out)
    require(not out.exists(),"refuse overwrite source audit")
    d=source_traits(source)
    client=boto3.client("s3",region_name=REGION,config=Config(signature_version=UNSIGNED,
        connect_timeout=15,read_timeout=90,retries={"max_attempts":3}))
    pages=client.get_paginator("list_objects_v2").paginate(Bucket=BUCKET,Prefix="occurrence/",
        Delimiter="/",PaginationConfig={"PageSize":1000,"MaxItems":10000})
    prefixes=[]
    for pg in pages:prefixes.extend(x["Prefix"] for x in pg.get("CommonPrefixes",[]))
    date=choose_recent_snapshot(prefixes,asof=asof)
    base=f"occurrence/{date}/"
    cit=client.get_object(Bucket=BUCKET,Key=base+"citation.txt")
    citation=cit["Body"].read(100000)
    require(citation and len(citation)<100000,"citation too large or empty")
    listing=client.list_objects_v2(Bucket=BUCKET,Prefix=base+"occurrence.parquet/",MaxKeys=50)
    objects=[{"key":x["Key"],"size":int(x["Size"]),"etag":str(x.get("ETag",""))}
      for x in listing.get("Contents",[]) if int(x.get("Size",0))>0 and not x["Key"].endswith("/")]
    require(objects,"GBIF public parquet objects not listable")
    one=sorted(objects,key=lambda x:(x["size"],x["key"]))[0]
    s3=fs.S3FileSystem(region=REGION,anonymous=True)
    with s3.open_input_file(f"{BUCKET}/{one['key']}") as opened:
        file=pq.ParquetFile(opened)
        schema=file.schema_arrow
        requested=["specieskey","species","scientificname","genus","kingdom",
           "decimallatitude","decimallongitude","countrycode","occurrencestatus",
           "taxonkey"]
        cols=[k for k in requested if k in schema.names]
        require("specieskey" in cols and "species" in cols and "kingdom" in cols,
                "expected GBIF Parquet taxonomy column absent")
        rg=file.read_row_group(0,columns=cols)
        nrows=file.metadata.num_rows
        ngroups=file.metadata.num_row_groups
    require(rg.num_rows>0,"first row group empty")
    z=rg.to_pandas()
    plants=z.loc[z.kingdom.astype(str).eq("Plantae")].copy()
    species=plants.species.fillna("").astype(str).map(canonical_pair)
    names=set(d.canonical)
    key_hits=key_overlap(plants.specieskey.dropna().astype(str),d.key)
    matched=species.isin(names)
    compatible=int(matched.sum())
    matched_names=species[matched].value_counts().head(30).to_dict()
    # Only 1 row group: a diagnostic example, NOT full source match completeness.
    report={
       "schema":SCHEMA,"date":date,"snapshot_bucket":BUCKET,
       "snapshot_source_prefix":base,
       "snapshot_citation_SHA256":sha(citation),
       "first_visible_object":one,
       "visible_object_count_first_page":len(objects),
       "source_file_schema":str(schema),
       "specieskey_arrow_type":str(schema.field("specieskey").type),
       "taxonomy_index_note":"COL XR accepted key may differ from GBIF snapshot's default checklist key",
       "first_object_parquet_rows":nrows,"first_object_row_groups":ngroups,
       "one_row_group_rows":int(len(z)),
       "one_row_group_Plantae_rows":int(len(plants)),
       "one_row_group_plant_matching_MOIS_binomial_rows":compatible,
       "one_row_group_plant_matching_MOIS_distinct_names":int(species[matched].nunique()),
       "one_row_group_plant_matching_MOIS_COLXR_specieskeys":int(key_hits),
       "join_diagnostic":classify_join(key_hits,compatible),
       "matching_name_examples":matched_names,
       "no_GBIF_username_or_password":True,"no_AWS_access_key_or_account":True,
       "full_snapshot_scanned":False,"37104_distributions_fully_recovered":False,
       "sample_is_not_global_geographic_distribution":True
    }
    out.mkdir(parents=True)
    (out/"citation.txt").write_bytes(citation)
    (out/"receipt.json").write_text(json.dumps(report,sort_keys=True,indent=2)+"\n")
    # Example only: source species in the first row group; no original occurrences redistributed.
    (out/"example_MOIS_species_matches.csv").write_text(
       "species,example_occurrences_in_one_parquet_rowgroup\n"+
       "".join(f"{name},{n}\n" for name,n in matched_names.items()))
    report["output_SHA256"]={p.name:sha(p.read_bytes()) for p in out.iterdir()}
    (out/"receipt.json").write_text(json.dumps(report,sort_keys=True,indent=2)+"\n")
    print("MOIS_PUBLIC_GBIF_NOAUTH_S3_PILOT",json.dumps({
        k:v for k,v in report.items() if k not in ("source_file_schema","matching_name_examples","output_SHA256")
    },sort_keys=True),flush=True)
    return report

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--source",required=True);p.add_argument("--output",required=True)
    p.add_argument("--asof",default="2026-11-15")
    a=p.parse_args();run(a.source,a.output,a.asof)
