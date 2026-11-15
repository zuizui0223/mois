#!/usr/bin/env python3
"""Fail-closed reduction of 16 FULL anonymous S3 GBIF snapshot partition scans.

Only final output is called GLOBAL when each of the 10,346 fixed Parquet
objects appears exactly once with its signed manifest key, ETag and length.
DuckDB externally aggregates species x EPSG:6933 100-km cells, never conflated
with the unrelated ISEA3H resolution 8 grid used by GBIF SQL downloads.
"""
from __future__ import annotations
import argparse,collections,hashlib,io,json
from pathlib import Path
from zipfile import ZipFile
import pandas as pd

INVENTORY_ZIP_SHA="473c9eff45abba0cafa589494973c48af918a0251fcc3d7fb857d8f0b090cc4c"
INVENTORY_CSV_SHA="79dd655c0b5acc0c99365198a91402f8f667fcf7edcd6d552feb7e30e553765d"
ACCEPTED_ZIP_SHA="abcd3aee766c85a92b86a8f93f8ed696bbe49d0f6cf2f64e5f75b457abf78e6a"
ACCEPTED_CSV_SHA="5a77f3c968cf221287cd98a855559b3c3e8a8cf1322c08add28878adeb293df8"

def require(ok,why):
    if not bool(ok):raise ValueError("MOIS_GBIF_ANON_FULL16_AGGREGATION: "+why)
def sha(data):return hashlib.sha256(data).hexdigest()
def file_sha(p):return sha(Path(p).read_bytes())

def signed_inventory(path):
    raw=Path(path).read_bytes()
    require(sha(raw)==INVENTORY_ZIP_SHA,"GBIF full S3 file roster ZIP SHA changed")
    with ZipFile(io.BytesIO(raw)) as z:
        content=z.read("GBIF_2026-11-01_all_anonymous_parquet_objects.csv.gz")
        source_receipt=json.loads(z.read("receipt.json"))
        citation=z.read("citation.txt")
    require(sha(content)==INVENTORY_CSV_SHA,"file roster CSV SHA changed")
    d=pd.read_csv(io.BytesIO(content),compression="gzip",dtype={"key":str})
    require(len(d)==10346 and d.key.is_unique
            and int(d.bytes.sum())==271777586222 and source_receipt["Parquet_file_count"]==10346,
            "full GBIF 2026 November snapshot denominator incorrect")
    require(set(d.partition)==set(range(16)),"source partition IDs incomplete")
    return d,citation

def accepted_species(path):
    raw=Path(path).read_bytes()
    require(sha(raw)==ACCEPTED_ZIP_SHA,"source GBIF COL XR accepted 37k ZIP changed")
    with ZipFile(io.BytesIO(raw)) as z:
        content=z.read("current_exact_accepted_COLXR_species.csv.gz")
    require(sha(content)==ACCEPTED_CSV_SHA,"source GBIF accepted CSV changed")
    d=pd.read_csv(io.BytesIO(content),compression="gzip",dtype=str,keep_default_na=False)
    require(len(d)==37104 and d.key.is_unique and d.inat_taxon_id.is_unique
            and d.status.eq("colxr_exact_accepted_species").all(),"source names failed")
    return d

def confirm_object_coverage(partition_root,source_manifest):
    root=Path(partition_root)
    expected=source_manifest.set_index("key")
    seen={}
    totals=collections.Counter()
    reports=[]
    paths=[]
    for part in range(16):
        # actions/download-artifact artifact name starts with 'mois-GBIF-anon-shard-'
        directory=root/f"mois-GBIF-anon-20261101-shard-{part:02d}"
        receipt_path=directory/"receipt.json"
        require(receipt_path.is_file(),f"missing completed shard receipt {part}")
        r=json.loads(receipt_path.read_text())
        require(r["partition"]==part and r["total_partitions"]==16
                and r["full_partition_scan"] is True
                and r["source_taxa"]==37104
                and r["source_accepted_species_SHA256"]==ACCEPTED_ZIP_SHA
                and r["source_inventory_SHA256"]==INVENTORY_ZIP_SHA
                and r["grid_projection"]=="EPSG:6933" and r["grid_width_m"]==100000,
                f"incompatible or incomplete shard {part}")
        prov_path=directory/"object_provenance.json"
        objs=json.loads(prov_path.read_text())
        require(len(objs)==r["files_actually_scanned"],f"partial shard {part}")
        for obj in objs:
            k=obj["key"]
            require(k in expected.index and int(expected.at[k,"partition"])==part,
                    "unexpected or wrong partition S3 file "+k)
            require(k not in seen,"S3 file listed by multiple shards "+k)
            require(obj["size"]==int(expected.at[k,"bytes"]) and
                    obj["etag"]==str(expected.at[k,"etag"]),
                    "S3 object bytes or ETag differs from pre-scan listing "+k)
            seen[k]=part
        for path in ("species_x_EPSG6933_100km_cells.csv.gz","species_x_countrycode.csv.gz"):
            p=directory/path
            require(p.is_file() and file_sha(p)==r["output_SHA256"][path],
                    "shard data SHA failure")
        for k,v in r["scan_counters"].items():totals[k]+=int(v)
        require(r["scan_counters"]["file_count_scanned"]==len(objs),
                "shard file count not exact")
        reports.append(r)
        paths.append(directory)
    require(set(seen)==set(expected.index) and len(seen)==10346,
            "full GBIF S3 file listing not completely scanned")
    require(totals["file_count_scanned"]==10346,"global shard file sums invalid")
    return reports,dict(totals),paths

def merge(*,inventory,accepted,shards,out):
    import duckdb
    out=Path(out)
    require(not out.exists(),"refuse overwriting global snapshot")
    expected,citation=signed_inventory(inventory)
    taxa=accepted_species(accepted)
    shard_receipts,stats,dirs=confirm_object_coverage(shards,expected)
    out.mkdir(parents=True)
    (out/"citation.txt").write_bytes(citation)
    con=duckdb.connect(str(out/"temporary_merge.duckdb"))
    con.execute("SET memory_limit='3GB'")
    con.execute("SET threads=2")
    con.execute("SET preserve_insertion_order=false")
    # SQL only reads the 16 output shards, not 271GB of cloud source a second time.
    g_glob=str(Path(shards)/"mois-GBIF-anon-20261101-shard-*" /
               "species_x_EPSG6933_100km_cells.csv.gz")
    c_glob=str(Path(shards)/"mois-GBIF-anon-20261101-shard-*" /
               "species_x_countrycode.csv.gz")
    g_out=str(out/"all_37104_species_x_EPSG6933_100km_cells.parquet")
    c_out=str(out/"all_37104_species_x_countrycode.parquet")
    con.execute("""COPY (
        SELECT specieskey,
            CAST(grid_x AS INTEGER) AS grid_x,
            CAST(grid_y AS INTEGER) AS grid_y,
            SUM(CAST(occurrence_count AS BIGINT))::BIGINT AS occurrence_count
        FROM read_csv_auto(?, all_varchar=true, header=true)
        GROUP BY specieskey, grid_x, grid_y
        ) TO ? (FORMAT parquet, COMPRESSION ZSTD)""", [g_glob,g_out])
    con.execute("""COPY (
        SELECT specieskey, countrycode,
            SUM(CAST(occurrence_count AS BIGINT))::BIGINT AS occurrence_count
        FROM read_csv_auto(?, all_varchar=true, header=true)
        GROUP BY specieskey, countrycode
        ) TO ? (FORMAT parquet, COMPRESSION ZSTD)""",[c_glob,c_out])
    summary=con.execute("""SELECT specieskey, COUNT(*) AS occupied_grid_cells,
        SUM(occurrence_count)::BIGINT AS records_passed_spatial_QA
        FROM read_parquet(?) GROUP BY specieskey""",[g_out]).df()
    country=con.execute("""SELECT specieskey,
        COUNT(*) FILTER (WHERE countrycode <> '') AS observed_countries
        FROM read_parquet(?) GROUP BY specieskey""",[c_out]).df()
    require(summary.specieskey.is_unique,"duplicate aggregated accepted taxa")
    require(set(summary.specieskey).issubset(set(taxa.key)),
            "non-accepted GBIF source taxa in global final grid")
    names=taxa[["key","inat_taxon_id","species","source_hue9","hue_status"]].copy()
    result=names.merge(summary,left_on="key",right_on="specieskey",how="left",validate="one_to_one")
    result=result.merge(country,left_on="key",right_on="specieskey",how="left",
                        validate="one_to_one",suffixes=("","_country"))
    require(len(result)==37104 and result.key.is_unique,
            "global taxon metadata disappeared during aggregation")
    for col in ("occupied_grid_cells","records_passed_spatial_QA","observed_countries"):
        result[col]=result[col].fillna(0).astype("int64")
    result["GBIF_record_status"]="source_registered_but_no_accepted_valid_geographic_occurrence"
    result.loc[result.occupied_grid_cells.gt(0),"GBIF_record_status"]="georeferenced_presence_confirmed"
    result["flower_color_at_occurrence_measured"]=False
    result["botanical_petal_hue_validated"]=False
    result["native_species_range_verified"]=False
    result["grid_projection"]="EPSG:6933"
    result["grid_cell_size_m"]=100000
    result.to_csv(out/"all_37104_accepted_species_grid_coverage_and_photo9.csv.gz",
        index=False,compression={"method":"gzip","mtime":0})
    require(int(result.records_passed_spatial_QA.sum())==
            int(stats["passed_strict_taxon_and_geography"]),
            "global filtered record denominator mismatch with shard receipts")
    g_count=int(con.execute("SELECT COUNT(*) FROM read_parquet(?)",[g_out]).fetchone()[0])
    c_count=int(con.execute("SELECT COUNT(*) FROM read_parquet(?)",[c_out]).fetchone()[0])
    con.close()
    try:(out/"temporary_merge.duckdb").unlink()
    except FileNotFoundError:pass
    receipt={"schema":"mois_GBIF_anonymous_S3_Nov2026_complete_37104_taxa_global_grid_v1",
      "GBIF_source_snapshot_date":"2026-11-01",
      "complete_public_S3_object_count":10346,
      "complete_public_S3_compressed_source_bytes":271777586222,
      "all_16_shards_source_verified":True,
      "all_37104_accepted_species_in_roster":True,
      "GBIF_COL_XR_accepted_species":37104,
      "GBIF_accepted_species_with_usable_occurrences":
         int(result.occupied_grid_cells.gt(0).sum()),
      "GBIF_accepted_species_without_usable_occurrences":
         int(result.occupied_grid_cells.eq(0).sum()),
      "GBIF_total_source_occurrence_rows_scanned":
         int(stats["parquet_rows_scanned"]),
      "GBIF_accepted_specieskey_and_binomial_agree":
         int(stats["accepted_key_and_binomial"]),
      "GBIF_rows_rejected_due_to_specieskey_binomial_disagreement":
         int(stats["key_name_disagree"]),
      "GBIF_occurrences_failed_spatial_QA":
         int(stats["failed_geospatial_QA"]),
      "GBIF_occurrences_passed_all_filters":
         int(stats["passed_strict_taxon_and_geography"]),
      "output_unique_species_100km_grid_cells":g_count,
      "output_species_countries":c_count,
      "coordinate_projection":"EPSG:6933",
      "grid_size_m":100000,
      "GBIF_SQL_ISEA3H8_compatible_grid":False,
      "global_distribution_not_native_only":True,
      "photo9_is_not_local_petal_color":True,
      "GBIF_download_API_authentication_used":False,
      "GBIF_AWS_credentials_used":False,
      "GBIF_authenticated_download_key":None,
      "source_snapshots_are_public_cloud_data_with_citation":True,
      "source_inventory_zip_sha256":INVENTORY_ZIP_SHA,
      "source_colxr_accepted_zip_sha256":ACCEPTED_ZIP_SHA,
      "source_shard_receipt_sha256":{str(r["partition"]):sha(
        (dirs[r["partition"]]/"receipt.json").read_bytes()) for r in shard_receipts},
      "source_filters":"COLXR specieskey AND binomial exact; kingdom Plantae; PRESENT; valid coords; exclude GBIF four severe geospatial flags; uncertainty <=50000m if known",
      "output_sha256":{p.name:file_sha(p) for p in out.iterdir() if p.is_file()}}
    (out/"receipt.json").write_text(json.dumps(receipt,sort_keys=True,indent=2)+"\n")
    print("MOIS_GBIF_NOAUTH_37104_GLOBAL_DISTRIBUTION_COMPLETE",json.dumps(
       {k:v for k,v in receipt.items() if k not in ("source_shard_receipt_sha256","output_sha256")},sort_keys=True),flush=True)
    return receipt

if __name__=="__main__":
    p=argparse.ArgumentParser()
    for name in ("inventory","accepted","shards","output"):
        p.add_argument("--"+name,required=True)
    a=p.parse_args()
    merge(inventory=a.inventory,accepted=a.accepted,shards=a.shards,out=a.output)
