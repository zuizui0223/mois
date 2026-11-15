#!/usr/bin/env python3
"""Streaming, anonymous GBIF S3 snapshot -> MOIS COL XR species x EPSG:6933 100km grid.

Only officially published unsigned GBIF monthly Parquet and source-signed
37,104 accepted species enter. No credential, no database server, no downloads
of entire raw snapshots to disk. Each partition owns a disjoint object set.
Snapshot specieskeys must also match accepted binomials for STRICT matches.
Do not report a subset scan as complete global distributions.
"""
from __future__ import annotations
import argparse,collections,hashlib,io,json,math,re
from pathlib import Path
from zipfile import ZipFile
import pandas as pd
import numpy as np

SOURCE_SHA="abcd3aee766c85a92b86a8f93f8ed696bbe49d0f6cf2f64e5f75b457abf78e6a"
SOURCE_CSV_SHA="5a77f3c968cf221287cd98a855559b3c3e8a8cf1322c08add28878adeb293df8"
MANIFEST_SHA="473c9eff45abba0cafa589494973c48af918a0251fcc3d7fb857d8f0b090cc4c"
MANIFEST_CSV_SHA="79dd655c0b5acc0c99365198a91402f8f667fcf7edcd6d552feb7e30e553765d"
BUCKET="gbif-open-data-us-east-1"
BASE="occurrence/2026-11-01/occurrence.parquet/"
BAD_ISSUES={"ZERO_COORDINATE","COUNTRY_COORDINATE_MISMATCH",
            "COORDINATE_INVALID","COORDINATE_OUT_OF_RANGE"}
GRID_METRES=100000

def require(x,why):
    if not bool(x):raise ValueError("MOIS_ANON_GBIF_GRID: "+why)
def sha(x):return hashlib.sha256(x).hexdigest()
def canonical_pair(s):
    parts=str(s).strip().split()
    if len(parts)<2:return None
    if not (parts[0][:1].isupper() and re.fullmatch(r'[a-z]+(?:-[a-z]+)*',parts[1])):return None
    return parts[0]+" "+parts[1]

def flags(x):
    # Parquet's issue is list of nested structs; choose conservative recursion,
    # not a substring guess based on Python repr.
    found=[]
    if x is None:return found
    if isinstance(x,str):return [x]
    if isinstance(x,dict):
        for v in x.values():found.extend(flags(v))
        return found
    if isinstance(x,(list,tuple,np.ndarray)):
        for y in x:found.extend(flags(y))
        return found
    return found

def allowed_coords(lon,lat,unc,issue):
    if not (np.isfinite(lon) and np.isfinite(lat) and -180<=lon<=180
            and -90<=lat<=90 and (lon!=0 or lat!=0)):return False
    if np.isfinite(unc) and (unc<0 or unc>50000):return False
    return not bool(set(flags(issue))&BAD_ISSUES)

def read_signed_inputs(source,manifest,partition,limit=None):
    require(0<=partition<16,"invalid source partition")
    require(limit is None or limit>0,"invalid bounded object count")
    b=Path(source).read_bytes();require(sha(b)==SOURCE_SHA,"accepted COL source ZIP changed")
    with ZipFile(io.BytesIO(b)) as z:csv=z.read("current_exact_accepted_COLXR_species.csv.gz")
    require(sha(csv)==SOURCE_CSV_SHA,"accepted COL source roster SHA changed")
    taxa=pd.read_csv(io.BytesIO(csv),compression="gzip",dtype=str,keep_default_na=False)
    require(len(taxa)==37104 and taxa.key.is_unique
            and taxa.status.eq("colxr_exact_accepted_species").all(),"37,104 COL XR taxa absent")
    b=Path(manifest).read_bytes();require(sha(b)==MANIFEST_SHA,"S3 inventory artifact SHA changed")
    with ZipFile(io.BytesIO(b)) as z:
        csv=z.read("GBIF_2026-11-01_all_anonymous_parquet_objects.csv.gz")
        receipt=json.loads(z.read("receipt.json"))
    require(sha(csv)==MANIFEST_CSV_SHA,"source object roster changed")
    objs=pd.read_csv(io.BytesIO(csv),compression="gzip",dtype={"key":str})
    require(len(objs)==receipt["Parquet_file_count"]==10346 and objs.key.is_unique
            and int(objs.bytes.sum())==receipt["total_Parquet_compressed_bytes"],
            "full snapshot file denominator changed")
    require(set(objs.partition)==set(range(16)),"source file partitions incomplete")
    selected=objs.loc[objs.partition.eq(partition)].sort_values("key")
    if limit is not None:selected=selected.head(limit)
    require(len(selected)>0,"empty shard")
    return taxa,selected,receipt

def accumulate_arrow_batch(batch, names_by_key, proj, stats, grid, country):
    import pyarrow as pa
    import pyarrow.compute as pc
    if batch.num_rows==0:return
    # Project source keys first; avoid materializing 100% of all GBIF rows in pandas.
    keyset=pa.array(list(names_by_key),type=pa.string())
    keymask=pc.is_in(batch.column("specieskey"),value_set=keyset)
    matches=batch.filter(keymask)
    stats["source_key_hits"]+=matches.num_rows
    if matches.num_rows==0:return
    frame=matches.to_pandas()
    # Missing kingdom/status is never equivalent to accepted source scope.
    mask=frame.kingdom.eq("Plantae") & frame.occurrencestatus.eq("PRESENT")
    frame=frame.loc[mask]
    stats["plant_present_key_hits"]+=len(frame)
    if frame.empty:return
    srcname=frame.specieskey.map(names_by_key)
    canon=frame.species.astype("string").map(canonical_pair)
    match=(srcname==canon)
    stats["key_name_disagree"]+=int((~match).sum())
    frame=frame.loc[match].copy()
    stats["accepted_key_and_binomial"]+=len(frame)
    if frame.empty:return
    lon=pd.to_numeric(frame.decimallongitude,errors="coerce").to_numpy(float)
    lat=pd.to_numeric(frame.decimallatitude,errors="coerce").to_numpy(float)
    unc=pd.to_numeric(frame.coordinateuncertaintyinmeters,errors="coerce").to_numpy(float)
    srcissue=frame.issue.to_list()
    good=np.array([allowed_coords(float(lo),float(la),float(un),iss)
                   for lo,la,un,iss in zip(lon,lat,unc,srcissue)],dtype=bool)
    stats["failed_geospatial_QA"]+=int((~good).sum())
    frame=frame.iloc[np.flatnonzero(good)]
    if frame.empty:return
    xs,ys=proj.transform(lon[good],lat[good])
    xyok=np.isfinite(xs)&np.isfinite(ys)
    stats["failed_projection"]+=int((~xyok).sum())
    frame=frame.iloc[np.flatnonzero(xyok)]
    gx=np.floor(np.asarray(xs)[xyok]/GRID_METRES).astype(np.int32)
    gy=np.floor(np.asarray(ys)[xyok]/GRID_METRES).astype(np.int32)
    stats["passed_strict_taxon_and_geography"]+=len(frame)
    mini=pd.DataFrame({"specieskey":frame.specieskey.to_numpy(),
                       "cell_x":gx,"cell_y":gy,"country":frame.countrycode.fillna("").astype(str).to_numpy()})
    group=mini.groupby(["specieskey","cell_x","cell_y"],sort=False).size()
    grid.update({k:int(n) for k,n in group.items()})
    country_group=mini.groupby(["specieskey","country"],sort=False).size()
    country.update({k:int(n) for k,n in country_group.items()})

def scan(source,manifest,partition,out,limit=None):
    import pyarrow.parquet as pq
    from pyarrow import fs
    from pyproj import Transformer
    out=Path(out);require(not out.exists(),"refuse to overwrite old scan")
    taxa,objects,full_receipt=read_signed_inputs(source,manifest,partition,limit)
    names_by_key={str(k):canonical_pair(name) for k,name in zip(taxa.key,taxa.species)}
    require(all(names_by_key.values()),"non-binomial source accepted taxa")
    proj=Transformer.from_crs("EPSG:4326","EPSG:6933",always_xy=True)
    client=fs.S3FileSystem(region="us-east-1",anonymous=True)
    wanted=["specieskey","species","kingdom","occurrencestatus",
            "decimallatitude","decimallongitude","coordinateuncertaintyinmeters",
            "issue","countrycode"]
    grid=collections.Counter();country=collections.Counter();stats=collections.Counter()
    objects_read=[]
    for j,row in enumerate(objects.itertuples(index=False),1):
        require(str(row.key).startswith(BASE),"source file not in pinned snapshot")
        with client.open_input_file(f"{BUCKET}/{row.key}") as f:
            pqfile=pq.ParquetFile(f)
            require(set(wanted).issubset(pqfile.schema_arrow.names),
                    "required GBIF Parquet field missing")
            stats["parquet_rows_scanned"]+=pqfile.metadata.num_rows
            for batch in pqfile.iter_batches(batch_size=131072,columns=wanted):
                accumulate_arrow_batch(batch,names_by_key,proj,stats,grid,country)
        objects_read.append({"key":row.key,"size":int(row.bytes),"etag":str(row.etag)})
        stats["file_count_scanned"]+=1
        if j%10==0 or j==len(objects):
            print("MOIS_GBIF_ANON_S3_SHARD",partition,"file",j,"of",len(objects),
                  "accepted_match_rows",stats["passed_strict_taxon_and_geography"],flush=True)
    out.mkdir(parents=True)
    a=pd.DataFrame([{"specieskey":k[0],"grid_x":k[1],"grid_y":k[2],"occurrence_count":n}
                    for k,n in grid.items()],columns=["specieskey","grid_x","grid_y","occurrence_count"])
    c=pd.DataFrame([{"specieskey":k[0],"countrycode":k[1],"occurrence_count":n}
                    for k,n in country.items()],columns=["specieskey","countrycode","occurrence_count"])
    a.to_csv(out/"species_x_EPSG6933_100km_cells.csv.gz",index=False,
       compression={"method":"gzip","mtime":0})
    c.to_csv(out/"species_x_countrycode.csv.gz",index=False,
       compression={"method":"gzip","mtime":0})
    (out/"object_provenance.json").write_text(json.dumps(objects_read,indent=2)+"\n")
    result={"schema":"mois_anonymous_GBIF_20261101_37k_species_equal_area_grid_shard_v1",
         "source_snapshot":"2026-11-01","source_storage":"S3 public anonymous",
         "partition":partition,"total_partitions":16,
         "files_in_this_partition":int((full_receipt["Parquet_file_count"]+15)//16) if limit is None else len(objects),
         "files_actually_scanned":len(objects),
         "full_partition_scan":limit is None,
         "entire_snapshot_objects":full_receipt["Parquet_file_count"],
         "entire_snapshot_compressed_bytes":full_receipt["total_Parquet_compressed_bytes"],
         "source_taxa":len(taxa),"source_taxonomy":"COL XR exact accepted specieskey+binomial cross-check",
         "grid_projection":"EPSG:6933","grid_width_m":GRID_METRES,
         "GBIF_SQL_ISEA3H_grid_equivalent":False,
         "grid_cell_pairs_in_this_partition":len(a),
         "country_pairs_in_this_partition":len(c),
         "scan_counters":dict(stats),
         "one_shard_not_global_distribution":True,
         "public_S3_signature_authentication_used":False,
         "GBIF_account_used":False,
         "GBIF_S3_object_ETags_and_lengths_retained":True,
         "rejected_geospatial_flags":sorted(BAD_ISSUES),
         "maximum_coordinate_uncertainty_m":50000,
         "source_inventory_SHA256":MANIFEST_SHA,
         "source_accepted_species_SHA256":SOURCE_SHA,
         "output_SHA256":{x.name:sha(x.read_bytes()) for x in out.iterdir()}}
    (out/"receipt.json").write_text(json.dumps(result,sort_keys=True,indent=2)+"\n")
    print("MOIS_GBIF_ANON_PARTITION_REAL_SPATIAL_ANALYSIS_SUCCESS",
          json.dumps({k:v for k,v in result.items() if k!="output_SHA256"},sort_keys=True),flush=True)
    return result

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--accepted",required=True);p.add_argument("--manifest",required=True)
    p.add_argument("--partition",type=int,required=True)
    p.add_argument("--max-files",type=int)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    scan(a.accepted,a.manifest,a.partition,a.output,a.max_files)
