#!/usr/bin/env python3
"""Restore real observation date/lat/lon of all 128 original-source photo shards.

The 1,181,146-row strict original ROI9 union was photo-ID / taxon accurate,
but the 131,072-photo bulk128 supplement dropped the original 128-shard
latitude, longitude and observed_on metadata during its numerical join.
Read exactly the already published SHA-verified 512 ORIGINAL worker CSVs,
NEVER photo images, and correct ONLY those three missing metadata columns.
New release (v2) must not overwrite source-frozen 9-colour baseline (v1).
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,tarfile
from pathlib import Path

import numpy as np
import pandas as pd

WORLD_SHA="3fbd0f0ec7eb347b3bbc7ab509ec9312ec6d797ed932dcf7bb425959898c3835"
SOURCE_MANIFEST_SHA="7e975e701281ba812451acd48be7e91a7d17f2bbe6f9f53ad15ef4f8db9e557e"
SOURCE_TAR_SHA="61c0eee377f6edead471436d3f75b92700e54687e6b9dd0a417819a472c23dcb"
SCHEMA="mois_signed1181146_trueROI9_bulk128_original_metadata_georestore_v2"
BULK="bulk128_131072"
N_BULK=131072
N_WORLD=1181146
N_STRICT=451643
REFERENCE="mois_bulk128_tier10_photo_reference_v1"
TARGETED="mois_bulk128_tier10_taxon_region_season_gap_v1"
FIELDS=("photo_id","observation_id","inat_taxon_id","observed_on",
        "latitude","longitude","source_arm")

def require(ok,msg):
    if not bool(ok):raise ValueError("MOIS_SOURCE_BULK128_GEORESTORE: "+msg)

def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda:f.read(4<<20),b""):h.update(b)
    return h.hexdigest()

def source_photo_roster(manifest,source_tar,*,unit_count=512):
    manifest=Path(manifest);source_tar=Path(source_tar)
    spec=json.loads(manifest.read_text())
    require(spec["schema"]=="mois_full_source_observation_sampling_cohort_v1"
            and spec["sampling_schema"]==
            "mois_original93m_tier10_bulk128_probability_reference_plus_gap_v1",
            "Not immutable original bulk128 signed photo metadata")
    require(spec["worker_ids"]==list(range(unit_count))
            and len(spec["photo_shard_sha256"])==unit_count,
            "Signed source photo worker IDs missing, unordered or duplicated")
    seen_photos=set();seen_obs=set();out=[]
    with tarfile.open(source_tar,"r:gz") as tar:
        members={m.name.removeprefix("./"):m for m in tar.getmembers()
                 if m.isfile()}
        for unit in range(unit_count):
            name=f"shards/shard-{unit:03d}.csv"
            require(name in members,"Missing original metadata photo shard: "+name)
            blob=tar.extractfile(members[name]).read()
            require(hashlib.sha256(blob).hexdigest()==
                    spec["photo_shard_sha256"][str(unit)],
                    "Original selected source 256-photo shard bytes changed")
            rows=list(csv.DictReader(io.StringIO(blob.decode("utf-8"))))
            require(0<len(rows)<=256,"Missing original photo identity rows")
            for r in rows:
                pid,obs,taxon=(r[k] for k in ("photo_id","observation_id","inat_taxon_id"))
                require(pid.isdecimal() and taxon.isdecimal() and obs
                        and pid not in seen_photos and obs not in seen_obs
                        and r["source_arm"] in (REFERENCE,TARGETED),
                        "Signed original photo/observation/taxon/arm invalid or duplicated")
                seen_photos.add(pid);seen_obs.add(obs)
                out.append({k:r[k] for k in FIELDS})
    require(len(out)==len(seen_photos)==len(seen_obs),
            "Original image source photo and observation IDs not unique")
    return pd.DataFrame(out),spec

def restore(frame,source,*,actual_n=None,strict_n=None):
    frame=frame.copy()
    source=source.copy()
    expected=len(source) if actual_n is None else actual_n
    bulk=frame.source_design.eq(BULK)
    require(int(bulk.sum())==expected and
            frame.photo_id.astype(str).is_unique and
            source.photo_id.astype(str).is_unique and
            source.observation_id.astype(str).is_unique,
            "Missing exact source selected bulk photo IDs or duplicates")
    if strict_n is not None:
        require(int(frame.strict9.sum())==strict_n,
                "Genuine original FCP nine-colour photo count changed")
    cols=["photo_id","observation_id","inat_taxon_id",
          "source_subarm","observed_on","latitude","longitude"]
    selected=frame.loc[bulk,cols].copy()
    selected["_row_index"]=selected.index
    for k in ("photo_id","observation_id","inat_taxon_id"):
        selected[k]=selected[k].fillna("").astype(str)
        source[k]=source[k].fillna("").astype(str)
    renames={
        "observation_id":"source_observation_id",
        "inat_taxon_id":"source_taxon_id",
        "observed_on":"source_observed_on",
        "latitude":"source_latitude",
        "longitude":"source_longitude",
        "source_arm":"source_source_arm",
    }
    linked=selected.merge(source.rename(columns=renames),
                          on="photo_id",how="left",validate="one_to_one",
                          sort=False)
    require(len(linked)==expected
            and linked.photo_id.is_unique
            and linked.observation_id.eq(linked.source_observation_id).all()
            and linked.inat_taxon_id.eq(linked.source_taxon_id).all(),
            "Original photo IDs have inconsistent source observation or plant taxon")
    wanted=linked.source_subarm.map({"reference":REFERENCE,"targeted":TARGETED})
    require(wanted.notna().all() and wanted.eq(linked.source_source_arm).all(),
            "Source probability-reference / gap-targeted role changed after image inference")
    a=pd.to_numeric(linked.source_latitude,errors="coerce")
    b=pd.to_numeric(linked.source_longitude,errors="coerce")
    correct_point=a.between(-90,90)&b.between(-180,180)
    new_lat=a.where(correct_point,np.nan).astype("float32")
    new_lon=b.where(correct_point,np.nan).astype("float32")
    old_lat=pd.to_numeric(linked.latitude,errors="coerce")
    old_lon=pd.to_numeric(linked.longitude,errors="coerce")
    old_point=old_lat.between(-90,90)&old_lon.between(-180,180)
    mismatch=old_point & (~correct_point|
               ((old_lat-new_lat).abs()>1e-3)|
               ((old_lon-new_lon).abs()>1e-3))
    require(not mismatch.any(),
            "Existing original geographic point differs from source photo metadata")
    old_dates=linked.observed_on.fillna("").astype(str).str.strip()
    source_dates=linked.source_observed_on.fillna("").astype(str).str.strip()
    require((old_dates.eq("")|old_dates.eq(source_dates)).all(),
            "Existing source observation date would be changed to conflicting date")
    idx=linked["_row_index"].to_numpy()
    frame.loc[idx,"latitude"]=new_lat.to_numpy(dtype="float32")
    frame.loc[idx,"longitude"]=new_lon.to_numpy(dtype="float32")
    frame.loc[idx,"observed_on"]=source_dates.to_numpy()
    report={
        "new_bulk_original_photo_IDs_restored":len(linked),
        "prior_bulk_source_geopoints_known":int(old_point.sum()),
        "source_original_bulk_photo_geopoints_available":int(correct_point.sum()),
        "bulk_photos_still_without_valid_coords":int((~correct_point).sum()),
        "bulk_original_photo_observation_days_nonempty":int(source_dates.ne("").sum()),
        "bulk_source_original_photo_IDs_matching_obs_and_taxon":len(linked),
        "bulk_source_reference_photos":int(
            linked.source_subarm.eq("reference").sum()),
        "bulk_source_targeted_photos":int(
            linked.source_subarm.eq("targeted").sum()),
        "old_existing_other_photo_rows_modified":0,
        "original_nine_palette_photos_modified":0,
        "original_photo_IDs_modified":0,
    }
    require(report["bulk_source_reference_photos"]+
            report["bulk_source_targeted_photos"]==len(linked),
            "Original bulk reference/targeted source counts changed")
    return frame,report

def run(world,manifest,source_tar,out,*,verify_sha=True,
        expected_world=N_WORLD,expected_bulk=N_BULK,expected_units=512):
    world,manifest,source_tar,out=map(Path,(world,manifest,source_tar,out))
    require(not out.exists(),"Never replace signed immutable photo/colour baseline")
    if verify_sha:
        require(sha(world)==WORLD_SHA and sha(manifest)==SOURCE_MANIFEST_SHA
                and sha(source_tar)==SOURCE_TAR_SHA,
                "Source original photo data or source metadata SHA changed")
    d=pd.read_parquet(world)
    source,spec=source_photo_roster(manifest,source_tar,unit_count=expected_units)
    require(len(d)==expected_world and len(source)==expected_bulk,
            "Original exact 1.18m world source or 131072 selected photos truncated")
    fixed,details=restore(
        d,source,actual_n=expected_bulk,strict_n=N_STRICT if verify_sha else None)
    require(len(fixed)==len(d) and
            int(fixed.strict9.sum())==int(d.strict9.sum()) and
            fixed.photo_id.astype(str).equals(d.photo_id.astype(str)),
            "Scientific photo-ID or nine-colour outcomes mutated during metadata-only restoration")
    out.mkdir(parents=True)
    target=out/"mois_1181146_true_original_strict9_geography_restored_v2.parquet"
    fixed.to_parquet(target,index=False,compression="zstd")
    receipt={
        "schema":SCHEMA,
        "source_original_all_photo_rows":len(d),
        "real_original_strict9_numerical_images":int(d.strict9.sum()),
        "photo_source_metadata_shard_count":expected_units,
        **details,
        "source_georestored_original_parquet_SHA256":sha(target),
        "input_world_original_parquet_SHA256":sha(world),
        "signed_tier10_source_worker_manifest_SHA256":sha(manifest),
        "signed_tier10_source_worker_archive_SHA256":sha(source_tar),
        "raw_original_image_downloads":0,
        "full_original_RGB_or_EfficientSAM_recomputed":False,
        "green_brown_unknown_or_flower_stage_relabelled":False,
        "same_original_image_ID_and_original_9palette_unmodified":True,
        "not_unbiased_wild_native_flora_sample":True
    }
    (out/"mois_1181146_true_original_source_georestoration_receipt.json").write_text(
        json.dumps(receipt,sort_keys=True,indent=2,allow_nan=False)+"\n")
    print("MOIS_SIGNED_1181146_ORIGINAL_131072_GEOGRAPHY_RESTORED_NO_PIXEL_FETCH",
          json.dumps({k:receipt[k] for k in (
            "new_bulk_original_photo_IDs_restored",
            "prior_bulk_source_geopoints_known",
            "source_original_bulk_photo_geopoints_available",
            "bulk_original_photo_observation_days_nonempty",
            "bulk_photos_still_without_valid_coords"
          )},sort_keys=True),flush=True)
    return receipt

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for arg in ("world-parquet","source-manifest","source-tar","output"):
        p.add_argument("--"+arg,type=Path,required=True)
    a=p.parse_args()
    run(a.world_parquet,a.source_manifest,a.source_tar,a.output)

if __name__=="__main__":
    main()
