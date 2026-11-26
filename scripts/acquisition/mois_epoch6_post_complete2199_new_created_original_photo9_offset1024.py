#!/usr/bin/env python3
"""MOIS epoch6 offset1024: immutable 2,199-taxon hash partition, after signed offset768.

After one previously zero-ROI9 taxon acquired numeric ROI at offset768, reconstruct
the original epoch6 taxon priority frame before slicing 1024:1152. This prevents
the newly numeric taxon from shifting source-search groups.
Historical full ~1.23M raw image SHA disjointness remains unproven.
"""
import argparse
import hashlib
import json
from pathlib import Path
import pandas as pd
from scripts.acquisition import mois_epoch6_post_complete2199_new_created_original_photo9_offset768 as predecessor
from scripts.acquisition import mois_20261115_created_epoch2_disjoint_original_photo9_after2048 as core

m = predecessor.m
NEW_DELTA = "old/epoch6_offset768/new_real_original_nine_photo9_delta.parquet"
NEW_DELTA_SHA = "6d7f7cb748c0ef349d777cf6c417ddffa262de5ee7234ddb7877c1cce2dabe01"
BASE_TAXON_FIRST_NUMERIC = "1249987"
m.SCHEMA = "mois_incremental_20261117_epoch6_offset1024_after_terminal768_1231349_original_photo9_v1"
m.SOURCE_ARM = "MOIS_EPOCH6_OFFSET1024_POST_SIGNED768_FULL_1231349_FROZEN2199_DISJOINT"
m.MIN_CREATED = "2026-11-15T21:27:57+00:00"
m.CREATED_D1 = "2026-11-15"
m.MAX_SOURCE_TAXA = 128
m.ATLAS_SHA = "124672d2073f490ceb26bd7132b234912221f278a7528c4063273743db03f9b4"
m.OLD_PHOTOS = 1231349
m.OLD_NUMERIC_PHOTOS = 467567
m.OLD_NUMERIC_TAXA = 39913
m.OLD_CLEAR_TAXA = 39173

def history_after_offset768():
    m.require(Path(NEW_DELTA).is_file() and m.sha256(NEW_DELTA) == NEW_DELTA_SHA,
              "offset768 signed independent strict9 pixel audit delta absent or changed")
    current = m.OLD_PHOTOS
    try:
        m.OLD_PHOTOS = 1231348
        ids, observations, raw = predecessor.history_after_offset640()
    finally:
        m.OLD_PHOTOS = current
    m.require(len(ids) == 1231348 and len(raw) == 104,
              "offset768 predecessor historical original/UUID/raw SHA counts drift")
    d = pd.read_parquet(NEW_DELTA)
    m.require(len(d) == 1 and {"photo_id","observation_id","image_sha256","strict9",
                                "inat_taxon_id"}.issubset(d.columns),
              "offset768 independent measured pixel delta lacks identity fields")
    row = d.iloc[0]
    photo, obs, img = str(row.photo_id), str(row.observation_id), str(row.image_sha256)
    m.require(photo == "764195408" and
              obs == "60616d1c-a999-453c-93b1-23483a0103dd" and
              img == "a0fa9f9b79189319d48ac27644834d350b7cb8672cd69241c9a4201ca686ca67" and
              str(row.inat_taxon_id) == BASE_TAXON_FIRST_NUMERIC and bool(row.strict9) and
              photo not in ids and obs not in observations and img not in raw,
              "offset768 measured original photo identity, UUID, strict9 or SHA replay")
    fractions = [row["strict9_"+s] for s in m.PALETTE]
    m.require(len(fractions) == 9 and all(pd.notna(x) and float(x) >= 0 for x in fractions)
              and abs(sum(float(x) for x in fractions)-1.0) <= .001,
              "offset768 signed original ROI nine-component numerical vector invalid")
    ids.add(photo); observations.add(obs); raw.add(img)
    m.require(len(ids) == 1231349 and len(raw) == 105,
              "offset1024 complete historical photo ID/UUID/recent raw SHA counts invalid")
    return ids, observations, raw

def select_taxa_epoch6(atlas, n=128, offset=1024):
    missing = atlas.loc[atlas.after_strict9.eq(0)].copy()
    new_numeric = atlas.loc[atlas.inat_taxon_id.eq(BASE_TAXON_FIRST_NUMERIC)].copy()
    m.require(len(atlas) == 42111 and atlas.inat_taxon_id.is_unique and
              len(missing) == 2198 and len(new_numeric) == 1 and
              int(new_numeric.iloc[0].after_strict9) == 1 and
              BASE_TAXON_FIRST_NUMERIC not in set(missing.inat_taxon_id) and
              n == 128 and offset == 1024,
              "epoch6 offset1024 must restore original 2,199 taxon frame before slicing")
    frozen = pd.concat([missing, new_numeric], ignore_index=True)
    m.require(frozen.inat_taxon_id.is_unique and
              frozen.inat_taxon_id.str.fullmatch(r"\d+").all() and len(frozen)==2199,
              "reconstructed epoch6 original taxa not unique or not numeric")
    frozen["priority"] = frozen.inat_taxon_id.map(lambda tid: hashlib.sha256(
        ("MOIS-20261117-EPOCH6-POST-2199-DISJOINT-NEW-CREATED|" + tid).encode()
    ).hexdigest())
    ordered = frozen.sort_values("priority", kind="stable").inat_taxon_id.tolist()
    m.require(BASE_TAXON_FIRST_NUMERIC in ordered[768:896],
              "numerically newly covered taxon was not in the already-searched offset768 cohort")
    selected = ordered[offset:offset+n]
    m.require(len(selected) == len(set(selected)) == 128 and
              not set(selected).intersection(ordered[:offset]) and
              BASE_TAXON_FIRST_NUMERIC not in selected,
              "offset1024 hash cohort overlaps previously queried taxa")
    return selected

m.history = history_after_offset768
m.select_taxa = select_taxa_epoch6
m.source = core.source

def main():
    p=argparse.ArgumentParser()
    stages=p.add_subparsers(dest="stage",required=True)
    for stage,keys in (
        ("source",("atlas","output")),
        ("audit",("source","worker","output")),
        ("append",("atlas","source","worker","audit","output")),
    ):
        q=stages.add_parser(stage)
        for key in keys:q.add_argument("--"+key,required=True)
        if stage=="source":q.add_argument("--offset",type=int,default=1024)
    args=p.parse_args()
    if args.stage=="source":
        m.source(args.atlas,args.output,quota=128,offset=args.offset)
        rpath=Path(args.output)/"receipt.json"
        r=json.loads(rpath.read_bytes())
        m.require((r["hash_disjoint_taxon_offset"],r["searched_zero_numeric_taxa"],
                   r["source_API_error_count"],r["old_original_photo_attempt_ID_count"],
                   r["old_original_photo_obs_ID_overlap"])==(1024,128,0,1231349,0)
                  and r["old_all_image_raw_SHA_disjointness_proven"] is False
                  and r["created_utc_lower_bound"]==m.MIN_CREATED
                  and r["original_pixels_fetched"]==r["strict9_numeric_photos"]==0,
                  "offset1024 original source-only provenance gate did not pass")
        r["recent_known_original_image_raw_SHA_count"]=105
        r["frozen_original_zero_numeric_taxon_frame"]=2199
        r["currently_zero_numeric_taxa"]=2198
        r["new_numeric_taxon_removed_from_sampling_only"]="1249987"
        r.pop("recent_97_nonempty_raw_image_SHA_count",None)
        rpath.write_text(json.dumps(r,sort_keys=True,indent=2)+"\n")
        print("MOIS_EPOCH6_OFFSET1024_SOURCE_SIGNED",json.dumps(r,sort_keys=True),flush=True)
    elif args.stage=="audit":
        m.audit(args.source,args.worker,args.output)
    else:
        m.append(args.atlas,args.source,args.worker,args.audit,args.output)

if __name__=="__main__":
    main()
