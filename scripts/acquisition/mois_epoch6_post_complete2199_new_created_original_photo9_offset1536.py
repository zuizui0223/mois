#!/usr/bin/env python3
"""MOIS epoch6 offset1536: immutable 2,199-taxon hash partition, after signed offset1152.

After one previously zero-ROI9 taxon acquired numeric ROI at offset768, reconstruct
the original epoch6 taxon priority frame before slicing 1536:1664. This prevents
the newly numeric taxon from shifting source-search groups.
Historical full ~1.23M raw image SHA disjointness remains unproven.
"""
import argparse
import hashlib
import json
from pathlib import Path
import pandas as pd
from scripts.acquisition import mois_epoch6_post_complete2199_new_created_original_photo9_offset1152 as predecessor
from scripts.acquisition import mois_20261115_created_epoch2_disjoint_original_photo9_after2048 as core

m = predecessor.m
NEW_DELTA = "old/epoch6_offset1152/new_real_original_nine_photo9_delta.parquet"
NEW_DELTA_SHA = "0e10a278f9879ed810beb8dcb7def86f732e43faaacb8282f691cd4bccce25fd"
BASE_TAXON_FIRST_NUMERIC = "1249987"
m.SCHEMA = "mois_incremental_20261117_epoch6_offset1536_after_terminal1152_1231350_original_photo9_v1"
m.SOURCE_ARM = "MOIS_EPOCH6_OFFSET1536_POST_SIGNED1152_FULL_1231350_FROZEN2199_DISJOINT"
m.MIN_CREATED = "2026-11-15T21:27:57+00:00"
m.CREATED_D1 = "2026-11-15"
m.MAX_SOURCE_TAXA = 128
m.ATLAS_SHA = "767b59274d66bd292584af1b338d1828ba6e4db8910c545b798e5b3bd9278d5b"
m.OLD_PHOTOS = 1231350
m.OLD_NUMERIC_PHOTOS = 467567
m.OLD_NUMERIC_TAXA = 39913
m.OLD_CLEAR_TAXA = 39173

def history_after_offset1152():
    m.require(Path(NEW_DELTA).is_file() and m.sha256(NEW_DELTA) == NEW_DELTA_SHA,
              "offset1152 signed independent strict9 pixel audit delta absent or changed")
    current = m.OLD_PHOTOS
    try:
        m.OLD_PHOTOS = 1231349
        ids, observations, raw = predecessor.history_after_offset768()
    finally:
        m.OLD_PHOTOS = current
    m.require(len(ids) == 1231349 and len(raw) == 105,
              "offset1152 predecessor historical original/UUID/raw SHA counts drift")
    d = pd.read_parquet(NEW_DELTA)
    m.require(len(d) == 1 and {"photo_id","observation_id","image_sha256","strict9",
                                "inat_taxon_id"}.issubset(d.columns),
              "offset768 independent measured pixel delta lacks identity fields")
    row = d.iloc[0]
    photo, obs, img = str(row.photo_id), str(row.observation_id), str(row.image_sha256)
    m.require(photo == "764423521" and
              obs == "b14f6504-d32a-485a-8c62-c0fd5cd52109" and
              img == "75dd574f0ca14f2767a38a86030f5fc22ec6f78157d0f6504cb28d1a05c8104d" and
              str(row.inat_taxon_id) == "708260" and not bool(row.strict9) and
              photo not in ids and obs not in observations and img not in raw,
              "offset1152 original photo identity, UUID, technical unknown or SHA replay")
    fractions = [row["strict9_"+s] for s in m.PALETTE]
    m.require(len(fractions) == 9 and all(pd.isna(x) for x in fractions),
              "offset1152 technical unknown must keep all nine fractions missing")
    ids.add(photo); observations.add(obs); raw.add(img)
    m.require(len(ids) == 1231350 and len(raw) == 106,
              "offset1536 complete historical photo ID/UUID/recent raw SHA counts invalid")
    return ids, observations, raw

def select_taxa_epoch6(atlas, n=128, offset=1536):
    missing = atlas.loc[atlas.after_strict9.eq(0)].copy()
    new_numeric = atlas.loc[atlas.inat_taxon_id.eq(BASE_TAXON_FIRST_NUMERIC)].copy()
    m.require(len(atlas) == 42111 and atlas.inat_taxon_id.is_unique and
              len(missing) == 2198 and len(new_numeric) == 1 and
              int(new_numeric.iloc[0].after_strict9) == 1 and
              BASE_TAXON_FIRST_NUMERIC not in set(missing.inat_taxon_id) and
              n == 128 and offset == 1536,
              "epoch6 offset1536 must restore original 2,199 taxon frame before slicing")
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
              "offset1536 hash cohort overlaps previously queried taxa")
    return selected

m.history = history_after_offset1152
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
        if stage=="source":q.add_argument("--offset",type=int,default=1536)
    args=p.parse_args()
    if args.stage=="source":
        m.source(args.atlas,args.output,quota=128,offset=args.offset)
        rpath=Path(args.output)/"receipt.json"
        r=json.loads(rpath.read_bytes())
        m.require((r["hash_disjoint_taxon_offset"],r["searched_zero_numeric_taxa"],
                   r["source_API_error_count"],r["old_original_photo_attempt_ID_count"],
                   r["old_original_photo_obs_ID_overlap"])==(1536,128,0,1231350,0)
                  and r["old_all_image_raw_SHA_disjointness_proven"] is False
                  and r["created_utc_lower_bound"]==m.MIN_CREATED
                  and r["original_pixels_fetched"]==r["strict9_numeric_photos"]==0,
                  "offset1152 original source-only provenance gate did not pass")
        r["recent_known_original_image_raw_SHA_count"]=106
        r["frozen_original_zero_numeric_taxon_frame"]=2199
        r["currently_zero_numeric_taxa"]=2198
        r["new_numeric_taxon_removed_from_sampling_only"]="1249987"
        r.pop("recent_97_nonempty_raw_image_SHA_count",None)
        rpath.write_text(json.dumps(r,sort_keys=True,indent=2)+"\n")
        print("MOIS_EPOCH6_OFFSET1536_SOURCE_SIGNED",json.dumps(r,sort_keys=True),flush=True)
    elif args.stage=="audit":
        m.audit(args.source,args.worker,args.output)
    else:
        m.append(args.atlas,args.source,args.worker,args.audit,args.output)

if __name__=="__main__":
    main()
