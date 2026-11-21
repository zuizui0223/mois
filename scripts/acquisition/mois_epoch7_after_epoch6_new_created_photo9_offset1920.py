#!/usr/bin/env python3
"""MOIS epoch7: new-created interval after signed complete epoch6.

This is a new temporal source interval, not a replay of any epoch6 source SHA.
The source selects 128 taxon IDs from the 2,197 still-zero-strict9 taxa in a
frozen epoch7 hash order; previously queried taxa are eligible only for
observations created after the epoch6 final source receipt. All older photo
IDs/observation UUIDs are excluded. Raw image SHA disjointness against the
full historical ~1.23 million images is NOT certified.
"""
import argparse
import hashlib
import json
from pathlib import Path
import pandas as pd

from scripts.acquisition import (
    mois_epoch6_post_complete2199_new_created_original_photo9_offset2176 as predecessor
)
from scripts.acquisition import mois_20261115_created_epoch2_disjoint_original_photo9_after2048 as core

from scripts.acquisition import mois_epoch7_after_epoch6_new_created_photo9_offset1536 as prior_epoch7
m = prior_epoch7.m
m.SCHEMA = "mois_incremental_20261117_epoch7_offset1920_new_time_window_1231354_original9_v1"
m.SOURCE_ARM = "MOIS_EPOCH7_OFFSET1920_AFTER_TERMINAL_EPOCH6_2199_NEW_CREATED_INTERVAL"
# Epoch6 final source artifact was frozen at 2026-11-17T11:47:05Z.
m.MIN_CREATED = "2026-11-17T11:47:06+00:00"
m.CREATED_D1 = "2026-11-17"
m.MAX_SOURCE_TAXA = 128
m.ATLAS_SHA = "22189c1746989e57f4ab6da1494d401c8a59feb841ed3242a31a8428bed93a41"
m.OLD_PHOTOS = 1231354
m.OLD_NUMERIC_PHOTOS = 467568
m.OLD_NUMERIC_TAXA = 39914
m.OLD_CLEAR_TAXA = 39174

def select_taxa_epoch7(atlas, n=128, offset=1920):
    m.require(len(atlas) == 42111 and atlas.inat_taxon_id.is_unique,
              "nominal atlas must contain exactly 42,111 distinct taxa")
    z = atlas.loc[atlas.after_strict9.eq(0)].copy()
    m.require(len(z) == 2197 and z.inat_taxon_id.str.fullmatch(r"\d+").all(),
              "epoch7 must select only the 2,197 presently unmeasured taxa")
    z["priority"] = z.inat_taxon_id.map(
        lambda tid: hashlib.sha256(
            ("MOIS-20261117-EPOCH7-AFTER-EPOCH6-TERMINAL-TEMPORAL-NEW-CREATED|" + tid).encode()
        ).hexdigest())
    ordered = z.sort_values("priority", kind="stable").inat_taxon_id.tolist()
    selected = ordered[offset:offset+n]
    m.require(offset == 1920 and n == 128
              and len(selected) == len(set(selected)) == n
              and not set(selected).intersection(ordered[:offset]),
              "epoch7 fifteenth disjoint taxon cohort changed")
    return selected

# Existing signed history includes the 1,231,352 photo-ID/UUID ledger, the
# independently audited epoch6 offset1920 technical unknown, and 108 recent
# actual raw-image hashes. The zero-source offset2176 added no photo IDs.

NEW_DELTA = "old/epoch7_offset1536/new_real_original_nine_photo9_delta.parquet"
NEW_DELTA_SHA = "c653d80677f10947a361ca61584450f68806392322cd9b3fcd13cc7dd00eada1"

def history_after_epoch7_offset1536():
    m.require(Path(NEW_DELTA).is_file() and m.sha256(NEW_DELTA) == NEW_DELTA_SHA,
              "epoch7 offset1536 signed independent delta absent or byte changed")
    current = m.OLD_PHOTOS
    try:
        m.OLD_PHOTOS = 1231353
        ids, observations, raw = prior_epoch7.history_after_epoch7_offset768()
    finally:
        m.OLD_PHOTOS = current
    m.require(len(ids) == 1231353 and len(raw) == 109,
              "prior 1,231,353-ID and 109 known-raw-SHA history drift")
    d = pd.read_parquet(NEW_DELTA)
    m.require(len(d) == 1 and {"photo_id", "observation_id", "image_sha256",
                               "strict9", "inat_taxon_id"}.issubset(d.columns),
              "epoch7 offset1536 independent audit missing columns")
    row = d.iloc[0]
    photo, obs, img = str(row.photo_id), str(row.observation_id), str(row.image_sha256)
    m.require(photo == "764630362" and
              obs == "c430799a-c8f3-43a7-a38c-a5a9c91781fa" and
              img == "dc0379d51440231ba01eb55b9799c130504956f3f6d3d679699e504510a538b9" and
              str(row.inat_taxon_id) == "1065359" and not bool(row.strict9) and
              photo not in ids and obs not in observations and img not in raw,
              "epoch7 offset1536 identity/UUID/recent raw SHA collision or unexpected numeric colour")
    m.require(all(pd.isna(row["strict9_"+p]) for p in m.PALETTE),
              "epoch7 offset1536 technical unknown must be NA in all nine fractions")
    ids.add(photo); observations.add(obs); raw.add(img)
    m.require(len(ids) == 1231354 and len(raw) == 110,
              "full photo-ID / recent-image-SHA census failed after offset1536")
    return ids, observations, raw

m.history = history_after_epoch7_offset1536
m.select_taxa = select_taxa_epoch7
m.source = core.source

def main():
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="stage", required=True)
    for stage, keys in (
        ("source", ("atlas", "output")),
        ("audit", ("source", "worker", "output")),
        ("append", ("atlas", "source", "worker", "audit", "output")),
    ):
        q = sub.add_parser(stage)
        for key in keys:
            q.add_argument("--" + key, required=True)
        if stage == "source":
            q.add_argument("--offset", type=int, default=1920)
    args = parser.parse_args()
    if args.stage == "source":
        m.source(args.atlas, args.output, quota=128, offset=args.offset)
        receipt_path = Path(args.output) / "receipt.json"
        r = json.loads(receipt_path.read_bytes())
        m.require(
            (r["hash_disjoint_taxon_offset"], r["searched_zero_numeric_taxa"],
             r["source_API_error_count"], r["old_original_photo_attempt_ID_count"],
             r["old_original_photo_obs_ID_overlap"]) == (1920, 128, 0, 1231354, 0)
            and r["old_all_image_raw_SHA_disjointness_proven"] is False
            and r["created_utc_lower_bound"] == m.MIN_CREATED
            and r["original_pixels_fetched"] == r["strict9_numeric_photos"] == 0,
            "new temporal source-only receipt gate failed")
        r["recent_known_original_image_raw_SHA_count"] = 110
        r["prior_signed_epoch6_taxon_frame_complete"] = 2199
        r["currently_zero_numeric_taxa"] = 2197
        r["previous_epoch6_terminal_run"] = 44076781021
        r["previous_epoch7_offset0_run"] = 44084420377
        r["previous_epoch7_offset0_source_ZIP_SHA256"] = "8eaee465a4d4b31c318a0abd26303cebcedaa23345fa32d57c16cbdd84ac28d2"
        r["previous_epoch7_offset0_source_candidate_count"] = 0
        r["previous_epoch7_offset128_run"] = 44088047232
        r["previous_epoch7_offset128_source_ZIP_SHA256"] = "a1fdd2b04e03a7d8455c7bd09c0959016a75427bd1c6d2a6c2c35dc63d5ddf46"
        r["previous_epoch7_offset128_source_candidate_count"] = 0
        r["previous_epoch7_offset256_run"] = 44105932217
        r["previous_epoch7_offset256_source_ZIP_SHA256"] = "7b9a159cc772dbd1b859139ff498cc900947971602a855617ace3d6322b0bf37"
        r["previous_epoch7_offset256_source_candidate_count"] = 0
        r["previous_epoch7_offset384_run"] = 44108050010
        r["previous_epoch7_offset384_source_ZIP_SHA256"] = "9b181cbed5ba5d33e5f5accc1036751c9aa120060319331601d7d96650b4788c"
        r["previous_epoch7_offset384_source_candidate_count"] = 0
        r["previous_epoch7_offset512_run"] = 44125684318
        r["previous_epoch7_offset512_source_ZIP_SHA256"] = "dc1f2e0171661170ed8b8d393b279169d8843f06f8def56778f436327dbaa9aa"
        r["previous_epoch7_offset512_source_candidate_count"] = 0
        r["previous_epoch7_offset640_run"] = 44138261989
        r["previous_epoch7_offset640_source_ZIP_SHA256"] = "02c9c1b60c375f1450fe4f299ce9c793f20924f939d3dd2d27b1d5334a3c5153"
        r["previous_epoch7_offset640_source_candidate_count"] = 0
        r["previous_epoch7_offset768_run"] = 44140649923
        r["previous_epoch7_offset768_source_ZIP_SHA256"] = "d797a36f97b6f26db6ebfeeb33ca4e8e9c7444ec5157580730af3a8be238ccd7"
        r["previous_epoch7_offset768_original_photo_ID"] = "764609465"
        r["previous_epoch7_offset768_audited_delta_SHA256"] = "2f5d7d4b905e071caf1198ecdbac9835b3793c0e0d37f2fc73aa5d41ab8229b8"
        r["previous_epoch7_offset768_metadata_bug_corrected_this_receipt_only"] = True
        r["previous_epoch7_offset896_run"] = 44159159412
        r["previous_epoch7_offset896_source_artifact"] = 14363370860
        r["previous_epoch7_offset896_source_ZIP_SHA256"] = "c2b2d53d678ff6590abaa7cc5d5bdabbb72d5a3829b2212a6439d02c6df1647a"
        r["previous_epoch7_offset896_source_candidate_count"] = 0
        r["previous_epoch7_offset1024_run"] = 44170221242
        r["previous_epoch7_offset1024_source_artifact"] = 14367587413
        r["previous_epoch7_offset1024_source_ZIP_SHA256"] = "0fcc9363f7598ec278e4eb8c5cbc621abac9d617b79591fd0f1f46135539d7a4"
        r["previous_epoch7_offset1024_source_candidate_count"] = 0
        r["previous_epoch7_offset1152_run"] = 44173253159
        r["previous_epoch7_offset1152_source_artifact"] = 14368388310
        r["previous_epoch7_offset1152_source_ZIP_SHA256"] = "e00fc9d0215f0fdb172b6fa7731db2a68307b7d4274b7e667c72a84bc2605520"
        r["previous_epoch7_offset1152_source_candidate_count"] = 0
        r["previous_epoch7_offset1280_run"] = 44192320978
        r["previous_epoch7_offset1280_source_artifact"] = 14376667948
        r["previous_epoch7_offset1280_source_ZIP_SHA256"] = "d28b743dbff4832fd85d6c2cb4b950ddbd88a23f962615643fecd56c46913664"
        r["previous_epoch7_offset1280_source_candidate_count"] = 0
        r["previous_epoch7_offset1408_run"] = 44201443264
        r["previous_epoch7_offset1408_source_artifact"] = 14381214117
        r["previous_epoch7_offset1408_source_ZIP_SHA256"] = "cdcba19d017db4420bddf874b7ebd5ee6e668474af0c5609f5b938013bb5f89b"
        r["previous_epoch7_offset1408_source_candidate_count"] = 0
        r["previous_epoch7_offset1536_run"] = 44208024637
        r["previous_epoch7_offset1536_source_artifact"] = 14385435693
        r["previous_epoch7_offset1536_source_ZIP_SHA256"] = "53e6f7e83be2b7f0ff68fd800c9d6675d695efad3ec4a4a8667cefd3ff1207fe"
        r["previous_epoch7_offset1536_audited_delta_SHA256"] = NEW_DELTA_SHA
        r["previous_epoch7_offset1536_photo_ID"] = "764630362"
        r["previous_epoch7_offset1664_run"] = 44215477957
        r["previous_epoch7_offset1664_source_artifact"] = 14387762866
        r["previous_epoch7_offset1664_source_ZIP_SHA256"] = "8675bd7989a38132c214d66980424d361cbb73646a15a9c0378a01b3688c0b5e"
        r["previous_epoch7_offset1664_receipt_SHA256"] = "684b182b0f4e5aec4a031ba9e165cb84d0bc3e713ac3cedbf29c2977120dd7f5"
        r["previous_epoch7_offset1664_source_candidate_count"] = 0
        r["previous_epoch7_offset1792_run"] = 44219184371
        r["previous_epoch7_offset1792_source_artifact"] = 14390790013
        r["previous_epoch7_offset1792_source_ZIP_SHA256"] = "6042d82ea0bca4894bd6b32ef9d48334c277ef1e73b8b06a527d42536181c6bb"
        r["previous_epoch7_offset1792_receipt_SHA256"] = "83614ad3a1c308fec0c13ab98aa4ae27437e8cc6c9aaff53977635dc1ceac76c"
        r["previous_epoch7_offset1792_source_candidate_count"] = 0
        r["previous_epoch6_terminal_source_ZIP_SHA256"] = (
            "3ce764b6d837348d02009173513539772d71de569088832a2f768d09c22fb87e"
        )
        r.pop("recent_97_nonempty_raw_image_SHA_count", None)
        receipt_path.write_text(json.dumps(r, sort_keys=True, indent=2) + "\n")
        print("MOIS_EPOCH7_OFFSET1920_NEW_TEMPORAL_SOURCE_SIGNED",
              json.dumps(r, sort_keys=True), flush=True)
    elif args.stage == "audit":
        m.audit(args.source, args.worker, args.output)
    else:
        m.append(args.atlas, args.source, args.worker, args.audit, args.output)

if __name__ == "__main__":
    main()
