#!/usr/bin/env python3
"""Fourth independent post-snapshot iNaturalist CC pixel wave (offset192).

Requires the 1,231,244 old ID ledger and all three signed prior new-image
source/audit pairs. Historical raw bytes are NOT available: never append old
atlas or misstate that raw SHA collision check as complete.
"""
from __future__ import annotations
import argparse, json, shutil
from pathlib import Path
from datetime import datetime, timezone
from urllib.parse import urlsplit
import pandas as pd
from scripts.acquisition.mois_new_inat_postsnapshot_verified_original_roi9 import (
    check, write_json, UUID, ARM, EXPECTED_ATTEMPTED, WORKER_SCHEMA,
    audit as base_audit,
)
from scripts.acquisition.mois_new_inat_postsnapshot_third_fullID_roi9 import (
    load_prior, HIST_NAME, AUD_NAME,
)
from scripts.acquisition.mois_zero9_recent_inat_new_created_photo_source import (
    ATLAS_SHA, source_sha, SCHEMA as API_SCHEMA, CC,
)
from scripts.acquisition.mois_full_index_real_original_photo_measure import sha as sha256

OFFSET = 192
PRIOR_COUNTS = (2, 1, 3)

def prior_chain(first_source, first_audit, second_source, second_audit,
                third_source, third_audit):
    first = load_prior(first_source, first_audit, PRIOR_COUNTS[0])
    second = load_prior(second_source, second_audit, PRIOR_COUNTS[1], offset=64)
    third = load_prior(third_source, third_audit, PRIOR_COUNTS[2], offset=128)
    parts = [first, second, third]
    full_hash = sha256(first[3])
    check(all(sha256(p[3]) == full_hash for p in parts),
          "Historical original-ID ledger changed between waves")
    ids = pd.concat([p[4] for p in parts], ignore_index=True)
    raw = {}
    for p in parts:
        check(not set(raw).intersection(p[5]),
              "Prior signed source photo ID repeated")
        check(not set(raw.values()).intersection(p[5].values()),
              "Prior signed raw original image SHA repeated")
        raw.update(p[5])
    check(len(ids) == sum(PRIOR_COUNTS) == 6
          and ids.photo_id.is_unique and ids.observation_id.is_unique
          and ids.inat_taxon_id.is_unique and len(raw) == 6,
          "Prior signed photo/observation/taxon lineage invalid")
    m3 = third[0]
    check(m3.get("source_prior_post_snapshot_photo_IDs") == ids.photo_id.iloc[:3].tolist()
          and m3.get("source_prior_post_snapshot_observation_UUIDs") ==
              ids.observation_id.iloc[:3].tolist()
          and m3.get("source_prior_post_snapshot_original_image_SHA256") ==
              {**first[5], **second[5]},
          "Third source not chained to original first+second")
    return first[3], ids, raw, [sha256(Path(p) / "receipt.json")
                              for p in (first_audit, second_audit, third_audit)]

def gate(meta, receipt, atlas, first_source, first_audit, second_source,
         second_audit, third_source, third_audit, output):
    meta, receipt, atlas, output = map(Path, (meta, receipt, atlas, output))
    check(not output.exists(), "Refuse source re-freeze or same-SHA re-download")
    r = json.loads(receipt.read_text())
    check(r.get("schema") == API_SCHEMA
          and r.get("source_deterministic_hash_offset") == OFFSET
          and r.get("sampled_zero_ROI9_source_taxa") == 64
          and r.get("original_source_zero_ROI9_taxa") == 2229
          and r.get("source_taxa_API_errors") == 0
          and r.get("old_photo_ID_1_231_244_ledger_disjointness_verified") is False
          and r.get("original_photo_pixels_downloaded") == 0
          and r.get("genuine_nine_fraction_ROI_measured") == 0
          and source_sha(meta) == r["output_sha256"][meta.name]
          and source_sha(atlas) == ATLAS_SHA,
          "Offset192 new-source receipt/atlas hash invalid")
    census = pd.read_csv(atlas, dtype={"inat_taxon_id": str})
    check(len(census) == 42111 and census.inat_taxon_id.is_unique
          and int(census.after_strict9.eq(0).sum()) == 2229,
          "Signed historical species frame altered")
    reg, prior_ids, prior_sha, prior_audits = prior_chain(
        first_source, first_audit, second_source, second_audit,
        third_source, third_audit)
    history = pd.read_csv(reg, dtype=str, keep_default_na=False,
                          usecols=["photo_id", "observation_id", "inat_taxon_id"])
    check(len(history) == EXPECTED_ATTEMPTED and history.photo_id.is_unique
          and history.photo_id.ne("").all(), "Historical 1,231,244 photo IDs incomplete")
    candidates = pd.read_csv(meta, dtype=str, keep_default_na=False)
    check(len(candidates) == r["source_taxa_with_new_post_snapshot_API_photo_opportunities"]
          == r["new_independent_original_iNaturalist_photo_IDs_seen"]
          and 0 <= len(candidates) <= 64
          and candidates.photo_id.is_unique and candidates.observation_id.is_unique
          and candidates.inat_taxon_id.is_unique
          and candidates.inat_taxon_id.isin(
            census.loc[census.after_strict9.eq(0), "inat_taxon_id"]).all(),
          "New metadata candidate count/taxon/photo UUID conflict")
    check(not (set(candidates.photo_id) & (set(history.photo_id) | set(prior_ids.photo_id)))
          and not (set(candidates.observation_id) &
                   (set(history.observation_id) | set(prior_ids.observation_id)))
          and not (set(candidates.inat_taxon_id) & set(prior_ids.inat_taxon_id)),
          "Reused original photo/observation/taxon ID across historic+three recent waves")
    for row in candidates.itertuples(index=False):
        u = urlsplit(row.photo_url_large)
        check(row.photo_id.isdigit()
              and UUID.fullmatch(row.observation_id) is not None
              and row.photo_license in CC
              and u.scheme == "https"
              and u.hostname == "inaturalist-open-data.s3.amazonaws.com"
              and u.path.startswith("/photos/" + row.photo_id + "/")
              and u.path.lower().endswith((".jpg", ".jpeg", ".png")),
              "CC/UUID/photo URL provenance invalid")
        when = datetime.fromisoformat(
            row.source_observation_created_at.replace("Z", "+00:00"))
        check(when.tzinfo is not None
              and when.astimezone(timezone.utc) >=
                datetime(2026, 10, 28, tzinfo=timezone.utc),
              "Existing old observation not independent post-snapshot")
    output.mkdir(parents=True)
    shutil.copyfile(reg, output / HIST_NAME)
    (output / "shards").mkdir()
    shard_sha = {}
    if len(candidates):
        shard = output / "shards" / "shard-000.csv"
        candidates = candidates.copy()
        candidates["measurement_id"] = (
            "MOIS-POST20261028-OFFSET192-" + candidates.photo_id)
        candidates["image_filename"] = candidates.measurement_id + ".jpg"
        candidates["shard"] = "0"
        candidates.to_csv(shard, index=False, lineterminator="\n")
        shard_sha["0"] = sha256(shard)
    manifest = {
        "schema": WORKER_SCHEMA, "source_arm": ARM,
        "colour_measured": False, "source_previous_attempts": EXPECTED_ATTEMPTED,
        "source_previous_ID_registry_SHA256": sha256(output / HIST_NAME),
        "source_prior_no_reused_original_photo_or_observation": True,
        "source_prior_post_snapshot_photo_IDs": prior_ids.photo_id.tolist(),
        "source_prior_post_snapshot_observation_UUIDs": prior_ids.observation_id.tolist(),
        "source_prior_post_snapshot_original_image_SHA256": prior_sha,
        "source_prior_three_audit_receipt_sha256": prior_audits,
        "source_post_snapshot_API_receipt_SHA256": source_sha(receipt),
        "source_post_snapshot_API_metadata_SHA256": source_sha(meta),
        "source_fixed_hash_offset": OFFSET,
        "new_photo_attempts_planned": len(candidates),
        "work_units": int(len(candidates) > 0),
        "worker_ids": [0] if len(candidates) else [],
        "photo_shard_sha256": shard_sha,
        "raw_image_pixels_downloaded": False,
        "new_nineROI_measured": False,
        "official_iNat_42111_original_frame_auto_appended": False,
        "species_wide_fixed_colour_proven": False,
    }
    write_json(output / "manifest.json", manifest)
    gate_receipt = {
        "schema": "mois_postsnapshot_offset192_hist_plus2plus1plus3_signed_gate_v1",
        "historical_verified_photo_ID_count": len(history),
        "previous_new_photo_source_count": len(prior_ids),
        "prior_raw_original_SHA_count": len(prior_sha),
        "new_CC_photo_opportunities": len(candidates),
        "new_photo_IDs_historical_collision_count": 0,
        "new_observation_UUID_historical_collision_count": 0,
        "source_manifest_sha256": sha256(output / "manifest.json"),
        "historical_full_registry_sha256": sha256(output / HIST_NAME),
        "source_only": True, "strict9_measured": 0,
        "new_photo_pixels_downloaded": 0,
    }
    write_json(output / "gate_receipt.json", gate_receipt)
    print("MOIS_OFFSET192_FULL_HIST_PLUS_SIX_NEW_SHA_GATE",
          json.dumps(gate_receipt, sort_keys=True), flush=True)

def audit_fourth(source, worker, output):
    source, worker, output = map(Path, (source, worker, output))
    m = json.loads((source / "manifest.json").read_text())
    g = json.loads((source / "gate_receipt.json").read_text())
    check(m.get("source_fixed_hash_offset") == OFFSET
          and g.get("schema") ==
              "mois_postsnapshot_offset192_hist_plus2plus1plus3_signed_gate_v1"
          and g.get("source_manifest_sha256") == sha256(source / "manifest.json")
          and g.get("historical_full_registry_sha256") == sha256(source / HIST_NAME)
          and len(m["source_prior_post_snapshot_photo_IDs"]) == 6
          and len(m["source_prior_post_snapshot_original_image_SHA256"]) == 6,
          "Unsigned fourth source lineage")
    now = pd.read_csv(source / "shards" / "shard-000.csv",
                      dtype=str, keep_default_na=False)
    b = pd.read_csv(worker / "shard_000" / "photos.csv",
                    dtype=str, keep_default_na=False)
    check(not (set(now.photo_id) & set(m["source_prior_post_snapshot_photo_IDs"]))
          and not (set(now.observation_id) &
                   set(m["source_prior_post_snapshot_observation_UUIDs"]))
          and not (set(b.image_sha256) - {""} &
                   set(m["source_prior_post_snapshot_original_image_SHA256"].values())),
          "Fourth original matches previous recent source/image")
    base_audit(source, worker, output)
    a = json.loads((output / "receipt.json").read_text())
    a.update({
        "additional_prior_new_photo_IDs_excluded": 6,
        "additional_prior_new_observation_UUID_collisions": 0,
        "additional_prior_new_raw_image_SHA_collisions": 0,
        "full_historical_raw_image_SHA_collision_check_possible": False,
        "official_iNat_42111_original_frame_auto_appended": False
    })
    write_json(output / "receipt.json", a)
    print("MOIS_OFFSET192_REAL_PIXEL_INDEPENDENT_AUDIT",
          json.dumps({k: a[k] for k in (
              "original_licensed_photo_attempts", "actual_image_downloaded",
              "actual_genuine_strict9_numeric_success",
              "actual_clear_photo_nine_palette_labels", "actual_technical_unknowns")},
              sort_keys=True), flush=True)

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="command", required=True)
    gate_parser = sub.add_parser("gate")
    for k in ("meta", "receipt", "atlas", "first-source", "first-audit",
              "second-source", "second-audit", "third-source", "third-audit",
              "output"):
        gate_parser.add_argument("--" + k, required=True)
    audit_parser = sub.add_parser("audit")
    for k in ("source", "worker", "output"):
        audit_parser.add_argument("--" + k, required=True)
    v = p.parse_args()
    if v.command == "gate":
        gate(v.meta, v.receipt, v.atlas, v.first_source, v.first_audit,
             v.second_source, v.second_audit, v.third_source, v.third_audit,
             v.output)
    else:
        audit_fourth(v.source, v.worker, v.output)
