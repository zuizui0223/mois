#!/usr/bin/env python3
"""MOIS QA180: verify frozen photo/observation identities and licence via current iNat API.

All 538 images remain blind to the 9-colour model labels.  One public API
query is made per batch of 30 observation IDs, not per photo and not for pixels.
Current public API metadata can change since the historical registry.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import time
from datetime import datetime, timezone
from pathlib import Path
import pandas as pd
import requests

from scripts.analysis.mois_qa180_original_photo_head_audit_20261117 import (
    audit, signed_csv, REGISTRY_ZIP_SHA, REGISTRY_MEMBER, REGISTRY_MEMBER_SHA,
    TRAITS_ZIP_SHA, TRAIT_MEMBER, TRAIT_MEMBER_SHA, ensure, sha,
)

API = "https://api.inaturalist.org/v1/observations"
BATCH_SIZE = 30
API_SPACING = 1.25


def get_batch(ids, session):
    ensure(1 <= len(ids) <= BATCH_SIZE and len(set(ids)) == len(ids),
           "invalid unique original observation ID batch")
    headers = {"User-Agent": "MOIS-blinded-180-source-metadata-QA/1.0 (noncommercial)"}
    last = None
    for attempt in range(3):
        try:
            # Official recommended bulk-by-ID path avoids incomplete search-index filters.
            # https://www.inaturalist.org/pages/api+recommended+practices
            r = session.get(API + "/" + ",".join(ids),
                headers=headers, timeout=(8, 30))
            if r.status_code in (429, 500, 502, 503, 504):
                if attempt != 2:
                    time.sleep(3 * (attempt + 1))
                    continue
            r.raise_for_status()
            obj = r.json()
            ensure(isinstance(obj, dict) and isinstance(obj.get("results"), list),
                   "public iNat observation API response schema drift")
            ensure(len(obj["results"]) <= len(ids)
                   and 0 <= int(obj.get("total_results", len(obj["results"]))) <= len(ids),
                   "unexpected pagination or ignored ID filter")
            found = {}
            for o in obj["results"]:
                k = str(o.get("id"))
                ensure(k in ids and k not in found, "unexpected or repeated public observation ID")
                found[k] = o
            digest = hashlib.sha256(json.dumps(
                obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")
            ).encode()).hexdigest()
            return found, digest
        except (requests.RequestException, ValueError) as exc:
            last = str(exc)[:150]
            if attempt != 2:
                time.sleep(2 * (attempt + 1))
    raise RuntimeError("public observations batch could not be verified: " + str(last))


def inspect_one(row, observation):
    pid = str(row.photo_id)
    oid = str(row.observation_id)
    out = {
        "blind_photo_review_id": str(row.blind_photo_review_id),
        "original_photo_id": pid,
        "original_observation_id": oid,
        "original_inat_taxon_id": str(row.inat_taxon_id),
        "current_observation_record_found": observation is not None,
        "same_original_photo_in_current_observation": False,
        "current_taxon_same_as_historical_source": "",
        "current_photo_license_code": "",
        "photo_attribution_metadata_present": False,
        "current_photo_url_metadata_present": False,
        "observed_on_public": "",
        "observation_level_flowering_annotation": False,
        "not_biologically_verified": True,
    }
    if observation is None:
        return out
    ensure(str(observation.get("id")) == oid, "API observation ID mismatch")
    current_taxon = (observation.get("taxon") or {}).get("id")
    if current_taxon is not None:
        out["current_taxon_same_as_historical_source"] = (
            str(current_taxon) == str(row.inat_taxon_id)
        )
    for p in observation.get("photos") or []:
        if str(p.get("id")) == pid:
            out["same_original_photo_in_current_observation"] = True
            out["current_photo_license_code"] = str(p.get("license_code") or "")
            out["photo_attribution_metadata_present"] = bool(p.get("attribution"))
            out["current_photo_url_metadata_present"] = bool(p.get("url"))
            break
    observed = observation.get("observed_on")
    if observed:
        out["observed_on_public"] = str(observed)[:10]
    labels = {str((a.get("controlled_value") or {}).get("label") or "")
              for a in observation.get("annotations") or []}
    out["observation_level_flowering_annotation"] = (
        "Flowers" in labels or "Flowers and Fruits" in labels
    )
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--blind", required=True)
    ap.add_argument("--registry", required=True)
    ap.add_argument("--traits", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    out = Path(args.out)
    ensure(not out.exists(), "never overwrite existing API audit")
    reg = signed_csv(args.registry, REGISTRY_ZIP_SHA,
                     REGISTRY_MEMBER, REGISTRY_MEMBER_SHA)
    traits = signed_csv(args.traits, TRAITS_ZIP_SHA, TRAIT_MEMBER, TRAIT_MEMBER_SHA)
    joined, source = audit(args.blind, reg, traits)
    obsids = sorted(joined.observation_id.astype(str).unique())
    ensure(len(obsids) == 538, "independent original observation count changed")
    found, batch_digests = {}, []
    with requests.Session() as session:
        for i in range(0, len(obsids), BATCH_SIZE):
            if i:
                time.sleep(API_SPACING)
            batch = obsids[i:i+BATCH_SIZE]
            part, digest = get_batch(batch, session)
            ensure(not set(part).intersection(found), "a public observation is repeated across API batches")
            found.update(part)
            batch_digests.append({
                "first_index": i, "request_ids": len(batch),
                "returned_public_observations": len(part),
                "public_response_sha256": digest,
            })
            print("MOIS_QA180_CURRENT_INAT_BATCH",i//BATCH_SIZE+1,
                  len(part),"/",len(batch),flush=True)
    df = pd.DataFrame([inspect_one(row, found.get(str(row.observation_id)))
                       for row in joined.itertuples(index=False)])
    ensure(len(df) == 538 and df.original_photo_id.is_unique,
           "original 538 photo mapping changed")
    out.mkdir(parents=True)
    target = out / "original_photo_current_iNaturalist_API_status.csv"
    df.to_csv(target, index=False)
    summary = {
        "n_frozen_original_photo_ids": 538,
        "n_original_accepted_species": 180,
        "returned_current_public_observation": int(df.current_observation_record_found.sum()),
        "original_photo_still_present_with_original_observation": int(df.same_original_photo_in_current_observation.sum()),
        "current_photo_license_code_present": int(df.current_photo_license_code.ne("").sum()),
        "attribution_metadata_present": int(df.photo_attribution_metadata_present.sum()),
        "photo_url_metadata_present": int(df.current_photo_url_metadata_present.sum()),
        "current_original_inat_taxon_id_matches": int(df.current_taxon_same_as_historical_source.eq(True).sum()),
        "current_original_inat_taxon_id_changed": int(df.current_taxon_same_as_historical_source.eq(False).sum()),
        "current_observation_level_flowering_annotation": int(df.observation_level_flowering_annotation.sum()),
        "photo_pixels_downloaded": 0,
        "individual_botanical_petal_colour_verified": 0,
        "current_public_API_record_is_not_original_historical_snapshot": True,
        "observation_level_flowering_annotation_is_not_photo_specific": True,
        "licence_metadata_is_not_permission_for_unrestricted_redistribution": True,
    }
    receipt = {
        "schema": "mois_QA180_original_538_current_photo_observation_and_licence_metadata_v1",
        "utc": datetime.now(timezone.utc).isoformat(),
        "current_api": API + "/{comma_separated_observation_ids}",
        "batch_size": BATCH_SIZE,
        "requests": len(batch_digests),
        "raw_api_responses_not_redistributed": True,
        "batch_input_response_sha": batch_digests,
        "frozen_reviewer_sheet_sha256_canonical": sha(Path(args.blind)),
        "historical_source_zip_sha": REGISTRY_ZIP_SHA,
        "trait_source_zip_sha": TRAITS_ZIP_SHA,
        "frozen_source_provenance": source,
        "current_metadata": summary,
        "numeric_output_sha": sha(target),
    }
    (out / "receipt.json").write_text(json.dumps(
        receipt, indent=2, ensure_ascii=False
    )+"\n")
    print("MOIS_QA180_CURRENT_PUBLIC_API_METADATA_DONE",
          json.dumps(summary),flush=True)


if __name__ == "__main__":
    main()
