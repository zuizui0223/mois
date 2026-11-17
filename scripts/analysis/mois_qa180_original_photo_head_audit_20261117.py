#!/usr/bin/env python3
"""MOIS QA180: authenticate the original blind 538-photo list, then HEAD photo URLs.

This code never downloads photo pixels, estimates petal colour, changes the
blinded sample, imputes missing photos, or opens the sealed colour key.
External URL availability and botanically correct floral-ROI measurement are
different scientific claims.
"""
from __future__ import annotations
import argparse
import hashlib
import io
import json
import re
import threading
import time
import zipfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import requests

REGISTRY_ZIP_SHA = "ac558f06e5a0f75cf0fab2d0a37e3f17b2fe3428f788cccb801b5a01dd2b81ac"
TRAITS_ZIP_SHA = "abcd3aee766c85a92b86a8f93f8ed696bbe49d0f6cf2f64e5f75b457abf78e6a"
REGISTRY_MEMBER = "historical_1231244_original_photo_observation_taxon_ID_registry.csv.gz"
REGISTRY_MEMBER_SHA = "f945779b760b46ace93bcfa10763d2ad56f251dbc9c8146d3e44c30b55dca710"
TRAIT_MEMBER = "current_exact_accepted_COLXR_numeric9_traits.csv.gz"
TRAIT_MEMBER_SHA = "b06d8980860ab11d1e9cff1af7dd1d92b3bfa9de4f2d227bc7e58809f931663a"
PALETTE = ("white", "yellow", "orange", "red", "pink", "magenta", "purple", "blue", "bronze")
PAGE_RE = re.compile(r"^https://www\.inaturalist\.org/photos/(\d+)$")
QA_RE = re.compile(r"^(MOIS-QA-\d{3})-P(0[1-3])$")
HEAD_MAX_WORKERS = 4
REQUEST_SPACING_S = 0.16
USER_AGENT = "MOIS-independent-flower-ROI-review-source-audit/1.0 (noncommercial research; HEAD only)"


def ensure(v, message):
    if not bool(v):
        raise ValueError("MOIS_QA180_SOURCE_AUDIT: " + message)


def sha(data):
    if isinstance(data, (Path, str)):
        data = Path(data).read_bytes()
    return hashlib.sha256(data).hexdigest()


def signed_csv(zip_path, outer_sha, member, inner_sha):
    ensure(sha(zip_path) == outer_sha, f"ZIP SHA mismatch {zip_path}")
    with zipfile.ZipFile(zip_path) as z:
        buf = z.read(member)
    ensure(sha(buf) == inner_sha, f"inner gzip SHA mismatch {member}")
    return pd.read_csv(io.BytesIO(buf), compression="gzip", dtype="string", low_memory=False)


def audit(sheet, registry, traits):
    r = pd.read_csv(sheet, dtype="string").fillna("")
    ensure(len(r) == 538, "sheet changed from 538 frozen photos")
    ensure({"blind_photo_review_id", "original_photo_page"}.issubset(r), "blind sheet schema drift")
    check = []
    for a, b in zip(r.blind_photo_review_id, r.original_photo_page):
        m = QA_RE.fullmatch(str(a))
        n = PAGE_RE.fullmatch(str(b))
        ensure(m is not None and n is not None, "unexpected blind or photo URL")
        check.append((m.group(1), m.group(2), n.group(1)))
    r["qa_species_id"] = [v[0] for v in check]
    r["photo_order"] = [v[1] for v in check]
    r["photo_id"] = [v[2] for v in check]
    ensure(r.photo_id.is_unique and r.blind_photo_review_id.is_unique,
           "repeated IDs in fixed blinded photo sheet")
    ensure(r.qa_species_id.nunique() == 180, "frozen blinded species count changed")
    counts = r.groupby("qa_species_id").photo_id.nunique()
    ensure((counts.eq(3).sum(), counts.eq(2).sum()) == (178, 2),
           "180 species: 178 x3 and 2 x2 integrity changed")
    ensure(set(counts[counts.eq(2)].index) == {"MOIS-QA-027", "MOIS-QA-046"},
           "two 2-photo review IDs changed")
    reviewer_fields = [
        "flower_visible_yes_no_uncertain", "flower_petal_visible_yes_no_uncertain",
        "visible_flower_colour_9class_or_unscorable", "ROI_mask_accessible_yes_no",
        "ROI_petals_yes_no_uncertain", "nonpetal_or_background_yes_no_uncertain",
        "observer_reviewer_id", "comment",
    ]
    ensure(set(reviewer_fields).issubset(r.columns), "missing original review form fields")
    ensure(r[reviewer_fields].eq("").all().all(), "previous reviewer annotations exist: do not overwrite")

    ensure(len(registry) == 1231244 and registry.photo_id.is_unique,
           "source registry ID-count changed")
    ensure(len(traits) == 37104
           and traits.inat_taxon_id.is_unique
           and traits.gbif_colxr_accepted_species_key.is_unique,
           "accepted trait universe changed")
    right = registry[["photo_id", "observation_id", "inat_taxon_id"]]
    z = r.merge(right, on="photo_id", how="left", validate="one_to_one", indicator=True)
    ensure(z._merge.eq("both").all() and z.observation_id.notna().all(),
           "a frozen photo was not found in authentic historical registry")
    z = z.drop(columns="_merge")
    ensure(z.groupby("qa_species_id").inat_taxon_id.nunique().eq(1).all(),
           "blind species group has multiple original iNaturalist taxa")
    ensure(z.groupby("qa_species_id").observation_id.nunique().equals(
        z.groupby("qa_species_id").photo_id.nunique()),
        "multiple photos from same observation selected for a blind species")
    ensure(z.observation_id.is_unique, "one observation appears for two frozen images")

    tr = traits[[
        "inat_taxon_id", "gbif_colxr_accepted_species_key", "photo_modal_palette9",
        "photo_evidence_status", "photo_source_species_name",
        *("observed_mean_" + c for c in PALETTE),
    ]].copy()
    z = z.merge(tr, on="inat_taxon_id", how="left", validate="many_to_one")
    ensure(z.gbif_colxr_accepted_species_key.notna().all(), "frozen taxa not accepted")
    ensure(z.groupby("qa_species_id").gbif_colxr_accepted_species_key.nunique().eq(1).all(),
           "blind species ID maps to more than one accepted species")
    ensure(z.gbif_colxr_accepted_species_key.nunique() == 180,
           "two blind species map to same accepted species key")

    # The strata are *checked only in aggregate* and never written per-row
    # to the reviewer sheet; analysis must not reveal the model-predicted hue.
    means = z[["observed_mean_" + c for c in PALETTE]].astype(float).to_numpy()
    observed_mean = np.asarray(PALETTE)[means.argmax(axis=1)]
    z["agree"] = z.photo_modal_palette9.to_numpy(str) == observed_mean
    one = z.drop_duplicates("qa_species_id")
    summary = one.groupby(["photo_modal_palette9", "agree"]).size()
    ensure(len(summary) == 18 and (summary == 10).all(),
           "frozen 18-colour strata × 10 species changed")
    return z, {
        "accepted_species": int(one.gbif_colxr_accepted_species_key.nunique()),
        "source_photo_ids": int(z.photo_id.nunique()),
        "source_observation_ids": int(z.observation_id.nunique()),
        "exact_strata_count": 18,
        "species_per_blinded_stratum": 10,
        "species_with_3_distinct_observations": int(counts.eq(3).sum()),
        "species_with_2_distinct_observations": int(counts.eq(2).sum()),
        "existing_biological_review_annotations": 0,
        "source_sample_not_changed": True,
    }


class Throttle:
    def __init__(self):
        self.lock = threading.Lock()
        self.next_time = 0.0

    def wait(self):
        with self.lock:
            now = time.monotonic()
            gap = max(0.0, self.next_time - now)
            self.next_time = max(now, self.next_time) + REQUEST_SPACING_S
        if gap:
            time.sleep(gap)


limiter = Throttle()


def safe_head(url):
    limiter.wait()
    # HEAD sends no photo pixels; do not automatically GET, download or archive images.
    try:
        with requests.Session() as s:
            s.trust_env = True
            resp = s.head(
                url,
                headers={"User-Agent": USER_AGENT},
                allow_redirects=True,
                timeout=(5, 10),
            )
            return {
                "http_status": int(resp.status_code),
                "mime_type": resp.headers.get("Content-Type", "").split(";")[0],
                "content_length": resp.headers.get("Content-Length", ""),
                "final_host": requests.utils.urlparse(resp.url).hostname or "",
                "error_class": "",
            }
    except requests.RequestException as e:
        return {
            "http_status": None, "mime_type": "", "content_length": "",
            "final_host": "", "error_class": type(e).__name__,
        }


def verify_one(row):
    page = str(row.original_photo_page)
    pid = str(row.photo_id)
    medium = f"https://inaturalist-open-data.s3.amazonaws.com/photos/{pid}/medium.jpg"
    page_status = safe_head(page)
    s3_status = safe_head(medium)
    return {
        "blind_photo_review_id": str(row.blind_photo_review_id),
        "photo_id": pid,
        "photo_page_http": page_status["http_status"],
        "photo_page_mime": page_status["mime_type"],
        "photo_page_error": page_status["error_class"],
        "provider_medium_http": s3_status["http_status"],
        "provider_medium_mime": s3_status["mime_type"],
        "provider_medium_length": s3_status["content_length"],
        "provider_medium_error": s3_status["error_class"],
        "photo_page_observation_and_license_confirmed": False,
        "image_pixels_downloaded": 0,
        "botanical_annotation_done": False,
    }


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--blind", required=True)
    parser.add_argument("--registry", required=True)
    parser.add_argument("--traits", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--probe", type=int, default=0,
                        help="0 = source integrity only, otherwise HEAD first N distinct photos")
    args = parser.parse_args(argv)
    out = Path(args.out)
    ensure(not out.exists(), "refuse output overwrite")
    registry = signed_csv(args.registry, REGISTRY_ZIP_SHA,
                          REGISTRY_MEMBER, REGISTRY_MEMBER_SHA)
    traits = signed_csv(args.traits, TRAITS_ZIP_SHA, TRAIT_MEMBER, TRAIT_MEMBER_SHA)
    joined, source_summary = audit(args.blind, registry, traits)
    out.mkdir(parents=True)
    print("MOIS_QA180_SOURCE_IDENTITY_PASS", json.dumps(source_summary), flush=True)

    probes = []
    count = min(max(0, args.probe), len(joined))
    if count:
        with ThreadPoolExecutor(max_workers=HEAD_MAX_WORKERS) as pool:
            future = {pool.submit(verify_one, row): str(row.blind_photo_review_id)
                      for row in joined.head(count).itertuples(index=False)}
            for f in as_completed(future):
                probes.append(f.result())
    if probes:
        df = pd.DataFrame(probes).sort_values("blind_photo_review_id")
        ensure(df.photo_id.is_unique, "duplicate probed URLs")
        df.to_csv(out / "original_photo_url_HEAD_only.csv", index=False)
        pstats = {
            "probed_photos": len(df),
            "s3_http_200_image": int(
                (df.provider_medium_http == 200)
                & df.provider_medium_mime.str.startswith("image/")
            ).sum(),
            "s3_http_403": int(df.provider_medium_http.eq(403).sum()),
            "s3_http_404": int(df.provider_medium_http.eq(404).sum()),
            "provider_network_or_timeout_errors": int(
                df.provider_medium_error.fillna("").ne("").sum()
            ),
            "photo_page_http_200_html": int(
                ((df.photo_page_http == 200)
                 & df.photo_page_mime.str.contains("html", case=False)).sum()
            ),
            "no_image_bytes_saved": True,
        }
    else:
        pstats = {"probed_photos": 0, "no_image_bytes_saved": True}
    receipt = {
        "schema": "mois_original_blind_QA180_HEAD_gate_v1",
        "date_utc": datetime.now(timezone.utc).isoformat(),
        "frozen_review_sheet_sha256_canonical": sha(Path(args.blind)),
        "original_historical_registry_zip_sha256": REGISTRY_ZIP_SHA,
        "COLXR_trait_zip_sha256": TRAITS_ZIP_SHA,
        "source": source_summary,
        "head_gate": pstats,
        "warning": "HTTP HEAD success does not verify photo ID-to-visible-image identity, licence, botanically valid flowers, ROI masks, or true petal colours",
    }
    receipt["files_sha256"] = {p.name: sha(p) for p in out.iterdir()}
    (out / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print("MOIS_QA180_REAL_SOURCE_AUDIT_DONE", json.dumps(receipt), flush=True)


if __name__ == "__main__":
    main()
