#!/usr/bin/env python3
"""Create a deterministic INTERNAL numerical-evidence ZIP, never original photo data.

This bundle is an editor/author technical snapshot of files already public in mois.
It is NOT a GEB anonymous reviewer archive, raw input data, or a DOI deposit.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import io
import json
import zipfile
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
FILES=(
    "results/fcp_macro_direction_20261013/direction_ledger.csv",
    "results/fcp_macro_direction_20261013/decision_receipt.json",
    "results/fcp_macro_direction_20261013/independent_verified.json",
    "results/fcp_m1_matched_spatial_support_20261014/matched_spatial_support_differences.csv",
    "results/fcp_m1_matched_spatial_support_20261014/decision_receipt.json",
    "results/fcp_m1_matched_spatial_support_20261014/independent_verified.json",
    "results/mois_original_source_independent_recovery_audit_20261014.json",
    "results/figure_source_data.json",
    "docs/MOIS_GEB_RESEARCH_ARTICLE_DRAFT_20261014.md",
    "docs/MOIS_GEB_SUPPORTING_INFORMATION_20261014.md",
    "docs/MOIS_GEB_ACTUAL_SOURCE_RECOVERY_AUDIT_20261014.md",
    "scripts/analysis/fcp_abiotic_hue_direction_macro_20261013.py",
    "scripts/analysis/fcp_macro_spatial_support_decomposition_20261014.py",
    "scripts/analysis/verify_mois_original_source_two_axis_wls_20261014.py",
    "figures/global_moisture_colour.svg",
    "figures/mois_photo_selection_and_matched_rainfall.svg",
    "figures/mois_environmental_transfer_and_model_form.svg",
    "figures/mois_original_four_hue_moisture.svg",
)
BAN_NAMES=("raw_inputs", "original_100543_photo_sites.zip", "original_100543_photo_sites_all36_environment.zip")
BANNED_CSV_COLUMNS=("photo_id", "observation_id", "inat_taxon_id", "latitude",
                    "longitude", "attribution", "photo_license", "image_sha256")
DEST=ROOT/"dist/mois_internal_numerical_evidence_20261014.zip"


def get_inputs():
    content={}
    for name in FILES:
        if any(x in name for x in BAN_NAMES):
            raise RuntimeError("Unsafe per-photo archive selected: "+name)
        path=ROOT/name
        if not path.is_file():
            raise RuntimeError("Required numerical evidence missing: "+name)
        content[name]=path.read_bytes()
        if name.endswith(".csv"):
            hdr=next(csv.reader(io.StringIO(content[name].decode("utf8").splitlines()[0])))
            assert not set(hdr).intersection(BANNED_CSV_COLUMNS),name
    assert sum(1 for p in content if p.startswith("figures/"))==4
    assert (len(content["results/fcp_macro_direction_20261013/direction_ledger.csv"].splitlines())-1)==1280
    assert (len(content["results/fcp_m1_matched_spatial_support_20261014/matched_spatial_support_differences.csv"].splitlines())-1)==576
    r=json.loads(content["results/mois_original_source_independent_recovery_audit_20261014.json"])
    assert r["source_photo_opportunities"]==100543
    assert r["source_rights_and_geoprivacy_human_review_required"]
    assert not r["per_photo_ids_coordinates_and_attribution_released"]
    return content


def build() -> bytes:
    content=get_inputs()
    manifest={
      "schema":"mois_internal_numerical_evidence_package_20261014_v1",
      "identity":"PUBLIC_NUMERICAL_RESULTS_INTERNAL_EDITORIAL_ARCHIVE_NOT_REVIEWER_DATA_ACCESS",
      "n_original_photo_ids":100543,
      "included_data_rows":{"original_world_hue_axes":1280,"source_matched_spatial":576},
      "original_full_source_archives_included":False,
      "individual_photographer_attribution_or_photo_site_records_included":False,
      "anonymized_reviewer_data_link_created":False,
      "permanent_archive_DOI_created":False,
      "licence_review_still_required":True,
      "files":{k:{"bytes":len(v),"sha256":hashlib.sha256(v).hexdigest()} for k,v in sorted(content.items())}
    }
    content["INTERNAL_SOURCE_MANIFEST.json"]=(json.dumps(manifest,sort_keys=True,indent=2)+"\n").encode()
    data=io.BytesIO()
    with zipfile.ZipFile(data,"w",compression=zipfile.ZIP_DEFLATED,compresslevel=9) as out:
        for name,raw in sorted(content.items()):
            info=zipfile.ZipInfo(name,date_time=(2026,10,14,0,0,0))
            info.compress_type=zipfile.ZIP_DEFLATED
            info.external_attr=0o644<<16
            out.writestr(info,raw)
    return data.getvalue()


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--write",action="store_true")
    p.add_argument("--check",action="store_true")
    args=p.parse_args()
    blob=build()
    if args.write:
        DEST.parent.mkdir(parents=True,exist_ok=True)
        DEST.write_bytes(blob)
    if args.check:
        assert DEST.exists(), "Run --write first"
        assert DEST.read_bytes()==blob, "Internal snapshot differs from frozen public source"
        with zipfile.ZipFile(io.BytesIO(blob)) as z:
            assert z.testzip() is None
            assert len(z.namelist())==len(FILES)+1
    print("PASS internal NUMERICAL snapshot",len(FILES),"source-visible files",
          len(blob),"bytes, original ZIPs and photo IDs EXCLUDED")


if __name__=="__main__":
    main()
