#!/usr/bin/env python3
"""REAL FCP photo opportunity depth constrains Mois seasonal sampling support.

Reads only the published FCP 42,111-taxon ORIGINAL SOURCE COUNTS receipt,
never original pixels, phenology annotations or their taxon/site identities.
The frequency thresholds bound the number of species eligible for N selected
distinct observations in two/four seasons. Since photo-count >= observation
count, these are valid NECESSARY upper bounds, not observed seasonal support.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path

SCHEMA="mois_real_original_fcp_photo_count_seasonal_design_support_bound_v1"
STATUS="REAL_ORIGINAL_SOURCE_COUNT_UPPER_BOUNDS_NO_SEASON_OR_COLOUR_MEASURED"
PINNED_COUNT_SHA="ef26c16252dd785f2d08bde3b18938c7624d73a54386fa3d3974dceb09ed97a9"
PINNED_ARCHIVE_SHA="bfc1a97c436bf5989fccd43f12392953128489b76b92cea24da335ad73da4d8e"
PINNED_FRAME_SHA="5d871bd1f1190c68fd513bb527e6c93de1606dc35ff8cfff5f6187850dc382cc"
SOURCE_RECEIPT_URL="https://github.com/zuizui0223/fcp/blob/analysis/fcp-fresh-uncapped-zero-20261015/results/fcp_real_42111_source_opportunity_20261016/verified_source_counts_capacity_receipt.json"


def derive(source:dict)->dict:
    if (source.get("schema")!="fcp_42111_source_photo_opportunity_scientific_receipt_v1"
        or source.get("status")!="ORIGINAL_PROVIDER_42111_SOURCE_PHOTO_COUNTS_INDEPENDENTLY_VERIFIED_NO_COLOUR_MEASURED"
        or source.get("source_snapshot_date")!="2026-09-27"
        or source.get("original_source_frame_sha256")!=PINNED_FRAME_SHA
        or source.get("original_bundle_sha256")!=PINNED_ARCHIVE_SHA
        or source.get("original_source_photo_count_ledger_sha256")!=PINNED_COUNT_SHA
        or source.get("source_taxa")!=42111
        or source.get("taxa_with_at_least_one_licensed_source_photo")!=42095
        or source.get("taxa_without_licensed_source_photos")!=16
        or source.get("source_photos_global_unique_licensed")!=93209351
        or source.get("source_photo_count_run")!=38749681793
        or source.get("independent_42111_identity_audit_run")!=38818718877
        or source.get("data_scope_limitations",{}).get("source_original_photo_ID_manifests_frozen") is not False
        or source.get("data_scope_limitations",{}).get("original_image_pixels_fully_downloaded_and_measured") is not False):
        raise ValueError("Not the real verified original 42,111-taxon FCP metadata opportunity audit")
    thresholds=source["species_by_minimum_source_photo_depth"]
    if (thresholds.get("10")!=39018 or thresholds.get("20")!=36379
        or thresholds.get("40")!=32794 or thresholds.get("100")!=26510):
        raise ValueError("Real photo opportunity depth/source version changed")
    source_taxa=42111
    n20=thresholds["20"]
    n40=thresholds["40"]
    if not 0<n40<=n20<=42095:
        raise ValueError("The real photo depth support hierarchy is invalid")
    return {
        "schema":SCHEMA,"status":STATUS,
        "original_real_FCP_source_receipt":SOURCE_RECEIPT_URL,
        "source_snapshot":"2026-09-27",
        "original_count_ledger_sha256":PINNED_COUNT_SHA,
        "original_archive_sha256":PINNED_ARCHIVE_SHA,
        "original_taxon_frame_sha256":PINNED_FRAME_SHA,
        "source_named_taxa":source_taxa,
        "source_photo_opportunities":93209351,
        "source_photo_count_gini_approx":0.85852981,
        "source_species_minimum_photo_opportunities":{
            "1":42095,"10":39018,"20":n20,"40":n40,"100":26510
        },
        "two_seasons_ten_distinct_observations_each":{
            "necessary_original_photo_opportunities":20,
            "MAX_species_not_ruled_out_by_photo_counts":n20,
            "at_least_species_ruled_out_by_insufficient_photo_counts":source_taxa-n20,
            "upper_fraction_of_source_taxa":n20/source_taxa
        },
        "four_seasons_ten_distinct_observations_each":{
            "necessary_original_photo_opportunities":40,
            "MAX_species_not_ruled_out_by_photo_counts":n40,
            "at_least_species_ruled_out_by_insufficient_photo_counts":source_taxa-n40,
            "upper_fraction_of_source_taxa":n40/source_taxa
        },
        "photographic_opportunities_are_not_unique_observations":True,
        "same_species_seasons_actually_observed":False,
        "ROI_four_colour_classified_photos_measured":0,
        "qualifying_different_observers_proven":False,
        "photo_sampling_inclusion_probs_known_only_after_full_ID_freeze":True,
        "not_evidence_of_seasonal_colour_or_weather_relation":True,
        "note":"Real complete source PHOTO counts upper-bound feasibility for 10 unique original observations in each phase; actual season/date/climate/observer support and image classifiability can only reduce support. Not a world flowering-population prevalence estimator."
    }


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fcp-real-count-receipt",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    a=ap.parse_args()
    if a.output.exists():
        raise ValueError("Refuse overwrite of an older seasonal design real-source count receipt")
    source=json.loads(a.fcp_real_count_receipt.read_text(encoding="utf8"))
    result=derive(source)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n",encoding="utf8")
    print(json.dumps({key:result[key] for key in (
        "status","source_named_taxa",
        "two_seasons_ten_distinct_observations_each",
        "four_seasons_ten_distinct_observations_each")},sort_keys=True))


if __name__=="__main__":
    main()
