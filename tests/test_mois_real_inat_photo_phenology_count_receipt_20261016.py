"""Lock 35 actual public iNaturalist photo-observation annotation count queries.

This source is 2026 API counts of year-2025 photographed Plantae, not
the dated 93,209,351 original FCP photo-ID frame or flower-colour pixels.
"""
from __future__ import annotations
import json
import math
from pathlib import Path
from urllib.parse import parse_qs,urlsplit

ROOT=Path(__file__).resolve().parents[1]
RECEIPT=ROOT/"results/mois_real_inat_photo_phenology_annotation_count_pilot_20261016.json"
FILTER_KEYS={
    "all_plant_photo_observations","any_flowers_fruits_annotation",
    "explicit_flowers_annotation","explicit_no_flowers_fruits_annotation",
    "without_flowers_fruits_annotation"
}

def verified():
    r=json.loads(RECEIPT.read_text(encoding="utf8"))
    assert r["schema"]=="mois_bounded_public_inat_photo_phenology_annotation_count_pilot_v1"
    assert r["status"]=="REAL_INAT_API_PUBLIC_PHOTO_OBSERVATION_PHENOLOGY_COUNTS_NOT_FCP_SOURCE"
    assert r["queries_total"]==35 and len(r["cohorts"])==7
    assert r["verified_source_execution_url"]=="https://github.com/zuizui0223/mois/actions/runs/38843210801"
    assert r["verified_saved_artifact_audit_execution_url"]=="https://github.com/zuizui0223/mois/actions/runs/38843566842"
    assert r["source_plant_taxon_id"]==47126
    assert r["no_iNaturalist_image_pixels_opened"] is True
    assert r["actual_original_source_FCP_42111_taxon_frame_sampled"] is False
    assert r["original_93209351_photo_ID_manifest_used"] is False
    assert r["annotations_assumed_frozen_as_of_2026_09_27"] is False
    assert r["whole_plant_population_or_all_flowering_species_inferred"] is False
    return r


def test_world_2025_public_plant_photo_observation_counts():
    r=verified()
    w=r["cohorts"]["world_2025_observed"]
    assert w["total_observed_plant_photo_records"]==28565884
    assert w["with_any_phenology_annotation"]==2550229
    assert w["explicit_flowers_observations"]==1615314
    assert w["explicit_no_flowers_fruits_observations"]==424670
    assert w["without_phenology_annotation_filter_records"]==26015038
    assert math.isclose(w["any_annotation_fraction_observed"],2550229/28565884)
    assert math.isclose(w["explicit_flowers_fraction_observed"],1615314/28565884)
    # The five serial API count queries are NOT an atomic snapshot; the
    # two union/exclusion filters differ by a small nonzero count.
    assert w["any_plus_without_minus_total"]==-617
    assert w["any_plus_without_partitions_exactly"] is False
    assert w["unannotated_is_not_nonflowering"] is True


def test_region_and_observation_season_are_heterogeneous_not_colours():
    r=verified()["cohorts"]
    expected={
        "northern_east_asia_spring_2025":(52403,6568,5351),
        "northern_east_asia_summer_2025":(38805,4481,3245),
        "southern_australia_spring_2025":(395531,84266,67646),
        "southern_australia_summer_2024_25":(143217,25636,16502),
        "tropical_south_america_calendar_q1_2025":(26576,2497,1519),
        "tropical_south_america_calendar_q3_2025":(28687,1648,1105),
    }
    for label,nums in expected.items():
        q=r[label]
        assert (q["total_observed_plant_photo_records"],
                q["with_any_phenology_annotation"],
                q["explicit_flowers_observations"])==nums
    assert (r["northern_east_asia_spring_2025"]["explicit_flowers_fraction_observed"] >
            r["northern_east_asia_summer_2025"]["explicit_flowers_fraction_observed"])
    assert (r["southern_australia_spring_2025"]["explicit_flowers_fraction_observed"] >
            r["southern_australia_summer_2024_25"]["explicit_flowers_fraction_observed"])
    assert (r["tropical_south_america_calendar_q1_2025"]["explicit_flowers_fraction_observed"] >
            r["tropical_south_america_calendar_q3_2025"]["explicit_flowers_fraction_observed"])
    assert r["tropical_south_america_calendar_q1_2025"]["filters"]["d1"]=="2025-01-01"
    assert r["southern_australia_summer_2024_25"]["filters"]["d1"]=="2024-12-01"


def test_reproducible_provider_query_urls_and_digest_only_not_individuals():
    r=verified()
    for label,row in r["cohorts"].items():
        assert set(row["source_query_provenance"])==FILTER_KEYS
        for key,entry in row["source_query_provenance"].items():
            assert set(entry)=={"request_url","original_api_response_sha256","n_observations"}
            url=urlsplit(entry["request_url"])
            assert url.scheme=="https" and url.netloc=="api.inaturalist.org"
            assert url.path=="/v1/observations"
            query=parse_qs(url.query)
            assert query["taxon_id"]==["47126"]
            assert query["photos"]==["true"] and query["per_page"]==["1"]
            assert query["d1"]==[row["filters"]["d1"]]
            assert query["d2"]==[row["filters"]["d2"]]
            assert entry["n_observations"]>=0
            assert len(entry["original_api_response_sha256"])==64
        assert not any(k in row for k in ("user_id","observer","latitude","longitude","photos","observations"))
    assert all("results" not in k for k in r)
