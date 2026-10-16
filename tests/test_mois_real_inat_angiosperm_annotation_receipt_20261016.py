"""Original public API Angiospermae (flowering-plant) counts, not FCP photo IDs."""
from __future__ import annotations
import json
import math
from pathlib import Path
from urllib.parse import parse_qs,urlsplit

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/"results/mois_real_inat_angiosperm_phenology_annotation_count_pilot_20261016.json"
PLANT=ROOT/"results/mois_real_inat_photo_phenology_annotation_count_pilot_20261016.json"

def source():
    r=json.loads(SOURCE.read_text(encoding="utf8"))
    assert r["schema"]=="mois_bounded_public_inat_photo_phenology_annotation_count_pilot_v1"
    assert r["status"]=="REAL_INAT_API_PUBLIC_PHOTO_OBSERVATION_PHENOLOGY_COUNTS_NOT_FCP_SOURCE"
    assert r["original_aggregate_artifact_file_sha256"]=="67c22183e6919be4ed5f7d466f6047c178ca64d18a34479cee71c51b95ccc03e"
    assert r["verified_original_execution_url"]=="https://github.com/zuizui0223/mois/actions/runs/38844098039"
    assert r["verified_artifact_recovery_url"]=="https://github.com/zuizui0223/mois/actions/runs/38844372307"
    assert r["source_plant_taxon_id"]==47125
    assert r["source_plant_group_name"]=="Angiospermae (flowering plants)"
    assert r["queries_total"]==35 and len(r["cohorts"])==7
    assert r["no_iNaturalist_image_pixels_opened"] is True
    assert r["actual_original_source_FCP_42111_taxon_frame_sampled"] is False
    assert r["original_93209351_photo_ID_manifest_used"] is False
    assert r["annotations_assumed_frozen_as_of_2026_09_27"] is False
    assert r["whole_plant_population_or_all_flowering_species_inferred"] is False
    return r

def test_angiosperm_real_public_photo_observation_universe_and_open_flower_tags():
    w=source()["cohorts"]["world_2025_observed"]
    assert w["total_observed_plant_photo_records"]==26069247
    assert w["with_any_phenology_annotation"]==2550231
    assert w["explicit_flowers_observations"]==1615314
    assert w["explicit_no_flowers_fruits_observations"]==424672
    assert w["without_phenology_annotation_filter_records"]==23518398
    assert math.isclose(w["any_annotation_fraction_observed"],2550231/26069247)
    assert math.isclose(w["explicit_flowers_fraction_observed"],1615314/26069247)
    assert round(100*w["explicit_flowers_fraction_observed"],2)==6.20
    assert w["any_plus_without_minus_total"]==-618
    assert w["any_plus_without_partitions_exactly"] is False
    assert w["unannotated_is_not_nonflowering"] is True

def test_angiosperm_uses_proper_flowering_taxa_not_fern_and_moss_plant_denominator():
    allp=json.loads(PLANT.read_text())
    ang=source()
    x=allp["cohorts"]["world_2025_observed"]
    y=ang["cohorts"]["world_2025_observed"]
    assert allp["source_plant_taxon_id"]==47126
    assert x["total_observed_plant_photo_records"]==28565884
    assert y["total_observed_plant_photo_records"]<x["total_observed_plant_photo_records"]
    assert y["explicit_flowers_fraction_observed"]>x["explicit_flowers_fraction_observed"]
    # Queries not atomic: do NOT treat small query-count changes as matched
    # observation pairs or evidence of source-photo FCP organism inclusion.
    assert y["with_any_phenology_annotation"]!=x["with_any_phenology_annotation"]

def test_fixed_hierarchical_source_phenology_difference_is_not_color_adaptation():
    r=source()["cohorts"]
    truth={
       "northern_east_asia_spring_2025":(47063,6568,5351),
       "northern_east_asia_summer_2025":(35072,4481,3245),
       "southern_australia_spring_2025":(379937,84266,67646),
       "southern_australia_summer_2024_25":(134793,25636,16502),
       "tropical_south_america_calendar_q1_2025":(24649,2497,1519),
       "tropical_south_america_calendar_q3_2025":(26648,1648,1105)
    }
    for name,t in truth.items():
        v=r[name]
        assert (v["total_observed_plant_photo_records"],
                v["with_any_phenology_annotation"],
                v["explicit_flowers_observations"])==t
        assert len(v["source_query_provenance"])==5
        for request in v["source_query_provenance"].values():
            assert parse_qs(urlsplit(request["request_url"]).query)["taxon_id"]==["47125"]
            assert len(request["original_api_response_sha256"])==64
    for a,b in (
       ("northern_east_asia_spring_2025","northern_east_asia_summer_2025"),
       ("southern_australia_spring_2025","southern_australia_summer_2024_25"),
       ("tropical_south_america_calendar_q1_2025","tropical_south_america_calendar_q3_2025")
    ):
        assert r[a]["explicit_flowers_fraction_observed"]>r[b]["explicit_flowers_fraction_observed"]
