#!/usr/bin/env python3
"""Bounded REAL iNaturalist public-API phenology annotation COUNTS pilot.

Counts only, NEVER image pixels, downloaded observation records, names, user
handles or individual coordinates. Five 1-record API queries per cohort,
strict global spacing 1.35 s and <=35 public requests total.

This is a 2026 API query about photographed Plantae observed during selected
2024–2025 date windows, NOT the frozen 2026-09-27 original FCP 42,111-taxon
93,209,351-photo frame. The whole-platform and fixed-geographic comparison
regions were selected before viewing output; their results cannot be passed
to Mois original-photo-ID source admission or paper's 13-variable regression.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

API = "https://api.inaturalist.org/v1/observations"
SCHEMA = "mois_bounded_public_inat_photo_phenology_annotation_count_pilot_v1"
STATUS = "REAL_INAT_API_PUBLIC_PHOTO_OBSERVATION_PHENOLOGY_COUNTS_NOT_FCP_SOURCE"
TAXON = 47126  # Original comparison: Plantae.
ANGIOSPERMS = 47125  # Angiospermae, proper flowering-plant frame.
MAX_QUERIES = 35
QUERY_MIN_DELAY_SEC = 1.35
FILTERS = {
    "all_plant_photo_observations": {},
    "any_flowers_fruits_annotation": {"term_id": 12},
    "explicit_flowers_annotation": {"term_id": 12, "term_value_id": 13},
    "explicit_no_flowers_fruits_annotation": {"term_id": 12, "term_value_id": 21},
    "without_flowers_fruits_annotation": {"without_term_id": 12},
}
# Two hemispheres each use LOCAL seasons, and a tropical monthly contrast
# deliberately does NOT equate Jan or Jul with rain/dry conditions.
COHORTS = [
    {"name":"world_2025_observed","d1":"2025-01-01","d2":"2025-12-31"},
    {"name":"northern_east_asia_spring_2025","d1":"2025-03-01","d2":"2025-05-31",
     "swlat":31,"swlng":129,"nelat":45,"nelng":145},
    {"name":"northern_east_asia_summer_2025","d1":"2025-06-01","d2":"2025-08-31",
     "swlat":31,"swlng":129,"nelat":45,"nelng":145},
    {"name":"southern_australia_spring_2025","d1":"2025-09-01","d2":"2025-11-30",
     "swlat":-38,"swlng":113,"nelat":-22,"nelng":153},
    {"name":"southern_australia_summer_2024_25","d1":"2024-12-01","d2":"2025-02-28",
     "swlat":-38,"swlng":113,"nelat":-22,"nelng":153},
    {"name":"tropical_south_america_calendar_q1_2025","d1":"2025-01-01","d2":"2025-03-31",
     "swlat":-10,"swlng":-75,"nelat":10,"nelng":-45},
    {"name":"tropical_south_america_calendar_q3_2025","d1":"2025-07-01","d2":"2025-09-30",
     "swlat":-10,"swlng":-75,"nelat":10,"nelng":-45},
]


def validate_cohorts(cohorts:list[dict]) -> None:
    if not 1<=len(cohorts)<=7 or len(cohorts)*len(FILTERS)>MAX_QUERIES:
        raise ValueError("Exceeded authorized small iNaturalist API pilot quota")
    if len({x["name"] for x in cohorts})!=len(cohorts):
        raise ValueError("Repeated exploratory photo region/date cohort")
    for row in cohorts:
        if set(row)-{"name","d1","d2","swlat","swlng","nelat","nelng"}:
            raise ValueError("Unreviewed iNaturalist observation search parameter")
        if not row.get("name") or not row.get("d1") or not row.get("d2"):
            raise ValueError("Missing original observed-on date window")
        if set(row)&{"swlat","swlng","nelat","nelng"}:
            if not {"swlat","swlng","nelat","nelng"}.issubset(row):
                raise ValueError("Partial geographic bounding box is not allowed")
            if not (-90<=row["swlat"]<row["nelat"]<=90
                    and -180<=row["swlng"]<row["nelng"]<=180):
                raise ValueError("Invalid original observation regional pilot bounding box")


def uri_query(cohort:dict, extra:dict, *, taxon_id:int=TAXON)->str:
    if taxon_id not in (TAXON,ANGIOSPERMS):
        raise ValueError("Unknown/unreviewed source Plant/Angiospermae taxon")
    params={"taxon_id":taxon_id,"photos":"true","per_page":1,"page":1,
            **{k:v for k,v in cohort.items() if k!="name"},**extra}
    return API+"?"+urllib.parse.urlencode(params)


def query_count(url:str, opener=None)->dict:
    if not url.startswith(API+"?"):
        raise ValueError("Never contact a non-iNaturalist provider in this pilot")
    req=urllib.request.Request(
        url,headers={
            "User-Agent":"Mois-global-flower-colour-research/1.0 (+https://github.com/zuizui0223/mois)",
            "Accept":"application/json",
        }
    )
    if opener is None:opener=urllib.request.urlopen
    with opener(req,timeout=35) as response:
        if getattr(response,"status",200)!=200:
            raise ValueError("Original iNaturalist response was not HTTP 200")
        raw=response.read(2*1024*1024+1)
    if len(raw)>2*1024*1024:
        raise ValueError("Unexpected oversized iNaturalist count-only response")
    payload=json.loads(raw)
    if (not isinstance(payload,dict)
        or "total_results" not in payload
        or not isinstance(payload["total_results"],int)
        or payload["total_results"]<0
        or not isinstance(payload.get("results"),list)
        or len(payload["results"])>1):
        raise ValueError("Could not independently verify iNaturalist paginated count-only response")
    return {
        "n_observations":payload["total_results"],
        "raw_response_sha256":hashlib.sha256(raw).hexdigest()
    }


def summarize(counts:dict)->dict:
    if set(counts)!=set(FILTERS) or any(
        not isinstance(n,int) or n<0 for n in counts.values()
    ):
        raise ValueError("Missing a prospective observation annotation category")
    total=counts["all_plant_photo_observations"]
    any_=counts["any_flowers_fruits_annotation"]
    yes=counts["explicit_flowers_annotation"]
    no=counts["explicit_no_flowers_fruits_annotation"]
    unannotated=counts["without_flowers_fruits_annotation"]
    if total<=0 or max(any_,yes,no,unannotated)>total or yes>any_ or no>any_:
        raise ValueError("API phenology classes cannot be larger than their observed-frame denominator")
    # Positive Flowers + positive No-Flowers could in principle overlap if
    # a user flags contradictory annotations. NEVER assume they partition the
    # observed photo opportunities. Same for no-annotation API filters: if
    # counts drift due query semantics, retain the mismatch as a diagnostic.
    return {
        "total_observed_plant_photo_records":total,
        "with_any_phenology_annotation":any_,
        "explicit_flowers_observations":yes,
        "explicit_no_flowers_fruits_observations":no,
        "without_phenology_annotation_filter_records":unannotated,
        "any_annotation_fraction_observed":any_/total,
        "explicit_flowers_fraction_observed":yes/total,
        "no_annotation_filter_fraction_observed":unannotated/total,
        "any_plus_without_minus_total":any_+unannotated-total,
        "any_plus_without_partitions_exactly":any_+unannotated==total,
        "unannotated_is_not_nonflowering":True,
    }


def run(out_path:Path,*,cohorts=COHORTS,get_count=query_count,sleep=time.sleep,
        taxon_id:int=TAXON)->dict:
    validate_cohorts(cohorts)
    if taxon_id not in (TAXON,ANGIOSPERMS):
        raise ValueError("Unreviewed public platform taxon group cannot replace original source frame")
    if out_path.exists() or out_path.with_name(out_path.name+".INCOMPLETE").exists():
        raise ValueError("No overwrite or interrupted iNaturalist count audit")
    out_path.parent.mkdir(parents=True,exist_ok=True)
    all_rows={}
    last_request_time=None
    for cohort in cohorts:
        raw={}
        receipt={}
        for label,term in FILTERS.items():
            if last_request_time is not None:
                # API guidance requests ~one per second. This pilot is more
                # conservative and NEVER starts concurrent requests.
                sleep(QUERY_MIN_DELAY_SEC)
            url=uri_query(cohort,term,taxon_id=taxon_id)
            source=get_count(url)
            last_request_time=True
            if set(source)!={"n_observations","raw_response_sha256"}:
                raise ValueError("No exact count and original response digest")
            raw[label]=source["n_observations"]
            receipt[label]={
                "request_url":url,
                "original_api_response_sha256":source["raw_response_sha256"],
                "n_observations":source["n_observations"],
            }
        all_rows[cohort["name"]]={
            "filters":{k:v for k,v in cohort.items() if k!="name"},
            **summarize(raw),
            "source_query_provenance":receipt,
        }
    result={
        "schema":SCHEMA,
        "status":STATUS,
        "provider":"iNaturalist current public API v1, observations total_results",
        "queries_total":len(cohorts)*len(FILTERS),
        "queries_rate_limited_and_sequential":True,
        "minimum_inter_request_spacing_seconds":QUERY_MIN_DELAY_SEC,
        "source_plant_taxon_id":taxon_id,
        "source_plant_group_name":("Angiospermae (flowering plants)" if
                                   taxon_id==ANGIOSPERMS else "Plantae (all plants)"),
        "source_photo_only":True,
        "observation_date_windows":True,
        "retrieved_on_utc":datetime.now(timezone.utc).isoformat(),
        "cohorts":all_rows,
        "actual_original_source_FCP_42111_taxon_frame_sampled":False,
        "original_93209351_photo_ID_manifest_used":False,
        "annotations_assumed_frozen_as_of_2026_09_27":False,
        "whole_plant_population_or_all_flowering_species_inferred":False,
        "flowers_observation_annotation_does_not_prove_selected_image_ROI":True,
        "no_iNaturalist_image_pixels_opened":True,
        "data_volume":"One observation JSON response per filtered aggregate query; raw individual responses discarded immediately; only counts and source digests persisted.",
        "limitation":"Current public iNaturalist annotation data in selected 2024-2025 observation-date windows; this is NOT a representative source-photo-ID pilot of the FCP 2026-09-27 42,111-taxon original archive, and API observations need not have licensed source photos. Geographic boxes, plant composition and observer use differ.",
    }
    candidate=out_path.with_name(out_path.name+".INCOMPLETE")
    candidate.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    candidate.replace(out_path)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-aggregate",type=Path,required=True)
    parser.add_argument("--taxon-id",type=int,choices=[TAXON,ANGIOSPERMS],default=TAXON)
    args=parser.parse_args()
    result=run(args.out_aggregate,taxon_id=args.taxon_id)
    print("MOIS_REAL_INAT_PHENOLOGY_COUNT_PILOT_ONLY_NO_FLOWER_COLOUR",
          json.dumps({name:{k:v[k] for k in (
              "total_observed_plant_photo_records",
              "with_any_phenology_annotation",
              "explicit_flowers_observations",
              "without_phenology_annotation_filter_records",
              "any_annotation_fraction_observed",
              "explicit_flowers_fraction_observed",
              "any_plus_without_partitions_exactly",
          )} for name,v in result["cohorts"].items()},sort_keys=True))


if __name__=="__main__":
    main()
