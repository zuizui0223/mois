#!/usr/bin/env python3
"""Metadata-only upper bounds on Mois repeated-species seasonal comparisons.

Reads the frozen, SHA-audited <=40 original observations per taxon selected
WITHOUT opening flower pixels or annotations. Reconciles all selected photo
identities and stratification before counting season-matched taxon opportunities.

A nominal taxon with >=3 sampled OBSERVATIONS in both local seasons is only a
possible future >=3 *four-hue-classifiable* taxon. None of these support counts
is a floral phenotype, iNaturalist flowering-annotation count, or evidence of
season-dependent pigmentation; duplicate photographers also limit inference.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0,str(ROOT))

from scripts.acquisition.mois_stratified_photo_selection_20261016 import (
    PRIVATE_FIELDS, SOURCE_SPECIES, SOURCE_PHOTOS, CAP,
    classify_site, sha256_file,
)
from scripts.analysis.mois_seasonal_flower_colour_opportunity_20261016 import (
    LOCAL, local_season,
)

SCHEMA="mois_seasonal_photo_sample_preclassification_support_v1"
STATUS="SELECTED_OBSERVATION_SEASON_SUPPORT_ONLY_NO_PHOTO_PIXELS"
LOCAL_PAIRS=tuple(combinations(LOCAL,2))
THRESHOLDS=(3,5,7,10)
MIN_OBSERVERS=2


def key_summary(source:dict)->dict:
    return {
        "source_snapshot_taxa":source["source_species"],
        "source_original_photo_opportunities":source["source_photo_ids_authenticated"],
        "selected_original_photo_ids":source["selected_original_photo_ids"],
        "photo_selection_seed_already_sealed":True,
    }


def _pair_count(groups:dict, minimum:int)->list[dict]:
    domain_pairs=defaultdict(lambda:{"n_nominal_comparisons":0,
                                     "n_with_shared_observer":0})
    for identity,seasons in groups.items():
        domain=identity[-1]
        if domain not in ("NORTH_EXTRATROPICAL","SOUTH_EXTRATROPICAL"):
            continue
        for first,second in LOCAL_PAIRS:
            a=seasons.get(first)
            b=seasons.get(second)
            if a is None or b is None:
                continue
            if (a["n"]<minimum or b["n"]<minimum
                or len(a["observers"])<MIN_OBSERVERS
                or len(b["observers"])<MIN_OBSERVERS):
                continue
            group=domain_pairs[(domain,first,second)]
            group["n_nominal_comparisons"]+=1
            group["n_with_shared_observer"]+=bool(a["observers"] & b["observers"])
    return [
        {"domain":domain,"season_1":first,"season_2":second,**values}
        for (domain,first,second),values in sorted(domain_pairs.items())
    ]


def audit(sample_csv:Path,receipt_json:Path,private_json:Path,public_json:Path,
          *,synthetic:bool=False)->dict:
    if (private_json.exists() or public_json.exists()
        or private_json.with_name(private_json.name+".partial").exists()
        or public_json.with_name(public_json.name+".partial").exists()):
        raise ValueError("Do not overwrite or use an interrupted seasonal support preflight")
    source=json.loads(receipt_json.read_text(encoding="utf8"))
    if (source.get("schema")!="mois_direct_fcp_source_stratified_photo_selection_v1"
        or source.get("status")!=
            "SOURCE_ADMITTED_SELECTED_PRIVATE_IDS_ONLY_NO_PHOTO_PIXELS_MEASURED"
        or source.get("selected_private_photo_csv_sha256")!=sha256_file(sample_csv)
        or source.get("selection_photo_cap_per_taxon")!=CAP
        or source.get("historical_colour_outcomes_seen") is not False
        or source.get("pixels_opened_or_colours_measured") is not False
        or source.get("original_13var_environment_join_complete") is not False):
        raise ValueError("Only an original-source-certified colour-blind photo sample is allowed")
    if not synthetic and (
        source.get("source_species")!=SOURCE_SPECIES
        or source.get("source_photo_ids_authenticated")!=SOURCE_PHOTOS):
        raise ValueError("Not the original real 42,111-taxon source population")

    species=defaultdict(lambda:defaultdict(lambda:{"n":0,"observers":set()}))
    geo=defaultdict(lambda:defaultdict(lambda:{"n":0,"observers":set()}))
    season_domains=Counter()
    missing_month=Counter()
    chosen_taxa=set()
    seen_photo=set()
    seen_observation=set()
    chosen=0
    taxa_season_support=Counter()
    with sample_csv.open(newline="",encoding="utf8") as handle:
        reader=csv.DictReader(handle)
        if reader.fieldnames!=PRIVATE_FIELDS:
            raise ValueError("Mois source-photo selection fields changed")
        for row in reader:
            tax,gen,photo,obs=(row[name].strip() for name in (
                "inat_taxon_id","genus","photo_id","observation_id"))
            if (not tax.isdecimal() or int(tax)<1 or not gen or not photo
                or photo in seen_photo or obs in seen_observation or not obs):
                raise ValueError("Duplicate/missing original selected photo or distinct observation")
            if int(row["technical_pilot_flag"]) not in (0,1):
                raise ValueError("Technical pilot was not outcome-blind and binary")
            if classify_site(row)!=(row["region_stratum"],row["precip_stratum"]):
                raise ValueError("Original source geographic/BIO12 selection stratum changed")
            domain,season,month,year=local_season(row["latitude"],row["observed_on"])
            pi=float(row["observation_inclusion_probability"])
            pw=float(row["photo_design_weight"])
            if (not math.isfinite(pi) or not 0<pi<=1 or
                not math.isfinite(pw) or pw<1):
                raise ValueError("Original photo selection probability drifted")
            observer=row["observer_id"].strip() or "UNKNOWN_OBSERVER"
            x=species[(tax,gen,domain)][season]
            y=geo[(tax,gen,row["region_stratum"],row["precip_stratum"],domain)][season]
            for item in (x,y):
                item["n"]+=1
                if observer!="UNKNOWN_OBSERVER":
                    item["observers"].add(observer)
            season_domains[(domain,season)]+=1
            if month=="UNKNOWN_MONTH":
                missing_month[domain]+=1
            taxa_season_support[(tax,domain,season)]+=1
            seen_observation.add(obs)
            seen_photo.add(photo)
            chosen_taxa.add(tax)
            chosen+=1
    if (chosen!=source["selected_original_photo_ids"]
        or len(chosen_taxa)!=source["source_species_with_photos"]
        or chosen>CAP*source["source_species_with_photos"]):
        raise ValueError("Sample species/photo counts differ from frozen source sampling receipt")
    by_tax=Counter()
    for (tax,_,domain),phases in species.items():
        by_tax[tax]+=sum(z["n"] for z in phases.values())
    if any(n>CAP for n in by_tax.values()):
        raise ValueError("Taxon exceeds original maximum 40 distinct observations")

    private={
        "schema":SCHEMA+"_private",
        "status":STATUS,
        "photo_sample_sha256":sha256_file(sample_csv),
        "frozen_source_receipt_sha256":sha256_file(receipt_json),
        "season_counts_unsuppressed":[
            {"domain":d,"season":s,"n":n}
            for (d,s),n in sorted(season_domains.items())],
        "matched_taxon_pairs_by_min_original_observations":{
            str(th):_pair_count(species,th) for th in THRESHOLDS},
        "matched_taxon_same_geographic_precip_stratum_pairs_by_min_original_observations":{
            str(th):_pair_count(geo,th) for th in THRESHOLDS},
        "source_taxon_ids_available_in_private_only":True,
        "source_photo_pixels_opened":False,
        "colour_outcomes_or_annotations_used":False,
    }
    public={
        "schema":SCHEMA,
        "status":STATUS,
        **key_summary(source),
        "preclassification_selected_unique_observations":chosen,
        "selected_distinct_taxa":len(chosen_taxa),
        "selected_missing_observation_month":sum(missing_month.values()),
        "selected_hemisphere_tropical_domain_support":[
            {"domain":d,"selected_observations":sum(n for (dd,s),n in season_domains.items() if dd==d)}
            for d in ("NORTH_EXTRATROPICAL","SOUTH_EXTRATROPICAL",
                      "TROPICAL","UNKNOWN_GEOGRAPHY")],
        "matched_taxon_season_support_upper_bounds":{
            str(th):_pair_count(species,th) for th in THRESHOLDS},
        "matched_taxon_same_geographic_precip_stratum_support_upper_bounds":{
            str(th):_pair_count(geo,th) for th in THRESHOLDS},
        "sample_original_photo_observation_cap_per_taxon":CAP,
        "at_least_two_source_observers_per_compared_local_season_required":True,
        "n_source_taxon_with_zero_photos":source["source_species"]-source["source_species_with_photos"],
        "real_classifiable_photos_measured":0,
        "annotation_positive_flowers_count_measured":False,
        "monthly_tropical_wet_dry_defined_here":False,
        "same_species_across_seasons_does_not_identify_genetic_colour_plasticity":True,
        "matched_counts_are_only_upper_bounds_before_flower_roi_classification":True,
        "synthetic_test_only":bool(synthetic),
        "note":"Counted only potentially matched selected source observations; >=3 selected per season is a necessary but NOT sufficient condition for >=3 actually classifiable flower photos. Does not identify the global flora or new 13-environment slopes."
    }
    for p in (private_json,public_json):
        p.parent.mkdir(parents=True,exist_ok=True)
    tmp_pr=private_json.with_name(private_json.name+".partial")
    tmp_pub=public_json.with_name(public_json.name+".partial")
    try:
        tmp_pr.write_text(json.dumps(private,indent=2,sort_keys=True)+"\n",encoding="utf8")
        tmp_pub.write_text(json.dumps(public,indent=2,sort_keys=True)+"\n",encoding="utf8")
        tmp_pr.replace(private_json)
        tmp_pub.replace(public_json)
    except Exception:
        for p in (tmp_pr,tmp_pub):
            if p.exists():p.unlink()
        raise
    return public


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--frozen-private-sample",type=Path,required=True)
    parser.add_argument("--frozen-sampling-receipt",type=Path,required=True)
    parser.add_argument("--out-private-support",type=Path,required=True)
    parser.add_argument("--out-public-aggregate",type=Path,required=True)
    a=parser.parse_args()
    result=audit(a.frozen_private_sample,a.frozen_sampling_receipt,
                 a.out_private_support,a.out_public_aggregate)
    print(json.dumps({key:result[key] for key in (
        "status","preclassification_selected_unique_observations",
        "selected_distinct_taxa","matched_taxon_season_support_upper_bounds")},
        sort_keys=True))


if __name__=="__main__":
    main()
