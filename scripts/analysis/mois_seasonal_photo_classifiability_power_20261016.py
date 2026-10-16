#!/usr/bin/env python3
"""Binomial PHOTO CLASSIFIABILITY planning, never flowering-colour inference.

Reuse only VERIFIED older Mois historical photographic denominators:
100,543 original unique opportunities; 45,953 classified 4-hue photographs.
For the original dry and wet opportunity deciles, the source support is
4,981/10,084 and 3,904/9,990 classifiable photos, respectively.

For a hypothetical independent-photo classifiability p, ask how often an
intact source selection of n observations PER season yields >=k=3 original
photos that ROI-v4 can classify IN BOTH seasons. p may differ by season,
taxon, observer, site, method or source. This is NOT a confidence interval,
and does NOT measure a new, independent real source population.
"""
from __future__ import annotations
import argparse
import json
import math
from math import comb
from pathlib import Path

SCHEMA="mois_prespecified_photo_classifiability_two_season_design_power_v1"
STATUS="HISTORICAL_PHOTO_NONRESPONSE_BINOMIAL_PLANNING_ONLY_NO_NEW_PHOTOS"
HISTORICAL={
    "all_source":(45953,100543),
    "dry_photo_opportunity_decile":(4981,10084),
    "wet_photo_opportunity_decile":(3904,9990),
}
GOAL=3
CAP=40
SIZES=(3,5,7,8,9,10,11,12,14,16,20)
TARGETS=(.8,.9,.95)


def classified_probability(n:int,p:float,k:int=GOAL)->float:
    if not isinstance(n,int) or n<0 or not 0<=p<=1 or not 1<=k<=CAP:
        raise ValueError("Invalid independent-photo binomial planning inputs")
    if n<k:
        return 0.
    return math.fsum(comb(n,j)*p**j*(1-p)**(n-j) for j in range(k,n+1))


def probability_two_seasons(n1:int,n2:int,p1:float,p2:float,
                            k:int=GOAL)->float:
    return classified_probability(n1,p1,k)*classified_probability(n2,p2,k)


def minimum_per_phase(p1:float,p2:float,threshold:float)->int|None:
    if not 0<threshold<1:
        raise ValueError("Two-season outcome gate must be a nontrivial probability")
    for n in range(GOAL,CAP//2+1):
        if probability_two_seasons(n,n,p1,p2)>=threshold:
            return n
    return None


def report()->dict:
    rates={k:n/d for k,(n,d) in HISTORICAL.items()}
    matrix={}
    for group,p in rates.items():
        matrix[group]={
            "previously_classifiable_photo_fraction":p,
            "scenarios":[{
                "independent_original_observations_per_season":n,
                "total_both_seasons":2*n,
                "probability_both_seasons_have_at_least_three_classifiable":
                    probability_two_seasons(n,n,p,p)
            } for n in SIZES],
            "minimum_equal_original_observations_per_season_for_pair_success":{
                f"{threshold:.2f}":minimum_per_phase(p,p,threshold)
                for threshold in TARGETS
            }
        }
    dry=rates["dry_photo_opportunity_decile"]
    wet=rates["wet_photo_opportunity_decile"]
    mixed={
        "dry_photo_classifiability_input":dry,
        "wet_photo_classifiability_input":wet,
        "minimum_equal_source_observations_per_dry_and_wet_phase":{
            f"{threshold:.2f}":minimum_per_phase(dry,wet,threshold)
            for threshold in TARGETS
        },
        "scenarios":[{
            "observations_dry":n,"observations_wet":n,
            "probability_both_have_at_least_three_classifiable":
                probability_two_seasons(n,n,dry,wet)
        } for n in SIZES]
    }
    return{
        "schema":SCHEMA,"status":STATUS,
        "source_historical_photographs_unique":100543,
        "historical_four_hue_classifiable_photographs":45953,
        "historical_photo_records_unclassifiable":100543-45953,
        "previous_original_source_dry_decile_classifiable_n_of_total":
            list(HISTORICAL["dry_photo_opportunity_decile"]),
        "previous_original_source_wet_decile_classifiable_n_of_total":
            list(HISTORICAL["wet_photo_opportunity_decile"]),
        "model":"two independent Binomial(n,p) four-hue classifiable counts",
        "compared_source_seasons":2,
        "target_classifiable_in_each_season":GOAL,
        "prespecified_original_observation_cap_per_species":CAP,
        "same_historical_classifiability_rate_scenarios":matrix,
        "unequal_original_dry_wet_classifiability_scenarios":mixed,
        "proposal_min_original_selected_observations_per_occupied_season_when_quota_allows":10,
        "independence_and_constant_p_are_hypothetical_not_empirical":True,
        "no_new_photo_or_phenology_annotation_labels_measured":True,
        "no_seasonal_white_climate_effect_computed":True,
        "no_genetic_or_wild_flora_frequency_claim":True,
        "note":(
            "Planning probabilities treat original observation-level ROI successes "
            "as independent Bernoulli trials with fixed p. Old SOURCE photo-level "
            "classifiability is not necessarily the new observation-level, "
            "seasonal or wet/dry classifiability. Seasonal selection and "
            "annotation may change p; these are sensitivity scenarios, not "
            "forecasts or confidence intervals."
        )
    }


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",required=True,type=Path)
    args=parser.parse_args()
    if args.output.exists():
        raise ValueError("Never silently replace an earlier design power receipt")
    args.output.parent.mkdir(parents=True,exist_ok=True)
    data=report()
    args.output.write_text(json.dumps(data,indent=2,sort_keys=True)+"\n",encoding="utf8")
    print(json.dumps({k:data[k] for k in (
        "status","source_historical_photographs_unique",
        "historical_four_hue_classifiable_photographs",
        "proposal_min_original_selected_observations_per_occupied_season_when_quota_allows"
    )},sort_keys=True))


if __name__=="__main__":
    main()
