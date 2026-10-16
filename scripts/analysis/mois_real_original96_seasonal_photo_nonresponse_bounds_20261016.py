#!/usr/bin/env python3
"""Actual 96 original FCP photographed-observation season × ROI missingness.

The pilot includes 96 actual measured photo pixels / four nonrandom taxa,
20 four-hue-classifiable photos / three classified white. Sharp all-photo
BINARY ASSIGNMENT bounds are only a mathematical missing-data sensitivity:
many nonclassified images have NO visible flower, so their counterfactual
white value may be undefined. Do not call the bounds identified flowering
community prevalence or climate-induced pigmentation.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path

SCHEMA="mois_original96_technical_photo_season_classifiability_mnar_bounds_v1"
STATUS="REAL_ORIGINAL96_PHOTOGRAPHIC_CLASSIFIABILITY_MISSINGNESS_NOT_FLORA_PREVALENCE"
SOURCE_SHA="911ec9edd10a7126bee28cfa67738cef0ff4d275b1f0ca1c06e3f26e7e290039"
CLIMATE_MANIFEST_SHA="9d0cf37a64534542d7788debef28c13f106a15c3bf7bf91be0653f3ce106baae"


def derive(src:dict)->dict:
    if (src.get("status")!="ACTUAL_96_IMAGE_ROI_TECHNICAL_PILOT_CLIMATE_JOIN_NOT_GLOBAL_SEASON_EFFECT"
        or src.get("source_original_FCP_96_pilot_csv_sha256")!=SOURCE_SHA
        or src.get("source_48_worldclim_monthly_manifest_sha256")!=CLIMATE_MANIFEST_SHA
        or src.get("n_original_source_photo_IDs")!=96
        or src.get("n_original_distinct_observations")!=96
        or src.get("n_original_taxa")!=4
        or src.get("n_roi_original_four_colour_classifiable")!=20
        or src.get("n_roi_original_classified_white")!=3
        or src.get("n_monthly_source_site_climate_joined")!=96
        or src.get("original_source_population_42111_taxa_photo_IDs_joined") is not False
        or src.get("seasonal_climate_white_association_or_fitness_estimated") is not False
        or src.get("synthetic_test_only") is not False):
        raise ValueError("Original 96 remeasured flower-photo source receipt not admitted")
    groups={}
    for entry in src["observed_season_rows"]:
        key=(entry["domain"],entry["local_season"])
        if key in groups:
            raise ValueError("Original local-season bin duplicated")
        n=entry["sampled_original_photo_opportunities"]
        c=entry["four_hue_classifiable"]
        white=entry["classified_white"]
        if not (0<=white<=c<=n) or n<1:
            raise ValueError("Observed image/photo denominator invalid")
        missing=n-c
        groups[key]={
            "source_photo_opportunities":n,
            "four_hue_classifiable_photographs":c,
            "classified_white_photographs":white,
            "source_unclassifiable_photographs":missing,
            "roi_four_colour_classifiability_fraction":c/n,
            "white_fraction_given_roi_classifiable":white/c if c else None,
            "mathematical_all_photo_binary_assignment_lower_white_fraction":white/n,
            "mathematical_all_photo_binary_assignment_upper_white_fraction":(white+missing)/n,
            "source_classified_climate_group_averages":entry["source_classified_monthly_climate_means"]
        }
    if (sum(d["source_photo_opportunities"] for d in groups.values())!=96
        or sum(d["four_hue_classifiable_photographs"] for d in groups.values())!=20
        or sum(d["classified_white_photographs"] for d in groups.values())!=3):
        raise ValueError("Seasonal photo opportunity assignment changed from complete real pilot")
    public=[]
    for (domain,season),d in sorted(groups.items()):
        # Keep unknown/unclassifiable origin even in small groups; however
        # do not promote 1-image observations into a stable seasonal white
        # probability or climate-comparison result.
        row={"domain":domain,"local_season":season,**d}
        row["sufficient_for_a_descriptive_season_classifiability_fraction"]=d[
            "source_photo_opportunities"]>=10
        row["sufficient_for_season_white_climate_inference"]=False
        public.append(row)
    spring=groups[("NORTH_EXTRATROPICAL","SPRING")]
    summer=groups[("NORTH_EXTRATROPICAL","SUMMER")]
    conditional=spring["white_fraction_given_roi_classifiable"]-summer[
        "white_fraction_given_roi_classifiable"]
    # Logical extreme assignments of ALL currently unclassifiable source
    # photographs; not a feasible plant-level white phenotype census.
    lower=(spring["mathematical_all_photo_binary_assignment_lower_white_fraction"]
           -summer["mathematical_all_photo_binary_assignment_upper_white_fraction"])
    upper=(spring["mathematical_all_photo_binary_assignment_upper_white_fraction"]
           -summer["mathematical_all_photo_binary_assignment_lower_white_fraction"])
    assert lower<=conditional<=upper
    return {
        "schema":SCHEMA,"status":STATUS,
        "real_original96_source_sha256":SOURCE_SHA,
        "real_48_climate_source_manifest_sha256":CLIMATE_MANIFEST_SHA,
        "real_source_photo_month_run":src["verified_original_FCP_photo_climate_execution_url"],
        "real_source_nonrepresentative_taxa":4,
        "all_original_photographs":96,
        "original_roi_classifiable":20,
        "original_roi_white":3,
        "photographs_without_classifiable_4hue":76,
        "source_season_counts_and_logical_missing_photo_bounds":public,
        "north_spring_minus_summer":{
            "source_white_given_classifiable_difference_percentage_points":100*conditional,
            "worst_case_all_photo_binary_assignment_difference_percentage_points":[
                100*lower,100*upper
            ],
            "both_signs_of_raw_all_photo_binary_assignment_difference_possible":
                lower<0<upper,
            "is_comparison_of_same_taxa_and_same_geographic_cell":False,
            "no_temperature_solar_rainfall_VPD_causal_effect_estimated":True
        },
        "missing_photo_binary_colour_assignments_not_valid_plant_flowering_traits":True,
        "overall_image_pilot_is_not_random_42111_taxon_photographic_census":True,
        "all_model_selection_and_geo_month_climate_fit_missingness_unresolved":True,
        "note":"A no-flower photo does not intrinsically have a white/coloured flower trait. Extreme-assignment bounds illustrate the inferential loss in photographic opportunities only; they must not be interpreted as estimated wild-plant floral colour frequencies or a statistical confidence interval."
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--actual-original96-monthly-pilot-receipt",required=True,type=Path)
    p.add_argument("--out-aggregate",required=True,type=Path)
    args=p.parse_args()
    if args.out_aggregate.exists():
        raise ValueError("Refuse to overwrite source-assured original-photo ascertainment summary")
    src=json.loads(args.actual_original96_monthly_pilot_receipt.read_text())
    result=derive(src)
    args.out_aggregate.parent.mkdir(parents=True,exist_ok=True)
    args.out_aggregate.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n",encoding="utf8")
    print(json.dumps(result,sort_keys=True))


if __name__=="__main__":main()
