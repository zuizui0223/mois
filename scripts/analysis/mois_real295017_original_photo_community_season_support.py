#!/usr/bin/env python3
"""Same-source photo-assemblage season contrast, after REAL 295017 photo union.

The original photographed taxa in summer and winter may differ. This diagnoses
coarse-cell source community colour, not within-plant pigment plasticity.
Importantly, the qualified observed dates are independently verified and the
original frozen colour labels and all original technical failures are retained.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

from scripts.analysis.mois_seasonal_267k_climate_normals import (
    HUES, labelled_calendar, require,
)
from scripts.analysis.mois_original_photo_community_season_geoscale import (
    matched_community, cell_cluster_interval,
)
from scripts.analysis.mois_merge_real_FCP180_and_MOIS_three4096_original_photos import (
    DATE_FILE, SOURCE_RECEIPT,
)

EXPECTED_ROWS=295017
EXPECTED_HUES=142538
EXPECTED_DATE_ROWS=289170
EXPECTED_DATED_HUES=140247
SOURCE="mois_real295017_original_photo_hues_plus_verified_observed_dates.csv.gz"
WIDTHS=(.25,.5,1.,2.)
N_MIN=5
TAX_MIN=2


def check_read(source:Path,receipt:Path):
    source_proof=json.loads(receipt.read_text())
    sha=hashlib.sha256(source.read_bytes()).hexdigest()
    require(sha==source_proof["files_SHA256"][SOURCE] and
            source_proof["original_measured_photographs"]==EXPECTED_ROWS and
            source_proof["original_measured_photo_four_hue_labels"]==EXPECTED_HUES,
            "Fully source-locked 295017 actual photographed hues and dates required")
    d=pd.read_csv(source,compression="gzip",
        usecols=["photo_id","observation_id","inat_taxon_id",
                 "morph","latitude","longitude","qualified_date"],
        dtype={"photo_id":"string","observation_id":"string",
               "inat_taxon_id":"string"},low_memory=False)
    require(len(d)==EXPECTED_ROWS and d.photo_id.is_unique and
            d.morph.isin(HUES).sum()==EXPECTED_HUES and
            d.qualified_date.notna().sum()==EXPECTED_DATE_ROWS and
            (d.qualified_date.notna()&d.morph.isin(HUES)).sum()==EXPECTED_DATED_HUES,
            "Original immutable photo/day/classification counts drift")
    d=d.rename(columns={"qualified_date":"observed_on"})
    return labelled_calendar(d),sha


def run(source:Path,receipt:Path,outdir:Path):
    require(not outdir.exists() or not any(outdir.iterdir()),
            "No silent overwrite of signed completed photo data")
    q,source_sha=check_read(source,receipt)
    results=[]
    for width in WIDTHS:
        matches=matched_community(q,width=width,min_obs=N_MIN,min_taxa=TAX_MIN)
        for zone in ("north_extratropical","south_extratropical"):
            sub=matches[matches.zone.eq(zone)]
            point=cell_cluster_interval(sub,draws=999,seed=20261017)
            results.append({
                "grid_width_degrees":width,
                "hemisphere":zone,
                "min_independent_source_observations_in_each_season":N_MIN,
                "min_distinct_named_photo_taxa_each_season":TAX_MIN,
                "matched_cell_season_year_groups":int(len(sub)),
                "different_original_geocells":int(sub[
                    ["latitude_cell","longitude_cell"]].drop_duplicates().shape[0]),
                "classifiable_photographs_in_matched_communities":int(
                    sub.summer_photos.sum()+sub.winter_photos.sum()) if len(sub) else 0,
                "summer_minus_winter_white_photo_pp":(
                    None if point is None else point["equal_cell_year_white_difference_pp"]),
                "bootstrap_geocell95_percentile":(
                    None if point is None else point["bootstrap_95pp"]),
                "inference_status":(
                    "NO_SOURCE_MATCHES" if point is None else point["uncertainty_status"]),
            })
    outdir.mkdir(parents=True,exist_ok=True)
    outfile=outdir/"mois_real295017_original_photo_community_local_year_white_scale.csv"
    pd.DataFrame(results).to_csv(outfile,index=False,lineterminator="\n")
    j={
        "schema":"mois_real295017_photographed_community_original_photo_year_spatial_scale_v1",
        "source_sha256":source_sha,
        "real_original_photo_attempts":EXPECTED_ROWS,
        "real_classified_photo_ROIs":EXPECTED_HUES,
        "source_or_independently_verified_public_dates":EXPECTED_DATE_ROWS,
        "real_classifiable_original_photo_ROIs_with_observed_dates":EXPECTED_DATED_HUES,
        "source_cell_and_observation_year_5observation_2taxon_results":results,
        "unadjusted_photo_hue_not_botanically_verified_open_flower":True,
        "original_technical_photo_hue_missing_not_coded_nonwhite":True,
        "source_community_different_taxa_both_seasons_allowed":True,
        "source_same_taxon_or_plant_pigment_plasticity_not_estimated":True,
        "0p25_geocell_not_biological_population":True,
        "source_observation_ID_not_independent_photographer_or_individual":True,
        "geocell_bootstrap_percentile95_exploratory_not_maxT":True,
        "no_actual_observation_year_month_climate_or_fitness_measured":True,
        "FCP_W_D_50km_intraspecific_estimator_not_imported":True,
        "community_source_season_CSV_sha256":
            hashlib.sha256(outfile.read_bytes()).hexdigest(),
    }
    (outdir/"mois_real295017_photo_community_season_grain_decision.json").write_text(
        json.dumps(j,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    print("MOIS_ACTUAL_295017_PHOTO_LOCALITY_YEAR_SAMPLING",
          json.dumps([{"zone":r["hemisphere"],"width":r["grid_width_degrees"],
                       "cell_years":r["matched_cell_season_year_groups"],
                       "local_summer_white_pp":r["summer_minus_winter_white_photo_pp"]}
                     for r in results if r["grid_width_degrees"]==.25],
                     sort_keys=True),flush=True)
    return j


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source",required=True,type=Path)
    p.add_argument("--source-receipt",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    a=p.parse_args()
    run(a.source,a.source_receipt,a.output)
