#!/usr/bin/env python3
"""After all 16 real photo ROI-v4 workers finish, add fourth 4096 to MOIS.

Verifies every source CSV, all 4096 true colour/unknown pixel terminal records,
unchanged original YOLO/EfficientSAM method SHA, and 309439 source photo/date
GZIP before adding anything to the real original-photo dataset.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

from scripts.analysis.mois_original_fourth4096_after_309439_source_freeze import (
    FREEZE_OUT, SOURCE_OUT, SALT, MEASURED_N, MEASURED_HUES,
    MEASURED_DATED_HUES, MEASURED_SHA,
)
from scripts.analysis.mois_reassemble_original_new4096_measurements import (
    BASE_COLS, HUES, METHOD, load_real_wave,
)
from scripts.analysis.mois_seasonal_267k_climate_normals import (
    require, labelled_calendar,
)
from scripts.analysis.mois_original_photo_community_season_geoscale import (
    matched_community, cell_cluster_interval,
)

ORIGINAL_SOURCE_FILENAME="mois_original309439_fullFCP256_actual_photo_hues_and_source_dates.csv.gz"
ORIGINAL_SOURCE_COHORT_SHA="546266e2bbaa340a2523864c28c5a877842eb4b0904a7e0dff90870c6519a47e"
RESULT_SOURCE="mois_full313535_original_pixel_ROIs_and_observation_dates.csv.gz"
RESULT_RECEIPT="mois_fourth_4096_actual_original_ROI_pixels_and_source_dates_decision.json"
RESULT_SEASONS="mois_fourth4096_new_real_photo_pixel_hemisphere_seasons.csv"
RESULT_COMMUNITY="mois_fourth4096_new_real_photo_pixel_siteyear_community.csv"
NEW_PIXEL_N=4096


def digest(binary:bytes)->str:
    return hashlib.sha256(binary).hexdigest()


def signed_original(archive:Path):
    with zipfile.ZipFile(archive) as zipped:
        raw=zipped.read(ORIGINAL_SOURCE_FILENAME)
    require(digest(raw)==ORIGINAL_SOURCE_COHORT_SHA==MEASURED_SHA,
            "Already measured original309439 photo/date source GZIP changed")
    p=pd.read_csv(io.BytesIO(raw),compression="gzip",dtype=str,
                  keep_default_na=False,low_memory=False)
    require(len(p)==MEASURED_N and p.photo_id.is_unique and
            p.morph.isin(HUES).sum()==MEASURED_HUES and
            (p.morph.isin(HUES)&p.qualified_date.str.len().eq(10)).sum()==MEASURED_DATED_HUES,
            "Original already measured photo, hue and source date denominators changed")
    require(tuple(p.columns)==(*BASE_COLS,"qualified_date",
                               "original_source_date_provenance"),
            "Original full256 photo identity and extra dated-source columns changed")
    return p


def append(old:pd.DataFrame,photos:list[dict]) -> tuple[pd.DataFrame,dict]:
    require(len(old)==MEASURED_N and old.photo_id.is_unique,
            "Previous source photo population changed")
    require(len(photos)==NEW_PIXEL_N and
            len({str(x["photo_id"]) for x in photos})==NEW_PIXEL_N,
            "Fourth original photo image pixel classifier not fully completed")
    genuine=pd.DataFrame(photos).fillna("").astype(str)
    require(set(BASE_COLS).issubset(genuine.columns) and
            genuine.morph.isin((*HUES,"mixed_uncertain")).all(),
            "Original frozen model photo ROI results were not source qualified")
    require(set(genuine.photo_id).isdisjoint(set(old.photo_id)),
            "Fourth new real photo pixels overlap existing original309439")
    require(genuine.observed_on.str.fullmatch(r"\d{4}-\d{2}-\d{2}").all(),
            "Fourth original photograph's source calendar day absent")
    to_add=genuine[list(BASE_COLS)].copy()
    to_add["qualified_date"]=to_add.observed_on
    to_add["original_source_date_provenance"] = (
        "fcp_full3430502_original_photo_observed_on_with_MOIS_wave4_unchanged_ROI_v4")
    output=pd.concat([old,to_add[list(old.columns)]],ignore_index=True)
    require(output.photo_id.is_unique and len(output)==MEASURED_N+NEW_PIXEL_N and
            output.iloc[:len(old)].reset_index(drop=True).equals(old.reset_index(drop=True)),
            "Original 309439 genuinely measured pixels, colour, source date were changed")
    nclass=int(to_add.morph.isin(HUES).sum())
    ndates=int((to_add.morph.isin(HUES)&
                to_add.qualified_date.str.len().eq(10)).sum())
    require(ndates==nclass and nclass>0 and
            int(output.morph.isin(HUES).sum())==MEASURED_HUES+nclass and
            int((output.morph.isin(HUES)&output.qualified_date.str.len().eq(10)).sum())
            ==MEASURED_DATED_HUES+nclass,
            "New genuine photographed pixel ROIs/date cardinality invalid")
    counts={
       "before_real_measured_photo_ids":MEASURED_N,
       "before_original_photo_four_hue_ROIs":MEASURED_HUES,
       "before_original_date_plus_four_hue_ROIs":MEASURED_DATED_HUES,
       "new_photo_ids_actually_processed":len(to_add),
       "new_photo_four_hue_original_ROI_classified":nclass,
       "new_photo_unknown_pixel_OR_ROI_color":int(len(to_add)-nclass),
       "new_photo_with_original_source_observed_on":NEW_PIXEL_N,
       "new_photo_both_source_observation_day_and_four_hue":ndates,
       "expanded_genuine_original_photo_ids":len(output),
       "expanded_original_four_hue_ROI_photo_ids":int(output.morph.isin(HUES).sum()),
       "expanded_actual_date_and_four_hue_ROI_photo_ids":int(
            (output.morph.isin(HUES)&output.qualified_date.str.len().eq(10)).sum()),
       "source_observation_IDs_also_seen_in_old_photos":int(
            to_add.observation_id.isin(old.observation_id).sum()),
       "no_old_source_hue_taxon_photo_or_observation_date_changed":True,
       "botanically_confirmed_open_petal_per_photo":False,
       "FCP_intraspecific_genetic_or_D_W_50km_models_imported":False,
    }
    return output,counts


def calendar_summary(allphotos:pd.DataFrame):
    q=allphotos[[
        "photo_id","observation_id","inat_taxon_id","morph",
        "latitude","longitude","old_mois_or_new_fcp","qualified_date"]].copy()
    q.rename(columns={"qualified_date":"observed_on"},inplace=True)
    for col in ("latitude","longitude"):
        q[col]=pd.to_numeric(q[col],errors="coerce")
    panel=labelled_calendar(q)
    require(int((panel.classified&panel.date.notna()).sum())==
            int((allphotos.morph.isin(HUES)&
                allphotos.qualified_date.str.len().eq(10)).sum()),
            "Source photo/season-year parser cannot account for original photo ID")
    seasonal=[]
    for zone in ("north_extratropical","south_extratropical"):
        for season in ("winter","spring","summer","autumn"):
            p=panel[(panel.zone==zone)&(panel.season==season)&
                    panel.date.notna()]
            c=p[p.classified]
            seasonal.append({
                "hemisphere":zone,"local_meteorological_season":season,
                "dated_original_photo_attempts":int(len(p)),
                "classifiable_original_photo_hues":int(len(c)),
                "white_original_photo_ROIs":int(c.white.sum()),
                "white_photo_share_among_four_hues":(
                   float(c.white.mean()) if len(c)>0 else None),
            })
    groups=matched_community(panel,width=.25,min_obs=5,min_taxa=2)
    spatial=[]
    for zone in ("north_extratropical","south_extratropical"):
        eligible=groups[groups.zone==zone]
        estimate=cell_cluster_interval(
            eligible,draws=999,seed=20261017+len(eligible))
        spatial.append({
            "hemisphere":zone,
            "same_025degree_source_cell_and_meteorological_seasonyear_groups":len(eligible),
            "actual_four_hue_source_photos_in_comparable_cell_years":int(
                eligible.summer_photos.sum()+eligible.winter_photos.sum()),
            "descriptive_white_summer_minus_winter_pp":None if estimate is None else
                estimate["equal_cell_year_white_difference_pp"],
            "geocell_cluster_exploratory95_pp":None if estimate is None else
                estimate["bootstrap_95pp"],
            "n_distinct_geocells":None if estimate is None else estimate["n_unique_geocells"],
        })
    return pd.DataFrame(seasonal),pd.DataFrame(spatial)


def run(archive:Path,plan:Path,terminal:Path,output:Path):
    require(not output.exists() or not any(output.iterdir()),
            "Never rewrite real original flower pixel or date source")
    base=signed_original(archive)
    records,worker=load_real_wave(plan,terminal,"wave4")
    original,result=append(base,records)
    require(worker["model_attestation"]==METHOD and
            worker["n_real_original_photo_IDs"]==4096 and
            worker["no_imported_FCP_intraspecific_D_W_or_50km_measures"],
            "Original FCP ROI-v4 pixel source or model attestation changed")
    source=json.loads(next(plan.rglob(FREEZE_OUT)).read_text())
    require(source["schema"]==
            "mois_original_full3430502_next4096_after309439_wave4_preoutcome_v1"
            and source["first_wave_original_4096_csv_sha256"]==
                worker["source_original_selection_sha256"]
            and source["source_measured_original_309439_SHA256"]==
                MEASURED_SHA
            and source["first_wave_original_source_photo_ids"]==4096,
            "Pre-image original source ID and photo date selection changed")
    output.mkdir(parents=True,exist_ok=True)
    full=output/RESULT_SOURCE
    original.to_csv(full,index=False,lineterminator="\n",
                    compression={"method":"gzip","mtime":0})
    summary,local=calendar_summary(original)
    season=output/RESULT_SEASONS
    cells=output/RESULT_COMMUNITY
    summary.to_csv(season,index=False,lineterminator="\n",float_format="%.8f")
    local.to_csv(cells,index=False,lineterminator="\n",float_format="%.6f")
    report={
        "schema":"mois_actual_fourth4096_genuine_original_RGB_ROI_309439plus4096_v1",
        "source_full_3430502_original_ID_population":source["source_photo_ids_original"],
        "initial_measured_original_photo_GZIP_SHA256":MEASURED_SHA,
        "blind_fourth4096_original_source_selection_SHA256":
            worker["source_original_selection_sha256"],
        "source_original_v4_photo_measurement_method_attestation":METHOD,
        "original_source_observer_unique_only_within_fourth_wave":True,
        "source_original_photo_pixels_actually_processed":True,
        "source_photo_hue_mixed_uncertain_not_imputed":True,
        "never_identify_original_photo_bloom_anthesis_from_source_date":True,
        "real_observed_year_month_temperature_solar_water_not_analyzed_here":True,
        "current_fourth4096_source_SHA256":digest(full.read_bytes()),
        "current_fourth4096_season_table_SHA256":digest(season.read_bytes()),
        "current_fourth4096_geocell_table_SHA256":digest(cells.read_bytes()),
        "local_photo_community_minimum_original_observations_each_season":5,
        "local_photo_community_minimum_named_taxa_each_season":2,
        "geocell_bootstrap_draws":999,
        "all_real_original_photographed_photo_season_source":summary.to_dict("records"),
        "same_cell_same_seasonyear_photo_community_uncertainty":local.where(
            pd.notna(local),None).to_dict("records"),
        **result,
    }
    (output/RESULT_RECEIPT).write_text(json.dumps(
        report,sort_keys=True,ensure_ascii=False,indent=2)+"\n")
    print("MOIS_TRUE_NEW_FOURTH4096_PIXEL_ROI_MEASUREMENTS",
          json.dumps({k:report[k] for k in (
           "new_photo_ids_actually_processed",
           "new_photo_four_hue_original_ROI_classified",
           "expanded_genuine_original_photo_ids",
           "expanded_actual_date_and_four_hue_ROI_photo_ids")},sort_keys=True),flush=True)
    return report


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--real-original309439-source-zip",type=Path,required=True)
    p.add_argument("--blind-original-photo-plan",type=Path,required=True)
    p.add_argument("--real-original4096-photo-ROI-workers",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    run(a.real_original309439_source_zip,a.blind_original_photo_plan,
        a.real_original4096_photo_ROI_workers,a.output)
