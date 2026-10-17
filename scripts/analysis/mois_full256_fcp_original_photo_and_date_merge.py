#!/usr/bin/env python3
"""MOIS source-exact finished FCP wave4 x 302242 date-qualified photograph union.

Only the *actually finished* 49,800 FCP source photograph pixels from all 256
original ROI worker partitions are admitted, and only 7,197 novel original
photo IDs extend the already completed MOIS 302,242. Any previous photo's
taxon, hue, observation or qualified date is held immutable. Automatic photo
hue is not botanically independently confirmed open petals.
"""
from __future__ import annotations
import argparse
import hashlib
import io
import json
import zipfile
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd

from scripts.analysis.mois_import_full_FCP_wave4_completed_original_photos import (
    read_fcp_wave4, BASE_COLS,
)
from scripts.analysis.mois_seasonal_267k_climate_normals import (
    HUES, labelled_calendar, require,
)
from scripts.analysis.mois_original_photo_community_season_geoscale import (
    matched_community, cell_cluster_interval,
)

BASE_SOURCE = "mois_original302242_exact_observed_on_full_uuid_numeric_verified.csv.gz"
BASE_RECEIPT = "mois_original302242_146052_photo_hue_full_observed_on_provenance_receipt.json"
BASE_SHA256 = "9c9692238f58cc61561270f7c07b1d782ce711adeade9c58550866f668fb4690"
FULL_FCP_SHA256 = "fc642d42b861bc79305dbd7e7e57d6b14561028e5552bc39dece93e59442f766"
BASE_N = 302242
BASE_ROI = 146061
BASE_DATED_N = 302223
BASE_DATED_ROI = 146052
FCP_N = 49800
NEW_N = 7197
NEW_ROI = 3601
EXPECTED_FINAL_N = 309439
EXPECTED_FINAL_ROI = 149662
EXPECTED_FINAL_DATED = 309420
EXPECTED_FINAL_DATED_ROI = 149653
OUTPUT_PHOTOS = "mois_original309439_fullFCP256_actual_photo_hues_and_source_dates.csv.gz"
OUTPUT_SEASON = "mois_309439_source_real_date_photo_white_season_numerators.csv"
OUTPUT_CELLS = "mois_309439_same_cell_year_summer_winter_photo_community.csv"
OUTPUT_RECEIPT = "mois_FCP_full256_309439_original_pixels_149653_dated_hues_source_receipt.json"


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def old_photographs(path: Path) -> tuple[pd.DataFrame, dict]:
    with zipfile.ZipFile(path) as zipped:
        blob = zipped.read(BASE_SOURCE)
        receipt = json.loads(zipped.read(BASE_RECEIPT))
    require(sha(blob) == BASE_SHA256 and
            receipt["final_source_data_csv_SHA256"] == BASE_SHA256 and
            receipt["original_302242_actual_measured_photo_pixels"] == BASE_N and
            receipt["original_146061_automatically_classified_four_hue_ROIs"] == BASE_ROI and
            receipt["final_source_observation_dated_photo_attempts"] == BASE_DATED_N and
            receipt["final_dated_and_original_four_hue_classifiable_ROIs"] == BASE_DATED_ROI,
            "Published original 302242 photo and historical date source checksum differs")
    photo = pd.read_csv(io.BytesIO(blob), compression="gzip",
                        dtype=str, keep_default_na=False)
    require(len(photo) == BASE_N and photo.photo_id.is_unique and
            tuple(photo.columns) == (*BASE_COLS, "qualified_date",
                                     "original_source_date_provenance") and
            photo.morph.isin(HUES).sum() == BASE_ROI and
            photo.qualified_date.str.len().eq(10).sum() == BASE_DATED_N and
            (photo.morph.isin(HUES)&photo.qualified_date.str.len().eq(10)).sum() == BASE_DATED_ROI,
            "Published originally classified photo identities or date metadata changed")
    return photo, receipt


def _date_real(day: str) -> bool:
    try:
        return date.fromisoformat(day).isoformat() == day and date.fromisoformat(day) <= date(2026,10,17)
    except (TypeError, ValueError):
        return False


def append_finished_pixels(old: pd.DataFrame, originals: pd.DataFrame,
                           *, strict:bool=True) -> tuple[pd.DataFrame, dict]:
    require(old.photo_id.is_unique and originals.photo_id.is_unique and
            old.photo_id.str.fullmatch(r"\d+").all() and
            originals.photo_id.str.fullmatch(r"\d+").all(),
            "Actual original photograph IDs duplicated or not numeric")
    seen = originals.photo_id.isin(old.photo_id)
    both = originals[seen].merge(old,on="photo_id",validate="one_to_one",
                                  suffixes=("_fcp","_mois"))
    addition = originals[~seen].copy()
    require(not addition.empty, "Full FCP 49,800 original photo pixels do not add any photos")
    for col in ("inat_taxon_id","morph","observation_id"):
        require(both[col+"_fcp"].eq(both[col+"_mois"]).all(),
                f"Original sampled photo duplicated but {col} disagrees")
    for col in ("latitude","longitude"):
        l = pd.to_numeric(both[col+"_fcp"], errors="coerce")
        r = pd.to_numeric(both[col+"_mois"], errors="coerce")
        require(l.notna().all() and r.notna().all() and
                np.isclose(l,r,atol=1e-8,rtol=0).all(),
                f"Original photo {col} changed between independently measured sources")
    # Older historical MOIS rows have blank raw date but a *verified* S3/API
    # source date. Exactly compare to qualified_date, never overwrite either.
    require(both.observed_on_fcp.eq(both.qualified_date).all(),
            "FCP original observed_on differs from the earlier verified original photo observation date")
    raw_day_conflicts=both.observed_on_mois.ne(both.observed_on_fcp)
    require((both.loc[raw_day_conflicts,"observed_on_mois"].eq("") &
             both.loc[raw_day_conflicts,"qualified_date"].eq(
                 both.loc[raw_day_conflicts,"observed_on_fcp"])).all(),
            "Source public original date cannot repair incompatible older raw dates")
    require(addition.inat_taxon_id.str.isdecimal().all() and
            addition.observation_id.str.isdecimal().all() and
            addition.observed_on.map(_date_real).all() and
            addition.morph.isin((*HUES, "mixed_uncertain")).all() and
            addition.photo_license.ne("").all(),
            "A new FCP image is unmeasured, has no original source date, or wrong palette")
    require(addition.morph.isin(HUES).eq(
            addition.measurement_status.eq("classified_four_state_morph")).all() and
            addition.global_classifiable.str.lower().eq(
                addition.morph.isin(HUES).astype(str).str.lower()).all(),
            "Completed FCP pixel classification status disagrees with ROI hue")
    if strict:
        require(len(old)==BASE_N and old.morph.isin(HUES).sum()==BASE_ROI and
                len(originals)==FCP_N and len(addition)==NEW_N and
                addition.morph.isin(HUES).sum()==NEW_ROI and len(both)==FCP_N-NEW_N,
                "Full original FCP 256-photo source or original overlapping IDs differ from signed source")
        require(int(raw_day_conflicts.sum())==960,
                "Prior historic missing raw date but qualified S3/API source-date overlap changed")
    newcols = pd.DataFrame({col:addition[col] if col in addition else ""
                            for col in BASE_COLS},index=addition.index)
    newcols["source_original_measured"] = "fcp_wave4_full256_original_ROI_v4"
    newcols["old_mois_or_new_fcp"] = "fcp_additional_ROI"
    newcols["qualified_date"] = addition.observed_on
    newcols["original_source_date_provenance"] = "fcp_frozen_full256_pixel_terminal_observed_on"
    combined = pd.concat([old,newcols[list(old.columns)]],ignore_index=True)
    require(combined.photo_id.is_unique and
            combined.iloc[:len(old)].reset_index(drop=True).equals(old.reset_index(drop=True)),
            "Previously complete original photographic source photo/date rows were altered")
    new_class = int(addition.morph.isin(HUES).sum())
    result = {
        "previous_original_measured_photo_IDs":len(old),
        "previous_original_classifiable_photo_ROIs":int(old.morph.isin(HUES).sum()),
        "previous_dated_and_classifiable_ROIs":int(
            (old.morph.isin(HUES)&old.qualified_date.str.len().eq(10)).sum()),
        "original_FCP_full256_photo_rows":len(originals),
        "original_FCP_pixels_overlapping_previous_source_photo_IDs":len(both),
        "source_original_observed_on_blank_but_independently_qualified_days_agree":
            int(raw_day_conflicts.sum()),
        "new_original_genuinely_measured_photo_pixels":len(addition),
        "new_original_four_hue_ROI_photo_labels":new_class,
        "new_original_photo_IDs_with_actual_observed_on":len(addition),
        "new_date_and_four_hue_original_photo_ROIs":new_class,
        "now_original_measured_photo_IDs":len(combined),
        "now_original_classifiable_four_hue_ROIs":int(combined.morph.isin(HUES).sum()),
        "now_original_observation_day_present_photo_IDs":int(
            combined.qualified_date.str.len().eq(10).sum()),
        "now_original_observation_day_and_four_hue_ROIs":int(
            (combined.morph.isin(HUES)&combined.qualified_date.str.len().eq(10)).sum()),
        "source_photo_taxon_hue_date_and_location_identities_unmodified":True,
        "historical_source_S3_API_date_priority_not_overridden":True,
        "technically_unclassifiable_photo_not_recoded_to_nonwhite":True,
        "FCP_within_species_genetic_D_W_50km_not_imported":True
    }
    require(result["now_original_measured_photo_IDs"]==len(old)+len(addition) and
            result["now_original_observation_day_and_four_hue_ROIs"]==
               result["previous_dated_and_classifiable_ROIs"]+new_class,
            "Newly measured photo and dated-hue source counts are not additive")
    if strict:
        require(result["now_original_measured_photo_IDs"]==EXPECTED_FINAL_N and
                result["now_original_classifiable_four_hue_ROIs"]==EXPECTED_FINAL_ROI and
                result["now_original_observation_day_present_photo_IDs"]==EXPECTED_FINAL_DATED and
                result["now_original_observation_day_and_four_hue_ROIs"]==EXPECTED_FINAL_DATED_ROI,
                "Actual completed full FCP photo ROI/date endpoint differs from independently audited counts")
    return combined,result


def source_photo_seasons(p:pd.DataFrame) -> tuple[pd.DataFrame,pd.DataFrame,dict]:
    d = p[["photo_id","observation_id","inat_taxon_id","morph",
           "latitude","longitude","old_mois_or_new_fcp","qualified_date"]].copy()
    d = d.rename(columns={"qualified_date":"observed_on"})
    for c in ("latitude","longitude"):
        d[c]=pd.to_numeric(d[c],errors="coerce")
    cal=labelled_calendar(d)
    require(int((cal.date.notna() & cal.classified).sum())==
            EXPECTED_FINAL_DATED_ROI, "New genuine source photo seasons no longer reconcile")
    rows=[]
    for zone in ("north_extratropical","south_extratropical"):
        for season in ("winter","spring","summer","autumn"):
            sub=cal[cal.zone.eq(zone)&cal.date.notna()&cal.season.eq(season)]
            valid=sub[sub.classified]
            rows.append({
                "source_hemisphere":zone,"meteorological_season":season,
                "original_photo_attempts":int(len(sub)),
                "classified_original_photo_ROIs":int(len(valid)),
                "original_white_label_ROIs":int(valid.white.sum()),
                "white_fraction_among_actual_classifiable_ROIs":
                    round(float(valid.white.mean()),8) if len(valid) else None})
    seasons=pd.DataFrame(rows)
    spatial=[]
    for width in (0.25,0.5,1.0,2.0):
        groups=matched_community(cal,width=width,min_obs=5,min_taxa=2)
        for zone in ("north_extratropical","south_extratropical"):
            sub=groups.loc[groups.zone.eq(zone)]
            stat=cell_cluster_interval(sub,draws=999,seed=20261017+int(width*100))
            spatial.append({
                "cell_degrees":width,"hemisphere":zone,
                "same_cell_same_season_year_photo_community_groups":int(len(sub)),
                "source_photo_ROIs_in_matched_communities":
                    int(sub.summer_photos.sum()+sub.winter_photos.sum()) if len(sub) else 0,
                "equal_local_cell_year_summer_minus_winter_white_pp":
                    None if stat is None else stat["equal_cell_year_white_difference_pp"],
                "exploratory_999_geocell_bootstrap_percentile95_pp":
                    None if stat is None else stat["bootstrap_95pp"]})
    result={
        "original_date_and_hue_nominal_source_taxa":int(
            cal.loc[cal.date.notna()&cal.classified,"inat_taxon_id"].nunique()),
        "actual_original_source_classifiability_denominator":int(
            cal.date.notna().sum()),
        "new_photographic_source_is_not_biologically_independent_of_FCP":True,
        "local_taxon_photo_presence_not_biological_bloom_presence":True,
        "open_flower_organ_not_independently_validated":True,
        "no_actual_source_observation_year_climate_or_pigment_fitness_fitted":True
    }
    return seasons,pd.DataFrame(spatial),result


def run(photo_archive:Path,fcp_directory:Path,output:Path):
    require(not output.exists() or not any(output.iterdir()),
            "Never silently replace a previously completed source photo export")
    base,old_receipt=old_photographs(photo_archive)
    candidate_paths=list(fcp_directory.rglob("wave4_original_method_measured_photos.csv"))
    require(len(candidate_paths)==1 and sha(candidate_paths[0].read_bytes())==FULL_FCP_SHA256,
            "Actual 49800 finished original FCP photo-pixel terminal SHA drifted")
    fcp,receipt=read_fcp_wave4(fcp_directory)
    source=pd.DataFrame(fcp,dtype=str).fillna("")
    combined,stats=append_finished_pixels(base,source)
    seasons,sites,ecology=source_photo_seasons(combined)
    output.mkdir(parents=True,exist_ok=True)
    target=output/OUTPUT_PHOTOS
    combined.to_csv(target,index=False,lineterminator="\n",
                    compression={"method":"gzip","mtime":0})
    season_file=output/OUTPUT_SEASON
    seasons.to_csv(season_file,index=False,lineterminator="\n")
    site_file=output/OUTPUT_CELLS
    sites.to_csv(site_file,index=False,lineterminator="\n")
    final={
        "schema":"mois_full_FCP256_309439_original_pixels_149653_dated_hues_v1",
        "original_MOIS_baseline_302242_photo_SHA256":BASE_SHA256,
        "genuine_FCP_full_49800_completed_original_photo_CSV_SHA256":FULL_FCP_SHA256,
        "FCP_full_original_photo_study_metadata_SHA256":receipt["source_metadata_sha256"],
        "FCP_full_original_256_worker_partitions":receipt["terminal_partitions"],
        "FCP_full_original_49800_photo_ids":receipt["wave4_new_photo_ids"],
        "FCP_full_49800_original_four_hue_classifications":receipt["wave4_classifiable_photo_rows"],
        "source_MOIS_302242_archive_original_run":
            "https://github.com/zuizui0223/mois/actions/runs/39020466312",
        "source_FCP_complete_256_ROI_original_run":
            "https://github.com/zuizui0223/fcp/actions/runs/38992543115",
        "source_completed_original_FCP_final_artifact_id":12097805951,
        "archive_exact_original_new_photo_source_SHA256":sha(target.read_bytes()),
        "hemisphere_photo_season_figure_source_SHA256":sha(season_file.read_bytes()),
        "same_site_year_seasonal_photo_community_source_SHA256":sha(site_file.read_bytes()),
        "summer_winter_photographic_community_geocell_results":sites.to_dict(orient="records"),
        "photographic_dates_are_original_public_observed_on_not_independently_image_EXIF":True,
        **stats,**ecology
    }
    (output/OUTPUT_RECEIPT).write_text(json.dumps(
        final,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    print("MOIS_FULL_FCP_ORIGINAL_256_FINISHED_PIXELS_SOURCE_EXPANSION",
          json.dumps({k:final[k] for k in (
            "new_original_genuinely_measured_photo_pixels",
            "new_original_four_hue_ROI_photo_labels",
            "now_original_measured_photo_IDs",
            "now_original_classifiable_four_hue_ROIs",
            "now_original_observation_day_and_four_hue_ROIs",
            "original_date_and_hue_nominal_source_taxa")},sort_keys=True),flush=True)
    return final


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--302242-original-photo-artifact",type=Path,required=True)
    ap.add_argument("--full256-original-FCP-artifact-dir",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    a=ap.parse_args()
    run(a.__dict__["302242_original_photo_artifact"],
        a.full256_original_FCP_artifact_dir,a.output)
