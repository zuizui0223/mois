#!/usr/bin/env python3
"""Final MOIS original photo/pixel/date union from THREE independently signed sources.

(1) 302242 authenticated actual original ROI-v4 photo pixels, including new
    finished FCP218 parts, with 143770 previously dated four-hue ROI photos.
(2) 1732 EXACT original photo+numeric observation ID current iNat dated repairs
    (706 previously already-classified ROI outcomes).
(3) 4096 original first-wave UUID-only observation photo IDs matched to the
    full 3430502 FCP authoritative Sept-2026 photo metadata with 1576 existing
    original automated four-hue ROI labels.

Never create new phenotypes from photo opportunities, never overwrite old
S3/public exact observed_on dates, never import FCP within-species W/D/50km
outcomes, and never identify photographed open petals from these source joins.
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

from scripts.analysis.mois_seasonal_267k_climate_normals import (
    HUES, require, labelled_calendar,
)
from scripts.analysis.mois_original_photo_community_season_geoscale import (
    matched_community, cell_cluster_interval,
)
from scripts.analysis.mois_final_original_photo_s3_api_seasonal_composition import (
    seasonal_taxon_decomposition,
)
from scripts.analysis.mois_seasonal_spatial_identifiability_20261017 import (
    eligible_site_photo_support,
)

PIXEL_N=302242
CLASSIFIABLE_N=146061
BASE_DATED_ORIG=296395
BASE_DATED_HUES=143770
OLD_N=295017
API_NEW_DATES=1732
API_NEW_HUES=706
UUID_NEW_DATES=4096
UUID_NEW_HUES=1576
FINAL_DATED_ORIG=302223
FINAL_DATED_HUES=146052
BASE_SHA="b83871c6cd69cb624693364742c2a1d61003253d084af0cb6d80ab9f898c3959"
API_SHA="65c072a97f7a59b1933dbc0541041a29069d468bf6fc889933d5a922c743a11a"
UUID_SHA="9ca764f7d3596fd3245a9ef079a74a2413ba0223622fd85227bcceed46ca9c81"
PHOTO_BASE="mois_next_finished_fcp_original_photo_date_hues_plus295017.csv.gz"
PHOTO_RECEIPT="mois_next_finished_fcp_original_photo_new_pixel_date_hue_receipt.json"
API_BASE="mois_real295017_photo_source_plus_exact_remaining_numeric_API_dates.csv.gz"
API_RECEIPT="mois_original295017_remaining1751_numeric_photo_date_recovery_QC_receipt.json"
UUID_BASE="mois_first4096_uuid_only_original_photo_source_date_exact3430502.csv"
UUID_RECEIPT="mois_first4096_uuid_observation_photo_day_original_source_QC_receipt.json"
OUTPUT_SOURCE="mois_real302242_all_verified_observed_on_four_hue_original_photos.csv.gz"
OUTPUT_SUMMARY="mois_real302242_146052_source_photo_season_and_classifiability.csv"
OUTPUT_LOCAL="mois_real302242_146052_original_photographic_community_geographic_scale.csv"
OUTPUT_RECEIPT="mois_real302242_146052_verified_date_hues_photo_source_season_decision.json"
SCALES=(0.25,0.5,1.,2.)
ZONES=("north_extratropical","south_extratropical")


def hashbytes(raw:bytes)->str:
    return hashlib.sha256(raw).hexdigest()


def signed_csv(archive:Path,member:str,receipt_member:str,
               expected_sha:str,source_key:str)->tuple[pd.DataFrame,dict]:
    with zipfile.ZipFile(archive) as z:
        raw=z.read(member)
        receipt=json.loads(z.read(receipt_member))
    require(hashbytes(raw)==expected_sha and
            receipt[source_key]==expected_sha,
            "Complete signed original photo date source SHA changed "+member)
    if member.endswith(".gz"):
        d=pd.read_csv(io.BytesIO(raw),compression="gzip",
                      dtype=str,keep_default_na=False)
    else:
        d=pd.read_csv(io.BytesIO(raw),dtype=str,keep_default_na=False)
    return d,receipt


def apply_all_verified_days(original:pd.DataFrame,numeric:pd.DataFrame,
                            firstuuid:pd.DataFrame)->tuple[pd.DataFrame,dict]:
    require(len(original)==PIXEL_N and original.photo_id.is_unique
            and original.morph.isin(HUES).sum()==CLASSIFIABLE_N and
            original.qualified_date.str.len().eq(10).sum()==BASE_DATED_ORIG and
            (original.morph.isin(HUES)&original.qualified_date.str.len().eq(10)).sum()==BASE_DATED_HUES,
            "Expanded source original full finished-pixel/photo-ID/metadata population drift")
    require(len(numeric)==OLD_N and numeric.photo_id.is_unique and
            original.photo_id.iloc[:OLD_N].tolist()==numeric.photo_id.tolist(),
            "API source original 295017 photo identities or order diverge from finished FCP218 source")
    original_columns=[c for c in original.columns if c not in (
        "qualified_date","original_source_date_provenance")]
    require(numeric[original_columns].equals(original.iloc[:OLD_N][original_columns].reset_index(drop=True)),
            "Original photo/taxon/location/raw-date/automatic-hue label changed in API date repair")

    # Only new, source-photo and source-numeric-ID-verified dates are eligible.
    new_numeric=numeric[numeric.original_source_date_provenance.eq(
        "current_iNat_exact_numeric_observation_original_photo_date_verified")]
    require(len(new_numeric)==API_NEW_DATES and
            new_numeric.photo_id.is_unique and
            new_numeric.morph.isin(HUES).sum()==API_NEW_HUES and
            new_numeric.qualified_date.str.fullmatch(r"\d{4}-\d{2}-\d{2}").all(),
            "Additional independent 1732 exact numeric photo/observation date proofs changed")
    require(len(firstuuid)==UUID_NEW_DATES and
            firstuuid.original_photo_id.is_unique and
            firstuuid.date_identity_status.eq(
              "photo_ID_original_taxon_and_snapshot_date_verified").all()
            and firstuuid.original_source_photo_hue.isin(HUES).sum()==UUID_NEW_HUES
            and firstuuid.verified_original_full_source_observed_on.str.fullmatch(
                r"\d{4}-\d{2}-\d{2}").all(),
            "Original full3430502 SHA-verified 4096 UUID-photo date join changed")
    require(set(new_numeric.photo_id).isdisjoint(set(firstuuid.original_photo_id)),
            "Public numeric original photo date and original UUID date double-count one photo")
    out=original.copy()
    index=pd.Series(np.arange(len(out)),index=out.photo_id)
    original_days=out.qualified_date.copy()
    for source in new_numeric.itertuples(index=False):
        pid=source.photo_id
        require(pid in index.index,"Public original numeric photo source missing in new measured photo census")
        j=int(index.loc[pid])
        require(out.at[j,"qualified_date"]=="" and
                out.at[j,"observation_id"]==source.observation_id and
                out.at[j,"inat_taxon_id"]==source.inat_taxon_id and
                out.at[j,"morph"]==source.morph,
                "Original numeric observed day would overwrite or change an original source photo")
        date.fromisoformat(source.qualified_date)
        out.at[j,"qualified_date"]=source.qualified_date
        out.at[j,"original_source_date_provenance"]=(
            "current_iNat_exact_numeric_observation_original_photo_date_verified")
    for source in firstuuid.itertuples(index=False):
        pid=source.original_photo_id
        require(pid in index.index,"UUID-source original photograph missing from real measured 302242 source")
        j=int(index.loc[pid])
        day=source.verified_original_full_source_observed_on
        date.fromisoformat(day)
        require(out.at[j,"qualified_date"]=="" and
                out.at[j,"inat_taxon_id"]==source.original_first4096_taxon_id and
                out.at[j,"morph"]==source.original_source_photo_hue and
                out.at[j,"observation_id"]==source.original_first4096_observation_source_UUID
                and out.at[j,"source_original_measured"]=="mois_wave1_new_original_ROI_v4",
                "Original UUID-only photograph/taxon/hue or old verified date conflicts")
        out.at[j,"qualified_date"]=day
        out.at[j,"original_source_date_provenance"]=(
            "FCP_full3430502_20260927_photo_ID_taxon_observation_day_verified")
    updated=out.qualified_date
    require(out.photo_id.tolist()==original.photo_id.tolist()
            and out[original_columns].equals(original[original_columns])
            and out.loc[original_days.ne(""),"qualified_date"].equals(
                original.loc[original_days.ne(""),"qualified_date"]),
            "Full original measured photos, colours or previously verified dates were changed")
    require(updated.str.len().eq(10).sum()==FINAL_DATED_ORIG and
            (out.morph.isin(HUES)&updated.str.len().eq(10)).sum()==FINAL_DATED_HUES and
            out.photo_id.is_unique and
            out.morph.isin(HUES).sum()==CLASSIFIABLE_N,
            "Full genuine source calendar/hue denominator does not add exactly")
    stats={
        "original_measured_photo_ids":len(out),
        "original_four_hue_roi_photo_ids":int(out.morph.isin(HUES).sum()),
        "initial_dated_actual_photo_attempts":int(original_days.str.len().eq(10).sum()),
        "initial_date_plus_four_hue_original_photos":int(
            (original.morph.isin(HUES)&original_days.str.len().eq(10)).sum()),
        "new_exact_numeric_observation_photo_ID_dates":len(new_numeric),
        "new_date_and_classified_four_hue_from_numeric_photo_source":
            int(new_numeric.morph.isin(HUES).sum()),
        "new_original_UUID_only_observation_photo_dates":len(firstuuid),
        "new_date_and_classified_four_hue_from_UUID_photo_source":
            int(firstuuid.original_source_photo_hue.isin(HUES).sum()),
        "final_real_photo_IDs_with_any_verified_observed_on":int(updated.str.len().eq(10).sum()),
        "final_real_photo_IDs_with_four_hue_and_verified_observed_on":int(
            (out.morph.isin(HUES)&updated.str.len().eq(10)).sum()),
        "missing_original_photo_observation_dates_after_recovery":
            len(out)-int(updated.str.len().eq(10).sum()),
        "automatic_four_hue_classifications_still_without_verified_observation_day":
            int((out.morph.isin(HUES)&updated.eq("")).sum()),
        "old_historical_99817_S3_API_source_dates_never_overwritten":True,
        "original_measured_photo_ROI_pixel_and_taxonomy_identifiers_not_changed":True,
        "observation_metadata_does_not_validate_independent_open_petal_stage":True,
        "FCP_within_species_D_W_50km_results_not_imported":True,
    }
    return out,stats


def source_analysis(enriched:pd.DataFrame):
    photo=enriched.copy()
    photo["observed_on"]=photo["qualified_date"]
    for col in ("latitude","longitude"):
        photo[col]=pd.to_numeric(photo[col],errors="coerce")
    p=labelled_calendar(photo)
    require(int((p.classified&p.date.notna()).sum())==FINAL_DATED_HUES,
            "Photographic season source does not match verified date+original hue denominator")
    region=[]
    for zone in ZONES:
        q=p[(p.zone==zone)&p.date.notna()&p.season.isin(("summer","winter"))]
        for season in ("summer","winter"):
            s=q[q.season.eq(season)]
            colored=s[s.classified]
            region.append({
                "hemisphere":zone,"season":season,
                "original_source_attempts_with_observation_day":len(s),
                "real_four_hue_ROI_with_observation_day":len(colored),
                "white_labeled_real_original_photos":int(colored.white.sum()),
                "conditional_white_photo_fraction":
                    round(float(colored.white.mean()),8) if len(colored) else None,
                "four_hue_classifiability_among_dated_original_photo_attempts":
                    round(float(s.classified.mean()),8) if len(s) else None,
            })
    desc={zone:seasonal_taxon_decomposition(p,zone) for zone in ZONES}
    geotab=[]
    for width in SCALES:
        matched=matched_community(p,width=width,min_obs=5,min_taxa=2)
        for zone in ZONES:
            sub=matched[matched.zone.eq(zone)] if not matched.empty else matched
            v=cell_cluster_interval(sub,draws=999,seed=20261017+int(width*100))
            geotab.append({
                "source_cell_degrees":width,"hemisphere":zone,
                "same_cell_season_year_original_photographed_communities":len(sub),
                "photographs_qualified_both_summer_winter":int(
                    sub.summer_photos.sum()+sub.winter_photos.sum()) if len(sub) else 0,
                "equal_cell_year_summer_minus_winter_photo_white_pp":
                    v["equal_cell_year_white_difference_pp"] if v else None,
                "geocell_cluster_bootstrap999_95_low_pp":
                    None if not v or v["bootstrap_95pp"] is None else v["bootstrap_95pp"][0],
                "geocell_cluster_bootstrap999_95_high_pp":
                    None if not v or v["bootstrap_95pp"] is None else v["bootstrap_95pp"][1],
                "uncertainty_label":v["uncertainty_status"] if v else "NO_PAIRED_GEOCELLS",
            })
    return region,desc,geotab


def run(pixel_zip:Path,api_zip:Path,uuid_zip:Path,output:Path):
    require(not output.exists() or not any(output.iterdir()),
            "Never rewrite an existing real original photo/date source")
    photos,photo_r=signed_csv(pixel_zip,PHOTO_BASE,PHOTO_RECEIPT,
                              BASE_SHA,"expanded_real_date_hue_photo_source_SHA256")
    api,api_r=signed_csv(api_zip,API_BASE,API_RECEIPT,
                        API_SHA,"full_295017_source_plus_verified_day_csv_SHA256")
    uuid,uuid_r=signed_csv(uuid_zip,UUID_BASE,UUID_RECEIPT,
                          UUID_SHA,"source_all_4096_original_photo_date_result_SHA256")
    require(photo_r["expanded_distinct_measured_source_photo_IDs"]==PIXEL_N and
            photo_r["expanded_dated_original_four_hue_ROIs"]==BASE_DATED_HUES and
            api_r["newly_verified_source_photo_numeric_original_obs_dates"]==API_NEW_DATES and
            uuid_r["original_photo_IDs_with_authenticated_full_source_calendar_day"]==UUID_NEW_DATES,
            "Original 3 independent source-receipt provenance counts changed")
    combined,counts=apply_all_verified_days(photos,api,uuid)
    output.mkdir(parents=True,exist_ok=True)
    full=output/OUTPUT_SOURCE
    combined.to_csv(full,index=False,lineterminator="\n",
                    compression={"method":"gzip","mtime":0})
    region,turnover,geo=source_analysis(combined)
    seasons=output/OUTPUT_SUMMARY
    locality=output/OUTPUT_LOCAL
    pd.DataFrame(region).to_csv(seasons,index=False,lineterminator="\n")
    pd.DataFrame(geo).to_csv(locality,index=False,lineterminator="\n")
    result={
        "schema":"mois_final_302242_actual_measured_original_photo_146052_dated_hues_v1",
        "original_FCP218_photo_source_SHA256":BASE_SHA,
        "original_1751_numeric_photo_date_provenance_SHA256":API_SHA,
        "original_4096_UUID_photo_date_source_SHA256":UUID_SHA,
        "full_original_photo_date_plus_automatic_four_hue_source_SHA256":hashbytes(full.read_bytes()),
        "actual_original_source_regional_season_white_CSV_SHA256":hashbytes(seasons.read_bytes()),
        "actual_original_source_geocell_community_year_CSV_SHA256":hashbytes(locality.read_bytes()),
        "source_hemispheric_season_counts":region,
        "original_source_taxon_season_exact_photo_frequency_decomposition":turnover,
        "original_photo_source_same_locality_year_community_white_scale":geo,
        "distinct_named_taxa_with_classified_photo_and_observed_day":
            int(combined.loc[combined.morph.isin(HUES)&
                combined.qualified_date.str.len().eq(10),"inat_taxon_id"].nunique()),
        "no_original_photo_RGB_or_four_hue_label_changed":True,
        "no_unclassified_photo_assigned_a_nonwhite_flowering_state":True,
        "no_FCP_D_W_50km_intraspecific_estimates_imported":True,
        "photo_specific_open_flowering_organ_unverified":True,
        "source_date_is_observation_metadata_not_camera_EXIF":True,
        "real_month_of_observation_TerraClimate_exposure_not_fitted":True,
        "geographically_coarse_source_cell_not_biological_population":True,
        **counts
    }
    (output/OUTPUT_RECEIPT).write_text(
        json.dumps(result,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    print("MOIS_FINAL_302242_MEASURED_ORIGINAL_PHOTOS_FULL_DATE",
         json.dumps({k:result[k] for k in (
            "final_real_photo_IDs_with_any_verified_observed_on",
            "final_real_photo_IDs_with_four_hue_and_verified_observed_on",
            "missing_original_photo_observation_dates_after_recovery",
            "automatic_four_hue_classifications_still_without_verified_observation_day",
            "distinct_named_taxa_with_classified_photo_and_observed_day"
         )},sort_keys=True),flush=True)
    return result


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--new302242-fcp218-source-archive",type=Path,required=True)
    p.add_argument("--numeric1751-observed-on-recovery-archive",type=Path,required=True)
    p.add_argument("--uuid4096-observed-on-recovery-archive",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    q=p.parse_args()
    run(q.new302242_fcp218_source_archive,
        q.numeric1751_observed_on_recovery_archive,
        q.uuid4096_observed_on_recovery_archive,q.output)
