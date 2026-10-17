#!/usr/bin/env python3
"""Reproduce the biological sampling-support gain after authentic date recovery.

A source-verified 9600 original-photo historical date cohort increases named
taxonomic breadth, not necessarily same-population within-season comparisons.
All denominators are actual preexisting photo IDs and original four-hue ROIs,
never unmeasured iNaturalist photo opportunities or plant anthesis assessments.
"""
from __future__ import annotations
import argparse
import hashlib
import io
import json
import zipfile
from pathlib import Path
import pandas as pd

from scripts.analysis.mois_seasonal_267k_climate_normals import (
    read_signed, labelled_calendar, HUES, MOIS_ROWS, MOIS_HUES,
    PHOTO_DATA, PHOTO_RECEIPT, require,
)
from scripts.analysis.mois_seasonal_spatial_identifiability_20261017 import eligible_site_photo_support

DATE_INPUT="mois_historical_9600_original_photo_dates_current_public.csv"
DATE_RECEIPT="mois_historical_9600_original_photo_dates_and_photo_seasons_receipt.json"
PHOTO_SHA="fbb54b58dc91e40e2e166081212588acc18adb19af4bdbcd1c7fe1778bee28d9"
DATE_SHA="dac4f2d9776a665c431d8ec251adffcf42e8c118c67626abd24d2e169db0064a"
SEASON_ORDER=("summer","winter")


def make_panel(original:pd.DataFrame, date_enrichment:pd.DataFrame,
               new_dates:bool) -> pd.DataFrame:
    require(original.photo_id.is_unique and len(original)==MOIS_ROWS,
            "Original immutable measured photographic population drift")
    all_dated=original.copy()
    if new_dates:
        allowed=date_enrichment.loc[
            date_enrichment.provenance_status.eq(
                "source_photo_and_observation_verified_with_current_public_date"),
            ["original_photo_id","retrieved_observed_on"]].rename(
                columns={"original_photo_id":"photo_id",
                         "retrieved_observed_on":"provider_verified_day"})
        all_dated=all_dated.merge(allowed,on="photo_id",how="left",validate="one_to_one")
        require(not (all_dated.observed_on.notna() &
                     all_dated.provider_verified_day.notna()).any(),
                "Do not overwrite previously recorded observation days")
        all_dated["observed_on"]=all_dated.observed_on.where(
            all_dated.observed_on.notna(),all_dated.provider_verified_day)
    return labelled_calendar(all_dated)


def spatial_support(p:pd.DataFrame, *, width:float, same_year:bool,
                    minimum:int, zone:str) -> dict:
    sub=p[p.classified&p.date.notna()&p.zone.eq(zone)&
          p.season.isin(SEASON_ORDER)&p.latitude.between(-90,90)&
          p.longitude.between(-180,180)]
    if not set(SEASON_ORDER).issubset(set(sub.season)):
        return {"matched_groups":0,"matched_taxa":0,"matched_photos":0,
                "matched_original_observations":0}
    return eligible_site_photo_support(sub,width,same_year,minimum)


def audit(original:pd.DataFrame, rec:pd.DataFrame) -> dict:
    require(len(rec)==9600 and rec.original_photo_id.is_unique,
            "Blind original 9600 sample changed")
    before=make_panel(original,rec,False)
    after=make_panel(original,rec,True)
    b=before[before.classified & before.date.notna()]
    a=after[after.classified & after.date.notna()]
    hist=a[a.old_mois_or_new_fcp.eq("historical")]
    old_tax=set(b.inat_taxon_id)
    hist_tax=set(hist.inat_taxon_id)
    after_tax=set(a.inat_taxon_id)
    sizes=hist.groupby("inat_taxon_id",observed=True).size()
    require(len(b)==81867 and len(a)==86172 and len(hist)==4305,
            "Source-locked historic observed-on and hue denominators changed")
    require(len(after_tax)==len(old_tax | hist_tax),
            "All-original photo IDs duplicated or source taxonomy changed")
    results={
        "before_actual_photo_dated_and_classifiable":int(len(b)),
        "after_actual_photo_dated_and_classifiable":int(len(a)),
        "additional_classified_original_historical_hues":int(len(hist)),
        "original_date_present_source_photo_attempts":int(before.date.notna().sum()),
        "date_enriched_source_photo_attempts":int(after.date.notna().sum()),
        "before_dated_photographed_taxa":int(len(old_tax)),
        "after_dated_photographed_taxa":int(len(after_tax)),
        "additional_dated_photographed_taxa":int(len(after_tax-old_tax)),
        "historical_sample_dated_hue_taxa":int(len(hist_tax)),
        "historical_sample_taxa_already_found_in_old_dated_source":int(len(old_tax & hist_tax)),
        "historical_sample_taxa_with_exactly_one_classified_dated_photo":int(sizes.eq(1).sum()),
        "historical_sample_singleton_taxon_fraction":round(float(sizes.eq(1).mean()),6),
        "support":[],
    }
    for zone in ("north_extratropical","south_extratropical"):
        for years in (False,True):
            for n in (1,2):
                for label,panel in (("before",before),("after",after),
                                    ("historical_date_recovery_only",hist)):
                    support=spatial_support(panel,width=.25,same_year=years,
                                            minimum=n,zone=zone)
                    results["support"].append({
                        "cohort":label,"hemisphere":zone,"site_cell_size_degrees":.25,
                        "matched_same_local_meteorological_season_year":years,
                        "minimum_original_observation_IDs_in_summer_and_winter":n,
                        **support,
                    })
    return results


def run(photozip:Path,dateszip:Path,output:Path) -> dict:
    require(not output.exists() or not any(output.iterdir()),
            "Never overwrite a frozen source-photo scientific conclusion")
    original,proof,digest=read_signed(photozip,PHOTO_DATA,PHOTO_RECEIPT,
        "source_photo_level_csv_sha256",["photo_id","observation_id",
        "inat_taxon_id","morph","observed_on","old_mois_or_new_fcp",
        "latitude","longitude"])
    require(digest==PHOTO_SHA and len(original)==MOIS_ROWS and
            int(original.morph.isin(HUES).sum())==MOIS_HUES,
            "Original archived measured photos or hues changed")
    with zipfile.ZipFile(dateszip) as z:
        receipt=json.loads(z.read(DATE_RECEIPT))
        raw=z.read(DATE_INPUT)
    require(receipt["new_original_photo_date_CSV_sha256"]==DATE_SHA and
            hashlib.sha256(raw).hexdigest()==DATE_SHA and
            receipt["exact_original_photo_observation_and_current_date_matches"]==9597 and
            receipt["original_full_source_photo_sha256"]==PHOTO_SHA,
            "Real 9597/9600 dated original photo/observation source drifted")
    dates=pd.read_csv(io.BytesIO(raw),dtype=str)
    result=audit(original,dates)
    result.update({
        "schema":"mois_actual_historical9600_taxonomic_breadth_vs_local_season_sample_support_v1",
        "originating_original_photo_source_sha256":PHOTO_SHA,
        "originating_verified_date_source_sha256":DATE_SHA,
        "original_full_genuine_measured_photo_count":MOIS_ROWS,
        "classification_not_recomputed":True,
        "same_taxon_site_source_year_not_population_identity":True,
        "original_photo_colour_not_necessarily_botanically_confirmed_open_petal":True,
        "observer_independence_is_not_guaranteed_by_distinct_observation_ID":True,
        "seasonal_source_siteyear_overlap_not_biological_disappearance":True,
        "photo_white_month_normals_not_observed_year_climate_response":True,
        "FCP_within_species_W_D_local_colour_inference_not_imported":True,
    })
    output.mkdir(parents=True,exist_ok=True)
    path=output/"mois_original9600_date_taxonomic_gain_vs_spatial_temporal_positivity.json"
    path.write_text(json.dumps(result,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    print("MOIS_PHOTOS_TAXONOMY_VS_LOCAL_TEMPORAL_SUPPORT",
          json.dumps({k:result[k] for k in (
              "additional_classified_original_historical_hues",
              "additional_dated_photographed_taxa",
              "historical_sample_singleton_taxon_fraction")},sort_keys=True),flush=True)
    return result


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--original-photo-archive",type=Path,required=True)
    p.add_argument("--recovered-historical9600-archive",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    run(a.original_photo_archive,a.recovered_historical9600_archive,a.output)
