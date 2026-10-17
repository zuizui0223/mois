#!/usr/bin/env python3
"""MOIS 146052 real original date+hue photos: EXACT summer/winter taxon account.

One original photo-id, observation date and frozen 4-hue ROI must be source
verified BEFORE any ecological inference. Each hemisphere's unpaired full
original photo white summer-winter difference is algebraically partitioned
into within-common-taxon conditional photo white rates, changed relative
photo frequencies of those common taxa and the support of taxa photographed
in only one local season. This does NOT identify biological flowering absence,
photographic bias mechanism, or seasonal pigment plasticity within plants.
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

from scripts.analysis.mois_seasonal_267k_climate_normals import (
    HUES, labelled_calendar, require,
)
from scripts.analysis.mois_taxon_season_composition_decomposition import (
    symmetric_photo_decomposition, bootstrap_taxa,
)

INPUT_SOURCE="mois_original302242_exact_observed_on_full_uuid_numeric_verified.csv.gz"
INPUT_RECEIPT="mois_original302242_146052_photo_hue_full_observed_on_provenance_receipt.json"
INPUT_SHA="9c9692238f58cc61561270f7c07b1d782ce711adeade9c58550866f668fb4690"
PHOTO_ATTEMPTS=302242
COLOUR_LABELS=146061
DATED_COLOUR_LABELS=146052
SOURCE_NAMED_TAXA=25638
OUT_JSON="mois_original146052_full_source_summer_winter_taxon_photo_mixture_decision.json"
OUT_CSV="mois_original146052_full_source_seasonal_white_decomposition_figure_data.csv"


def load_signed(artifact_zip:Path)->pd.DataFrame:
    with zipfile.ZipFile(artifact_zip) as z:
        raw=z.read(INPUT_SOURCE)
        receipt=json.loads(z.read(INPUT_RECEIPT))
    require(hashlib.sha256(raw).hexdigest()==INPUT_SHA and
            receipt["final_source_data_csv_SHA256"]==INPUT_SHA and
            receipt["original_302242_actual_measured_photo_pixels"]==PHOTO_ATTEMPTS and
            receipt["final_dated_and_original_four_hue_classifiable_ROIs"]==DATED_COLOUR_LABELS and
            receipt["named_taxa_with_original_day_and_four_hue_ROI"]==SOURCE_NAMED_TAXA,
            "Actual immutable 302242 original photos or verified source dates changed")
    d=pd.read_csv(io.BytesIO(raw),compression="gzip",dtype=str,keep_default_na=False)
    require(len(d)==PHOTO_ATTEMPTS and d.photo_id.is_unique and
            d.morph.isin(HUES).sum()==COLOUR_LABELS and
            ((d.qualified_date.str.len()==10)&d.morph.isin(HUES)).sum()==DATED_COLOUR_LABELS,
            "Cannot extend source 146052 actually dated original photo four-hue census")
    original=d[["photo_id","observation_id","inat_taxon_id","morph",
                "latitude","longitude","old_mois_or_new_fcp","qualified_date"]].rename(
                columns={"qualified_date":"observed_on"}).copy()
    original["latitude"]=pd.to_numeric(original.latitude,errors="coerce")
    original["longitude"]=pd.to_numeric(original.longitude,errors="coerce")
    calendar=labelled_calendar(original)
    return calendar[calendar.classified&calendar.date.notna()&
                    calendar.zone.isin(("north_extratropical","south_extratropical"))&
                    calendar.season.isin(("summer","winter"))].copy()


def decompose_photographic_taxa(sub:pd.DataFrame,*,draws:int=999,seed:int=20261017)->dict:
    """Source photo ID frequencies not unbiased source taxon floral abundances."""
    summer=sub[sub.season.eq("summer")]
    winter=sub[sub.season.eq("winter")]
    require(len(summer)>0 and len(winter)>0 and
            summer.photo_id.is_unique and winter.photo_id.is_unique,
            "Missing actual original classified/date-bearing summer or winter photos")
    taxaS=set(summer.inat_taxon_id)
    taxaW=set(winter.inat_taxon_id)
    shared=taxaS&taxaW
    require(len(shared)>=2, "Insufficient genuine common photographed taxa")
    commonS=summer[summer.inat_taxon_id.isin(shared)]
    commonW=winter[winter.inat_taxon_id.isin(shared)]
    exclusiveS=summer[~summer.inat_taxon_id.isin(shared)]
    exclusiveW=winter[~winter.inat_taxon_id.isin(shared)]
    tabS=commonS.groupby("inat_taxon_id",observed=True).white.agg(
        ["size","sum"]).rename(columns={"size":"ns","sum":"ws"})
    tabW=commonW.groupby("inat_taxon_id",observed=True).white.agg(
        ["size","sum"]).rename(columns={"size":"nw","sum":"ww"})
    tab=tabS.join(tabW,how="inner")
    require(len(tab)==len(shared), "Changed original common-taxon photo source support")
    counts=tab[["ns","nw","ws","ww"]].to_numpy(dtype=float)
    point,limits=bootstrap_taxa(counts,n_draws=draws,seed=seed)
    # The common photographed taxa are renormalized separately within each
    # season before the symmetric Kitagawa mixture/within photo composition.
    within,mix,common,equal=(float(x) for x in point)
    exclusive_S=100*(len(exclusiveS)/len(summer))*(
        float(exclusiveS.white.mean())-float(commonS.white.mean()))
    exclusive_W=-100*(len(exclusiveW)/len(winter))*(
        float(exclusiveW.white.mean())-float(commonW.white.mean()))
    raw=100*(float(summer.white.mean())-float(winter.white.mean()))
    require(abs(raw-within-mix-exclusive_S-exclusive_W)<1e-8 and
            abs(within+mix-common)<1e-8,
            "Original full date/photo census photo mixture composition not additive")
    summary={
       "summer_classified_original_photos":len(summer),
       "winter_classified_original_photos":len(winter),
       "summer_white_labeled_photo_ROIs":int(summer.white.sum()),
       "winter_white_labeled_photo_ROIs":int(winter.white.sum()),
       "summer_original_named_taxa_with_photo_hues":len(taxaS),
       "winter_original_named_taxa_with_photo_hues":len(taxaW),
       "source_named_taxa_common_both_summer_winter":len(shared),
       "original_named_taxa_union_both_seasons":len(taxaS|taxaW),
       "shared_taxon_photographic_Jaccard":round(len(shared)/len(taxaS|taxaW),6),
       "fraction_summer_photos_from_taxa_without_classified_winter_photos":
           round(len(exclusiveS)/len(summer),6),
       "raw_photo_white_summer_minus_winter_pp":round(raw,6),
       "within_common_taxon_photographic_white_rate_component_pp":round(within,6),
       "changed_common_taxon_relative_photo_frequency_component_pp":round(mix,6),
       "summer_photo_exclusive_taxon_support_component_pp":round(exclusive_S,6),
       "winter_photo_exclusive_taxon_support_component_pp":round(exclusive_W,6),
       "common_taxon_original_photographic_white_rate_difference_pp":round(common,6),
       "equal_taxon_mean_common_source_white_photo_difference_pp":round(equal,6),
       "common_taxon_within_white_photo_rate_999_bootstrap95_pp":
           [round(float(x),6) for x in limits[0]],
       "common_taxon_changed_photo_frequency_999_bootstrap95_pp":
           [round(float(x),6) for x in limits[1]],
       "taxon_bootstrap_draws":draws,
       "all_observation_photo_IDs_source_not_independent_photographers":True,
       "source_photos_not_botanically_verified_open_petal_census":True,
       "exact_full_census_mixture_4terms_sum_to_raw":True,
    }
    return summary


def run(artifact_zip:Path,output:Path):
    require(not output.exists() or not any(output.iterdir()),
            "Do not overwrite finished original-photo seasonal mixture result")
    q=load_signed(artifact_zip)
    results={}
    fig=[]
    for zone in ("north_extratropical","south_extratropical"):
        r=decompose_photographic_taxa(q[q.zone.eq(zone)],
            draws=999,seed=20261017)
        results[zone]=r
        for k,label in [
           ("within_common_taxon_photographic_white_rate_component_pp",
             "within_common_taxon_photo_white_share"),
           ("changed_common_taxon_relative_photo_frequency_component_pp",
             "photo_frequency_mixture_of_shared_taxa"),
           ("summer_photo_exclusive_taxon_support_component_pp",
             "summer_only_photographed_taxon_support"),
           ("winter_photo_exclusive_taxon_support_component_pp",
             "winter_only_photographed_taxon_support")]:
            fig.append({"hemisphere":zone,"photo_composition_term":label,
                        "symmetric_decomposition_effect_pp":r[k],
                        "original_photo_white_summer_minus_winter_pp":
                            r["raw_photo_white_summer_minus_winter_pp"]})
    output.mkdir(parents=True,exist_ok=True)
    data_path=output/OUT_CSV
    pd.DataFrame(fig).to_csv(data_path,index=False,lineterminator="\n")
    receipt={
       "schema":"mois_original146052_genuine_observed_on_photo_colour_taxon_turnover_v1",
       "original_real_photo_measured_pixels":PHOTO_ATTEMPTS,
       "original_actual_four_hue_photo_ROIs":COLOUR_LABELS,
       "original_real_date_and_hue_photo_ROIs":DATED_COLOUR_LABELS,
       "source_photo_full_SHA256":INPUT_SHA,
       "source_observation_day_ID_and_original_hue_not_reassigned":True,
       "exact_full_census_source_composition_not_causal_selection":True,
       "zero_original_photo_opportunities_without_measured_colour_counted":True,
       "hemispheres":results,
       "figure_source_data_SHA256":hashlib.sha256(data_path.read_bytes()).hexdigest(),
       "limits":[
         "Photo taxa absent from one photographed season are not biologically absent or nonflowering",
         "Taxa present both seasons may occur in different regions/years/observers, not identical plant populations",
         "Automatic ROI colour may be flower/petal/bract/fruit misclassification; no independent flowering organ audit",
         "Photographed taxon frequency is not population or flowering abundance",
         "Whole-taxon bootstrap is exploratory percentile interval, not geographically independent or familywise corrected",
         "No observed-year solar/weather, pollinator, pigment genetic cost, or FCP within-species W/D fitness tested"
       ],
    }
    (output/OUT_JSON).write_text(json.dumps(receipt,indent=2,
                               ensure_ascii=False,sort_keys=True)+"\n")
    print("MOIS_REAL146052_PHOTO_TAXON_SEASON_COMPOSITION",
          json.dumps({h:{k:v[k] for k in (
             "raw_photo_white_summer_minus_winter_pp",
             "source_named_taxa_common_both_summer_winter",
             "changed_common_taxon_relative_photo_frequency_component_pp",
             "within_common_taxon_photographic_white_rate_component_pp")}
             for h,v in results.items()},sort_keys=True),flush=True)
    return receipt


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source-verified302242-day-photo-artifact",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    a=ap.parse_args()
    run(a.source_verified302242_day_photo_artifact,a.output)
