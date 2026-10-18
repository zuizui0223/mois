#!/usr/bin/env python3
"""Audit the FIRST 256 genuine new-genus original iNaturalist RGB photo outcomes.

This is an immutable PHOTO-ATTEMPT/YIELD checkpoint, not the main worldwide
floral hue prevalence. Count real new species and genera only when original
four-hue ROI is classifiable, preserve technical unknown and acquisition
failure statuses distinctly, compare exact source/photo/model receipts.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

import pandas as pd

from scripts.acquisition import mois_global_angiosperm_breadth_source as broad
from scripts.analysis import mois_integrate_broad_angiosperm_real4096 as verified

OUT_CSV="mois_broad_first256_real_photo_ROI_hue_taxon_quality.csv"
OUT_JSON="mois_broad_first256_photo_RBG_novel_taxa_genera_QC_receipt.json"
PILOT=256


def require(ok,why):
    if not ok:
        raise ValueError(why)


def checksum(path:Path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for item in iter(lambda:f.read(1<<20),b""):
            h.update(item)
    return h.hexdigest()


def pilot_metrics(old,planned,measured):
    """Use source-matched original photo IDs, not colour selected for sampling."""
    require(len(planned)==len(measured)==PILOT and
            planned.measurement_id.is_unique and
            measured.measurement_id.is_unique and
            planned.photo_id.is_unique and
            measured.photo_id.is_unique,
            "Pilot is not exactly 256 unique signed original photo attempts")
    known=(measured.morph.isin(verified.HUES))
    require(known.eq(measured.classified_four_hue.eq("true")).all(),
            "An unknown original ROI has been recoded as a coloured flower")
    joined=planned.merge(
        measured[["measurement_id","photo_id","morph","measurement_status",
                  "classified_four_hue"]],
        on="measurement_id",how="inner",validate="one_to_one",
        suffixes=("_source","_measured"))
    require(len(joined)==PILOT and
            joined.photo_id_source.eq(joined.photo_id_measured).all() and
            not joined.photo_id_source.isin(old.photo_id).any(),
            "Pilot photo IDs do not match blinded original world source")
    genuine=joined[joined.morph.isin(verified.HUES)].copy()
    oldtaxa=set(old.inat_taxon_id.astype(str))
    oldgenera=set(old.species.astype(str).str.split().str[0])
    newtaxa=set(genuine.inat_taxon_id.astype(str))-oldtaxa
    newgenera=set(genuine.genus.astype(str))-oldgenera
    hue_counts={c:int(genuine.morph.eq(c).sum()) for c in verified.HUES}
    unknown=int(PILOT-len(genuine))
    require(sum(hue_counts.values())+unknown==PILOT,
            "Colour-unknown photographed attempts were counted as chromatic")
    status={str(k):int(v) for k,v in
            joined.measurement_status.value_counts().items()}
    metric={
        "actual_original_photo_attempts":PILOT,
        "new_photo_hue_classifications":int(len(genuine)),
        "technical_unknown_or_acquisition_failed_attempts":unknown,
        "automatic_four_hue_classification_success_rate":
            float(len(genuine)/PILOT),
        "visible_hue_photo_counts":hue_counts,
        "newly_classified_named_taxa_not_in_main_original_source":
            len(newtaxa),
        "newly_classified_genera_not_in_main_original_source":
            len(newgenera),
        "new_taxa_per_1000_original_RGB_photo_attempts":
            float(1000*len(newtaxa)/PILOT),
        "new_genera_per_1000_original_RGB_photo_attempts":
            float(1000*len(newgenera)/PILOT),
        "source_original_photo_species":int(planned.inat_taxon_id.nunique()),
        "source_original_photo_genera":int(planned.genus.nunique()),
        "original_source_photo_10degree_cells":int(planned.grid10.nunique()),
        "new_colour_classified_10degree_cells":int(genuine.grid10.nunique()),
        "source_prioritisation_stage_counts":{
            str(k):int(v) for k,v in
            planned.source_sampling_priority.value_counts().items()},
        "original_ROI_measurement_status_counts":status,
        "genuinely_classifiable_four_hues_not_same_as_all_source_species":True,
        "photo_acquisition_failure_is_distinct_from_unknown_visible_hue":True,
        "all_candidate_source_picks_before_RBG_photo_colour_extraction":True,
        "not_estimated_adaptive_white_chromatic_equilibrium_or_global_fitness":True,
        "not_botanically_verified_open_flower_petals":True,
        "not_independent_of_iNaturalist_taxonomic_photo_selection":True,
    }
    return joined,metric


def run(prior:Path,plan:Path,pilot:Path,out:Path):
    require(not out.exists() or not any(out.iterdir()),
            "Refuse to overwrite previous source-verified pilot result")
    old,base=verified.verify_old(prior)
    source,freezerec=verified.verified_blind_plan(plan)
    planned=source.loc[source.shard.eq("0")].reset_index(drop=True)
    require(len(planned)==PILOT,
            "Pilot source is not frozen original first256")
    measured,counts,terminal=verified.one_worker(0,planned,pilot,freezerec)
    assert sum(counts.values())==PILOT
    joined,metrics=pilot_metrics(old,planned,measured)
    require(metrics["new_photo_hue_classifications"]==
            terminal["four_hue_classified"] and
            metrics["technical_unknown_or_acquisition_failed_attempts"]==
            terminal["unclassified"],
            "Image source ROI terminal counts differ from inspected pilot")
    out.mkdir(parents=True,exist_ok=True)
    filename=out/OUT_CSV
    joined.to_csv(filename,index=False,lineterminator="\n")
    receipt={
       "schema":"mois_broad_angiosperm_first256_actual_photo_pixel_ROI_taxonomic_yield_v1",
       "main_prior_processed_original_photo_attempts":len(old),
       "source_snapshot_date":freezerec["source_snapshot_date"],
       "source_original_first4096_photo_metadata_SHA256":
           freezerec["first_wave_original_4096_csv_sha256"],
       "source_frozen_pilot256_original_photo_SHA256":
           freezerec["per_shard_sha256"]["0"],
       "source_pilot256_actual_original_ROI_terminal_SHA256":
           terminal["original_new_photo_hue_terminal_SHA256"],
       "source_prior_main_photo_table_SHA256":
           base["files_sha256"][broad.CORE_PHOTO],
       "actual_pilot_metrics":metrics,
       "original_photo_ID_taxon_licence_site_and_model_attestation_verified":True,
       "photo_color_unknowns_not_assumed_chromatic":True,
       "all_four_hues_model_inference_NOT_chemical_anthocyanin":True,
       "not_updating_main_313535_floral_white_fraction":True,
       "photo_source_vs_colour_selection_checked_for_taxon_gain":True,
       "source_outcomes_csv_SHA256":checksum(filename),
    }
    (out/OUT_JSON).write_text(json.dumps(
        receipt,indent=2,sort_keys=True,ensure_ascii=False,allow_nan=False)+"\n")
    print("MOIS_BROAD_ORIGINAL_FIRST256_ACTUAL_ROI_NEW_SPECIES_GENUS_YIELD",
          json.dumps(metrics,sort_keys=True),flush=True)
    return receipt


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--original313535-verified-Parquet",type=Path,required=True)
    ap.add_argument("--frozen-broad4096-source",type=Path,required=True)
    ap.add_argument("--signed-pilot256-original-ROI",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    q=ap.parse_args()
    run(q.original313535_verified_Parquet,
        q.frozen_broad4096_source,q.signed_pilot256_original_ROI,q.output)
