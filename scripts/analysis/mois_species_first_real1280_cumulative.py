#!/usr/bin/env python3
"""Append actual second-wave 1024 original-photo strict9 results and select new species.

Reads authenticated original-photo worker terminals, never downloads or
re-infers any old image. Prioritizes never-targeted zero-ROI9 taxa over deeper
repeat sampling. Output is a source-backed photographed-colour proxy only.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from scripts.analysis.mois_species_first_palette9 import (
    PALETTE, COLUMNS, build_catalogue, class_labels, sha256,
)
from scripts.analysis.mois_species_first_incremental_first256 import (
    require, parse_boolean, SOURCE_PHOTO_SHA, SOURCE_COVERAGE_SHA,
    SOURCE_CENSUS_SHA, FIRST_ATTEMPTS, BASE_PHOTOS, BASE_STRICT,
)

SCHEMA = "mois_species_first_real1280_exact_original_photo9_append_v1"
FIRST_DELTA_SHA = "4eb884b12322573143277858f8f774def3eeecc7c4ecf212d8d2a67a93d0199a"
WAVE2_SOURCE_RUN = 42913236517
WAVE2_TOTAL = 1024
WAVE2_STRICT = 363
HUES = {"white", "yellow_orange", "red_pink", "blue_purple"}


def _read_path(base, name):
    path = Path(base) / name
    if path.exists():
        return path
    alternate = Path(base) / "measured" / name
    require(alternate.exists(), "Missing original wave2 file: "+name)
    return alternate


def read_wave2(source_dir, actual_dir):
    """Verify signed 4×256 original worker receipts and derived numeric overlay."""
    source_dir, actual_dir = Path(source_dir), Path(actual_dir)
    manifest = json.loads((source_dir/"manifest.json").read_text())
    require(manifest.get("schema")=="mois_full_source_observation_sampling_cohort_v1"
            and manifest.get("colour_measured") is False
            and manifest.get("work_units")==4
            and manifest.get("new_photo_attempts_planned")==WAVE2_TOTAL
            and manifest.get("one_independent_observation_per_taxon") is True,
            "Second wave photo-source manifest is not authentic")
    summary_path = _read_path(actual_dir, "mois_next1024_actual_roi9_receipt.json")
    summary = json.loads(summary_path.read_text())
    require(summary.get("schema")=="mois_species_first_real1024_wave2_actual_photo_roi9_v1"
            and summary.get("new_distinct_original_photo_attempts")==WAVE2_TOTAL
            and summary.get("genuine_strict9_numeric_success")==WAVE2_STRICT
            and summary.get("new_technical_unknowns")==WAVE2_TOTAL-WAVE2_STRICT
            and summary.get("duplicate_original_photo_or_observation_or_taxon") is False,
            "Completed authentic original pixel wave2 run signature drift")
    original = []
    planned = []
    for shard in range(4):
        src = source_dir/"shards"/f"shard-{shard:03d}.csv"
        receipt = _read_path(actual_dir, f"shard_{shard:03d}/receipt.json")
        photos = _read_path(actual_dir, f"shard_{shard:03d}/photos.csv")
        require(sha256(src)==manifest["photo_shard_sha256"][str(shard)]
                and sha256(receipt)==summary["all_1024_real_terminal_receipts"][str(shard)],
                "Original frozen source SHA or genuine wave2 receipt drift")
        rec=json.loads(receipt.read_text())
        require(rec.get("schema")=="mois_full_original93m_real_RGB_photo_unit_v1"
                and rec.get("attempted")==256
                and rec.get("shard")==shard
                and rec.get("source_sha256")==sha256(src)
                and rec.get("photos_sha256")==sha256(photos)
                and rec.get("photo_colours_measured") is True,
                "Actual original pixel worker receipt mismatched the photo source")
        p = pd.read_csv(src,dtype=str,keep_default_na=False)
        d = pd.read_csv(photos,dtype=str,keep_default_na=False)
        require(len(p)==len(d)==256 and
                p.photo_id.is_unique and p.inat_taxon_id.is_unique and
                p.observation_id.is_unique and d.photo_id.is_unique
                and set(p.photo_id)==set(d.photo_id),
                "Original group second wave returned duplicate/missing photos")
        for key in ("photo_id","inat_taxon_id","observation_id"):
            a=p.set_index("photo_id")[key] if key!="photo_id" else pd.Index(p.photo_id)
            b=d.set_index("photo_id")[key] if key!="photo_id" else pd.Index(d.photo_id)
            if key!="photo_id":
                require(a.sort_index().equals(b.sort_index()),
                        "Original source taxon or observation ID changed")
        planned.append(p)
        original.append(d)

    src_all=pd.concat(planned,ignore_index=True)
    raw=pd.concat(original,ignore_index=True)
    require(len(src_all)==len(raw)==WAVE2_TOTAL and
            raw.photo_id.is_unique and raw.observation_id.is_unique and
            raw.inat_taxon_id.is_unique and src_all.photo_id.is_unique and
            src_all.observation_id.is_unique and src_all.inat_taxon_id.is_unique,
            "1024 source photo/observation/taxon independence violated")
    data_file = _read_path(actual_dir,"mois_next1024_genuine_original9_and_unknowns.csv.gz")
    derivative=pd.read_csv(data_file,dtype=str,keep_default_na=False)
    require(len(derivative)==WAVE2_TOTAL and derivative.photo_id.is_unique and
            set(raw.photo_id)==set(derivative.photo_id),
            "Wave2 derived strict9 photo set differs from signed 4 workers")
    orig=raw.set_index("photo_id").sort_index()
    der=derivative.set_index("photo_id").sort_index()
    untouched=["inat_taxon_id","observation_id","morph",
               "measurement_status","image_sha256"]
    require(orig[untouched].equals(der[untouched]),
            "Derived 9-colour labels changed original observation/morph/status/photo hash")
    pix=["flower_fraction_"+x for x in PALETTE]
    numeric=orig[pix].apply(pd.to_numeric,errors="coerce").to_numpy(float)
    valid = (np.isfinite(numeric).all(axis=1) &
             (numeric >= 0).all(axis=1) &
             (numeric <= 1.00001).all(axis=1) &
             (np.abs(np.nan_to_num(numeric).sum(axis=1)-1) <= 0.001))
    strict=orig.morph.isin(HUES).to_numpy(bool) & valid
    flags=parse_boolean(der.strict9_original_photo_numeric).to_numpy(bool)
    require(np.array_equal(flags,strict) and int(strict.sum())==WAVE2_STRICT,
            "Derived strict9 status disagrees with real original image ROI pixels")
    require((orig.loc[der.index[flags],pix].to_numpy() ==
             der.loc[flags,pix].to_numpy()).all(),
            "Numerical nine-palette values were modified after signed real ROI")
    require(der.loc[~flags,pix].eq("").all().all(),
            "Technical nonresponse must not be imputed by an apparent colour")
    labels=np.full(WAVE2_TOTAL,"",dtype=object)
    labels[strict]=class_labels(numeric[strict])
    require((labels==der.nine_palette_dominant_label.to_numpy()).all()
            and int((labels!="").sum())==
               int(summary["nine_colour_dominant_label_success"]),
            "Dominant palette label not reproducible from genuine source fractions")

    frame=pd.DataFrame({"photo_id":orig.index.to_numpy(),
                        "inat_taxon_id":orig.inat_taxon_id.to_numpy(),
                        "observation_id":orig.observation_id.to_numpy(),
                        "strict9":strict,
                        "measured_palette_dominant_label":labels})
    for i,col in enumerate(COLUMNS):
        frame[col]=np.where(strict,numeric[:,i],np.nan)
    return frame,src_all,summary


def assemble_catalogue(baseline,census,coverage,first,wave2):
    """One immutable photo per ID, observation-level taxon aggregation."""
    photo_cols=["photo_id","inat_taxon_id","observation_id","strict9",*COLUMNS]
    require(len(baseline)==BASE_PHOTOS and baseline.photo_id.is_unique
            and int(baseline.strict9.sum())==BASE_STRICT,
            "True historical original-photo 1.18m numeric strict9 source changed")
    require(len(first)==FIRST_ATTEMPTS and len(wave2)==WAVE2_TOTAL
            and int(first.strict9.sum())==91 and int(wave2.strict9.sum())==WAVE2_STRICT,
            "Verified exact real photo-wave sizes or numeric ROI totals changed")
    three=pd.concat([baseline[photo_cols],first[photo_cols],wave2[photo_cols]],
                    ignore_index=True)
    require(three.photo_id.is_unique and len(three)==BASE_PHOTOS+FIRST_ATTEMPTS+WAVE2_TOTAL,
            "Original photo ID reused across real existing measured waves")
    new_obs=set(first.observation_id)|set(wave2.observation_id)
    require(len(new_obs)==FIRST_ATTEMPTS+WAVE2_TOTAL
            and not new_obs.intersection(
                set(baseline.observation_id.dropna().astype(str))-{""}),
            "A new source observation UUID was previously measured")
    source=coverage.copy()
    source.inat_taxon_id=source.inat_taxon_id.astype(str)
    source=source.set_index("inat_taxon_id")
    novel=pd.concat([first[photo_cols],wave2[photo_cols]],ignore_index=True)
    counts=novel.groupby("inat_taxon_id",observed=True).agg(
        attempts=("photo_id","size"),success=("strict9","sum"))
    require(counts.index.isin(source.index).all(),"New photo taxon absent from signed census")
    source.loc[counts.index,"after_photo_attempts"] += counts.attempts
    source.loc[counts.index,"after_strict9"] += counts.success
    new_coverage=source.reset_index()
    atlas=build_catalogue(three,census,new_coverage)
    require(int(atlas.after_strict9.sum())<=int(three.strict9.sum()),
            "Original photo strict9 totals inconsistent with nominal source-taxon counts")
    return atlas,new_coverage,three


def choose_third_wave(catalogue,attempts,global_options,*,quota=1024):
    """Breadth-first: 1st new zero9 taxon photos before any unresolved retries."""
    require(0<quota<=4096,"Invalid new wave size")
    candidates=global_options.copy()
    for key in ("photo_id","observation_id","inat_taxon_id"):
        candidates[key]=candidates[key].fillna("").astype(str).str.strip()
    require(candidates.photo_id.is_unique and
            candidates.photo_id.str.fullmatch(r"\d+").all(),
            "Original prospective photo IDs are duplicate or malformed")
    # Some historical source IDs are not UUIDs; exclude all nonempty previous
    # observation identifiers to avoid falsely counting repeat photographs.
    used_photo=set(attempts.photo_id.astype(str))
    used_obs=set(attempts.observation_id.fillna("").astype(str))-{""}
    d=candidates.loc[
        ~candidates.photo_id.isin(used_photo) &
        candidates.observation_id.ne("") &
        ~candidates.observation_id.isin(used_obs)].copy()
    d=d.merge(catalogue[["inat_taxon_id","hue_status","dominant_hue9",
                         "after_strict9","after_photo_attempts"]],
              on="inat_taxon_id",how="inner",validate="many_to_one")
    # Stop all species with a clear 9-category label, even if single observation.
    # This is a resource-allocation stopping rule, not fixation evidence.
    d=d.loc[d.dominant_hue9.eq("") & d.hue_status.isin(
        ["no_strict9","mixed_image_unresolved"])].copy()
    original_focus=set(attempts.iloc[BASE_PHOTOS:].inat_taxon_id)
    d["acquisition_priority"]=np.where(
        d.inat_taxon_id.isin(original_focus),"retry_unclear","first_unmeasured")
    d["fallback_order"]=pd.to_numeric(d.source_observation_rank,errors="raise")
    d["hash_species"]=d.inat_taxon_id.map(
        lambda x:hashlib.sha256(("MOIS-SPECIES-FIRST-THIRD-WAVE|"+x).encode()).hexdigest())
    d=d.sort_values(["acquisition_priority","hash_species","fallback_order"],
                    kind="stable")
    d=d.drop_duplicates("inat_taxon_id").head(quota).copy()
    require(len(d)==quota and d.photo_id.is_unique and
            d.observation_id.is_unique and d.inat_taxon_id.is_unique,
            "Unable to select complete independent third-wave first-photo taxa")
    d["measurement_id"]="MOIS-ZERO9-W3-"+d.photo_id
    d["image_filename"]=d.measurement_id+".jpg"
    d["shard"]=(np.arange(len(d))//256).astype(str)
    return d


def run(photo_path,census_path,coverage_path,first_delta_path,wave2_source_dir,
        wave2_actual_dir,candidate_path,output,*,quota=1024,verify_sha=True):
    output=Path(output)
    require(not output.exists(),"Do not overwrite completed data-only source join")
    if verify_sha:
        require(sha256(photo_path)==SOURCE_PHOTO_SHA and
                sha256(census_path)==SOURCE_CENSUS_SHA and
                sha256(coverage_path)==SOURCE_COVERAGE_SHA and
                sha256(first_delta_path)==FIRST_DELTA_SHA,
                "Frozen original photographic RGB9 source or first delta SHA changed")
    original=pd.read_parquet(photo_path,columns=[
        "photo_id","inat_taxon_id","observation_id","strict9",*COLUMNS])
    first=pd.read_parquet(first_delta_path)
    census=pd.read_csv(census_path,dtype={"inat_taxon_id":str})
    coverage=pd.read_csv(coverage_path,compression="infer",dtype={"inat_taxon_id":str})
    second,source,measurement=read_wave2(wave2_source_dir,wave2_actual_dir)
    atlas,next_coverage,combined=assemble_catalogue(
        original,census,coverage,first,second)
    opts=pd.read_csv(candidate_path,compression="infer",dtype=str)
    queue=choose_third_wave(atlas,combined,opts,quota=quota)
    output.mkdir(parents=True)
    delta_file=output/"mois_real1024_wave2_original9_photo_delta.parquet"
    atlas_file=output/"mois_global_42111_species_first_after1280.csv.gz"
    queue_file=output/"mois_wave3_original1024_new_species_first_queue.csv"
    coverage_file=output/"mois_42111_original_colour_coverage_after1280.csv.gz"
    second.to_parquet(delta_file,index=False)
    atlas.to_csv(atlas_file,index=False,compression="gzip",float_format="%.6f")
    queue.to_csv(queue_file,index=False)
    next_coverage.to_csv(coverage_file,index=False,compression="gzip")
    gbif=atlas[["inat_taxon_id","species","hue_status","dominant_hue9",
                "observed_hue9_set"]].copy()
    gbif["gbif_accepted_species_key"]=""
    gbif["gbif_name_match_status"]="not_yet_matched"
    gbif.to_csv(output/"mois_42111_species_palette9_GBIF_name_roster.csv.gz",
                index=False,compression="gzip")
    old_taxa=set(catalogue_id for catalogue_id in
                 catalogue_row_ids(first,original))
    new_strict_taxa=set(second.loc[second.strict9,"inat_taxon_id"])
    fresh_taxa=new_strict_taxa-old_taxa
    now_hue=atlas.set_index("inat_taxon_id").loc[second.inat_taxon_id]
    reported={"schema":SCHEMA,
              "prior_signed_wave2_run":WAVE2_SOURCE_RUN,
              "original_photo_attempts_cumulative":len(combined),
              "strict9_original_image_ROIs_cumulative":int(combined.strict9.sum()),
              "source_taxa_total":len(atlas),
              "source_taxa_with_strict9":int(atlas.after_strict9.gt(0).sum()),
              "source_taxa_without_strict9":int(atlas.after_strict9.eq(0).sum()),
              "second_wave_attempted_photos":len(second),
              "second_wave_true_numeric_strict9":int(second.strict9.sum()),
              "second_wave_nine_palette_labelled_photos":int(second.measured_palette_dominant_label.ne("").sum()),
              "second_wave_new_taxa_with_first_strict9":len(fresh_taxa),
              "second_wave_strict9_on_previously_colour_supported_taxa":int(second.strict9.sum())-len(fresh_taxa),
              "second_wave_labelled_taxa_now_with_colour":int(now_hue.dominant_hue9.ne("").sum()),
              "second_wave_technical_unknown":int((~second.strict9).sum()),
              "all_species_with_dominant_nine_palette_label":int(atlas.dominant_hue9.ne("").sum()),
              "next_wave_true_photo_opportunities":len(queue),
              "next_wave_first_unmeasured":int(queue.acquisition_priority.eq("first_unmeasured").sum()),
              "next_wave_retry_unclear":int(queue.acquisition_priority.eq("retry_unclear").sum()),
              "GBIF_occurrences_downloaded":0,
              "species_wide_fixed_colour_proven":False,
              "fresh_original_photo_or_observation_collisions":0}
    reported["outputs_sha256"]={
        p.name:sha256(p) for p in output.iterdir() if p.is_file()}
    (output/"receipt.json").write_text(
        json.dumps(reported,indent=2,sort_keys=True)+"\n")
    print("MOIS_SPECIES_FIRST_REAL1280_CUMULATIVE_AND_WAVE3",
          json.dumps(reported,sort_keys=True),flush=True)
    return reported


def catalogue_row_ids(first,original):
    """IDs already with genuine strict9 before wave2."""
    return pd.concat([original.loc[original.strict9,"inat_taxon_id"],
                      first.loc[first.strict9,"inat_taxon_id"]],ignore_index=True).astype(str).unique()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("photos","census","coverage","first-delta","wave2-source",
                 "wave2-actual","global-options","output"):
        p.add_argument("--"+name,required=True)
    p.add_argument("--quota",type=int,default=1024)
    a=p.parse_args()
    run(a.photos,a.census,a.coverage,a.first_delta,
        a.wave2_source,a.wave2_actual,a.global_options,a.output,quota=a.quota)


if __name__=="__main__":
    main()
