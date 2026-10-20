#!/usr/bin/env python3
"""Audit *measured* original photos in the COMPLETE p0008 wave.

Each row is an original photograph, not an independent biological individual.
Unknown/failed photo colours stay unknown. Geographic/seasonal results are
conditional on the original iNaturalist opportunity frame and frozen ROI-v4.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

HUES=("white","yellow_orange","red_pink","blue_purple")
SOURCE_SHA="0038024ee40adfbc71aa427a112948ec5ce53cb53a7aae9879c1b32bc06eabde"
SOURCE_ARM="global93m_onepass_target_tier_02"
EXPECTED_N=36_375
PROBABILITY=0.0008


def require(ok, message):
    if not bool(ok):
        raise ValueError(message)


def sha256(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for part in iter(lambda:f.read(1<<20),b""):
            h.update(part)
    return h.hexdigest()


def latitude_band(value):
    try: latitude=float(value)
    except (TypeError,ValueError): return "unknown"
    if not -90<=latitude<=90: return "unknown"
    if abs(latitude)<=23.44: return "tropical"
    return "south_extratropical" if latitude<0 else "north_extratropical"


def describe(frame):
    attempted=len(frame)
    classified=int(frame["classified"].sum())
    failures=int(frame["acquisition_failed"].sum())
    white=int(frame["white"].sum())
    return {
        "photo_attempts":attempted,
        "acquired_images":attempted-failures,
        "acquisition_failed":failures,
        "classified_four_hue":classified,
        "acquired_unclassifiable":attempted-failures-classified,
        "white_classified":white,
        "classifiable_per_attempt":classified/attempted if attempted else None,
        "white_fraction_conditional_classified":white/classified if classified else None,
    }


def summarize(photos):
    expected={"photo_id","observation_id","inat_taxon_id","observed_on","latitude",
              "longitude","morph","measurement_status","classified_four_hue",
              "source_arm","selection_probability"}
    require(expected<=set(photos.columns),"Actual original photo fields missing")
    ids=photos.photo_id.astype(str)
    require(ids.str.fullmatch(r"\d+").all() and ids.is_unique,
            "Duplicate or malformed original photo ID")
    obs=photos.observation_id.fillna("").astype(str)
    taxa=photos.inat_taxon_id.fillna("").astype(str)
    require(obs.str.len().gt(0).all() and taxa.str.len().gt(0).all(),
            "Missing original observation or taxon identity")
    require(photos.source_arm.eq(SOURCE_ARM).all(),
            "Unexpected source arm: probability and enrichment arms must not mix")
    probability=pd.to_numeric(photos.selection_probability,errors="coerce")
    require(probability.notna().all() and probability.eq(PROBABILITY).all(),
            "Original observation sampling probability drift")
    require(photos.morph.isin((*HUES,"mixed_uncertain")).all(),
            "Unrecognized frozen FCP ROI flower hue")
    classified=photos.morph.isin(HUES)
    failed=photos.measurement_status.eq("image_acquisition_failed")
    require(classified.eq(photos.measurement_status.eq("classified_four_state_morph")).all(),
            "Original ROI classification status is inconsistent with photo hue")
    recorded=photos.classified_four_hue.astype(str).str.lower().map({
        "true":True,"false":False})
    require(recorded.notna().all() and classified.eq(recorded).all(),
            "Original ROI photo classification flag differs from hue")
    require(not (classified & failed).any(),
            "An original image acquisition failure cannot be a classified flower")
    day=pd.to_datetime(photos.observed_on,errors="coerce")
    require(day.notna().all(),"Original observation date missing")
    lat=pd.to_numeric(photos.latitude,errors="coerce")
    lon=pd.to_numeric(photos.longitude,errors="coerce")
    require(lat.between(-90,90).all() and lon.between(-180,180).all(),
            "Original geographic photo point invalid or missing")
    frame=pd.DataFrame({
        "photo_id":ids,
        "observation_id":obs,
        "inat_taxon_id":taxa,
        "month":day.dt.month.astype(int),
        "latitude_band":lat.map(latitude_band),
        "morph":photos.morph.astype(str),
        "classified":classified.to_numpy(),
        "acquisition_failed":failed.to_numpy(),
        "white":photos.morph.eq("white").to_numpy(),
    })
    overall=describe(frame)
    region=[dict(latitude_band=x,**describe(g))
            for x,g in frame.groupby("latitude_band",sort=True)]
    months=[dict(month=int(x),**describe(g))
            for x,g in frame.groupby("month",sort=True)]
    cells=[dict(latitude_band=str(k[0]),month=int(k[1]),**describe(g))
           for k,g in frame.groupby(["latitude_band","month"],sort=True)]
    region_taxa=[{
        "latitude_band":str(x),
        "distinct_attempted_taxon_ids":int(g.inat_taxon_id.nunique()),
        "taxon_ids_with_classified_photos":int(
            g.loc[g.classified,"inat_taxon_id"].nunique()),
        "distinct_observations":int(g.observation_id.nunique()),
        "classified_observations":int(
            g.loc[g.classified,"observation_id"].nunique()),
    } for x,g in frame.groupby("latitude_band",sort=True)]
    photo_per_observation=frame.groupby("observation_id",sort=False).size()
    classified_taxon_counts=frame.loc[frame.classified].groupby("inat_taxon_id").size()
    obs_class=frame.groupby("observation_id",sort=False).agg(
        classified=("classified","sum"),
        white=("white","sum"))
    viable_obs=obs_class.loc[obs_class.classified.gt(0)]
    observation_equal_white=(
        float((viable_obs.white/viable_obs.classified).mean())
        if len(viable_obs) else None)
    require(sum(x["photo_attempts"] for x in region)==len(frame) and
            sum(x["classified_four_hue"] for x in region)==overall["classified_four_hue"] and
            sum(x["photo_attempts"] for x in months)==len(frame) and
            sum(x["photo_attempts"] for x in cells)==len(frame),
            "Original region/month/photo output accounting changed")
    return {
        "schema":"mois_p0008_complete_original_photo_pixel_coverage_v1",
        "complete_photo_attempts":len(frame),
        "overall":overall,
        "hue_terminals":frame.morph.value_counts().sort_index().astype(int).to_dict(),
        "distinct_observations":int(photo_per_observation.size),
        "multi_photo_observations":int(photo_per_observation.gt(1).sum()),
        "photos_in_multi_photo_observations":int(
            photo_per_observation.loc[photo_per_observation.gt(1)].sum()),
        "distinct_photo_taxon_ids":int(frame.inat_taxon_id.nunique()),
        "taxon_ids_with_classified_photos":int(classified_taxon_counts.size),
        "taxon_ids_with_at_least_5_classifications":int(
            classified_taxon_counts.ge(5).sum()),
        "taxon_ids_with_at_least_30_classifications":int(
            classified_taxon_counts.ge(30).sum()),
        "classified_observations":int(len(viable_obs)),
        "white_fraction_equal_weight_over_classifiable_observations":observation_equal_white,
        "date_min":day.min().strftime("%Y-%m-%d"),
        "date_max":day.max().strftime("%Y-%m-%d"),
        "latitude_bands":region,
        "months":months,
        "latitude_by_month":cells,
        "taxon_observation_coverage_by_latitude":region_taxa,
        "region_is_original_latitude_not_botanical_bioregion":True,
        "original_photo_pixels_reacquired_or_reclassified":False,
        "globally_representative_white_fraction_claimed":False,
        "original_unknowns_not_assumed_nonflowering":True,
        "native_vs_introduced_individual_provenance_not_inferred":True,
        "geographic_enrichment_arm_included":False,
    }


def run(photos,receipt,source,output):
    photos,receipt,source,output=map(Path,(photos,receipt,source,output))
    require(not output.exists() or not any(output.iterdir()),
            "Refuse to overwrite a published original-photo audit")
    r=json.loads(receipt.read_text())
    s=json.loads(source.read_text())
    require(r.get("schema")=="mois_global93m_observation_baseline_real_RGB_photo_results_v1"
            and r.get("attempted")==EXPECTED_N
            and r.get("source_plan_sha256")==SOURCE_SHA
            and r.get("technical_unknowns_not_assumed_nonwhite") is True
            and r.get("ecological_analysis_run") is False
            and sha256(photos)==r.get("measured_photo_parquet_sha256")
            and sha256(source)==SOURCE_SHA,
            "Complete genuine p0008 original photo source SHA/receipt mismatch")
    require(s.get("schema")=="mois_full_source_observation_sampling_cohort_v1"
            and s.get("full_frame_photo_ids")==93_209_351
            and s.get("new_photo_attempts_planned")==EXPECTED_N
            and s.get("colour_measured") is False,
            "Unexpected signed original observation photo roster")
    frame=pd.read_parquet(photos)
    require(len(frame)==EXPECTED_N,"Full p0008 original photos not present")
    report=summarize(frame)
    q=report["overall"]
    require(q["photo_attempts"]==r["attempted"] and
            q["acquired_images"]==r["downloaded"] and
            q["classified_four_hue"]==r["classifiable"] and
            q["acquisition_failed"]==r["acquisition_failed"] and
            q["acquired_unclassifiable"]==r["unclassifiable_acquired"],
            "Complete original photo actual ROI counts changed from signed receipt")
    output.mkdir(parents=True,exist_ok=True)
    for key,name in [
        ("latitude_bands","original_photo_latitude_band.csv"),
        ("months","original_photo_observation_month.csv"),
        ("latitude_by_month","original_photo_latitude_by_month.csv"),
        ("taxon_observation_coverage_by_latitude","original_taxon_coverage_by_latitude.csv"),
    ]:
        pd.DataFrame(report[key]).to_csv(output/name,index=False)
    report["source_photo_data_sha256"]=sha256(photos)
    report["original_receipt_sha256"]=sha256(receipt)
    report["source_manifest_sha256"]=sha256(source)
    report["output_sha256"]={
        f.name:sha256(f) for f in sorted(output.glob("*.csv"))}
    (output/"p0008_real_photo_coverage_receipt.json").write_text(
        json.dumps(report,sort_keys=True,indent=2,allow_nan=False)+"\n")
    print("MOIS_P0008_COMPLETE_ACTUAL_PHOTO_LATITUDE_SEASON_RESULTS",
          json.dumps({
              "photo_attempts":q["photo_attempts"],
              "classified":q["classified_four_hue"],
              "observations":report["distinct_observations"],
              "taxa_with_classified_photos":report["taxon_ids_with_classified_photos"],
              "latitude_cells":len(report["latitude_bands"]),
          },sort_keys=True),flush=True)
    return report


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ("photos","receipt","source","output"):
        parser.add_argument("--"+key,required=True,type=Path)
    args=parser.parse_args()
    run(args.photos,args.receipt,args.source,args.output)
