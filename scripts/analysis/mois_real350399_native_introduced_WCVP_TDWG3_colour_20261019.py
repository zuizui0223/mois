#!/usr/bin/env python3
"""Photographic WHITE-labelled native-vs-introduced *regional-range* source analysis.

MOIS 350,399 original RGB-photo attempts were previously classified by the
same frozen original pixel ROI model. This script NEVER determines whether a
photographed individual is wild/native on an island. It compares species that
Kew WCVP calls native/introduced *within the broad TDWG level-3 region* of
the original photographed site, using genuinely classifiable original ROI
hues, unique observation IDs, and accepted WCVP species IDs.

The geographically/genus-matched source support is a selected DIFFERENT-species
comparison, not evidence that introduction changes phenotype or phenotype
fitness. Presence/classifiability and photography effort are not standardized.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

from scripts.analysis.mois_actual_photo_cohort_season_audit_20261019 import (
    HUES, require, sha256,
)
# Keep the frozen source filenames local: no GeoPandas dependency is needed
# to analyse the already authenticated TDWG3/photo-location Parquet.
PHOTO_OUT="mois350399_original_photo_TDWG3_WCVP_region_source_status.parquet"
AUDIT_OUT="mois350399_TDWG3_WCVP_geography_source_receipt.json"
WCVP_DIST="mois350399_wcvp_accepted_native_distribution_TDWG3.csv.gz"
WCVP_RECEIPT="mois350399_wcvp_TDWG3_source_only_taxon_native_receipt.json"

OUT_JSON="mois350399_region_native_introduced_original_photo_white_diagnostic.json"
OUT_RAW="mois350399_region_range_photo_classifiability_and_white.csv"
OUT_MATCH="mois350399_native_vs_introduced_matched_genus_TDWG3.csv"
OUT_HUES="mois350399_native_introduced_original_photo_four_visible_hues.csv"
FULL=350399
KNOWN=164297
REGIONAL=("native_at_TDWG3","introduced_at_TDWG3")
COHORTS=("core","new_species_first","all_source_arms")
GENUS=re.compile(r"^[A-Z][A-Za-z-]+$")


def source(photo_parquet,photo_receipt,wcvp_distribution,wcvp_receipt):
    photo_rec=json.loads(Path(photo_receipt).read_text())
    ranges_rec=json.loads(Path(wcvp_receipt).read_text())
    require(photo_rec.get("schema")==
            "mois_original350399_photo_taxon_WCVP_regional_TDWG3_geoposition_v1"
            and photo_rec.get("original_actual_photo_attempts")==FULL
            and photo_rec.get("four_hue_classifiable_original_ROIs")==KNOWN
            and photo_rec.get("WCVP_native_distributions_are_only_TDWG_level3") is True
            and photo_rec.get("H3_native_island_adaptation_tested") is False
            and sha256(photo_parquet)==photo_rec["files_sha256"][PHOTO_OUT],
            "Actual source photo and region status Parquet are not authenticated")
    require(ranges_rec.get("schema")==
            "mois350399_kew_wcvp_official_accepted_species_TDWG3_native_source_v1"
            and sha256(wcvp_distribution)==ranges_rec["data_files_sha256"][WCVP_DIST]
            and sha256(wcvp_receipt)==photo_rec["WCVP_source_receipt_sha256"],
            "Original official Kew WCVP accepted species and region source SHA mismatch")
    data=pd.read_parquet(photo_parquet)
    require(len(data)==FULL and data.photo_id.is_unique
            and int(data.classified.sum())==KNOWN and
            data.morph.isin((*HUES,"mixed_uncertain","")).all(),
            "ROI original four-hue photo/source denominator drifted")
    require(data.classified.eq(data.morph.isin(HUES)).all() and
            data.white.eq(data.morph.eq("white")).all(),
            "WCVP regional status must never relabel unknown photograph hues")
    ranges=pd.read_csv(wcvp_distribution,dtype=str,keep_default_na=False)
    require(not ranges.duplicated(["wcvp_accepted_plant_name_id","tdwg_level3_code"]).any(),
            "Regional native WCVP source contains unresolved duplicate status rows")
    allowed=ranges.loc[ranges.native_introduced_state.isin(REGIONAL),
                       ["wcvp_accepted_plant_name_id","tdwg_level3_code",
                        "wcvp_accepted_name","native_introduced_state"]].copy()
    allowed["accepted_genus"]=allowed.wcvp_accepted_name.str.split().str[0]
    # Accepted names with hybrid signs or unsupported genus spellings cannot
    # be forced into a nominal genus cluster; retain original source photo
    # evidence and leave that genus missing in the *matched* comparison.
    accepted_genus_ok=allowed.accepted_genus.map(
        lambda v:isinstance(v,str) and bool(GENUS.fullmatch(v)))
    allowed.loc[~accepted_genus_ok,"accepted_genus"]=""
    accepted=data.merge(
        allowed,
        on=["wcvp_accepted_plant_name_id","tdwg_level3_code",
            "native_introduced_state"],
        how="inner",sort=False,validate="many_to_one"
    )
    require(accepted.taxon_TDWG3_range_annotation.isin(
        ("WCVP_native_at_TDWG3","WCVP_introduced_at_TDWG3")).all(),
        "Native/introduced source status is not aligned to original TDWG3 annotation")
    require(len(accepted)==
            int(data.taxon_TDWG3_range_annotation.isin(
                ("WCVP_native_at_TDWG3","WCVP_introduced_at_TDWG3")).sum()),
            "A valid accepted/native or introduced photographed photo was lost")
    accepted["source_cohort"]=np.where(accepted.arm.eq("core"),"core","new_species_first")
    return data,accepted,photo_rec,ranges_rec


def photo_status_table(full):
    rows=[]
    for category, frame in full.groupby("taxon_TDWG3_range_annotation",
                                        observed=True,dropna=False,sort=True):
        classified=frame.loc[frame.classified]
        row={
            "region_species_status":str(category),
            "original_source_photo_attempts":len(frame),
            "classifiable_photo_ROIs":len(classified),
            "unknown_or_unclassifiable_photo_ROIs":len(frame)-len(classified),
            "classification_yield":len(classified)/len(frame) if len(frame) else None,
            "white_classified_ROIs":int(classified.white.sum()),
            "white_fraction_among_classified":
                float(classified.white.mean()) if len(classified) else None,
            "distinct_source_observations":int(frame.observation_id.nunique()),
            "source_named_taxa":int(frame.inat_taxon_id.nunique()),
            "source_cohort_core_attempts":int(frame.arm.eq("core").sum())
        }
        row.update({h+"_photographic_ROI":int(classified.morph.eq(h).sum()) for h in HUES})
        rows.append(row)
    result=pd.DataFrame(rows)
    require(int(result.original_source_photo_attempts.sum())==FULL
            and int(result.classifiable_photo_ROIs.sum())==KNOWN,
            "Photo WCVP region status partition must retain full original image denominator")
    return result


def distinct_original_observations(frame):
    """Observation equality precedes species within-place composition."""
    known=frame.loc[frame.classified & frame.observation_id.astype(str).ne("")].copy()
    keys=["observation_id","wcvp_accepted_plant_name_id","tdwg_level3_code",
          "native_introduced_state","accepted_genus"]
    # Two images of the same observation are one encounter, while an
    # observation occurring across source cohorts is quarantined as mixed.
    obs=known.groupby(keys,observed=True,sort=False).agg(
        white_label_fraction=("white","mean"),
        photos=("photo_id","size"),
        n_arms=("source_cohort","nunique"),
        cohort_first=("source_cohort","first")
    ).reset_index()
    obs["cohort"]=np.where(obs.n_arms.eq(1),obs.cohort_first,"mixed_source_arm")
    require(int(obs.photos.sum())==len(known) and
            obs.white_label_fraction.between(0,1).all(),
            "Original classified photo ROIs were lost or have nonprobabilistic colour")
    # A photographic observation may include two original taxon assignments;
    # such conflicting original metadata cannot count as two independent
    # encounters in the taxon/species comparison.
    conflicts=obs.observation_id.value_counts()
    repeated=set(conflicts.index[conflicts.gt(1)])
    obs["conflicting_observation_identity"]=obs.observation_id.isin(repeated)
    return obs


def species_region_table(obs,cohort):
    z=obs.loc[~obs.conflicting_observation_identity &
              obs.accepted_genus.map(
                  lambda v: isinstance(v,str) and bool(GENUS.fullmatch(v)))].copy()
    if cohort!="all_source_arms":
        z=z.loc[z.cohort.eq(cohort)]
    keys=["tdwg_level3_code","accepted_genus",
          "native_introduced_state","wcvp_accepted_plant_name_id"]
    result=z.groupby(keys,observed=True,sort=False).agg(
        p_white=("white_label_fraction","mean"),
        observation_count=("observation_id","size"),
        photo_count=("photos","sum")
    ).reset_index()
    require(result.p_white.between(0,1).all()
            and not result.duplicated(keys).any(),
            "One accepted species × region contributes more than one taxon weight")
    return result


def matched_native_introduced(species,*,min_species=1,draws=999,seed=20261019):
    """Within accepted GENUS and same TDWG3, but DIFFERENT accepted species."""
    require(min_species in (1,2) and draws>0,"Unsupported original-photo matching design")
    gkeys=["tdwg_level3_code","accepted_genus"]
    group=species.groupby([*gkeys,"native_introduced_state"],observed=True).agg(
        mean_white=("p_white","mean"),
        species=("wcvp_accepted_plant_name_id","size"),
        observations=("observation_count","sum"),
        actual_photos=("photo_count","sum"),
    ).reset_index()
    native=group.loc[group.native_introduced_state.eq(REGIONAL[0])].set_index(gkeys)
    intro=group.loc[group.native_introduced_state.eq(REGIONAL[1])].set_index(gkeys)
    shared=native.join(intro,how="inner",lsuffix="_native",rsuffix="_introduced")
    shared=shared.loc[shared.species_native.ge(min_species)&
                      shared.species_introduced.ge(min_species)].copy()
    denoms={state:{
        "accepted_species_TDWG3_rows":int(len(species.loc[
            species.native_introduced_state.eq(state)])),
        "original_observations":int(species.loc[
            species.native_introduced_state.eq(state),"observation_count"].sum())
    } for state in REGIONAL}
    n=len(shared)
    base={
        "minimum_accepted_species_each_status_same_genus_TDWG3":min_species,
        "matched_genus_region_strata":int(n),
        "matched_named_genera":int(shared.reset_index().accepted_genus.nunique())
        if n else 0,
        "matched_TDWG3_regions":int(shared.reset_index().tdwg_level3_code.nunique())
        if n else 0,
        "source_eligible_by_status":denoms,
        "matched_accepted_species_by_status":{
            REGIONAL[0]:int(shared.species_native.sum()) if n else 0,
            REGIONAL[1]:int(shared.species_introduced.sum()) if n else 0
        },
        "matched_original_photo_ROIs_by_status":{
            REGIONAL[0]:int(shared.actual_photos_native.sum()) if n else 0,
            REGIONAL[1]:int(shared.actual_photos_introduced.sum()) if n else 0
        },
        "accepted_taxon_equal_not_photo_equal":True,
        "same_genus_same_TDWG3_different_species":True,
        "wild_native_locality_or_island_identity_known":False,
    }
    if n==0:
        return {**base,
                "native_minus_introduced_white_pp_equal_matched_strata":None,
                "native_minus_introduced_white_pp_equal_genus":None,
                "whole_named_genus_bootstrap95_pp":None}
    delta=100*(shared.mean_white_native-shared.mean_white_introduced)
    pergenus=shared.assign(delta_pp=delta).reset_index().groupby(
        "accepted_genus",observed=True).delta_pp.mean()
    n_genus=len(pergenus)
    if n_genus>=10:
        rng=np.random.default_rng(seed+min_species)
        values=pergenus.to_numpy(float)
        samples=rng.integers(0,len(values),(draws,len(values)))
        bootstrap=values[samples].mean(axis=1)
        ci=[float(s) for s in np.quantile(bootstrap,[.025,.975])]
    else:
        ci=None
    return {
        **base,
        "native_minus_introduced_white_pp_equal_matched_strata":float(delta.mean()),
        "native_minus_introduced_white_pp_equal_genus":float(pergenus.mean()),
        "whole_named_genus_bootstrap95_pp":ci,
    }


def run(photo_parquet,photo_receipt,wcvp_distribution,wcvp_receipt,outdir):
    out=Path(outdir)
    require(not out.exists() or not any(out.iterdir()),
            "Never overwrite a completed and source-authenticated WCVP ecological diagnostic")
    all_photos,regional,photo_source,wcvp_source=source(
        photo_parquet,photo_receipt,wcvp_distribution,wcvp_receipt)
    qc=photo_status_table(all_photos)
    obs=distinct_original_observations(regional)
    models=[]
    species_summary=[]
    for cohort in COHORTS:
        tab=species_region_table(obs,cohort)
        for status in REGIONAL:
            z=tab.loc[tab.native_introduced_state.eq(status)]
            species_summary.append({
                "source_cohort":cohort,
                "regional_species_status":status,
                "accepted_species_region_rows":len(z),
                "species_equal_photo_white_fraction":float(z.p_white.mean()) if len(z) else None,
                "classified_original_source_observations":int(z.observation_count.sum()),
                "classifiable_original_source_photos":int(z.photo_count.sum())
            })
        for gate in (1,2):
            models.append({"cohort":cohort,**matched_native_introduced(
                tab,min_species=gate,draws=999)})
    out.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(models).to_csv(out/OUT_MATCH,index=False)
    qc.to_csv(out/OUT_HUES,index=False)
    pd.DataFrame(species_summary).to_csv(out/OUT_RAW,index=False)
    report={
        "schema":"mois_real350399_native_introduced_WCVP_TDWG3_photographic_hue_diagnostic_v1",
        "region_level_only_not_point_native":True,
        "island_specific_native_or_endemism_verified":False,
        "causal_effect_of_introduction_or_pollination_tested":False,
        "original_photo_attempts":len(all_photos),
        "all_original_four_hue_classified_photo_ROIs":int(all_photos.classified.sum()),
        "source_signatures":{
            "original_photo_geo_parquet_sha256":sha256(photo_parquet),
            "original_photo_geo_receipt_sha256":sha256(photo_receipt),
            "official_kew_wcvp_tdwg3_ranges_sha256":sha256(wcvp_distribution),
            "official_kew_wcvp_tdwg3_receipt_sha256":sha256(wcvp_receipt)
        },
        "regional_WCVP_record_photo_attempts":len(regional),
        "regional_WCVP_record_classified_ROIs":int(regional.classified.sum()),
        "regional_classified_original_ROIs_without_verified_accepted_genus":int(
            (regional.classified & regional.accepted_genus.eq("")).sum()),
        "unique_classified_source_observations_aggregated":len(obs),
        "conflicting_original_observation_identity_rows":
            int(obs.conflicting_observation_identity.sum()),
        "observations_with_mixed_photo_acquisition_arms":int(
            obs.cohort.eq("mixed_source_arm").sum()),
        "five_acquisition_arms_pooled_only_for_exploratory_descriptive_analysis":True,
        "species_equal_by_source_arm":species_summary,
        "matched_contrasts":models,
        "limitations":[
            "Native and introduced are accepted WCVP species regional distributions, NOT nativeness of individual photographed plants, local island endemic status, or a causal introduction effect.",
            "A native species photographed in a broad native-range region may be cultivated or introduced at a smaller island/photo locality.",
            "Different source species in a named genus×broad TDWG3 region do not imply matched biological habitats, equal observer effort, equal plant abundance or independent evolutionary transitions.",
            "All true photo RGB machine ROIs have unverified open-flower organs, and unclassifiable original image attempts remain missing, never chromatic.",
            "Photo species-first recruitment, WCVP distribution completeness and visibility/brightness can differ strongly among regions, taxa and source arms.",
            "Source original observation IDs are collapsed prior to taxon comparisons; repeated photographs are not independent observation replicates.",
            "Exploratory whole-nominal-genus bootstrap intervals are NOT spatial- or phylogeny-independent and not corrected for multiple tests.",
            "H3 native-island hypothesis not identified from broad TDWG level-3 metadata; biogeographic context and climate mechanisms remain untested."
        ],
        "files_sha256":{
            file:sha256(out/file) for file in (OUT_MATCH,OUT_HUES,OUT_RAW)
        }
    }
    (out/OUT_JSON).write_text(json.dumps(
        report,sort_keys=True,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("MOIS_REAL350399_NATIVE_TDWG3_INTRODUCED_PHOTO_HUE_EXPLORATORY",
          json.dumps({
              "regional_photo_attempts":len(regional),
              "regional_classifiable_photos":int(regional.classified.sum()),
              "accepted_species_regional_hue":species_summary,
              "matched_genus_region":models
          },sort_keys=True),flush=True)
    return report


def main():
    a=argparse.ArgumentParser(description=__doc__)
    for name in ("photo-parquet","photo-receipt","wcvp-distribution",
                 "wcvp-receipt","output"):
        a.add_argument("--"+name,type=Path,required=True)
    x=a.parse_args()
    run(x.photo_parquet,x.photo_receipt,x.wcvp_distribution,
        x.wcvp_receipt,x.output)


if __name__=="__main__":
    main()
