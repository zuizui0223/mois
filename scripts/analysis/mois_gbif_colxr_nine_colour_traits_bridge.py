#!/usr/bin/env python3
"""Merge GBIF COL XR accepted species keys with original nine-colour MOIS traits.

This is a species-level trait overlay for GBIF occurrence grids. Distinguishes
single-hue photographed support from observed multicolour taxa, does NOT
pretend colour was measured at occurrence positions or genetically fixed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd

PALETTE=("white","yellow","orange","red","pink","magenta","purple","blue","bronze")
SCHEMA="mois_26852_colxr_exact_species_original9_continuous_trait_link_v1"
ATLAS_SHA="e1438342218275078620cdd503d3a4b5e774eaa7c6cce0c160e8969eab875c52"
MATCH_STATUS="colxr_exact_accepted_species"
OBSERVED_SINGLE=("provisional_single","sampling_consensus_single")
OBSERVED_MULTI=("multicolour_candidate","sampled_multicolour")


def require(ok, reason):
    if not bool(ok):raise ValueError("MOIS_COLXR_NINE_TRAITS: "+reason)

def digest(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda:stream.read(1<<20),b""):h.update(block)
    return h.hexdigest()


def build_traits(atlas,matches,*,verify_counts=True):
    require({"inat_taxon_id","species","hue_status","dominant_hue9",
             "observed_hue9_set","after_strict9","after_photo_attempts",
             "independent_roi9_observations",
             *("observed_mean_"+h for h in PALETTE)}<=set(atlas.columns),
            "Authentic 9-colour photo trait fields absent")
    require({"inat_taxon_id","species","source_hue9","hue_status",
             "key","COL_XR_checklistKey","status"}<=set(matches.columns),
            "COL XR original source accepted species identity fields absent")
    atlas=atlas.copy();matches=matches.copy()
    for d in (atlas,matches):
        d["inat_taxon_id"]=d.inat_taxon_id.astype(str).str.strip()
        require(d.inat_taxon_id.is_unique and d.inat_taxon_id.ne("").all(),
                "Repeated or missing iNaturalist source-taxon ID")
    require(len(atlas)==42111,"42,111 nominal source-taxon original hue register truncated")
    require(matches.status.eq(MATCH_STATUS).all() and matches.key.fillna("").astype(str).str.match(r"^[A-Z0-9]+$").all()
            and matches.key.is_unique and matches.COL_XR_checklistKey.nunique()==1
            and matches.COL_XR_checklistKey.iloc[0]=="7ddf754f-d193-4cc9-b351-99906754a03b",
            "Source GBIF keys are not exact accepted nonduplicated COL XR species")
    columns=["inat_taxon_id","species","hue_status","dominant_hue9",
             "observed_hue9_set","after_strict9","after_photo_attempts",
             "independent_roi9_observations","distinct_palette_votes",
             "ambiguous_observation_rois",
             *("observed_mean_"+h for h in PALETTE)]
    combined=matches[["inat_taxon_id","species","source_hue9","hue_status",
                      "key","COL_XR_checklistKey"]].merge(
            atlas[columns],on="inat_taxon_id",validate="one_to_one",
            how="inner",suffixes=("_match","_source"))
    require(len(combined)==len(matches),"GBIF exact-match species dropped from original atlas")
    require(combined.species_match.equals(combined.species_source)
            and combined.hue_status_match.equals(combined.hue_status_source)
            and combined.source_hue9.equals(combined.dominant_hue9),
            "Original named species or hue changed between GBIF and real ROI9 atlas")
    fraccols=["observed_mean_"+h for h in PALETTE]
    fractions=combined[fraccols].apply(pd.to_numeric,errors="coerce").to_numpy(float)
    valid=np.isfinite(fractions).all(axis=1)&(fractions>=0).all(axis=1)&(
        np.abs(fractions.sum(axis=1)-1)<0.002)
    require(valid.all(),"COL XR accepted taxa cannot have missing or unnormalised nine-palette ROI")
    status=combined.hue_status_source
    require(status.isin((*OBSERVED_SINGLE,*OBSERVED_MULTI)).all(),
            "Unexpected observed hue evidence state in named accepted species")
    combined=combined.rename(columns={
        "key":"gbif_colxr_accepted_species_key",
        "COL_XR_checklistKey":"gbif_colxr_checklist_key",
        "species_source":"photo_source_species_name",
        "hue_status_source":"photo_evidence_status",
        "dominant_hue9":"photo_modal_palette9",
        "observed_hue9_set":"photo_observed_palette9_states",
        "source_hue9":"gbif_linked_photo_palette9"
    })
    combined=combined.drop(columns=["species_match","hue_status_match"])
    combined["hue_evidence_class"]=np.where(combined.photo_evidence_status.isin(OBSERVED_SINGLE),
                                             "sampled_single_hue","sampled_multicolour")
    combined["all_range_hue_fixed_proven"]=False
    combined["hue_measured_at_gbif_occurrence"]=False
    combined["photo_fractions_are_physical_pigment_spectra"]=False
    if verify_counts:
        require(len(combined)==26852 and
                (combined.hue_evidence_class=="sampled_single_hue").sum()==14435 and
                (combined.hue_evidence_class=="sampled_multicolour").sum()==12417,
                "Already authenticated global 26852 COLXR species/observed multi profile drift")
    return combined.sort_values("gbif_colxr_accepted_species_key").reset_index(drop=True)


def run(atlas_file,accepted_file,output,*,verify_sha=True):
    output=Path(output)
    require(not output.exists(),"Do not overwrite immutable COL XR nine-colour linked traits")
    if verify_sha:
        require(digest(atlas_file)==ATLAS_SHA,"Real 42111 global 9-colour atlas SHA changed")
    atlas=pd.read_csv(atlas_file,compression="infer",dtype={"inat_taxon_id":str})
    matches=pd.read_csv(accepted_file,compression="infer",dtype=str)
    traits=build_traits(atlas,matches)
    output.mkdir(parents=True)
    paths={}
    for label,sub in (
            ("all",traits),
            ("observed_single",traits.loc[traits.hue_evidence_class.eq("sampled_single_hue")]),
            ("observed_multicolour",traits.loc[traits.hue_evidence_class.eq("sampled_multicolour")])):
        dest=output/f"mois_colxr_species9_{label}_photo_traits.csv.gz"
        sub.to_csv(dest,index=False,compression="gzip",float_format="%.7f")
        paths[dest.name]=digest(dest)
    grouped=traits.groupby(["hue_evidence_class","photo_modal_palette9"],observed=True).agg(
        accepted_species=("gbif_colxr_accepted_species_key","size"),
        authentic_original_ROI9_count=("after_strict9","sum")
    ).reset_index()
    groups=output/"mois_colxr_species9_trait_evidence_by_hue.csv"
    grouped.to_csv(groups,index=False)
    paths[groups.name]=digest(groups)
    report={"schema":SCHEMA,
            "original_nominal_42111_taxa":len(atlas),
            "GBIF_COLXR_exact_accepted_species_with_numeric9":len(traits),
            "sampled_single_hue_species":int(traits.hue_evidence_class.eq("sampled_single_hue").sum()),
            "sampled_multicolour_species":int(traits.hue_evidence_class.eq("sampled_multicolour").sum()),
            "nine_colour_components_preserved":list(PALETTE),
            "nine_colour_classes_present":int(traits.photo_modal_palette9.nunique()),
            "source_atlas_sha256":digest(atlas_file),
            "source_colxr_accepted_sha256":digest(accepted_file),
            "output_sha256":paths,
            "GBIF_actual_occurrences_joined_in_this_step":0,
            "photographic_hue_species_wide_fixed_proven":False,
            "photographic_roi9_true_pigment_proven":False}
    (output/"receipt.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print("MOIS_COLXR_GBIF_SPECIES_NINE_TRAIT_LINK",json.dumps(report,sort_keys=True),flush=True)
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--atlas",required=True);p.add_argument("--accepted",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    run(a.atlas,a.accepted,a.output)


if __name__=="__main__":
    main()
