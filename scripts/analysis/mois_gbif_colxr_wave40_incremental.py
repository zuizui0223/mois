#!/usr/bin/env python3
"""GBIF COL XR exact name match for ONLY newly nine-hue-labelled MOIS wave40 taxa.

The 37,057 prior GBIF accepted species keys and their exact nine observed
original-photo fractions are reused from an immutable receipt. This is a
species-level trait/accepted-name relation, NOT a site-colour measurement.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from scripts.analysis.mois_gbif_colxr_new_hue_incremental import (
    classify_new_rows, combine_names
)
from scripts.analysis.mois_gbif_colxr_batch_match import source_sha
from scripts.analysis.mois_gbif_colxr_nine_palette_name_pilot import (
    PALETTE, NAME_RULE, COL_XR, require
)
from scripts.analysis.mois_gbif_colxr_nine_colour_traits_bridge import build_traits
from scripts.analysis.mois_gbif_colxr_sql_isea3h_requests import make_requests

SCHEMA="mois_species_first_real_wave40_gbif_colxr_new_names_v1"
PRIOR_ACCEPTED=37057
PRIOR_HUE_LABELS=39124
AFTER_WAVE40_HUE_LABELS=39127
SOURCE_TAXA=42111


def select_incremental_hue_names(photo_atlas,old_crosswalk):
    new=photo_atlas.copy()
    for field in ("inat_taxon_id","species","dominant_hue9","hue_status"):
        new[field]=new[field].fillna("").astype(str).str.strip()
    require(len(new)==SOURCE_TAXA and new.inat_taxon_id.is_unique,
            "Updated original-source 42,111 nominal taxa not unique")
    require(int(new.dominant_hue9.isin(PALETTE).sum())==AFTER_WAVE40_HUE_LABELS,
            "Genuine wave40 3 newly labelled nine-hue taxa absent")
    old=old_crosswalk.copy()
    for field in ("inat_taxon_id","species","source_hue9","hue_status",
                  "status","key","COL_XR_checklistKey"):
        old[field]=old[field].fillna("").astype(str).str.strip()
    require(old.inat_taxon_id.is_unique and
            int(old.status.eq("colxr_exact_accepted_species").sum())==PRIOR_ACCEPTED and
            old.COL_XR_checklistKey.eq(COL_XR).all(),
            "Prior exact GBIF species bridge changed")
    require(old.loc[old.status.eq("colxr_exact_accepted_species"),"key"].is_unique,
            "Prior accepted GBIF keys contain duplicate species")
    lookup=new.set_index("inat_taxon_id").reindex(old.inat_taxon_id)
    require(lookup.species.notna().all() and
            (lookup.species.to_numpy()==old.species.to_numpy()).all() and
            (lookup.dominant_hue9.to_numpy()==old.source_hue9.to_numpy()).all() and
            (lookup.hue_status.to_numpy()==old.hue_status.to_numpy()).all(),
            "GBIF old name/photographed hue labels drifted with wave40 append")
    newnames=new.loc[new.dominant_hue9.isin(PALETTE) &
                     new.species.str.match(NAME_RULE) &
                     ~new.inat_taxon_id.isin(old.inat_taxon_id)].copy()
    require(0<len(newnames)<=AFTER_WAVE40_HUE_LABELS-PRIOR_HUE_LABELS,
            "Unexpectedly many source names newly labelled from wave40")
    return newnames.sort_values("inat_taxon_id").reset_index(drop=True)


def run(atlas,prior_names,prior_receipt,append_receipt,output,*,batch_size=200,pause=1.0):
    output=Path(output)
    require(not output.exists(),"Do not overwrite prior signed wave40 COL XR name bridge")
    require(source_sha(atlas)==json.loads(Path(append_receipt).read_text())["output_sha256"]["atlas"],"Signed wave40 original-photo 9-fraction atlas SHA changed")
    appended=json.loads(Path(append_receipt).read_text())
    require(appended.get("schema")=="mois_species_first_fortieth4096_append_fortyfirst4096_v1"
            and appended.get("fortieth_actual_independent_original_photo_attempts")==27
            and appended.get("fortieth_actual_strict9_photo_successes")==7
            and appended.get("fortieth_new_source_taxa_with_clear_dominant_hue9")==3
            and appended.get("source_taxa_with_clear_sampled_nine_palette")==AFTER_WAVE40_HUE_LABELS
            and appended.get("prior_photo_or_observation_ID_reused") is False
            and appended.get("output_sha256",{}).get("atlas")==source_sha(atlas),
            "Wave40 original nine-fraction image atlas SHA not authenticated by signed append")
    old=pd.read_csv(prior_names,compression="infer",dtype=str,keep_default_na=False)
    receipt=json.loads(Path(prior_receipt).read_text())
    require(receipt.get("schema")=="mois_species_first_real_wave39_gbif_colxr_new_names_v1"
            and receipt.get("current_exact_colxr_accepted_species")==PRIOR_ACCEPTED
            and receipt.get("signed_file_SHA256",{}).get(Path(prior_names).name)==
                source_sha(prior_names),
            "Prior GBIF wave37 accepted-name SHA is not authenticated")
    df=pd.read_csv(atlas,compression="infer",dtype={"inat_taxon_id":str})
    to_match=select_incremental_hue_names(df,old)
    response=classify_new_rows(to_match,batch_size=batch_size,pause=pause)
    combined,accepted=combine_names(old,response)
    traits=build_traits(df,accepted,verify_counts=False)
    require(len(accepted)==len(traits)>=PRIOR_ACCEPTED and
            len(combined)==len(old)+len(to_match),
            "Genuine nine-fraction photo traits or old accepted key bridge lost")
    output.mkdir(parents=True)
    files={}
    for name,frame in [
        ("wave40_only_new_label_gbif_colxr_matches.csv.gz",response),
        ("all_current_photo9_gbif_colxr_names.csv.gz",combined),
        ("all_current_accepted_gbif_colxr_species.csv.gz",accepted),
        ("all_current_accepted_gbif_colxr_nine_fractions.csv.gz",traits),
    ]:
        path=output/name
        frame.to_csv(path,index=False,compression="gzip",float_format="%.7f")
        files[name]=source_sha(path)
    requests=make_requests(combined,chunk_size=250)
    sql=output/"gbif_isea3h8_sql_unsubmitted";sql.mkdir()
    for i,r in enumerate(requests):
        filename=sql/f"colxr_photo9_r8_{i:03d}.json"
        filename.write_text(json.dumps(r["GBIF_download_request"],indent=2)+"\n")
        files[str(filename.relative_to(output))]=source_sha(filename)
    info={
        "schema":SCHEMA,
        "photo_source_taxa":SOURCE_TAXA,
        "photo_source_taxa_dominant_hue":AFTER_WAVE40_HUE_LABELS,
        "previous_gbif_accepted_species":PRIOR_ACCEPTED,
        "new_nine_hue_nominal_taxa_queried":len(to_match),
        "new_exact_colxr_accepted_species":int(
            response.status.eq("colxr_exact_accepted_species").sum()),
        "current_exact_colxr_accepted_species":len(accepted),
        "current_GBIF_species_numeric9_fractional_traits":len(traits),
        "unresolved_or_nonaccepted_names":len(combined)-len(accepted),
        "new_match_status_counts":{
            str(k):int(v) for k,v in response.status.value_counts().items()},
        "GBIF_ISEA3H_SQL_request_count":len(requests),
        "GBIF_authenticated_bulk_occurrence_download_submitted":False,
        "GBIF_occurrence_records_downloaded":0,
        "photo_hue_fixed_at_GBIF_occurrences_proven":False,
        "source_photo_atlas_sha256":source_sha(atlas),
        "prior_matching_sha256":source_sha(prior_names),
        "signed_file_SHA256":files,
    }
    (output/"receipt.json").write_text(json.dumps(info,sort_keys=True,indent=2)+"\n")
    print("MOIS_WAVE40_NEW_GBIF_COLXR_TRUE9_NAME_MATCH",json.dumps(info,sort_keys=True),flush=True)
    return info


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("atlas","prior-names","prior-receipt","append-receipt","output"):
        p.add_argument("--"+k,required=True)
    p.add_argument("--batch-size",type=int,default=200)
    p.add_argument("--pause",type=float,default=1.0)
    a=p.parse_args()
    run(a.atlas,a.prior_names,a.prior_receipt,a.append_receipt,a.output,
        batch_size=a.batch_size,pause=a.pause)


if __name__=="__main__":
    main()
